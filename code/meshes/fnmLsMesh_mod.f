#ifndef SRO_VEL_SCALING
#define SRO_VEL_SCALING 1
#endif

#ifndef SLASH
#ifdef ABQ_WIN86_64
#define SLASH "\"
#else
#define SLASH "/"
#endif
#endif

!DEC$ FREEFORM

module fnmLsMesh_mod
    use iso_fortran_env
    use component_mod
    use floatingNodeLevelSetHex_mod
    use levelsetOptimization_mod
    use vtkLegacy_mod
    use fnmLsAnalysis_mod
    use adjointHex_mod
    implicit none
    
    type, extends(Component) :: FnmLsMesh

        integer                     :: nnodel   = 18
        integer                     :: nstr     = -1

        integer                     :: status   = -1
        integer                     :: iter     = -1
        integer                     :: currentStep = -1
        integer                     :: currentIncrement = -1

        real(real64)                :: dt       = 0.0d0
        real(real64)                :: h        = 0.0d0
        real(real64)                :: vol      = 0.0d0

        type(FnmLsHex), allocatable :: els(:)
        type(AdjHex),   allocatable :: adjels(:)
        type(INPFile)               :: inp
        type(LSOptimizer)           :: opt
        type(FnmLsSettings)         :: settings

        integer,        allocatable :: boundaryNodes(:)

        integer,        allocatable :: files(:)
        integer                     :: resFile = -1
        integer                     :: logFile = -1

        real(real64),   allocatable :: phi(:)
        ! real(real64),   allocatable :: dof(:,:)

        integer,           allocatable :: statusList(:)
        character(len=50), allocatable :: statusLabels(:)
        character(len=15)              :: statusTypes(5) =  &
            [                                               &
                character(len=15) ::                        &
                "sync",                                     &
                "displacement",                             &
                "adjoint",                                  &
                "velocity",                                 &
                "levelset"                                  &
            ]
        character(len=:), allocatable :: dir
        character(len=:), allocatable :: jobname

        
    contains

        procedure :: assembleLS                 => fnmlsmesh_assembLS
        procedure :: assembleVelocity           => fnmlsmesh_assembVel
        procedure :: assembleMeanStrain         => fnmlsmesh_assembMeanStrain
        procedure :: assembleMeanStress         => fnmlsmesh_assembMeanStress
        procedure :: assembleMeanStrainAdjoint  => fnmlsmesh_assembMeanStrainAdjoint
        procedure :: assembleMeanStressAdjoint  => fnmlsmesh_assembMeanStressAdjoint
        procedure :: assembleBoundaryNodes      => fnmlsmesh_assembBoundaryNds
        procedure :: assembleVtkDataArrays      => fnmlsmesh_assembVtkDataArrays

        procedure :: computeLsVelocity          => fnmlsmesh_computeLsVelocity
        procedure :: extrapolateLsVelocity      => fnmlsmesh_extrapolLsVelocity

        procedure :: initialiseLs               => fnmlsmesh_initls
        procedure :: reinitialiseLs             => fnmlsmesh_reinit
        procedure :: updateElmCoords            => fnmlsmesh_updateElmCoords
        procedure :: enforceLsConstraint        => fnmlsmesh_enforceLsConstraint
        procedure :: partition                  => fnmlsmesh_partition
        procedure :: update                     => fnmlsmesh_update
        procedure :: updateElement              => fnmlsmesh_updateElm
        procedure :: updateAnalysisStatus       => fnmlsmesh_updateAnalysisStatus
        procedure :: registerSteps              => fnmlsmesh_registerSteps
        procedure :: getStepType                => fnmlsmesh_getStepType
        
        procedure :: computeDomainVolume        => fnmlsmesh_computeDomainVolume
        procedure :: computeCompliance          => fnmlsmesh_computeCompliance
        procedure :: computeEngRelRate          => fnmlsmesh_computeEngRelRate
        
        procedure :: print                      => fnmlsmesh_print
        procedure :: handleFileOpening          => fnmlsmesh_handleFileOpening
        procedure :: handleFileClosing          => fnmlsmesh_handleFileClosing
        procedure :: writeOutput                => fnmlsmesh_writeOutput

        procedure :: getAnalysisSettings        => fnmlsmesh_getAnalysisSettings
        procedure :: setAnalysisSettings        => fnmlsmesh_setAnalysisSettings

        final :: fnmlsmesh_dtor
    end type
    
    interface FnmLsMesh
        procedure :: fnmlsmesh_ctor
    end interface
    
contains

    
    function fnmlsmesh_ctor(jobname, outdir, logunit) result(this)
        implicit none 

        character(len=*)    :: jobname 
        character(len=*)    :: outdir
        integer             :: logunit
        type(FnmLsMesh)     :: this 

        type(FnmLsHex), allocatable :: auxEls(:)
        type(AdjHex), allocatable :: auxAdjels(:)
        type(Transform3D) :: tls, tshell

        real(real64) :: ogls(3), vol, area
        integer :: iel, auxAdj(1), auxAdjId, nels
        integer, allocatable :: conel(:)

        ! vars init
        ogls        = 0.0d0
        vol         = 0.0d0
        area        = 0.0d0
        iel         = 0
        auxAdj      = 0
        auxAdjId    = 0
        nels        = 0

        call sroLog(" >> Sro::FnmLsMesh << | ctor")

        this % name = "FnmLsMesh"
        this % jobname = jobname

        this % inp = INPFile(jobname, outdir)
        this % dir = outdir

        tls = Transform3D(this % inp % transform_ls, Imat)
        tshell = Transform3D(this % inp % transform_shell, Imat)
        ogls = this % inp % lsorigin

        this % h = this % inp % meshsize

        nels = 0

        allocate(conel(size(this % inp % con, 2)), source = 0)
        allocate(auxEls(this % inp % nuels))
        do iel = 1, this % inp % nuels 
            call assign(conel, this % inp % con(iel,:))
            auxEls(iel) = FnmLsHex(                             &
                this % inp % uels(iel),                         &
                this % inp % x(conel, :),                       &
                conel,                                          & 
                this % inp % props,                             &
                this % inp % orientations,                      &
                this % inp % thicknesses,                       &
                tls,                                            &
                tshell,                                         &
                ogls                                            &
            )
        end do

        allocate(this % els(this % inp % nuels), source = auxEls)
        sro = SroGlobalData()
        call sro % allocateDof(this % inp % nuels, 2*auxEls(1)%ndof + 2*auxEls(1)%nfieldvars)
        call sro % allocateElmap(this % inp % nuels)
        call assign(sro % elmap, this % inp % uels)
        deallocate(auxEls)
        deallocate(conel)


        nels = nels + this % inp % nuels

        this % nstr = this % els(1) % nstr
        this % iter = 0


        if (this % inp % hasDisbond) then
            allocate(conel(this % inp % nndeladjuels), source=0)
            allocate(auxAdjels(this % inp % nadjuels))
            do iel = 1, this % inp % nadjuels
                call assign(conel, this % inp % adjcon(iel,:))
                auxAdjels(iel) = AdjHex(            &
                    this % inp % adjuels(iel),      &
                    this % inp % adjx(conel, :),    &
                    conel,                          &
                    this % inp % k(iel, :, :),      &
                    this % inp % dkda(iel, :, :)    &
                ) 
            end do
            allocate(this % adjels(this % inp % nadjuels), source = auxAdjels)
            deallocate(auxAdjels)
            deallocate(conel)

            call sro % allocateAdjdof(this % inp % nadjuels, 48)
            call sro % allocateAdjelmap(this % inp % nadjuels)
            call assign(sro % adjelmap, this % inp % adjuels)

            nels = nels + this % inp % nadjuels
        end if

        allocate(this % phi(this % inp % nnodesls), source = 0.0d0)

        vol = 0.0d0
        area = 0.0d0
        do iel = 1, size(this % els)
            vol = vol + this % els(iel) % volume 
            area = area + this % els(iel) % getArea()
        end do
        this % opt = LSOptimizer(this % inp % alpha, this % inp % vfrac, area)

        this % logFile = logunit
        open(newunit=this % resFile, file=outdir // "_log" // SLASH // "results.csv")

        write(this % resFile, fmt100) &
            "iter, comp, enrr, vol, vfrac, lambda, lagrMult, penalty, lagr, alpha, volFrac, volInit"

        close(this % resFile)

    end function 

    
    subroutine fnmlsmesh_dtor(this)
        implicit none
    
        type(FnmLsMesh), intent(inout) :: this

        integer :: i

        ! vars init
        i = 0

        call sroLog(" >> Sro::FnmLsMesh << | dtor")
        
        call this % inp % finalize 
        do i = 1, this % inp % nuels 
            call this % els(i) % finalize
        end do
        if (allocated(this % els)) deallocate(this % els)
        if (this % inp % hasDisbond) then
            do i = 1, size(this % adjels)
                call this % adjels(i) % finalize
            end do
        end if
        if (allocated(this % adjels)) deallocate(this % adjels)
        if (allocated(this % files)) then
            call this % handleFileClosing
        end if
        if (allocated(this % phi)) deallocate(this % phi)
        if (allocated(this % boundaryNodes)) deallocate(this % boundaryNodes)
        if (allocated(this % statusList)) deallocate(this % statusList)
        if (allocated(this % statusLabels)) deallocate(this % statusLabels)
    end subroutine

    
    subroutine fnmlsmesh_print(this, varName, unit)
        implicit none
        
        class(FnmLsMesh), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid, iel
        character(len=3) :: cs

        ! vars init
        uid = 0
        iel = 0
        cs = "   "

        call sroLog(" >> Sro::FnmLsMesh << | print")

        uid = 6
        if (present(unit)) then
            uid = unit 
        
            write(uid, fmt100) ""
            write(uid, fmt100) sep2
            write(uid, fmt100) varName // " <" // this % name // ">"
            write(uid, fmt100) ""

            write(uid, fmt101) "status", this % status 
            write(uid, fmt101) "iter", this % iter
            write(uid, fmt100) "" 
            write(uid, fmt102) "dt", this % dt
            write(uid, fmt100) "" 
            call this % opt % print(varName // ".opt", uid)
            if (this % inp % hasDisbond) then
                ! TODO print
                ! do iel = 1, size(this % enrrEls)
                !     write(uid, fmt100) ""
                !     write(cs, "(i3)") iel 
                !     call this % enrrEls(iel) % print(varName // ".enrrEls(" // cs // ")", uid)
                ! end do
            end if
            do iel = 1, size(this % els)
                write(uid, fmt100) ""
                write(cs, "(i3)") iel 
                call this % els(iel) % print(varName // ".els(" // cs // ")", uid)
            end do

        else 
            call this % handleFileOpening
            uid = this % files (1)

            write(uid, fmt100) ""
            write(uid, fmt100) sep2
            write(uid, fmt100) varName // " <" // this % name // ">"
            write(uid, fmt100) ""

            write(uid, fmt101) "status", this % status 
            write(uid, fmt101) "iter", this % iter
            write(uid, fmt100) "" 
            write(uid, fmt102) "dt", this % dt
            write(uid, fmt100) "" 
            write(uid, fmt101) "nuels", this % inp % nuels
            write(uid, fmt101) "nnodes", this % inp % nnodes
            write(uid, fmt100) "" 

            call this % opt % print(varName // ".opt", uid)

            do iel = 1, size(this % els)
                write(uid, fmt100) ""
                write(cs, "(i3)") iel 
                uid = this % files(1 + iel)
                call this % els(iel) % print(varName // ".els(" // cs // ")", uid)
            end do

            if (this % inp % hasDisbond) then
                ! TODO print
                ! do iel = 1, size(this % enrrEls)
                !     write(uid, fmt100) ""
                !     write(cs, "(i3)") iel 
                !     uid = this % files(1 + size(this % els) + iel)
                !     call this % enrrEls(iel) % print(varName // ".enrrEls(" // cs // ")", uid)
                ! end do
            end if
            call this % handleFileClosing
        end if
    end subroutine

    
    subroutine fnmlsmesh_handleFileOpening(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
    
        integer :: i, n, fileid, fileoffset
        character(len=5) :: cs, it
        character(len=:), allocatable :: its, dir

        ! vars init
        i           = 0
        n           = 0
        fileid      = 0
        fileoffset  = 0
        cs          = "     "
        it          = "     "
        its         = ""
        dir         = ""

        call sroLog(" >> Sro::FnmLsMesh << | handleFileOpening")

        if (this % inp % hasDisbond) then
            n = size(this % els) + size(this % adjels) + 1
        else 
            n = size(this % els) + 1
        end if

        write(it, "(i5)") this % iter
        its = trim(adjustl(it))
        call execute_command_line(                                                                      &
            trim(adjustl("mkdir " // achar(34) // this % dir // "_log" // SLASH // its // SLASH // achar(34)))    &
        )

        allocate(this % files(n), source = 0)
        this % files(:) = (/ ( 1000+i, i=1, n ) /)

        fileoffset = 1000
        open(this % files(1), file=this % dir // "_log" // SLASH // its // SLASH // "mesh.txt", status="new")
        fileid = 2
        do i = 1, size(this % els)
            this % files(fileid) = fileoffset + fileid
            write(cs, "(i5)") i
            open(this % files(fileid), file=this % dir // "_log" // SLASH // its // SLASH // "el" // trim(adjustl(cs)) // ".txt", status="new")
            fileid = fileid + 1
        end do
        if (this % inp % hasDisbond) then
            ! TODO open files for print
            ! do i = 1, size(this % enrrEls)
            !     this % files(fileid) = fileoffset + fileid
            !     write(cs, "(i5)") i
            !     open(this % files(fileid), &
            !         file=this % dir // "_log" // SLASH // its // SLASH // "enrrEl" // trim(adjustl(cs)) // ".txt", &
            !         status="new")
            !     fileid = fileid + 1
            ! end do
        end if
    end subroutine

    
    subroutine fnmlsmesh_handleFileClosing(this)
        implicit none 

        class(FnmLsMesh), intent(inout) :: this 

        integer :: ifile 

        ! vars init
        ifile = 0

        call sroLog(" >> Sro::FnmLsMesh << | handleFileClosing")

        if (allocated(this % files)) then
            do ifile = 1, size(this % files)
                close(this % files(ifile))
            end do
            deallocate(this % files)
        end if
    end subroutine

    
    function fnmlsmesh_assembLS(this) result(ls)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
        real(real64) :: ls(this % inp % nnodesls)
    
        integer :: iel, jnd, con(4), lcon(4), aux(1)

        ! vars init
        iel = 0
        con = 0
        lcon = 0
        aux = 0
        ls  = 0.0d0

        call sroLog(" >> Sro::FnmLsMesh << | assembLS")
        
        do iel = 1, size(this % els)
            con = this % els(iel) % getLsCon()
            do jnd = 1, size(con)
                aux = findloc(this % inp % lsnodes, con(jnd))
                if (.not. aux(1) == 0) then
                    if (aux(1) > size(ls)) then
                        call sroXIT(" >>>>>> fnmlsmesh_assembLS <<<<<< Bounds violation")
                    end if
                    lcon(jnd) = aux(1)
                else
                    call sroXIT(" >>>>>> fnmlsmesh_assembLS <<<<<< LS node index not found")
                end if
            end do
            ls(lcon) = this % els(iel) % getLevelset()
        end do
    end function

    
    function fnmlsmesh_assembVel(this, evo) result(v)
        implicit none 

        class(FnmLsMesh), intent(inout) :: this 
        logical, intent(in), optional :: evo
        real(real64) :: v(this % inp % nnodesls)
        logical :: use_evo

        integer :: iel, ind, idx(1), con(4)
        real(real64) :: vel(4), absval

        ! vars init
        iel = 0
        con = 0
        v   = 0.0d0
        absval = 0.0d0

        if (.not. present(evo)) then
            use_evo = .false.
        else
            use_evo = evo
        end if

        call sroLog(" >> Sro::FnmLsMesh << | assembVel")

        do iel=1, size(this % els)
            con = this % els(iel) % getLsCon()
            vel = this % els(iel) % getVelocity()
            do ind = 1, 4
                idx = findloc(this % inp % lsnodes, con(ind))
                if (.not. idx(1) == 0) then
                    if (idx(1) > size(v)) then
                        call sroXIT(" >>>>>> fnmlsmesh_assembVel <<<<<< Bounds violation")
                    end if
                    if (vel(ind) /= vel(ind)) then
                        call sroXIT(" >>>>>> fnmlsmesh_assembVel <<<<<< NaN detected")
                    end if
                    v(idx(1)) = vel(ind)
                end if
            end do
        end do

        if (use_evo .eqv. .true.) then
            this % dt = 1.0d0 * (this % h / maxval(abs(v)))
        else
            absval = maxval(abs(v))
            if (absval > 0.00001d0) then
                v = v / absval
                do iel = 1, size(this % els)
                    call assign(this % els(iel) % levelsetElm % v, this % els(iel) % levelsetElm % v / absval)
                end do
            end if
            this % dt = 0.5d0 * this % h
        end if

        call sroLog(" >> Sro::FnmLsMesh << | assembVel | dt = " // float2str(this % dt))

    end function

    
    function fnmlsmesh_assembMeanStrain(this) result(e)
        implicit none
    
        class(FnmLsMesh) :: this
        real(real64) :: e(this % inp % nnodes, this % nstr)

        integer :: iel, istr, ind
        integer :: con(this % nnodel)
        real(real64) :: ee(this % inp % nnodes, this % nstr)
        integer :: cnt(this % inp % nnodes)

        ! vars init
        e       = 0.0d0
        iel     = 0
        istr    = 0
        ind     = 0
        con     = 0
        ee      = 0.0d0
        cnt     = 0

        call sroLog(" >> Sro::FnmLsMesh << | assembMeanStrain")

        do iel = 1, size(this % els)
            con = this % els(iel) % con
            ee(con, :) = ee(con, :) + this % els(iel) % getMeanStrain()
            cnt(con) = cnt(con) + 1
        end do
        
        do istr = 1, this % nstr
            do ind = 1, this % inp % nnodes
                if (cnt(ind) > 0) then
                    e(ind, istr) = ee(ind, istr) / cnt(ind)
                end if
            end do        
        end do
    end function

    
    function fnmlsmesh_assembMeanStress(this) result(s)
        implicit none
    
        class(FnmLsMesh) :: this
        real(real64) :: s(this % inp % nnodes, this % nstr)
        
        integer :: iel, istr, ind
        integer :: con(this % nnodel)
        real(real64) :: ss(this % inp % nnodes, this % nstr)
        integer :: cnt(this % inp % nnodes)

        ! vars init
        s       = 0.0d0
        iel     = 0
        istr    = 0
        ind     = 0
        con     = 0
        ss      = 0.0d0
        cnt     = 0

        call sroLog(" >> Sro::FnmLsMesh << | assembMeanStress")

        do iel = 1, size(this % els)
            con = this % els(iel) % con
            ss(con,:) = ss(con, :) + this % els(iel) % getMeanStress()
            cnt(con) = cnt(con) + 1
        end do

        do istr = 1, this % nstr 
            do ind = 1, this % inp % nnodes
                if (cnt(ind) > 0) then
                    s(ind, istr) = ss(ind, istr) / cnt(ind)
                end if
            end do 
        end do
    end function

    
    function fnmlsmesh_assembMeanStrainAdjoint(this) result(e)
        implicit none
    
        class(FnmLsMesh) :: this
        real(real64) :: e(this % inp % nnodes, this % nstr)

        integer :: iel, istr, ind
        integer :: con(this % nnodel)
        real(real64) :: ee(this % inp % nnodes, this % nstr)
        integer :: cnt(this % inp % nnodes)

        ! vars init
        e       = 0.0d0
        iel     = 0
        istr    = 0
        ind     = 0
        con     = 0
        ee      = 0.0d0
        cnt     = 0

        call sroLog(" >> Sro::FnmLsMesh << | assembMeanStrainAdjoint")

        do iel = 1, size(this % els)
            con = this % els(iel) % con
            ee(con, :) = ee(con, :) + this % els(iel) % getMeanStrainAdjoint()
            cnt(con) = cnt(con) + 1
        end do

        do istr = 1, this % nstr 
            do ind = 1, this % inp % nnodes
                if (cnt(ind) > 0) then
                    e(ind, istr) = ee(ind, istr) / cnt(ind)
                end if
            end do 
        end do
    end function

    
    function fnmlsmesh_assembMeanStressAdjoint(this) result(s)
        implicit none
    
        class(FnmLsMesh) :: this
        real(real64) :: s(this % inp % nnodes, this % nstr)
        
        integer :: iel, istr, ind
        integer :: con(this % nnodel)
        real(real64) :: ss(this % inp % nnodes, this % nstr)
        integer :: cnt(this % inp % nnodes)

        ! vars init
        s       = 0.0d0
        iel     = 0
        istr    = 0
        ind     = 0
        con     = 0
        ss      = 0.0d0
        cnt     = 0

        call sroLog(" >> Sro::FnmLsMesh << | assembMeanStressAdjoint")

        do iel = 1, size(this % els)
            con = this % els(iel) % con
            ss(con,:) = ss(con, :) + this % els(iel) % getMeanStressAdjoint()
            cnt(con) = cnt(con) + 1
        end do

        do istr = 1, this % nstr 
            do ind = 1, this % inp % nnodes
                if (cnt(ind) > 0) then
                    s(ind, istr) = ss(ind, istr) / cnt(ind)
                end if
            end do 
        end do
    end function

    
    subroutine fnmlsmesh_assembBoundaryNds(this)
        implicit none 

        class(FnmLsMesh), intent(inout) :: this 
        
        integer :: iel, ind, auxId, nbnd, bnd(this % inp % nnodes)

        ! vars init
        iel     = 0
        ind     = 0
        auxId   = 0
        nbnd    = 0
        bnd     = 0

        call sroLog(" >> Sro::FnmLsMesh << | assembBoundaryNds")

        if (allocated(this % boundaryNodes)) deallocate(this % boundaryNodes)

        bnd = -1
        auxId = 1
        nbnd = 0
        do iel = 1, size(this % els)
            if (allocated(this % els(iel) % boundaryNds)) then
                do ind = 1, size(this % els(iel) % boundaryNds)
                    if (.not. any(bnd == this % els(iel) % boundaryNds(ind))) then 
                        bnd(auxId) = this % els(iel) % boundaryNds(ind)
                        auxId = auxId + 1 
                        nbnd = nbnd + 1
                    end if
                end do 
            end if
        end do

        if (nbnd > 0) then
            allocate(this % boundaryNodes(nbnd), source=bnd(1:nbnd))
        end if
    end subroutine

    
    subroutine fnmlsmesh_initls(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this

        integer :: iel, con(4), ind, indl(4), aux(1)

        ! vars init
        iel     = 0
        con     = 0
        ind     = 0
        indl    = 0
        aux     = 0

        call sroLog(" >> Sro::FnmLsMesh << | initls")
    
        do iel = 1, size(this % els)
            con = this % els(iel) % getLsCon()
            do ind = 1, size(con)
                aux = findloc(this % inp % lsnodes, con(ind))
                if (.not. aux(1) == 0) then
                    indl(ind) = aux(1)
                else 
                    call sroXIT(" >>>>>> fnmlsmesh_initls <<<<<< Unable to find node index within lsnodes array ")
                end if
            end do
            call this % els(iel) % updateLs(this % inp % lsinit(indl))
        end do
    end subroutine

    
    subroutine fnmlsmesh_reinit(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
    
        integer :: iel, ind, ibnd, con(4)
        real(real64) :: xel(4, 3), xnd(3), xnd2(3), dmin, daux, phi(4)
        real(real64), allocatable :: xbnd(:,:)

        ! vars init
        iel     = 0
        ind     = 0
        ibnd    = 0
        con     = 0
        xel     = 0.0d0
        xnd     = 0.0d0
        xnd2    = 0.0d0
        dmin    = 0.0d0
        daux    = 0.0d0
        phi     = 0.0d0

        call sroLog(" >> Sro::FnmLsMesh << | reinit")

        if (allocated(this % boundaryNodes)) then
            allocate(xbnd(size(this % boundaryNodes), 3), source=0.0d0)
            call assign(xbnd, this % inp % x(this % boundaryNodes, :))
            do iel = 1, size(this % els)
                con = this % els(iel) % getLsCon() 
                xel = this % inp % x(con, :)
                phi = this % els(iel) % getLevelset()

                do ind = 1, size(con)
                    xnd = xel(ind, :)

                    dmin = 100000.0d0
                    do ibnd = 1, size(xbnd, 1)
                        xnd2 = xbnd(ibnd, :)
                        daux = vec3distance(xnd, xnd2)
                        if (daux < dmin) dmin = daux
                    end do

                    phi(ind) = sign(dmin, phi(ind))
                end do

                call this % els(iel) % updateLs(phi)
            end do
    
            deallocate(xbnd)
        end if
    end subroutine

    
    subroutine fnmlsmesh_updateElmCoords(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
        
        integer :: iel 

        ! vars init 
        iel = 0

        call sroLog(" >> Sro::FnmLsMesh << | updateElmCoords")

        do iel = 1, size(this % els)
            this % inp % x(this % els(iel) % con, :) = this % els(iel) % x
        end do    
    end subroutine

    
    subroutine fnmlsmesh_enforceLsConstraint(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
    
        integer :: iel, ind, aux(1), con(4), lcon(4), elind, bnd(size(this % inp % boundary))
        real(real64) :: ls(4), phi(this % inp % nnodesls)

        ! vars init
        iel     = 0
        elind   = 0
        bnd     = 0
        ls      = 0.0d0
        con  = 0
        lcon = 0
        aux =0
        ind = 0

        call sroLog(" >> Sro::FnmLsMesh << | enforceLsConstraint")

        bnd = this % inp % boundary
        ls = 0.5d0
        phi = this % assembleLS()

        do iel = 1, size(bnd)
            elind = bnd(iel)

            con = this % els(elind) % getLsCon()
            do ind = 1, 4
                aux = findloc(this % inp % lsnodes, con(ind))
                if (.not. aux(1) == 0) then
                    if (aux(1) > size(phi)) then
                        call sroXIT(" >>>>>> fnmlsmesh_enforceLsConstraint <<<<<< Bounds violation")
                    end if
                    phi(aux(1)) = 0.5d0
                else
                    call sroXIT(" >>>>>> fnmlsmesh_enforceLsConstraint <<<<<< LS node index not found")
                end if
            end do
        end do

        do iel = 1, size(this % els)
            con = this % els(iel) % getLsCon()
            do ind = 1, size(con)
                aux = findloc(this % inp % lsnodes, con(ind))
                if (.not. aux(1) == 0) then
                    if (aux(1) > size(phi)) then
                        call sroXIT(" >>>>>> fnmlsmesh_enforceLsConstraint <<<<<< Bounds violation")
                    end if
                    lcon(ind) = aux(1)
                else
                    call sroXIT(" >>>>>> fnmlsmesh_enforceLsConstraint <<<<<< LS node index not found")
                end if
            end do
            call this % els(iel) % updateLs(phi(lcon))
        end do
    
    end subroutine

    
    subroutine fnmlsmesh_partition(this, forceStatus)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
        integer, optional, intent(in)   :: forceStatus
    
        integer :: iel 

        ! vars init
        iel = 0
        
        call sroLog(" >> Sro::FnmLsMesh << | partition")

        if (this % inp % hasBoundary) then
            call this % enforceLsConstraint
        end if

        do iel = 1, size(this % els)
            if (present(forceStatus)) then 
                call this % els(iel) % partition(forceStatus)
            else
                call this % els(iel) % partition
            end if
        end do

        call this % updateElmCoords
        if (.not. present(forceStatus)) then
            call this % assembleBoundaryNodes
            call this % reinitialiseLs
        end if
    end subroutine

    
    subroutine fnmlsmesh_extrapolLsVelocity(this, velocity)
        implicit none 
        class(FnmLsMesh), intent(inout) :: this 
        real(real64), intent(inout) :: velocity(:)

        integer         :: iel, ind, con(this % nnodel), rcon(8), bfcon(5), ufcon(5), cnt(this % inp % nnodes)
        real(real64)    :: xel(this % nnodel, 3), vect(8,3), xx(this % nnodel), vecti(3), vecti2(3)
        real(real64)    :: vel(this % nnodel), v(4)
        real(real64)    :: x1, x3, x9, x10, x11, x12, x21, x41, x22, x42

        ! vars init
        iel     = 0
        ind     = 0
        con     = 0
        rcon    = 0
        bfcon   = 0
        ufcon   = 0
        cnt     = 0
        xel     = 0.0d0
        vect    = 0.0d0
        xx      = 0.0d0
        vecti   = 0.0d0
        vecti2  = 0.0d0
        vel     = 0.0d0
        v       = 0.0d0
        x1      = 0.0d0
        x3      = 0.0d0
        x9      = 0.0d0
        x10     = 0.0d0
        x11     = 0.0d0
        x12     = 0.0d0
        x21     = 0.0d0
        x41     = 0.0d0
        x22     = 0.0d0
        x42     = 0.0d0

        call sroLog(" >> Sro::FnmLsMesh << | extrapolLsVelocity")

        do iel = 1, size(this % els)
            xel = this % els(iel) % x
            con = this % els(iel) % con
            vel = velocity(con)
            rcon = con(this % els(iel) % realNds)
            ufcon = con((/13, 14, 15, 16, 18/))
            bfcon = con((/9, 10, 11, 12, 17/))

            select case(this % els(iel) % status)
            case(12)
                vect(1:5,:) = xel((/10,1,2,3,4/), :) - xel((/9,9,9,9,9/),:)
                vecti = vect(1,:) / vecNorm(vect(1,:))

                xx(9)  = 0.0d0
                xx(1)  = xx(9) + dot_product(vecti, vect(2,:))
                xx(2)  = xx(9) + dot_product(vecti, vect(3,:))
                xx(3)  = xx(9) + dot_product(vecti, vect(4,:))
                xx(4)  = xx(9) + dot_product(vecti, vect(5,:))
                xx(10) = xx(9) + vecNorm(vect(1,:))

                v(1) = ( (vel(9) - vel(10)) / (xx(9) - xx(10)) ) * (xx(1) - xx(9)) + vel(9)
                v(2) = ( (vel(9) - vel(10)) / (xx(9) - xx(10)) ) * (xx(2) - xx(9)) + vel(9)
                v(3) = ( (vel(9) - vel(10)) / (xx(9) - xx(10)) ) * (xx(3) - xx(9)) + vel(9)
                v(4) = ( (vel(9) - vel(10)) / (xx(9) - xx(10)) ) * (xx(4) - xx(9)) + vel(9)

                velocity(rcon(1:4)) = velocity(rcon(1:4)) + v
                velocity(rcon(5:8)) = velocity(rcon(5:8)) + v
                velocity(ufcon) = velocity(bfcon)
                cnt(rcon) = cnt(rcon) + 1
            case(13)
                vect(1:5,:) = xel((/11,1,2,3,4/), :) - xel((/9,9,9,9,9/),:)
                vecti = vect(1,:) / vecNorm(vect(1,:))

                xx(9)  = 0.0d0
                xx(1)  = xx(9) + dot_product(vecti, vect(2,:))
                xx(2)  = xx(9) + dot_product(vecti, vect(3,:))
                xx(3)  = xx(9) + dot_product(vecti, vect(4,:))
                xx(4)  = xx(9) + dot_product(vecti, vect(5,:))
                xx(11) = xx(9) + vecNorm(vect(1,:))

                v(1) = ( (vel(9) - vel(11)) / (xx(9) - xx(11)) ) * (xx(1) - xx(9)) + vel(9)
                v(2) = ( (vel(9) - vel(11)) / (xx(9) - xx(11)) ) * (xx(2) - xx(9)) + vel(9)
                v(3) = ( (vel(9) - vel(11)) / (xx(9) - xx(11)) ) * (xx(3) - xx(9)) + vel(9)
                v(4) = ( (vel(9) - vel(11)) / (xx(9) - xx(11)) ) * (xx(4) - xx(9)) + vel(9)

                velocity(rcon(1:4)) = velocity(rcon(1:4)) + v
                velocity(rcon(5:8)) = velocity(rcon(5:8)) + v
                velocity(ufcon) = velocity(bfcon)
                cnt(rcon) = cnt(rcon) + 1
            case(14)
                vect(1:5,:) = xel((/12,1,2,3,4/), :) - xel((/9,9,9,9,9/),:)
                vecti = vect(1,:) / vecNorm(vect(1,:))

                xx(9)  = 0.0d0
                xx(1)  = xx(9) + dot_product(vecti, vect(2,:))
                xx(2)  = xx(9) + dot_product(vecti, vect(3,:))
                xx(3)  = xx(9) + dot_product(vecti, vect(4,:))
                xx(4)  = xx(9) + dot_product(vecti, vect(5,:))
                xx(12) = xx(9) + vecNorm(vect(1,:))

                v(1) = ( (vel(9) - vel(12)) / (xx(9) - xx(12)) ) * (xx(1) - xx(9)) + vel(9)
                v(2) = ( (vel(9) - vel(12)) / (xx(9) - xx(12)) ) * (xx(2) - xx(9)) + vel(9)
                v(3) = ( (vel(9) - vel(12)) / (xx(9) - xx(12)) ) * (xx(3) - xx(9)) + vel(9)
                v(4) = ( (vel(9) - vel(12)) / (xx(9) - xx(12)) ) * (xx(4) - xx(9)) + vel(9)

                velocity(rcon(1:4)) = velocity(rcon(1:4)) + v
                velocity(rcon(5:8)) = velocity(rcon(5:8)) + v
                velocity(ufcon) = velocity(bfcon)
                cnt(rcon) = cnt(rcon) + 1
            case(23)
                vect(1:5,:) = xel((/11,1,2,3,4/), :) - xel((/10,10,10,10,10/),:)
                vecti = vect(1,:) / vecNorm(vect(1,:))

                xx(10)  = 0.0d0
                xx(1)  = xx(10) + dot_product(vecti, vect(2,:))
                xx(2)  = xx(10) + dot_product(vecti, vect(3,:))
                xx(3)  = xx(10) + dot_product(vecti, vect(4,:))
                xx(4)  = xx(10) + dot_product(vecti, vect(5,:))
                xx(11) = xx(10) + vecNorm(vect(1,:))

                v(1) = ( (vel(10) - vel(11)) / (xx(10) - xx(11)) ) * (xx(1) - xx(10)) + vel(10)
                v(2) = ( (vel(10) - vel(11)) / (xx(10) - xx(11)) ) * (xx(2) - xx(10)) + vel(10)
                v(3) = ( (vel(10) - vel(11)) / (xx(10) - xx(11)) ) * (xx(3) - xx(10)) + vel(10)
                v(4) = ( (vel(10) - vel(11)) / (xx(10) - xx(11)) ) * (xx(4) - xx(10)) + vel(10)

                velocity(rcon(1:4)) = velocity(rcon(1:4)) + v
                velocity(rcon(5:8)) = velocity(rcon(5:8)) + v
                velocity(ufcon) = velocity(bfcon)
                cnt(rcon) = cnt(rcon) + 1
            case(24)
                vect(1:5,:) = xel((/12,1,2,3,4/), :) - xel((/10,10,10,10,10/),:)
                vecti = vect(1,:) / vecNorm(vect(1,:))

                xx(10)  = 0.0d0
                xx(1)  = xx(10) + dot_product(vecti, vect(2,:))
                xx(2)  = xx(10) + dot_product(vecti, vect(3,:))
                xx(3)  = xx(10) + dot_product(vecti, vect(4,:))
                xx(4)  = xx(10) + dot_product(vecti, vect(5,:))
                xx(12) = xx(10) + vecNorm(vect(1,:))

                v(1) = ( (vel(10) - vel(12)) / (xx(10) - xx(12)) ) * (xx(1) - xx(10)) + vel(10)
                v(2) = ( (vel(10) - vel(12)) / (xx(10) - xx(12)) ) * (xx(2) - xx(10)) + vel(10)
                v(3) = ( (vel(10) - vel(12)) / (xx(10) - xx(12)) ) * (xx(3) - xx(10)) + vel(10)
                v(4) = ( (vel(10) - vel(12)) / (xx(10) - xx(12)) ) * (xx(4) - xx(10)) + vel(10)

                velocity(rcon(1:4)) = velocity(rcon(1:4)) + v
                velocity(rcon(5:8)) = velocity(rcon(5:8)) + v
                velocity(ufcon) = velocity(bfcon)
                cnt(rcon) = cnt(rcon) + 1
            case(34)
                vect(1:5,:) = xel((/12,1,2,3,4/), :) - xel((/11,11,11,11,11/),:)
                vecti = vect(1,:) / vecNorm(vect(1,:))

                xx(11)  = 0.0d0
                xx(1)  = xx(11) + dot_product(vecti, vect(2,:))
                xx(2)  = xx(11) + dot_product(vecti, vect(3,:))
                xx(3)  = xx(11) + dot_product(vecti, vect(4,:))
                xx(4)  = xx(11) + dot_product(vecti, vect(5,:))
                xx(12) = xx(11) + vecNorm(vect(1,:))

                v(1) = ( (vel(11) - vel(12)) / (xx(11) - xx(12)) ) * (xx(1) - xx(11)) + vel(11)
                v(2) = ( (vel(11) - vel(12)) / (xx(11) - xx(12)) ) * (xx(2) - xx(11)) + vel(11)
                v(3) = ( (vel(11) - vel(12)) / (xx(11) - xx(12)) ) * (xx(3) - xx(11)) + vel(11)
                v(4) = ( (vel(11) - vel(12)) / (xx(11) - xx(12)) ) * (xx(4) - xx(11)) + vel(11)

                velocity(rcon(1:4)) = velocity(rcon(1:4)) + v
                velocity(rcon(5:8)) = velocity(rcon(5:8)) + v
                velocity(ufcon) = velocity(bfcon)
                cnt(rcon) = cnt(rcon) + 1
            case(1234)
                vect =                                              &
                    xel((/ 12,  1,  2,  4, 11,  2,  3,  4/), :) -   &
                    xel((/  9,  9,  9,  9, 10, 10, 10, 10/), :)
                    
                vecti = vect(1,:) / vecNorm(vect(1,:))
                vecti2 = vect(5, :) / vecNorm(vect(5, :))

                x9  = 0.0d0
                x1  = x9 + dot_product(vecti, vect(2,:))
                x21 = x9 + dot_product(vecti, vect(3,:))
                x41 = x9 + dot_product(vecti, vect(4,:))
                x12 = x9 + vecNorm(vect(1,:))
                x10 = 0.0d0
                x22 = x10 + dot_product(vecti2, vect(6,:))
                x3  = x10 + dot_product(vecti2, vect(7,:))
                x42 = x10 + dot_product(vecti2, vect(8,:))
                x11 = x10 + vecNorm(vect(5,:))

                v(1) = ( (vel(9) - vel(12)) / (x9 - x12) ) * ( x1 - x9 ) + vel(9)
                v(2) = 0.5d0 * ( ( (vel( 9) - vel(12)) / ( x9 - x12) ) * ( x21 -  x9 ) + vel( 9) ) + &
                       0.5d0 * ( ( (vel(10) - vel(11)) / (x10 - x11) ) * ( x22 - x10 ) + vel(10) )
                v(3) = ( (vel(10) - vel(11)) / (x10 - x11) ) * ( x3 - x10 ) + vel(10)
                v(4) = 0.5 * ( ( (vel( 9) - vel(12)) / ( x9 - x12) ) * ( x41 -  x9 ) + vel( 9) ) + &
                       0.5 * ( ( (vel(10) - vel(11)) / (x10 - x11) ) * ( x42 - x10 ) + vel(10) ) 

                velocity(rcon(1:4)) = velocity(rcon(1:4)) + v
                velocity(rcon(5:8)) = velocity(rcon(5:8)) + v
                velocity(ufcon) = velocity(bfcon)
                cnt(rcon) = cnt(rcon) + 1
            end select
        end do
        do ind = 1, size(velocity)
            if (cnt(ind) > 0) velocity(ind) = velocity(ind) / cnt(ind)
        end do
    end subroutine

    
    subroutine fnmlsmesh_computeLsVelocity(this, constantVal)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
        real(real64), optional, intent(in) :: constantVal
    
        real(real64) :: ls(this % inp % nnodesls)
        real(real64) :: e(this % inp % nnodes, this % nstr)
        real(real64) :: s(this % inp % nnodes, this % nstr)
        real(real64) :: eadj(this % inp % nnodes, this % nstr)
        real(real64) :: sadj(this % inp % nnodes, this % nstr)
        real(real64) :: eps(3)
        real(real64) :: sig(3)
        real(real64) :: sigadj(3)
        real(real64) :: constls(4)

        integer :: iel, ind, nbnd, ndid
        real(real64) :: vg, vc, v(this % inp % nnodes), vscale
        real(real64) :: vvg(this % inp % nnodes), vvc(this % inp % nnodes)
        
        ! vars init
        ls      = 0.0d0
        e       = 0.0d0
        s       = 0.0d0
        eadj    = 0.0d0
        sadj    = 0.0d0
        eps     = 0.0d0
        sig     = 0.0d0
        sigadj  = 0.0d0
        iel     = 0
        ind     = 0
        ndid    = 0
        nbnd    = 0
        vg      = 0.0d0
        vc      = 0.0d0
        v       = 0.0d0
        vvg     = 0.0d0
        vvc     = 0.0d0
        constls = 0.0d0
        vscale  = 0.0d0

        call sroLog(" >> Sro::FnmLsMesh << | computeLsVelocity")

        if (present(constantVal)) then 
            constls = (/ constantVal, constantVal, constantVal, constantVal /)
            
            do iel = 1, size(this % els)
                call this % els(iel) % updateLsVelocity( constls )
            end do
            
            v = constantVal
        else 

            ls = this % assembleLs()
            e = this % assembleMeanStrain()
            s = this % assembleMeanStress()
            if (this % inp % hasDisbond) then
                eadj = this % assembleMeanStrainAdjoint()
                sadj = this % assembleMeanStressAdjoint()
            end if

            call this % computeCompliance
            call this % computeEngRelRate
            call this % computeDomainVolume
            call this % opt % update

            nbnd = size(this % boundaryNodes)
            
            do ind = 1, nbnd
                ndid = this % boundaryNodes(ind)
                if (this % inp % hasDisbond) then
                    eps = e(ndid, (/1, 2, 4/))
                    sig = s(ndid, (/1, 2, 4/))
                    sigadj = sadj(ndid, (/1, 2, 4/))
                    
                    vc = dot_product(sig, eps)
                    vg = -dot_product(sigadj, eps)

                    ! v(ndid) = (this % opt % alpha * vg) + ((1.0d0 - this % opt % alpha) * vc)
                    vvg(ndid) = vg
                    vvc(ndid) = vc
                else 
                    eps = e(ndid, (/1, 2, 4/))
                    sig = s(ndid, (/1, 2, 4/))

                    vc = dot_product(sig, eps)

                    v(ndid) = vc
                end if
            end do 

            call printVar((/maxval(abs(vvc))/), "vc_maxabs", errunit)
            call printVar((/maxval(abs(vvg))/), "vg_maxabs", errunit)

            vvg = vvg / maxval(abs(vvg))
            vvg = vvg * maxval(abs(vvc))

            v = (this % opt % alpha * vvg) + ((1.0d0 - this % opt % alpha) * vvc)

            call printVar((/maxval(abs(v))/), "v_maxabs", errunit)

            v(this % boundaryNodes) = (v(this % boundaryNodes) / SRO_VEL_SCALING ) - this % opt % lambda

            call this % extrapolateLsVelocity(v)

            do iel = 1, size(this % els)
                call this % els(iel) % updateLsVelocity(v(this % els(iel) % getLsCon()))
            end do
        end if

    end subroutine

    
    subroutine fnmlsmesh_computeDomainVolume(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
    
        integer :: iel, isub

        ! vars init
        iel     = 0
        isub    = 0

        call sroLog(" >> Sro::FnmLsMesh << | computeDomainVolume")

        this % opt % vol = 0.0d0
        do iel = 1, this % inp % nuels 
            call this % els(iel) % computeVolume
            do isub = 1, this % els(iel) % nsubs
                if (this % els(iel) % getOnOffStatus(isub) /= 0) then 
                    this % vol = this % vol + this % els(iel) % subs(isub) % volume
                end if
            end do
            this % opt % vol = this % opt % vol + this % els(iel) % getArea()
        end do
    end subroutine

    
    subroutine fnmlsmesh_computeCompliance(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
    
        integer :: iel, isub
        real(real64) :: comp

        ! vars init
        iel = 0
        isub = 0
        comp = 0.0d0

        call sroLog(" >> Sro::FnmLsMesh << | computeCompliance")

        this % opt % comp = 0.0d0
        do iel = 1, this % inp % nuels
            comp = 0.0d0
            do isub = 1, this % els(iel) % nsubs
                if (this % els(iel) % getOnOffStatus(isub) /= 0) then
                    comp = comp + &
                        dot_product(this % els(iel) % subs(isub) % u, &
                            matmul(this % els(iel) % subs(isub) % k, this % els(iel) % subs(isub) % u))
                end if
            end do
            this % els(iel) % compliance = comp
            this % opt % comp = this % opt % comp + comp
        end do
    end subroutine

    
    subroutine fnmlsmesh_computeEngRelRate(this)
        implicit none
    
        class(FnmLsMesh), intent(inout) :: this
    
        integer :: iel 

        ! vars init
        iel = 0

        call sroLog(" >> Sro::FnmLsMesh << | computeEngRelRate")

        this % opt % enrr = 0.0d0 
        do iel = 1, this % inp % nadjuels 
            call this % adjels(iel) % computeG
            this % opt % enrr = this % opt % enrr + this % adjels(iel) % G
        end do
    end subroutine

    
    subroutine fnmlsmesh_assembVtkDataArrays(this, u, phi, v, s, e, uu, adj, aadj, eadj, sadj)
        implicit none
    
        class(FnmLsMesh),   intent(inout) :: this
        real(real64),       intent(inout) :: u(this % inp % nnodes, 3)
        real(real64),       intent(inout) :: phi(this % inp % nnodes)
        real(real64),       intent(inout) :: v(this % inp % nnodes)
        real(real64),       intent(inout) :: s(this % inp % nnodes, this % nstr)
        real(real64),       intent(inout) :: e(this % inp % nnodes, this % nstr)
        real(real64),       intent(inout) :: uu(:, :)
        real(real64),       intent(inout) :: adj(:, :)
        real(real64),       intent(inout) :: aadj(:, :)
        real(real64),       intent(inout) :: eadj(:, :)
        real(real64),       intent(inout) :: sadj(:, :)

        integer :: iel, con(18), lscon(4), vcon(8), vu(24)

        ! vars init
        u   = 0.0d0
        phi = 0.0d0 
        v   = 0.0d0
        s   = 0.0d0
        e   = 0.0d0

        call sroLog(" >> Sro::FnmLsMesh << | assembVtkDataArrays")
        
        e = this % assembleMeanStrain()
        s = this % assembleMeanStress()
        v(this % inp % lsnodes) = this % assembleVelocity()

        eadj = this % assembleMeanStrainAdjoint()
        sadj = this % assembleMeanStressAdjoint()

        do iel = 1, size(this % els)
            con = this % els(iel) % con
            u(con,:) = this % els(iel) % getDisplacement()
            phi(con) = this % els(iel) % phi 
            adj(con,:) = this % els(iel) % getAdjoint()
        end do

        do iel = 1, size(u, 1)
            if (abs(u(iel, 1)) < 1e-20) u(iel,1) = 0.0d0
            if (abs(u(iel, 2)) < 1e-20) u(iel,1) = 0.0d0
            if (abs(u(iel, 3)) < 1e-20) u(iel,1) = 0.0d0
            if (abs(adj(iel, 1)) < 1e-20) adj(iel,1) = 0.0d0
            if (abs(adj(iel, 2)) < 1e-20) adj(iel,1) = 0.0d0
            if (abs(adj(iel, 3)) < 1e-20) adj(iel,1) = 0.0d0
            if (abs(phi(iel)) < 1e-20) phi(iel) = 0.0d0
            if (abs(v(iel)) < 1e-20) v(iel) = 0.0d0
        end do

        if (this % inp % hasDisbond) then
            uu = 0.0d0
            do iel = 1, size(this % adjels)
                vcon = this % adjels(iel) % con 
                uu(vcon, :) = this % adjels(iel) % getDisplacement()
                aadj(vcon, :) = this % adjels(iel) % getAdjoint()
            end do

            do iel = 1, size(uu, 1)
                if (abs(uu(iel, 1)) < 1e-20) uu(iel,1) = 0.0d0
                if (abs(uu(iel, 2)) < 1e-20) uu(iel,2) = 0.0d0
                if (abs(uu(iel, 3)) < 1e-20) uu(iel,3) = 0.0d0
                if (abs(adj(iel, 1)) < 1e-20) adj(iel,1) = 0.0d0
                if (abs(adj(iel, 2)) < 1e-20) adj(iel,2) = 0.0d0
                if (abs(adj(iel, 3)) < 1e-20) adj(iel,3) = 0.0d0
            end do
        end if
    end subroutine


    subroutine fnmlsmesh_writeOutput(this, filename)
        implicit none 

        class(FnmLsMesh), intent(inout) :: this 
        character(len=*), intent(in)    :: filename 

        character(len=5) :: iter 
        character(len=:), allocatable :: its
        integer :: n, ii, jj, nc, nn, idxc, idxn, lo, s(2)
        integer, allocatable :: vc(:)
        integer, allocatable :: vo(:)
        integer, allocatable :: vt(:)
        integer, parameter:: hex=12, wedge=13
        real(real64), allocatable :: onoff(:)
        type(VTKFile) :: vtkf
        
        real(real64) :: u(this % inp % nnodes, 3)
        real(real64) :: phi(this % inp % nnodes)
        real(real64) :: v(this % inp % nnodes)
        real(real64) :: sig(this % inp % nnodes, this % nstr)
        real(real64) :: eps(this % inp % nnodes, this % nstr)
        real(real64), allocatable :: cx(:,:), cy(:,:), cz(:,:)
        real(real64), allocatable :: uu(:,:)
        
        integer :: vnc, vnn
        integer, allocatable :: vvc(:), vvo(:)
        integer, allocatable :: vvt(:)

        real(real64), allocatable :: bnd(:)
        integer :: bndaux(1)
        logical :: bndbool = .false.

        ! set vars
        integer, allocatable :: nd_id(:)
        integer, allocatable :: lsnds(:)
        integer, allocatable :: nd_surf(:)

        ! adjoint vars
        real(real64) :: adj(this % inp % nnodes, 3)
        real(real64), allocatable :: aadj(:,:)
        real(real64) :: sigadj(this % inp % nnodes, this % nstr)
        real(real64) :: epsadj(this % inp % nnodes, this % nstr)

        ! vars init
        iter    = "     "
        its     = ""
        n       = 0
        ii      = 0
        jj      = 0
        nc      = 0
        nn      = 0
        idxc    = 0
        idxn    = 0
        lo      = 0
        s       = 0
        u       = 0.0d0
        phi     = 0.0d0
        v       = 0.0d0
        sig     = 0.0d0
        eps     = 0.0d0
        vnc     = 0
        vnn     = 0
        bndaux  = 0
        
        call sroLog(" >> Sro::FnmLsMesh << | writeOutput")

        write(iter, "(i5)") this % iter
        its = trim(adjustl(iter))
        
        nc = 0
        nn = 0
        do ii = 1, this % inp % nuels
            s = shape(this % els(ii) % activeNds)
            nc = nc + s(1)
            nn = nn + s(1)*s(2)
        
            this % inp % x(this % els(ii) % con, :) = &
                this % els(ii) % x
        end do
        
        allocate(vc(nn), source=0)
        allocate(vo(nc), source=0)
        allocate(vt(nc), source=hex)
        allocate(onoff(nc), source=1.0d0)
        allocate(bnd(nc), source = 0.0d0)

        allocate(cx(nc,3), source = 0.0d0)
        allocate(cy(nc,3), source = 0.0d0)
        allocate(cz(nc,3), source = 0.0d0)
        
        idxc = 1
        idxn = 1
        do ii = 1, this % inp % nuels 
            s = shape(this % els(ii) % activeNds)
            
            bndbool = .false.
            bndaux = findloc(this % inp % boundary, ii)
            if (.not. bndaux(1) == 0) bndbool = .true.

            do jj = 1, s(1)
                if (idxc > 1) then 
                    lo = vo(idxc - 1)
                else 
                    lo = 0 
                end if
                vo(idxc) = lo + s(2)
        
                select case(s(2))
                case(6)
                    vt(idxc) = wedge 
                case(8)
                    vt(idxc) = hex 
                end select
        
                vc(idxn:(idxn-1+s(2))) = this % els(ii) % activeNds(jj,:) - 1
                
                onoff(idxc) = this % els(ii) % getOnOffStatus(jj)

                cx(idxc,:) = this % els(ii) % localCS(:,1)
                cy(idxc,:) = this % els(ii) % localCS(:,2)
                cz(idxc,:) = this % els(ii) % localCS(:,3)

                if (bndbool .eqv. .true.) bnd(idxc) = 1.0d0
                
                idxc = idxc + 1 
                idxn = idxn + s(2)
            end do
        end do

        ! set vars
        allocate(nd_id(this % inp % nnodes), source = 0)
        allocate(lsnds(this % inp % nnodes), source = 0)
        allocate(nd_surf(this % inp % nnodes), source = 0)

        lsnds(this % inp % lsnodes) = 1
        do ii = 1, this % inp % nnodes
            nd_id(ii) = this % inp % nds_global(ii)
        end do
        do ii = 1, this % inp % nuels
            nd_surf(this % els(ii) % con((/1, 2, 3, 4, 9, 10, 11, 12, 17/))) = 0
            nd_surf(this % els(ii) % con((/5, 6, 7, 8, 13, 14, 15, 16, 18/))) = 1
        end do

        if (this % inp % hasDisbond) then
            allocate(uu(this % inp % nnodesadjuel, 3), source = 0.0d0)
            allocate(aadj(this % inp % nnodesadjuel, 3), source = 0.0d0)

            vnc = size(this % adjels)
            vnn = 8 * vnc

            allocate(vvc(vnn), source = 0)
            allocate(vvo(vnc), source = 0)
            allocate(vvt(vnc), source = hex)

            idxc = 1
            idxn = 1
            do ii = 1, vnc
                if (idxc > 1) then
                    lo = vvo(idxc - 1)
                else
                    lo = 0
                end if

                vvo(idxc) = lo + 8
                vvc(idxn:(idxn-1+8)) = this % adjels(ii) % con - 1

                idxc = idxc + 1
                idxn = idxn + 8
            end do
        end if

        call this % assembleVtkDataArrays(u, phi, v, sig, eps, uu, adj, aadj, sigadj, epsadj)

        vtkf = VTKFile(                                                         &
            filename    = this % dir // "_vtk" // SLASH // filename // its // ".vtu",   &
            description = "StringerRunoutOptimisation"                          &
        )

        call vtkf % open
        call vtkf % writeUnstructuredGrid(  &
            points  = this % inp % x,       &
            cells   = vc,                   &
            offsets = vo,                   &
            types   = vt                    &
        )
        call vtkf % openCellData(nc)
        call vtkf % writeScalar(v = onoff,      name = "Density")
        call vtkf % writeScalar(v = bnd,        name = "Boundary")
        call vtkf % writeVector(v = cx,         name = "CSx")
        call vtkf % writeVector(v = cy,         name = "CSy")
        call vtkf % writeVector(v = cz,         name = "CSz")
        call vtkf % openPointData(this % inp % nnodes)
        call vtkf % writeVector(v = u,          name = "Displacement")
        if (this % inp % hasDisbond) then
            call vtkf % writeVector(v = adj,    name = "Adjoint")
        end if
        call vtkf % writeScalar(v = phi,        name = "LevelSet")
        call vtkf % writeScalar(v = v,          name = "Velocity")
        call vtkf % writeScalar(v = sig(:,1),   name = "S11")
        call vtkf % writeScalar(v = sig(:,2),   name = "S22")
        call vtkf % writeScalar(v = sig(:,3),   name = "S33")
        call vtkf % writeScalar(v = sig(:,4),   name = "S12")
        call vtkf % writeScalar(v = sig(:,5),   name = "S13")
        call vtkf % writeScalar(v = sig(:,6),   name = "S23")
        call vtkf % writeScalar(v = eps(:,1),   name = "E11")
        call vtkf % writeScalar(v = eps(:,2),   name = "E22")
        call vtkf % writeScalar(v = eps(:,3),   name = "E33")
        call vtkf % writeScalar(v = eps(:,4),   name = "E12")
        call vtkf % writeScalar(v = eps(:,5),   name = "E13")
        call vtkf % writeScalar(v = eps(:,6),   name = "E23")
        call vtkf % writeScalar(v = sigadj(:,1),   name = "SADJ11")
        call vtkf % writeScalar(v = sigadj(:,2),   name = "SADJ22")
        call vtkf % writeScalar(v = sigadj(:,3),   name = "SADJ33")
        call vtkf % writeScalar(v = sigadj(:,4),   name = "SADJ12")
        call vtkf % writeScalar(v = sigadj(:,5),   name = "SADJ13")
        call vtkf % writeScalar(v = sigadj(:,6),   name = "SADJ23")
        call vtkf % writeScalar(v = epsadj(:,1),   name = "EADJ11")
        call vtkf % writeScalar(v = epsadj(:,2),   name = "EADJ22")
        call vtkf % writeScalar(v = epsadj(:,3),   name = "EADJ33")
        call vtkf % writeScalar(v = epsadj(:,4),   name = "EADJ12")
        call vtkf % writeScalar(v = epsadj(:,5),   name = "EADJ13")
        call vtkf % writeScalar(v = epsadj(:,6),   name = "EADJ23")
        call vtkf % writeScalarInt(v = nd_id,   name = "GlobalId")
        call vtkf % writeScalarInt(v = lsnds,   name = "LsNodes")
        call vtkf % writeScalarInt(v = nd_surf, name = "NodeSurfaceId") 
        call vtkf % close

        if (this % inp % hasDisbond) then
            vtkf = vtkfile(&
                filename = this % dir // "_vtk" // SLASH // filename // "base" // its // ".vtu", &
                description = "sro" & 
            )

            call vtkf % open
            call vtkf % writeunstructuredgrid( &
                points = this % inp % adjx, &
                cells = vvc, &
                offsets = vvo, &
                types = vvt & 
            )

            call vtkf % openpointdata(this % inp % nnodesadjuel)
            call vtkf % writeVector(v = uu, name = "Displacement")
            call vtkf % writeVector(v = aadj, name = "Adjoint")
            call vtkf % close
        end if

        open(newunit=this % resFile, file=this % dir // "_log" // SLASH // "results.csv", status="OLD", access="APPEND")
        write(this % resFile, '((i5), a, 10(es16.9, a), (es16.9))') &
            this % iter,            ", ", &
            this % opt % comp,      ", ", &
            this % opt % enrr,      ", ", &
            this % opt % vol,       ", ", &
            this % opt % currentVolFrac,    ", ", &
            this % opt % lambda,    ", ", &
            this % opt % lagrMult,  ", ", &
            this % opt % penalty,   ", ", &
            this % opt % lagr,      ", ", &
            this % opt % alpha,     ", ", &
            this % opt % volFrac,   ", ", &
            this % opt % volInit
        close(this % resFile)
        
        if (allocated(vc)) deallocate(vc)
        if (allocated(vo)) deallocate(vo)
        if (allocated(vt)) deallocate(vt)
        if (allocated(onoff)) deallocate(onoff)
        if (allocated(uu)) deallocate(uu)
        if (allocated(vvc)) deallocate(vvc)
        if (allocated(vvo)) deallocate(vvo)
        if (allocated(vvt)) deallocate(vvt)
        if (allocated(cx)) deallocate(cx)
        if (allocated(cy)) deallocate(cy)
        if (allocated(cz)) deallocate(cz)
        if (allocated(nd_id)) deallocate(nd_id)
        if (allocated(lsnds)) deallocate(lsnds)
        if (allocated(nd_surf)) deallocate(nd_surf)
    end subroutine
    
    function fnmlsmesh_getAnalysisSettings(this) result(set)
        implicit none
        class(FnmLsMesh) :: this 
        type(FnmLsSettings) :: set

        call sroLog(" >> Sro::FnmLsMesh << | getAnalysisSettings")

        set = this % settings 
    end function

    subroutine fnmlsmesh_setAnalysisSettings(this, set)
        implicit none 
        class(FnmLsMesh), intent(inout) :: this 
        class(FnmLsSettings), intent(in) :: set

        call sroLog(" >> Sro::FnmLsMesh << | setAnalysisSettings")

        this % settings = set 
    end subroutine
    
    subroutine fnmlsmesh_registerSteps(this, statusNames)
        implicit none 
        class(FnmLsMesh), intent(inout) :: this
        character(len=50), intent(in) :: statusNames(:)

        integer :: i, j, idx
        integer :: statusIds(size(statusNames))
        character(len=50) :: stages(size(statusNames))
        ! integer :: statusIds(2 * size(statusNames))
        ! character(len=50) :: stages(2 * size(statusNames))
        character(len=:), allocatable :: aux1, aux2, stepType

        ! vars init 
        i           = 0
        j           = 0
        statusIds   = 0
        idx         = 0
        aux1        = ""
        aux2        = ""
        stepType    = ""

        call sroLog(" >> Sro::FnmLsMesh << | registerSteps")

        ! j = 1
        ! do i = 1, size(stages)
        !     if (mod(i, 2) == 0) then
        !         stages(i) = "sync"
        !     else
        !         stages(i) = statusNames(j)
        !         j = j + 1
        !     end if
        ! end do
        stages = statusNames

        do i = 1, size(stages)
            aux1 = trim(adjustl(stages(i)))
            do j = 1, size(this % statusTypes)
                aux2 = trim(adjustl(this % statusTypes(i)))
                idx = index(aux1, aux2)
                if (idx /= 0) then
                    statusIds(i) = idx
                end if
            end do
        end do
        
        if (allocated(this % statusList)) deallocate(this % statusList)
        if (allocated(this % statusLabels)) deallocate(this % statusLabels)

        allocate(this % statusList(size(statusIds)), source=statusIds)
        allocate(this % statusLabels(size(stages)), source=stages)
    end subroutine

    function fnmlsmesh_getStepType(this) result(type)
        implicit none

        character(len=:), allocatable :: type 
        class(FnmLsMesh) :: this

        integer :: i
        character(len=:), allocatable :: aux1, aux2

        ! vars init
        type = "initialisation"
        aux1 = ""
        aux2 = ""
        i    = 0

        aux1 = trim(adjustl(this % statusLabels(this % status)))
        do i = 1, size(this % statusTypes)
            aux2 = trim(adjustl(this % statusTypes(i)))
            if (index(aux1, aux2) /= 0) then
                type = aux2
            end if
        end do

    end function
    
    subroutine fnmlsmesh_updateAnalysisStatus(this)
        implicit none 
        class(FnmLsMesh), intent(inout) :: this

        integer :: curr, new

        ! vars init
        curr = 0
        new  = 0

        if (.not. allocated(this % statusList)) then 
            call sroXIT(" >>>>>> fnmlsmesh_updateAnalysisStatus <<<<<< Status list is not initialised ") 
        end if
        if (.not. allocated(this % statusLabels)) then 
            call sroXIT(" >>>>>> fnmlsmesh_updateAnalysisStatus <<<<<< Status labels are not initialised ") 
        end if 

        curr = this % status 
        new = this % status + 1
        if (new > size(this % statusList)) new = 1

        if (curr == -1) then 
            call sroLog(" >> Sro::FnmLsMesh << | updateAnalysisStatus | Current = beginning" // char(10))
        else if (curr == 0) then 
            call sroLog(" >> Sro::FnmLsMesh << | updateAnalysisStatus | Current = initialisation" // char(10))
        else
            call sroLog(" >> Sro::FnmLsMesh << | updateAnalysisStatus | Current = " &
               // trim(adjustl(this % statusLabels(curr))) // char(10))
        end if
        if (new == 0) then
            call sroLog(" >> Sro::FnmLsMesh << | updateAnalysisStatus | New     = initialisation")
        else 
            call sroLog(" >> Sro::FnmLsMesh << | updateAnalysisStatus | New     = " &
                // trim(adjustl(this % statusLabels(new))))
        end if

        this % status = new
    end subroutine

    
    subroutine fnmlsmesh_updateElm(this, elid, eltype, step, inc, u, k, rhs, enrg)
        implicit none
    
        class(FnmLsMesh),   intent(inout)   :: this
        integer,            intent(in)      :: elid, eltype, step, inc
        real(real64),       intent(in)      :: u(:)
        real(real64),       intent(inout)   :: k(:,:)
        real(real64),       intent(inout)   :: rhs(:)
        real(real64),       intent(inout)   :: enrg(:)
    
        integer :: elindex(1), i, iel
        character(len=:), allocatable :: aux1, aux2, stepType
        real(real64) :: v(this % inp % nnodesls)

        ! vars init 
        elindex     = 0
        i           = 0
        iel         = 0
        aux1        = ""
        aux2        = ""
        stepType    = ""
        v           = 0.0d0

        if (this % status == -1) then 
            call sroXIT(" >>>>>> fnmlsmesh_updateElm <<<<<< Analysis status was not updated ")
        else if (this % status == 0) then ! initialisation
            stepType = "initialisation"
            select case(eltype)
            case(810)
                elindex = findloc(this % inp % uels_global, elid)
                iel = this % inp % uels(elindex(1))
                
                call this % els(iel) % update(u, k, rhs, enrg, stepType, this % dt, this % logFile)
            end select
        else
            stepType = this % getStepType()

            if ((stepType == "adjoint") .and. (.not. this % inp % hasDisbond)) then 
                call sroXIT(" >>>>>> fnmlsmesh_updateElm <<<<<< Cannot run <adjoint> analysis on model without disbond ")
            end if
    
            select case(eltype)
            case(810) ! FnmLsHex
                elindex = findloc(this % inp % uels_global, elid)
                if (elindex(1) == 0) call sroXIT(" >>>>>> fnmlsmesh_updateElm <<<<<< Unable to match id FnmLsHex")
                iel = this % inp % uels(elindex(1))

                call this % els(iel) % update(u, k, rhs, enrg, stepType, this % dt, this % logFile)
                
            case(8) ! AdjHex
                elindex = findloc(this % inp % adjuels_global, elid)
                if (elindex(1) == 0) call sroXIT(" >>>>>> fnmlsmesh_updateElm <<<<<< Unable to match id AdjHex")
                iel = this % inp % adjuels(elindex(1))
    
                call this % adjels(iel) % update(u, k, rhs, stepType, this % logFile)
    
            case default
                call sroXIT(" >>>>>> fnmlsmesh_updateElm <<<<<< unknown element type ") 
            end select
        end if

    end subroutine

    subroutine fnmlsmesh_update(this, stepType, stepIdx, incIdx)
        implicit none
    
        class(FnmLsMesh),               intent(inout)   :: this
        character(len=:), allocatable,  intent(in)      :: stepType
        integer,                        intent(in)      :: stepIdx
        integer,                        intent(in)      :: incIdx

        integer                         :: iel
        character(len=:), allocatable   :: type
        real(real64)                    :: v(this % inp % nnodesls)
    
        ! vars init
        iel = 0
        type = ""
        v = 0.0d0

        type = this % getStepType()
        call sroLog(" >> Sro::FnmLsMesh << | update | (" // &
         int2str(stepIdx) // "," // int2str(incIdx) // ") --- " // type // " --- " // stepType)

        select case(stepType)
        case("analysisStart")

        case("analysisEnd")
        case("stepStart")

            call sroLog(" >> Sro::FnmLsMesh << | update | Syncing with global dof")
            do iel = 1, this % inp % nuels
                sro % dof(iel, :) = this % els(iel) % getDof()
            end do
            
            if (this % inp % hasDisbond) then
                do iel = 1, this % inp % nadjuels
                    sro % adjdof(iel, :) = this % adjels(iel) % getDof()
                end do
            end if

            select case(type)
            case("initialisation")

                if (this % settings % useConstantVelocity .eqv. .true.) then
                    call sroLog(" >> Sro::FnmLsMesh << | update | Assigning constant velocity")
                    call this % computeLsVelocity(this % settings % constantVelocity)
                end if

                call sroLog(" >> Sro::FnmLsMesh << | update | Initialising global dof")
                do iel = 1, this % inp % nuels
                    sro % dof(iel, :) = this % els(iel) % getDof()
                end do

            case("sync")
                
            case("displacement")

                call sroLog(" >> Sro::FnmLsMesh << | update | Partitioning mesh")
                if (this % settings % useForcedPartition .eqv. .true.) then
                    call this % partition(this % settings % partitionMode)
                else 
                    call this % partition()
                end if

                call sroLog(" >> Sro::FnmLsMesh << | update | Computing stiffness for displacement")
                do iel = 1, this % inp % nuels
                    call this % els(iel) % computeStiffness 
                end do

            case("adjoint")

                call sroLog(" >> Sro::FnmLsMesh << | update | Computing stiffness for adjoint")
                do iel = 1, this % inp % nuels
                    call this % els(iel) % computeStiffness 
                end do

            case("velocity")

                if (this % settings % velocityIter == 0) then
                    if (this % settings % useConstantVelocity .eqv. .true.) then
                        call sroLog(" >> Sro::FnmLsMesh << | update | Assigning constant velocity")
                        call this % computeLsVelocity(this % settings % constantVelocity)
                    else 
                        call sroLog(" >> Sro::FnmLsMesh << | update | Computing velocity")
                        call this % computeLsVelocity()
                    end if

                    call sroLog(" >> Sro::FnmLsMesh << | update | Assembling velocity - compute dt")
                    v = this % assembleVelocity(.true.)
                end if
                this % settings % velocityIter = this % settings % velocityIter + 1

                call sroLog(" >> Sro::FnmLsMesh << | update | Computing velocity extension load")
                do iel = 1, this % inp % nuels
                    call this % els(iel) % levelsetElm % computeExtensForce(this % dt)
                end do

            case("levelset")

                this % settings % velocityIter = 0

                call sroLog(" >> Sro::FnmLsMesh << | update | Writing vtk output for current step")
                call this % writeOutput(this % jobname)
                this % iter = this % iter + 1

                call sroLog(" >> Sro::FnmLsMesh << | update | Computing levelset evolution load")
                do iel = 1, this % inp % nuels
                    call this % els(iel) % levelsetElm % computeForce(this % dt)
                end do

            end select

        case("stepEnd")

            call sroLog(" >> Sro::FnmLsMesh << | update | Syncing with elements dof")
            do iel = 1, this % inp % nuels
                call this % els(iel) % setDof(sro % dof(iel, :))
            end do

            if (this % inp % hasDisbond) then
                do iel = 1, this % inp % nadjuels
                    call this % adjels(iel) % setDof(sro % adjdof(iel,:))
                end do
            end if

            select case(type)
            case("initialisation")

            case("sync")

            case("displacement")

                call sroLog(" >> Sro::FnmLsMesh << | update | Computing stress and strain")
                do iel = 1, this % inp % nuels
                    call this % els(iel) % computeStrainStress
                end do

                call sroLog(" >> Sro::FnmLsMesh << | update | Computing mean compliance")
                call this % computeCompliance

                if (this % inp % hasDisbond) then
                    call sroLog(" >> Sro::FnmLsMesh << | update | Computing energy release rate")
                    call this % computeEngRelRate
                end if

                ! call sroLog(" >> Sro::FnmLsMesh << | update | Writing vtk output for current step")
                ! call this % writeOutput(this % jobname)
                ! this % iter = this % iter + 1

            case("adjoint")

                call sroLog(" >> Sro::FnmLsMesh << | update | Computing stress and strain for adjoint")
                do iel = 1, this % inp % nuels
                    call this % els(iel) % computeStrainStressAdjoint
                end do

            case("velocity")

            case("levelset")

            end select

            call this % updateAnalysisStatus

        case("incrementStart")

        case("incrementEnd")

        end select

        if (.not. (stepIdx == this % currentStep .and. incIdx == this % currentIncrement)) then

            this % currentStep = stepIdx
            this % currentIncrement = incIdx 

            select case(stepType)
            case("analysisStart")
            case("analysisEnd")
            case("stepStart")
            case("stepEnd")
            case("incrementStart")
            case("incrementEnd")
            end select
            
        end if

    end subroutine
   
end module
