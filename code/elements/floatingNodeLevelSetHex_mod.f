!DEC$ FREEFORM

module floatingNodeLevelSetHex_mod
    use iso_fortran_env
    use dof_mod
    use element_mod
    use layered3D_mod
    use solidShellHex_mod
    use solidShellWedge_mod
    use levelsetQuad_mod
    implicit none

    type, extends(Element) :: FnmLsHex

        integer                             :: nrealnodes   = -1
        integer                             :: nfloatnodes  = -1
        integer                             :: nrealdof     = -1
        integer                             :: nfloatdof    = -1
        integer                             :: npartitions  = -1
        integer                             :: nfieldvars   = -1
        integer                             :: nactivenodes = -1
        integer                             :: nsubs        = -1
        integer                             :: status       = -1

        real(real64)                        :: volume       = 0.0d0
        real(real64)                        :: compliance   = 0.0d0
        real(real64)                        :: energRelRate = 0.0d0

        integer                             :: realNds(8)   = 0
        integer                             :: floatNds(10) = 0
        integer                             :: upperNds(4)  = 0
        integer                             :: bottomNds(4) = 0
        integer                             :: lsNds(4)     = 0
        integer,                allocatable :: activeNds(:,:)
        integer,                allocatable :: boundaryNds(:)
        integer,                allocatable :: lsCon(:)

        real(real64),           allocatable :: orientations(:)
        real(real64),           allocatable :: thicknesses(:)
        real(real64),           allocatable :: phi(:)
        real(real64),           allocatable :: adj(:)
        real(real64),           allocatable :: dkda(:,:)
        real(real64),           allocatable :: eadj(:,:,:,:)
        real(real64),           allocatable :: sadj(:,:,:,:)

        integer,                allocatable :: dofMap(:)

        type(Q4LS)                          :: levelsetElm
        class(LayeredElm3D),    allocatable :: subs(:)
        class(PartitionCase),   allocatable :: cases(:)

        type(Transform3D) :: transformLs
        type(Transform3D) :: transformShell

        real(real64)                        :: localCS(3,3) = 0.0d0
        type(Transform3D)                   :: Tgl
        type(Transform3D)                   :: Tlg

    contains

        procedure :: print                  => fnmlshex_print

        procedure :: initPartitionTable     => fnmlshex_initPartitionTable
        procedure :: partition              => fnmlshex_partition
        procedure :: partitionByMode        => fnmlshex_partitionByMode
        procedure :: getPartitionById       => fnmlshex_getPartitionById
        procedure :: updateFloatingNodes    => fnmlshex_updateFloatingNodes

        procedure :: finalize               => fnmlshex_final

        procedure :: isInit                 => fnmlshex_isInit
        procedure :: areSubsInit            => fnmlshex_areSubsInit

        procedure :: computeLocalCS         => fnmlshex_computeLocalCS
        procedure :: updateLocalCS          => fnmlshex_updateLocalCS

        procedure :: computeVolume          => fnmlshex_computeVolume
        procedure :: getArea                => fnmlshex_getArea
        procedure :: computeStiffness       => fnmlshex_computeStiffness

        procedure :: computeStrainStress            => fnmlshex_computeStrainStress
        procedure :: computeStrainStressAdjoint     => fnmlshex_computeStrainStressAdjoint

        procedure :: getCon                 => fnmlshex_getCon
        procedure :: getLsCon               => fnmlshex_getLsCon
        procedure :: getLevelset            => fnmlshex_getlevelset
        procedure :: getVelocity            => fnmlshex_getVelocity
        procedure :: getOnOffStatus         => fnmlshex_getOnOffStatus
        procedure :: getDisplacement        => fnmlshex_getDisplacement
        procedure :: getAdjoint             => fnmlshex_getAdjoint
        procedure :: getDof                 => fnmlshex_getDof
        procedure :: setDof                 => fnmlshex_setDof

        procedure :: getMeanStrain          => fnmlshex_getMeanStrain
        procedure :: getMeanStress          => fnmlshex_getMeanStress
        procedure :: getMeanStrainAdjoint   => fnmlshex_getMeanStrainAdjoint
        procedure :: getMeanStressAdjoint   => fnmlshex_getMeanStressAdjoint

        procedure :: updateDisplacement     => fnmlshex_updateDisp
        procedure :: updateAdjoint          => fnmlshex_updateAdjoint
        procedure :: updateLs               => fnmlshex_updateLs
        procedure :: updateLsVelocity       => fnmlshex_updateLsVelocity
        procedure :: update                 => fnmlshex_update

        final :: fnmlshex_dtor
    end type

    interface FnmLsHex
        procedure :: fnmlshex_ctor
    end interface

contains


    function fnmlshex_ctor(id, coords, con, props, theta, t, transformLs, transformShell, originLs) result(this)
        implicit none

        type(FnmLsHex)                      :: this
        integer,            intent(in)      :: id
        integer,            intent(in)      :: con(:)
        real(real64),       intent(in)      :: coords(:,:)
        real(real64),       intent(in)      :: props(:)
        real(real64),       intent(in)      :: theta(:)
        real(real64),       intent(in)      :: t(:)
        type(Transform3D),  intent(in)      :: transformLs
        type(Transform3D),  intent(in)      :: transformShell
        real(real64),       intent(in)      :: originLs(:)

        integer :: i, auxdofmap(116)

        ! vars init
        i           = 0
        auxdofmap   = 0

        this % name = "FnmLsHex"

        this % nnodes       = 18
        this % ndof         = 54
        this % ndofnode     = 3
        this % ndim         = 3
        this % nstr         = 6
        this % npts         = 5
        this % nprops       = size(props)

        this % nrealnodes   = 8
        this % nfloatnodes  = 10
        this % nrealdof     = 24
        this % nfloatdof    = 30
        this % npartitions  = 10
        this % nfieldvars   = 4
        this % nactivenodes = 8
        this % nsubs        = 1
        this % status       = 0

        this % realNds      = (/ (i, i=1,  8) /)
        this % upperNds     = (/ (i, i=1,  4) /)
        this % bottomNds    = (/ (i, i=5,  8) /)
        this % floatNds     = (/ (i, i=9, 18) /)
        this % lsNds        = (/ (i, i=1,  4) /)

        this % transformLs      = transformLs
        this % transformShell   = transformShell

        call this % initialize(id, coords, con, props)

        allocate(this % orientations(size(theta)), source=theta)
        allocate(this % thicknesses(size(t)),      source=t)

        allocate(this % activeNds(1, this % nactivenodes), source = 0)
        allocate(this % phi(this % nnodes), source = 1.0d0)
        allocate(this % adj(this % ndof), source = 0.0d0)
        allocate(this % dkda(this % ndof, this % ndof), source=0.0d0)

        auxdofmap = (/                              &
              1,   2,   3,   9,  10,  11,  17,  18, &
             19,  25,  26,  27,  33,  34,  35,  39, &
             40,  41,  45,  46,  47,  51,  52,  53, &
             57,  58,  59,  63,  64,  65,  69,  70, &
             71,  75,  76,  77,  81,  82,  83,  87, &
             88,  89,  93,  94,  95,  99, 100, 101, &
            105, 106, 107, 111, 112, 113,   4,   5, &
              6,  12,  13,  14,  20,  21,  22,  28, &
             29,  30,  36,  37,  38,  42,  43,  44, &
             48,  49,  50,  54,  55,  56,  60,  61, &
             62,  66,  67,  68,  72,  73,  74,  78, &
             79,  80,  84,  85,  86,  90,  91,  92, &
             96,  97,  98, 102, 103, 104, 108, 109, &
            110, 114, 115, 116,   7,  15,  23,  31, &
              8,  16,  24,  32                      &
        /)
        allocate(this % dofMap(116), source = auxdofmap)

        call this % initPartitionTable

        this % levelsetElm = Q4LS(      &
            1,                          &
            this % x((/1,2,3,4/),:),    &
            con((/1,2,3,4/)),           &
            (/1.0d0/),                  &
            transformLs,                &
            originLs                    &
        )

        call this % computeLocalCS
        ! placeholder dummy subelement for initialization purposes only
        allocate(this % subs(1),                        &
            source=S8PS(                                &
                    1,                                  &
                    this % x(1:this % nrealnodes, :),   &
                    con(1:this % nrealnodes),           &
                    props,                              &
                    theta,                              &
                    t,                                  &
                    this % localCS,                     &
                    this % Tlg,                         &
                    this % Tgl                          &
                )                                       &
        )

        call this % partition

        call this % computeVolume
    end function


    subroutine fnmlshex_dtor(this)
        implicit none

        type(FnmLsHex), intent(inout) :: this

        call this % finalize
    end subroutine


    function fnmlshex_isInit(this) result(b)
        implicit none

        logical :: b
        class(FnmLsHex) :: this

        ! vars init
        b = .false.

        if ( this % nrealnodes     /= -1    .and. &
             this % nfloatnodes    /= -1    .and. &
             this % nrealdof       /= -1    .and. &
             this % nfloatdof      /= -1    .and. &
             this % npartitions    /= -1    .and. &
             this % nfieldvars     /= -1    .and. &
             this % nactivenodes   /= -1    .and. &
             this % nsubs          /= -1    .and. &
             allocated(this % activeNds)    .and. &
             allocated(this % orientations) .and. &
             allocated(this % thicknesses)  .and. &
             allocated(this % phi)          .and. &
             allocated(this % adj)          .and. &
             allocated(this % subs)         .and. &
             allocated(this % cases)        .and. &
             allocated(this % dkda)         .and. &
             allocated(this % dofMap)       .and. &
             this % areSubsInit()           .and. &
             this % levelsetElm % isInit()  .and. &
             this % Element % isInit() ) then
            b = .true.
        else
            call sroLog(" >> Sro::FnmLsHex << | isInit | Attempted use before initialisation ")
            b = .false.
        end if

    end function


    function fnmlshex_areSubsInit(this) result(b)
        implicit none

        logical :: b
        class(FnmLsHex) :: this

        integer :: isub

        ! vars init
        b = .true.
        isub = 0

        do isub = 1, this % nsubs
            if (.not. this % subs(isub) % isInit()) b = .false.
        end do
    end function


    subroutine fnmlshex_final(this)
        implicit none

        class(FnmLsHex), intent(inout) :: this

        integer :: isub

        ! vars init
        isub = 0

        call this % levelsetElm % finalize

        do isub = 1, this % nsubs
            call this % subs(isub) % finalize
        end do

        if (allocated(this % activeNds)) deallocate(this % activeNds)
        if (allocated(this % boundaryNds)) deallocate(this % boundaryNds)
        if (allocated(this % lsCon)) deallocate(this % lsCon)
        if (allocated(this % orientations)) deallocate(this % orientations)
        if (allocated(this % thicknesses)) deallocate(this % thicknesses)
        if (allocated(this % phi)) deallocate(this % phi)
        if (allocated(this % adj)) deallocate(this % adj)
        if (allocated(this % eadj)) deallocate(this % eadj)
        if (allocated(this % sadj)) deallocate(this % sadj)
        if (allocated(this % subs)) deallocate(this % subs)
        if (allocated(this % cases)) deallocate(this % cases)
        if (allocated(this % dkda)) deallocate(this % dkda)
        if (allocated(this % dofMap)) deallocate(this % dofMap)

        call this % levelsetElm % finalize
        call this % Element % finalize
    end subroutine


    subroutine fnmlshex_computeLocalCS(this)
        implicit none

        class(FnmLsHex), intent(inout) :: this

        real(real64) :: r21(this % ndim), r31(this % ndim), r42(this % ndim)
        real(real64) :: s1(this % ndim), s3(this % ndim)
        real(real64) :: e1(this % ndim), e2(this % ndim), e3(this % ndim)
        real(real64) :: lcs(this % ndim, this % ndim)
        real(real64) :: c(this % nrealnodes, this % ndim)
        real(real64) :: xdir(this % ndim), xprojdir(this % ndim)
        real(real64) :: csys(this % ndim, this % ndim)

        ! vars init
        r21         = 0.0d0
        r31         = 0.0d0
        r42         = 0.0d0
        s1          = 0.0d0
        s3          = 0.0d0
        e1          = 0.0d0
        e2          = 0.0d0
        e3          = 0.0d0
        lcs         = 0.0d0
        c           = 0.0d0
        xdir        = 0.0d0
        xprojdir    = 0.0d0
        csys        = 0.0d0

        csys = this % transformShell % toCS

        if (this % Element % isInit()) then
            c(:,1) = this%x(this % realNds,1)! + this % u((/1,4,7,10,13,16,19,22/))
            c(:,2) = this%x(this % realNds,2)! + this % u((/2,5,8,11,14,17,20,23/))
            c(:,3) = this%x(this % realNds,3)! + this % u((/3,6,9,12,15,18,21,24/))

            r31 = c(3, :) - c(1, :)
            r42 = c(4, :) - c(2, :)

            s3 = vec3cross(r31, r42)
            e3 = s3 / vecNorm(s3)

            xdir = csys(:,1)
            xprojdir = xdir - (dot_product(xdir, e3) * e3)
            e1 = xprojdir / vecNorm(xprojdir)

            ! r21 = c(2, :) - c(1, :)
            ! s1 = r21 - (dot_product(r21, e3) * e3)
            ! e1 = s1 / vecNorm(s1)

            e2 = vec3cross(e3, e1)

            lcs(:,:) = 0.0d0
            lcs(:,1) = e1
            lcs(:,2) = e2
            lcs(:,3) = e3

            this % localCS = lcs
            this % Tgl = Transform3D(this % localCS, Imat)
            this % Tlg = Transform3D(this % Tgl % T3, Imat)
        else
            call sroXIT(" >>>>>> fnmlshex_computeLocalCS <<<<<< Base object not initialised ")
        end if

    end subroutine


    subroutine fnmlshex_updateLocalCS(this)
        implicit none

        class(FnmLsHex), intent(inout) :: this

        integer :: isub

        ! vars init
        isub = 0

        call this % computeLocalCS
        do isub = 1, this % nsubs
            call this % subs(isub) % updateLocalCS(this % localCS, this % Tlg, this % Tgl)
        end do
    end subroutine


    subroutine fnmlshex_print(this, varName, unit)
        implicit none

        class(FnmLsHex), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit

        integer :: uid, isub, icas
        character(len=3) :: cs

        ! vars init
        uid     = 0
        isub    = 0
        icas    = 0
        cs      = "   "

        uid = 6
        if (present(unit)) uid = unit

        if (this % isInit()) then
            call this % Element % print(varName, uid)
            write(uid, fmt100) ""
            write(uid, fmt101) "nrealnodes", this % nrealnodes
            write(uid, fmt101) "nfloatnodes", this % nfloatnodes
            write(uid, fmt101) "nrealdof", this % nrealdof
            write(uid, fmt101) "nfloatdof", this % nfloatdof
            write(uid, fmt101) "npartitions", this % npartitions
            write(uid, fmt101) "nfieldvars", this % nfieldvars
            write(uid, fmt101) "nactivenodes", this % nactivenodes
            write(uid, fmt101) "nsubs", this % nsubs
            write(uid, fmt101) "status", this % status
            write(uid, fmt100) ""
            write(uid, fmt102) "volume", this % volume
            write(uid, fmt102) "compliance", this % compliance
            write(uid, fmt102) "energRelRate", this % energRelRate
            write(uid, fmt100) ""
            call printVar(this % activeNds, "activeNds", uid)
            call printVar(this % orientations, "orientations", uid)
            call printVar(this % thicknesses, "thicknesses", uid)
            call printVar(this % phi, "phi", uid)
            call printVar(this % adj, "adj", uid)
            write(uid, fmt100) ""
            call this % levelsetElm % print(varName // ".levelsetElm", uid)
            do isub = 1, this % nsubs
                write(uid, fmt100) ""
                write(cs, "(i3)") isub
                call this % subs(isub) % print(varName // ".subs("// cs // ")", uid)
            end do
            do icas = 1, this % npartitions
                write(uid, fmt100) ""
                write(cs, "(i3)") icas
                call this % cases(icas) % print(varName // ".cases("// cs // ")", uid)
            end do
            write(uid, fmt100) ""

        else
            call sroXIT(" >>>>>> fnmlshex_print <<<<<< Object not initialised ")
        end if

    end subroutine


    subroutine fnmlshex_initPartitionTable(this)
        implicit none

        class(FnmLsHex), intent(inout) :: this

        type(PartitionCase) :: pt(this % npartitions)
        integer :: con0(1,8)
        integer :: con12(6,6), con13(2,8), con14(6,6)
        integer :: con23(6,6), con24(2,8)
        integer :: con34(6,6)
        integer :: con1234(8,6)
        integer :: con99(2,6)
        integer :: con98(4,6)

        ! vars init
        con1234     = 0
        con0        = 0
        con12       = 0
        con13       = 0
        con14       = 0
        con23       = 0
        con24       = 0
        con34       = 0
        con99       = 0

        if (.not. allocated(this % cases)) then
            con0(1,:) = (/1, 2, 3, 4, 5, 6, 7, 8/)
            pt(1) = PartitionCase(0, con0, (/0/))

            con12(1,:) = (/  1,  9, 17,  5, 13, 18/)
            con12(2,:) = (/  9, 10, 17, 13, 14, 18/)
            con12(3,:) = (/  9,  2, 10, 13,  6, 14/)
            con12(4,:) = (/  1, 17,  4,  5, 18,  8/)
            con12(5,:) = (/ 17,  3,  4, 18,  7,  8/)
            con12(6,:) = (/ 17, 10,  3, 18, 14,  7/)
            pt(2) = PartitionCase(12, con12, (/9, 10/))

            con13(1,:) = (/  1,  9, 11,  4,  5, 13, 15,  8/)
            con13(2,:) = (/  9,  2,  3, 11, 13,  6,  7, 15/)
            pt(3) = PartitionCase(13, con13, (/9, 11/))

            con14(1,:) = (/  1,  9, 12,  5, 13, 16/)
            con14(2,:) = (/  9, 17, 12, 13, 18, 16/)
            con14(3,:) = (/  9,  2, 17, 13,  6, 18/)
            con14(4,:) = (/  2,  3, 17,  6,  7, 18/)
            con14(5,:) = (/ 12, 17,  4, 16, 18,  8/)
            con14(6,:) = (/ 17,  3,  4, 18,  7,  8/)
            pt(4) = PartitionCase(14, con14, (/9, 12/))

            con23(1,:) = (/  1,  2, 17,  5,  6, 18/)
            con23(2,:) = (/  1, 17,  4,  5, 18,  8/)
            con23(3,:) = (/  2, 10, 17,  6, 14, 18/)
            con23(4,:) = (/ 17, 11,  4, 18, 15,  8/)
            con23(5,:) = (/ 17, 10, 11, 18, 14, 15/)
            con23(6,:) = (/ 10,  3, 11, 14,  7, 15/)
            pt(5) = PartitionCase(23, con23, (/10, 11/))

            con24(1,:) = (/  1,  2, 10, 12,  5,  6, 14, 16/)
            con24(2,:) = (/ 12, 10,  3,  4, 16, 14,  7,  8/)
            pt(6) = PartitionCase(24, con24, (/10, 12/))

            con34(1,:) = (/  1,  2, 17,  5,  6, 18/)
            con34(2,:) = (/  1, 17, 12,  5, 18, 16/)
            con34(3,:) = (/  2,  3, 17,  6,  7, 18/)
            con34(4,:) = (/ 12, 17, 11, 16, 18, 15/)
            con34(5,:) = (/ 12, 11,  4, 16, 15,  8/)
            con34(6,:) = (/ 17,  3, 11, 18,  7, 15/)
            pt(7) = PartitionCase(34, con34, (/11, 12/))

            con1234(1,:) = (/  1,  9, 12,  5, 13, 16/)
            con1234(2,:) = (/  9, 17, 12, 13, 18, 16/)
            con1234(3,:) = (/  9, 10, 17, 13, 14, 18/)
            con1234(4,:) = (/  9,  2, 10, 13,  6, 14/)
            con1234(5,:) = (/ 12, 11,  4, 16, 15,  8/)
            con1234(6,:) = (/ 12, 17, 11, 16, 18, 15/)
            con1234(7,:) = (/ 17, 10, 11, 18, 14, 15/)
            con1234(8,:) = (/ 10,  3, 11, 14,  7, 15/)
            pt(8) = PartitionCase(1234, con1234, (/9, 10, 11, 12/))

            con99(1,:) = (/ 1, 2, 3, 5, 6, 7 /)
            con99(2,:) = (/ 1, 3, 4, 5, 7, 8 /)
            pt(9) = PartitionCase(99, con99, (/0/))

            con98(1,:) = (/ 1, 2, 17, 5, 6, 18 /)
            con98(2,:) = (/ 2, 3, 17, 6, 7, 18 /)
            con98(3,:) = (/ 3, 4, 17, 7, 8, 18 /)
            con98(4,:) = (/ 4, 1, 17, 8, 5, 18 /)
            pt(10) = PartitionCase(98, con98, (/0/))

            allocate(this % cases(this % npartitions), source=pt)
        end if
    end subroutine


    subroutine fnmlshex_partition(this, forceStatus)
        implicit none
        class(FnmLsHex), intent(inout) :: this
        integer, optional, intent(in)  :: forceStatus

        real(real64) :: field(this % nfieldvars)
        real(real64) :: aux1(this % nfieldvars), aux2(this % nfieldvars)
        real(real64) :: xloc(this % nfieldvars,2), xintloc(this % nfieldvars,2)
        real(real64) :: threshold, aux3
        logical :: ints(this % nfieldvars)
        logical :: infs(this % nfieldvars)
        integer :: i, j, ids(this % nfieldvars), nids, status, power

        ! vars init
        field       = 0.0d0
        aux1        = 0.0d0
        aux2        = 0.0d0
        xloc        = 0.0d0
        xintloc     = 0.0d0
        threshold   = 0.0d0
        aux3        = 0.0d0
        ints        = .false.
        infs        = .false.
        i           = 0
        j           = 0
        ids         = 0
        nids        = 0
        status      = 0
        power       = 0

        if (this % isInit()) then
            if (present(forceStatus)) then
                if (forceStatus == 0 .or. forceStatus == 99) then
                    this % nactivenodes = 8
                else if (forceStatus == 98) then
                    this % nactivenodes = 10
                else
                    call sroXIT(" >>>>>> fnmlshex_partition <<<<<< forceStatus only possible if status is 0, 99, 98 ")
                end if

                xloc(1,:) = (/-1.0d0, -1.0d0/)
                xloc(2,:) = (/ 1.0d0, -1.0d0/)
                xloc(3,:) = (/ 1.0d0,  1.0d0/)
                xloc(4,:) = (/-1.0d0,  1.0d0/)

                ids = (/-1, -1, -1, -1/)
                nids = 0

                call this % partitionByMode(forceStatus, xloc, ids(1:nids))
            else
                field = this % phi(1:(this % nfieldvars))

                threshold = 0.9d0

                xloc(1,:) = (/-1.0d0, -1.0d0/)
                xloc(2,:) = (/ 1.0d0, -1.0d0/)
                xloc(3,:) = (/ 1.0d0,  1.0d0/)
                xloc(4,:) = (/-1.0d0,  1.0d0/)

                ids = (/-1, -1, -1, -1/)
                nids = 0

                aux1 = field * field((/2,3,4,1/))
                aux2 = -field / (field((/2,3,4,1/)) - field)

                ints = aux1 < 0.0d0
                infs = (aux2 - 1.0d0) == aux2
                xintloc = xloc

                do i = 1, size(ints)
                    if (ints(i)) then
                        nids = nids + 1
                        ids(nids) = i
                        if (i == 1) then
                            aux3 = -1.0d0 + 2.0d0 * aux2(i)
                            xintloc(i,1) = max(min(aux3, threshold), -threshold)
                        else if (i == 2) then
                            aux3 = -1.0d0 + 2.0d0 * aux2(i)
                            xintloc(i,2) = max(min(aux3, threshold), -threshold)
                        else if (i == 3) then
                            aux3 =  1.0d0 - 2.0d0 * aux2(i)
                            xintloc(i,1) = max(min(aux3, threshold), -threshold)
                        else
                            aux3 =  1.0d0 - 2.0d0 * aux2(i)
                            xintloc(i,2) = max(min(aux3, threshold), -threshold)
                        end if
                    end if
                end do

                status = 0
                power = nids - 1
                do i = 1, nids
                    status = status + ids(i) * ( 10 ** power )
                    power = power - 1
                end do

                select case(status)
                case(0)
                    this % nactivenodes = 8
                case(12,14,23,34)
                    this % nactivenodes = 14
                case(13,24)
                    this % nactivenodes = 12
                case(1234)
                    this % nactivenodes = 18
                end select

                call this % partitionByMode(status, xintloc, ids(1:nids))
            end if
        else
            call sroXIT(" >>>>>> fnmlshex_partition <<<<<< Object not initialised ")
        end if
    end subroutine


    subroutine fnmlshex_partitionByMode(this, mode, ratios, edges)
        implicit none
        class(FnmLsHex),    intent(inout)   :: this
        integer,            intent(in)      :: mode
        real(real64),       intent(in)      :: ratios(:,:)
        integer,            intent(in)      :: edges(:)

        type(PartitionCase) :: pc
        integer             :: isub

        type(S6PS), allocatable :: s6list(:)
        type(S8PS), allocatable :: s8list(:)

        integer, allocatable :: con(:)

        ! vars init
        isub = 0

        pc = this % getPartitionById(mode)
        if (.not. (pc % id .eq. -1)) then

            this % status = mode
            this % nsubs = pc % nsubs
            if (allocated(this % subs)) then
                do isub = 1, size(this % subs)
                    call this % subs(isub) % finalize
                end do
                deallocate(this % subs)
            end if
            if (allocated(this % activeNds)) deallocate(this % activeNds)
            if (allocated(this % lsCon)) deallocate(this % lsCon)

            if (allocated(this % boundaryNds)) deallocate(this % boundaryNds)
            if (.not. any(pc % boundaryNds == (/0/))) then
                allocate(this % boundaryNds(size(pc % boundaryNds)), source=this % con(pc % boundaryNds))
                allocate(this % lsCon(this % nfieldvars + size(pc % boundaryNds) / 2), source=0)
                this % lsCon(1:(this % nfieldvars)) = this % lsNds
                this % lsCon((this % nfieldvars + 1):) = this % con(pc % boundaryNds)
            end if

            call this % updateFloatingNodes(mode, ratios, edges)
            if (pc % nnodes .eq. 6) then

                allocate(this % activeNds(this % nsubs, pc % nnodes), source = 0)
                allocate(s6list(this % nsubs))
                allocate(con(6), source = 0)
                do isub = 1, this % nsubs
                    call assign(con, pc % connectivity(isub, :))
                    s6list(isub) = &
                        S6PS( &
                            isub, &
                            this % x(con,:), &
                            con, &
                            this % props, &
                            this % orientations, &
                            this % thicknesses, &
                            this % localCS,                     &
                            this % Tlg,                         &
                            this % Tgl                          &
                    )
                    this % activeNds(isub,:) = this % con(pc % connectivity(isub,:))
                end do
                allocate(this % subs(this % nsubs), source = s6list)

            else if (pc % nnodes .eq. 8) then

                allocate(this % activeNds(this % nsubs, pc % nnodes), source=0)
                allocate(s8list(this % nsubs))
                allocate(con(8), source = 0)
                do isub = 1, this % nsubs
                    call assign(con, pc % connectivity(isub, :))
                    s8list(isub) = &
                        S8PS( &
                            isub, &
                            this % x(con,:), &
                            con, &
                            this % props, &
                            this % orientations, &
                            this % thicknesses, &
                            this % localCS,                     &
                            this % Tlg,                         &
                            this % Tgl                          &
                        )
                    this % activeNds(isub,:) = this % con(pc % connectivity(isub,:))
                end do
                allocate(this % subs(this % nsubs), source = s8list)

            else
                call sroXIT(" >>>>>> fnmlshex_partitionByMode <<<<<< Ill defined partition ")
            end if

            if (allocated(s8list)) deallocate(s8list)
            if (allocated(s6list)) deallocate(s6list)
            if (allocated(con)) deallocate(con)

        else
            call sroXIT(" >>>>>> fnmlshex_partitionByMode <<<<<< Unknown partition ")
        end if
    end subroutine


    subroutine fnmlshex_updateFloatingNodes(this, mode, ratios, edges)
        implicit none
        class(FnmLsHex), intent(inout) :: this
        integer,                     intent(in)    :: mode
        real(real64),                intent(in)    :: ratios(:, :)
        integer,                     intent(in)    :: edges(:)

        real(real64) :: nmat(1, this % nrealnodes)
        real(real64) :: intloc(2)
        real(real64) :: c(this % nnodes, this % ndim)
        real(real64) :: cr(this % nrealnodes, this % ndim)
        integer :: ii, eid, fnid(4,2)

        type(HexIsoSF) :: dummySF

        ! vars init
        nmat    = 0.0d0
        intloc  = 0.0d0
        c       = 0.0d0
        cr      = 0.0d0
        ii      = 0
        eid     = 0
        fnid    = 0

        fnid(1,:) = (/ 9, 13/)
        fnid(2,:) = (/10, 14/)
        fnid(3,:) = (/11, 15/)
        fnid(4,:) = (/12, 16/)

        cr = this % x(1:(this % nrealnodes),:)
        c = 0.0d0
        c(1:this % nrealnodes, :) = cr

        if (this % Element % isInit()) then
            do ii = 1, size(edges)
                eid = edges(ii)
                intloc = ratios(eid,:)
                nmat = dummySF % N_at(intloc(1), intloc(2), -1.0d0, 8)
                c(fnid(eid,1),:) = getPointCoords(cr, nmat, this % ndim)
                nmat = dummySF % N_at(intloc(1), intloc(2),  1.0d0, 8)
                c(fnid(eid,2),:) = getPointCoords(cr, nmat, this % ndim)
            end do
            if (mode == 12 .or. &
                mode == 14 .or. &
                mode == 23 .or. &
                mode == 34 .or. &
                mode == 1234 .or. &
                mode == 98) then

                nmat = dummySF % N_at(0.0d0, 0.0d0, -1.0d0, 8)
                c(17,:) = getPointCoords(cr, nmat, this % ndim)
                nmat = dummySF % N_at(0.0d0, 0.0d0,  1.0d0, 8)
                c(18,:) = getPointCoords(cr, nmat, this % ndim)

                select case (mode)
                case (12)
                    this % phi((/17,18/)) = &
                        sum(this % phi((/1,3,4/)))/3.0d0
                case (14)
                    this % phi((/17,18/)) = &
                        sum(this % phi((/2,3,4/)))/3.0d0
                case (23)
                    this % phi((/17,18/)) = &
                        sum(this % phi((/1,2,4/)))/3.0d0
                case (34)
                    this % phi((/17,18/)) = &
                        sum(this % phi((/1,2,3/)))/3.0d0
                case (1234)
                    this % phi((/17,18/)) = &
                        sum(this % phi((/1,3/)))/2.0d0
                case (98)
                    this % phi((/17,18/)) = &
                        sum(this % phi((/1,2,3,4/)))/4.0d0
                end select

            end if
            call assign(this % x, c)
        end if
    end subroutine


    function fnmlshex_getPartitionById(this, pId) result(p)
        implicit none

        type(PartitionCase) :: p

        class(FnmLsHex),    intent(in) :: this
        integer,            intent(in) :: pId

        integer :: i, index

        ! vars init
        i       = 0
        index   = 0

        p = PartitionCase(-1, 1, 1, reshape((/-1, -1/), (/1,1/)))
        do i = 1, size(this % cases)
            if (this % cases(i) % id .eq. pId) then
                p = this % cases(i)
            end if
        end do
    end function


    subroutine fnmlshex_computeVolume(this)
        implicit none
        class(FnmLsHex), intent(inout) :: this

        integer :: ii

        ! vars init
        ii = 0
        this % volume = 0.0d0

        if (this % isInit()) then
            do ii = 1, this % nsubs
                call this % subs(ii) % computeVolume
                this % volume = this % volume + this % subs(ii) % volume
            end do
        else
            call sroXIT(" >>>>>> fnmlshex_computeStrainStress <<<<<< Object not initialized ")
        end if
    end subroutine


    function fnmlshex_getArea(this) result(area)
        implicit none 

        class(FnmLsHex) :: this
        real(real64) :: area

        integer :: isub
        real(real64) :: a(3), b(3), c(3), d(3)
        real(real64) :: l

        area    = 0.0d0
        a = 0.0d0
        b = 0.0d0
        c = 0.0d0
        d = 0.0d0

        do isub = 1, this % nsubs
            if (this % getOnOffStatus(isub) /= 0) then 
                if (this % subs(isub) % nnodes == 6) then
                    a = this % x(this % subs(isub) % con(1), :)
                    b = this % x(this % subs(isub) % con(2), :)
                    c = this % x(this % subs(isub) % con(3), :)
                    area = area + areaTriangle3D(a, b, c)
                else 
                    a = this % x(this % subs(isub) % con(1), :)
                    b = this % x(this % subs(isub) % con(2), :)
                    c = this % x(this % subs(isub) % con(3), :)
                    d = this % x(this % subs(isub) % con(4), :)
                    area = area + areaQuadrilateral3D(a, b, c, d)
                end if
            end if
        end do

        ! l = 0.0d0
        ! l = l + vecNorm(this % x(1,:) - this % x(5,:))
        ! l = l + vecNorm(this % x(2,:) - this % x(6,:))
        ! l = l + vecNorm(this % x(3,:) - this % x(7,:))
        ! l = l + vecNorm(this % x(4,:) - this % x(8,:))
        ! l = l / 4.0d0

        ! area = this % volume / l
    end function


    subroutine fnmlshex_computeStiffness(this)
        implicit none
        class(FnmLsHex), intent(inout) :: this

        integer :: ii, dofid(this % ndof)
        real(real64) :: onoff

        ! vars init
        ii = 0
        dofid = 0
        onoff = 0.0d0


        if (this % isInit()) then
            call this % updateLocalCS
            this % k = 0.0d0
            do ii = 1, this % nsubs
                onoff = this % getOnOffStatus(ii)
                if (onoff > 0.0d0) then
                    call this % subs(ii) % computeStiffness
                    this % k(this % subs(ii) % dof, this % subs(ii) % dof) = &
                        this % k(this % subs(ii) % dof, this % subs(ii) % dof) + this % subs(ii) % k
                    dofid(this % subs(ii) % dof) = 1
                end if
            end do

            ! Disabling this to get all zeros on unwanted dof
            ! do ii = 1, this % ndof
            !     if (dofid(ii) == 0) then
            !         this % k(ii, ii) = 1.0d12
            !     end if
            ! end do
        else
            call sroXIT(" >>>>>> fnmlshex_computeStrainStress <<<<<< Object not initialized ")
        end if
    end subroutine


    subroutine fnmlshex_computeStrainStress(this)
        implicit none
        class(FnmLsHex), intent(inout) :: this

        integer :: ii, jj

        ! vars init
        ii = 0
        jj = 0

        if (this % isInit()) then
            call this % updateLocalCS
            do ii = 1, this % nsubs
                if (this % getOnOffStatus(ii) /= 0) then
                    call this % subs(ii) % computeStrainStress
                else
                    do jj = 1, this % subs(ii) % nlayers 
                        this % subs(ii) % layers(jj) % e = 0.0d0
                        this % subs(ii) % layers(jj) % s = 0.0d0
                    end do
                end if
            end do
        else
            call sroXIT(" >>>>>> fnmlshex_computeStrainStress <<<<<< Object not initialized ")
        end if
    end subroutine


    subroutine fnmlshex_computeStrainStressAdjoint(this)
        implicit none
        class(FnmLsHex), intent(inout) :: this

        integer :: isub, ilay, ipt, istr
        integer :: nsub, nlay, npts, nstr
        real(real64) :: dmat(this % nstr, this % nstr)
        real(real64) :: eps(this % nstr), sig(this % nstr)

        ! vars init
        isub    = 0
        ilay    = 0
        ipt     = 0
        istr    = 0
        nsub    = 0
        nlay    = 0
        npts    = 0
        nstr    = 0
        dmat    = 0.0d0
        eps     = 0.0d0
        sig     = 0.0d0

        if (allocated(this % eadj)) deallocate(this % eadj)
        if (allocated(this % sadj)) deallocate(this % sadj)

        nsub = this % nsubs ! this one can change with partitioning
        nlay = size(this % orientations)
        npts = this % npts
        nstr = this % nstr

        allocate(this % eadj(nsub, nlay, npts, nstr), source = 0.0d0)
        allocate(this % sadj(nsub, nlay, npts, nstr), source = 0.0d0)

        if (this % isInit()) then
            do isub = 1, nsub
                if (this % getOnOffStatus(isub) /= 0) then
                    call this % subs(isub) % computeStrainStressAdjoint
                else
                    do ilay = 1, this % subs(isub) % nlayers
                        this % subs(isub) % layers(ilay) % eadj = 0.0d0
                        this % subs(isub) % layers(ilay) % sadj = 0.0d0
                    end do
                end if
            end do
        else
            call sroXIT(" >>>>>> fnmlshex_computeStrainStressAdjoint <<<<<< Object not initialized ")
        end if

    end subroutine


    function fnmlshex_getCon(this) result(con)
        implicit none

        class(FnmLsHex) :: this
        integer         :: con(this % nnodes)

        ! vars init
        con = 0

        con = this % con
    end function


    function fnmlshex_getLsCon(this) result(con)
        implicit none

        class(FnmLsHex) :: this
        integer         :: con(this % nfieldvars)

        ! vars init
        con = 0

        con = this % levelsetElm % con
    end function


    function fnmlshex_getlevelset(this) result(ls)
        implicit none

        class(FnmLsHex) :: this
        real(real64)    :: ls(this % nfieldvars)

        ! vars init
        ls = 0.0d0

        ls = this % levelsetElm % u
    end function


    function fnmlshex_getVelocity(this) result(v)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: v(this % nfieldvars)

        ! vars init
        v = 0.0d0

        v = this % levelsetElm % v
    end function


    function fnmlshex_getMeanStrain(this) result(e)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: e(this % nnodes, this % nstr)

        real(real64) :: ee(this % nnodes, this % nstr), auxE(this % nstr)

        integer :: cnt(this % nnodes)
        integer :: isub, istr, ind, nnodessub, isubnd, subnd

        ! vars init
        e           = 0.0d0
        ee          = 0.0d0
        auxE        = 0.0d0
        cnt         = 0
        isub        = 0
        istr        = 0
        ind         = 0
        nnodessub   = 0
        isubnd      = 0
        subnd       = 0

        do isub = 1, this % nsubs
            if (this % getOnOffStatus(isub) /= 0) then
                auxE = this % subs(isub) % computeMeanStrain()

                nnodessub = this % subs(isub) % nnodes
                do isubnd = 1, nnodessub
                    subnd = this % subs(isub) % con(isubnd)
                    ee(subnd,:) = ee(subnd,:) + auxE
                    cnt(subnd)  = cnt(subnd)  + 1
                end do
            end if
        end do

        do istr = 1, this % nstr
            do ind = 1, this % nnodes
                if (cnt(ind) > 0) then
                    e(ind,istr) = ee(ind,istr) / cnt(ind)
                end if
            end do
        end do
    end function


    function fnmlshex_getMeanStress(this) result(s)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: s(this % nnodes, this % nstr)

        real(real64) :: ss(this % nnodes, this % nstr), auxS(this % nstr)

        integer :: cnt(this % nnodes)
        integer :: isub, istr, ind, nnodessub, isubnd, subnd

        ! vars init
        s           = 0.0d0
        ss          = 0.0d0
        auxS        = 0.0d0
        cnt         = 0
        isub        = 0
        istr        = 0
        ind         = 0
        nnodessub   = 0
        isubnd      = 0
        subnd       = 0

        do isub = 1, this % nsubs
            if (this % getOnOffStatus(isub) /= 0) then
                auxS = this % subs(isub) % computeMeanStress()

                nnodessub = this % subs(isub) % nnodes
                do isubnd = 1, nnodessub
                    subnd = this % subs(isub) % con(isubnd)
                    ss(subnd,:) = ss(subnd,:) + auxS
                    cnt(subnd)  = cnt(subnd)  + 1
                end do
            end if
        end do

        do istr = 1, this % nstr
            do ind = 1, this % nnodes
                if (cnt(ind) > 0) then
                    s(ind,istr) = ss(ind,istr) / cnt(ind)
                end if
            end do
        end do
    end function


    function fnmlshex_getMeanStrainAdjoint(this) result(e)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: e(this % nnodes, this % nstr)

        real(real64) :: ee(this % nnodes, this % nstr), auxE(this % nstr)

        integer :: cnt(this % nnodes)
        integer :: isub, ilay, istr, ind, nnodessub, isubnd, subnd
        integer :: nsub, nlay, npts, nstr

        ! vars init
        e           = 0.0d0
        ee          = 0.0d0
        auxE        = 0.0d0
        cnt         = 0
        isub        = 0
        ilay        = 0
        istr        = 0
        ind         = 0
        nnodessub   = 0
        isubnd      = 0
        subnd       = 0
        nsub        = 0
        nlay        = 0
        npts        = 0
        nstr        = 0

        nsub = this % nsubs
        nlay = size(this % orientations)
        npts = this % npts
        nstr = this % nstr

        do isub = 1, nsub
            if (this % getOnOffStatus(isub) /= 0) then
                ! auxE = 0.0d0
                ! do ilay = 1, nlay
                !     auxE = auxE + (sum(this % eadj(isub, ilay, :, :), 1) / npts)
                ! end do
                auxE = this % subs(isub) % computeMeanStrainAdjoint()

                nnodessub = this % subs(isub) % nnodes
                do isubnd = 1, nnodessub
                    subnd = this % subs(isub) % con(isubnd)
                    ee(subnd,:) = ee(subnd,:) + auxE
                    cnt(subnd)  = cnt(subnd)  + 1
                end do
            end if
        end do

        do istr = 1, this % nstr
            do ind = 1, this % nnodes
                if (cnt(ind) > 0) then
                    e(ind,istr) = ee(ind,istr) / cnt(ind)
                end if
            end do
        end do
    end function


    function fnmlshex_getMeanStressAdjoint(this) result(s)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: s(this % nnodes, this % nstr)

        real(real64) :: ss(this % nnodes, this % nstr), auxS(this % nstr)

        integer :: cnt(this % nnodes)
        integer :: isub, ilay, istr, ind, nnodessub, isubnd, subnd
        integer :: nsub, nlay, npts, nstr

        ! vars init
        s           = 0.0d0
        ss          = 0.0d0
        auxS        = 0.0d0
        cnt         = 0
        isub        = 0
        ilay        = 0
        istr        = 0
        ind         = 0
        nnodessub   = 0
        isubnd      = 0
        subnd       = 0
        nsub        = 0
        nlay        = 0
        npts        = 0
        nstr        = 0

        nsub = this % nsubs
        nlay = size(this % orientations)
        npts = this % npts
        nstr = this % nstr

        do isub = 1, this % nsubs
            if (this % getOnOffStatus(isub) /= 0) then
                ! auxS = 0.0d0
                ! do ilay = 1, nlay
                !     auxS = auxS + (sum(this % sadj(isub, ilay, :, :), 1) / npts)
                ! end do
                auxS = this % subs(isub) % computeMeanStressAdjoint()

                nnodessub = this % subs(isub) % nnodes
                do isubnd = 1, nnodessub
                    subnd = this % subs(isub) % con(isubnd)
                    ss(subnd,:) = ss(subnd,:) + auxS
                    cnt(subnd)  = cnt(subnd)  + 1
                end do
            end if
        end do

        do istr = 1, this % nstr
            do ind = 1, this % nnodes
                if (cnt(ind) > 0) then
                    s(ind,istr) = ss(ind,istr) / cnt(ind)
                end if
            end do
        end do
    end function


    function fnmlshex_getDisplacement(this) result(u)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: u(this % nnodes, this % ndim)

        integer :: i, ids(this % nnodes)

        ! vars init
        i   = 0
        ids = 0
        u   = 0.0d0

        ids = (/ (this % ndim * (i - 1) + 1, i = 1, this % nnodes) /)

        u(:,1) = this % u(ids)
        u(:,2) = this % u(ids + 1)
        u(:,3) = this % u(ids + 2)

    end function

    function fnmlshex_getAdjoint(this) result(adj)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: adj(this % nnodes, this % ndim)

        integer :: i, ids(this % nnodes)

        ! vars init
        i   = 0
        ids = 0
        adj = 0.0d0

        ids = (/ (this % ndim * (i - 1) + 1, i = 1, this % nnodes) /)

        adj(:,1) = this % adj(ids)
        adj(:,2) = this % adj(ids + 1)
        adj(:,3) = this % adj(ids + 2)

    end function


    function fnmlshex_getOnOffStatus(this, isub) result(onoff)
        implicit none

        class(FnmLsHex) :: this
        integer :: isub
        real(real64) :: onoff

        real(real64) :: sum1

        ! vars init
        onoff   = 1.0d0
        sum1    = 0.0d0

        sum1 = sum(this % phi(this % subs(isub) % con))
        if (sum1 <= 0.00001d0) then
            onoff = 0.0d0
        end if
    end function

    function fnmlshex_getDof(this) result(q)
        implicit none

        class(FnmLsHex) :: this
        real(real64) :: q(116)
        integer :: i, udof(54), adof(54), vdof(4), pdof(4)

        !vars init
        i = 0
        q = 0.0d0
        udof    = 0
        adof    = 0
        vdof    = 0
        pdof    = 0

        udof = this % dofMap((/ (i, i=1        , 54       ) /))
        adof = this % dofMap((/ (i, i=1+54     , 54+54    ) /))
        vdof = this % dofMap((/ (i, i=1+54+54  , 54+54+4  ) /))
        pdof = this % dofMap((/ (i, i=1+54+54+4, 54+54+4+4) /))

        q(udof) = this % u 
        q(adof) = this % adj 
        q(vdof) = this % getVelocity()
        q(pdof) = this % getLevelset()

    end function

    subroutine fnmlshex_setDof(this, q)
        implicit none

        class(FnmLsHex), intent(inout) :: this
        real(real64), intent(in) :: q(:)

        integer :: i, udof(54), adof(54), vdof(4), pdof(4)

        !vars init
        i = 0
        udof    = 0
        adof    = 0
        vdof    = 0
        pdof    = 0

        udof = this % dofMap((/ (i, i=1        , 54       ) /))
        adof = this % dofMap((/ (i, i=1+54     , 54+54    ) /))
        vdof = this % dofMap((/ (i, i=1+54+54  , 54+54+4  ) /))
        pdof = this % dofMap((/ (i, i=1+54+54+4, 54+54+4+4) /))

        call this % updateDisplacement(q(udof))
        call this % updateAdjoint(q(adof))
        call this % updateLsVelocity(q(vdof))
        call this % updateLs(q(pdof))

    end subroutine

    subroutine fnmlshex_updateDisp(this, u)
        implicit none

        class(FnmLsHex), intent(inout) :: this
        real(real64), intent(in) :: u(this % ndof)

        integer :: isub

        ! vars init
        isub = 0

        call assign(this % u, u)
        do isub = 1, this % nsubs
            call assign(this % subs(isub) % u, this % u(this % subs(isub) % dof))
        end do
    end subroutine 

    subroutine fnmlshex_updateAdjoint(this, p)
        implicit none

        class(FnmLsHex), intent(inout) :: this
        real(real64), intent(in) :: p(this % ndof)

        integer :: isub

        isub = 0

        call assign(this % adj, p)
        do isub = 1, this % nsubs
            call assign(this % subs(isub) % adj, this % adj(this % subs(isub) % dof)) 
        end do
    end subroutine

    subroutine fnmlshex_updateLs(this, ls)
        implicit none

        class(FnmLsHex),    intent(inout)   :: this
        real(real64),       intent(in)      :: ls(this % nfieldvars)

        integer :: i
        real(real64) :: aux

        ! vars init
        i           = 0
        aux         = 0.0d0
        this % phi  = 0.0d0

        this % phi( (/ (i, i=1,  4) /) ) = ls
        this % phi( (/ (i, i=5,  8) /) ) = ls
        this % phi( (/ (i, i=9, 16) /) ) = 0.0d0

        aux = 0.0d0
        select case(this % status)
        case(12)
            aux = sum(ls((/1, 3, 4/))) / 3
        case(14)
            aux = sum(ls((/2, 3, 4/))) / 3
        case(23)
            aux = sum(ls((/1, 2, 4/))) / 3
        case(34)
            aux = sum(ls((/1, 2, 3/))) / 3
        case(1234)
            aux = sum(ls((/1, 3/))) / 2
        end select

        this % phi( (/17, 18/) ) = aux

        call assign(this % levelsetElm % u, ls)
    end subroutine


    subroutine fnmlshex_updateLsVelocity(this, velocity)
        implicit none
        class(FnmLsHex), intent(inout) :: this
        real(real64), intent(in) :: velocity(:)

        call assign(this % levelsetElm % v, velocity)
    end subroutine


    subroutine fnmlshex_update(this, u, k, rhs, enrg, analysisType, dt, unit)
        implicit none

        class(FnmLsHex),    intent(inout)   :: this
        real(real64),       intent(in)      :: u(:)
        real(real64),       intent(inout)   :: k(:,:)
        real(real64),       intent(inout)   :: rhs(:)
        real(real64),       intent(inout)   :: enrg(:)
        character(len=*),   intent(in)      :: analysisType
        real(real64),       intent(in)      :: dt
        integer, optional,  intent(in)      :: unit

        real(real64) :: kel(116, 116), fel(116), fex(116)
        real(real64) :: eye54(54,54), eye4(4,4)

        integer :: i, j, isub, idx(1)
        integer :: udof(54), adof(54), vdof(4), pdof(4)

        ! vars init
        kel     = 0.0d0
        fel     = 0.0d0
        fex     = 0.0d0
        i       = 0
        j       = 0
        isub    = 0
        udof    = 0
        adof    = 0
        vdof    = 0
        pdof    = 0
        eye54   = 0.0d0
        eye4    = 0.0d0

        udof = this % dofMap((/ (i, i=1        , 54       ) /))
        adof = this % dofMap((/ (i, i=1+54     , 54+54    ) /))
        vdof = this % dofMap((/ (i, i=1+54+54  , 54+54+4  ) /))
        pdof = this % dofMap((/ (i, i=1+54+54+4, 54+54+4+4) /))

        ! create identity matrix
        do i = 1, 116
            kel(i,i) = 1.0d0
        end do
        ! do i = 1, 54
        !     eye54(i,i) = 1.0d0
        ! end do
        ! do i = 1, 4
        !     eye4(i,i) = 1.0d0
        ! end do

        idx = findloc(sro % elmap, this % id)
        sro % dof(idx(1), :) = u 

        select case(analysisType)
        case("initialisation")

            fex(udof) = this % u
            fex(adof) = this % adj
            fex(vdof) = this % getVelocity()
            fex(pdof) = this % getLevelset()
            fel       = fex - matmul(kel, u)

        case("sync")

            fex(udof) = this % u
            fex(adof) = this % adj
            fex(vdof) = this % getVelocity()
            fex(pdof) = this % getLevelset()
            fel       = fex - matmul(kel, u)

        case("displacement")

            if (.not. idx(1) == 0) then 
                sro % dof(idx(1), adof) = this % adj
            end if

            kel(udof, udof) = this % k

            fex(adof)       = this % adj
            fex(vdof)       = this % getVelocity()
            fex(pdof)       = this % getLevelset()
            fel             = fex - matmul(kel, u)

        case("adjoint")

            if (.not. idx(1) == 0) then 
                sro % dof(idx(1), udof) = this % u
            end if

            kel(adof, adof) = this % k

            if (any(abs(u(adof)) >= 0.1d0)) then
                kel(adof, adof) = 0.0d0
            end if

            fex(udof)       = this % u
            fex(adof)       = 0.0d0 !matmul(this % dkda, this % u)
            fex(vdof)       = this % getVelocity()
            fex(pdof)       = this % getLevelset()
            fel             = fex - matmul(kel, u)

        case("velocity")
            
            if (.not. idx(1) == 0) then 
                sro % dof(idx(1), udof) = this % u
                sro % dof(idx(1), adof) = this % adj
            end if

            kel(vdof, vdof) = this % levelsetElm % k

            fex(udof)       = this % u
            fex(adof)       = this % adj
            fex(vdof)       = this % levelsetElm % fv
            fex(pdof)       = this % getLevelset()
            fel             = fex - matmul(kel, u)

        case("levelset")

            if (.not. idx(1) == 0) then 
                sro % dof(idx(1), udof) = this % u
                sro % dof(idx(1), adof) = this % adj
            end if

            kel(pdof, pdof) = this % levelsetElm % k

            fex(udof)       = this % u
            fex(adof)       = this % adj
            fex(vdof)       = this % getVelocity()
            fex(pdof)       = this % levelsetElm % f
            fel             = fex - matmul(kel, u)

        end select

        k       = kel
        rhs     = fel
        enrg    = 0.0d0
    end subroutine

end module