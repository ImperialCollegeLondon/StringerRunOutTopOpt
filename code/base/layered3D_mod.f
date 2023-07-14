!DEC$ FREEFORM

module layered3D_mod
    use iso_fortran_env
    use element_mod
    use layer_mod
    implicit none
    
    type, extends(Element), abstract :: LayeredElm3D

        integer :: nlayers = -1

        real(real64) :: volume = 0.0d0

        real(real64), allocatable :: localCS(:,:)

        class(Layer), allocatable :: layers(:)

        real(real64), allocatable :: adj(:)

        type(Transform3D) :: Tgl 
        type(Transform3D) :: Tlg

    contains 
    
        procedure(layeredelm3D_computeVolume),  deferred :: computeVolume
        procedure(layeredelm3D_computeLocalCS), deferred :: computeLocalCS  
        procedure(layeredelm3D_setupLayers),    deferred :: setupLayers 

        procedure :: print                      => layeredelm3D_print
        
        procedure :: updateLocalCS              => layeredelm3D_updateLocalCS

        procedure :: computeLayerDisplacement   => layeredelm3D_computeLayerField

        procedure :: computeStiffness           => layeredelm3D_computeStiffness
        procedure :: computeStrainStress        => layeredelm3D_computeStrainStress
        procedure :: computeMeanStrain          => layeredelm3D_computeMeanStrain
        procedure :: computeMeanStress          => layeredelm3D_computeMeanStress

        procedure :: computeLayerAdjoint        => layeredelm3D_computeLayerAdjoint

        procedure :: computeStrainStressAdjoint => layeredelm3D_computeStrainStressAdjoint
        procedure :: computeMeanStrainAdjoint   => layeredelm3D_computeMeanStrainAdjoint
        procedure :: computeMeanStressAdjoint   => layeredelm3D_computeMeanStressAdjoint

        procedure :: finalize                   => layeredelm3D_final
        procedure :: initialize                 => layeredelm3D_init

        procedure :: isInit                     => layeredelm3D_isInit
        procedure :: areLayersInit              => layeredelm3D_areLayersInit

    end type

    interface
        subroutine layeredelm3D_computeVolume(this)
            import LayeredElm3D
            implicit none
        
            class(LayeredElm3D), intent(inout) :: this
        end subroutine

        subroutine layeredelm3D_computeLocalCS(this, csys)
            use iso_fortran_env
            import LayeredElm3D
            implicit none
        
            class(LayeredElm3D),    intent(inout)   :: this
            real(real64),           intent(in)      :: csys(:,:)
        end subroutine

        subroutine layeredelm3D_setupLayers(this, orientations, thicknesses)
            use iso_fortran_env
            import LayeredElm3D
            implicit none 

            class(LayeredElm3D),    intent(inout)   :: this 
            real(real64),           intent(in)      :: orientations(:)
            real(real64),           intent(in)      :: thicknesses(:)
        end subroutine
    end interface
    
contains

    subroutine layeredelm3D_print(this, varName, unit)
        implicit none
        
        class(LayeredElm3D),    intent(inout)   :: this
        character(len=*),       intent(in)      :: varName
        integer, optional,      intent(in)      :: unit
        
        integer :: uid, ilay
        character(len=3) :: cs

        ! vars init
        uid     = 0
        ilay    = 0
        cs      = "   "
        
        uid = 6
        if (present(unit)) uid = unit 

        if (this % isInit()) then 
            call this % Element % print(varName, uid)
            write(uid, fmt100) ""
            write(uid, fmt102) "volume", this % volume 
            write(uid, fmt100) ""
            call printVar(this % localCS, "localCS", uid)
            write(uid, fmt100) ""
            call this % Tgl % print(varName // ".Tgl", uid)
            call this % Tlg % print(varName // ".Tlg", uid)
            write(uid, fmt101) "nlayers", this % nlayers 
            write(uid, fmt100) ""
            do ilay = 1, this % nlayers 
                write(uid, fmt100) ""
                write(cs, "(i3)") ilay
                call this % layers(ilay) % print(varName // ".layers(" // cs // ")", uid)
            end do
        else 
            call sroXIT(" >>>>>> layeredelm3D_print <<<<<< Object not initialised ")
        end if
        
    end subroutine

    subroutine layeredelm3D_init(this, id, coords, con, props)
        implicit none
    
        class(LayeredElm3D),    intent(inout)   :: this
        integer,                intent(in)      :: id, con(:)
        real(real64),           intent(in)      :: coords(:,:)
        real(real64),           intent(in)      :: props(:)
    
        if (this % nlayers > 0) then 

            this % name = "LayeredElm3D | " // this % name
            call this % element % initialize(id, coords, con, props)

            allocate(this % localCS(this % ndim, this % ndim), source=0.0d0)

            allocate(this % adj(this % ndof), source = 0.0d0)
        else 
            call sroXIT(" >>>>>> layeredelm3D_init <<<<<< Bad input data ")
        end if

    end subroutine

    subroutine layeredelm3D_updateLocalCS(this, localCS, Tlg, Tgl)
        implicit none
    
        class(LayeredElm3D),    intent(inout)   :: this
        real(real64),           intent(in)      :: localCS(:,:)
        type(Transform3D),      intent(in)      :: Tlg, Tgl
    
        integer :: ilay

        ! vars init
        ilay = 0

        do ilay = 1, this % nlayers
            call assign(this % localCS, localCS)
            this % Tlg = Tlg 
            this % Tgl = Tgl
            call this % layers(ilay) % updateLocalCS(Tlg, Tgl)
        end do 
    
    end subroutine

    function layeredelm3D_computeLayerField(this, ilay) result(u)
        implicit none

        class(LayeredElm3D)         :: this
        integer,        intent(in)  :: ilay

        real(real64)                :: u(this % ndof)

        real(real64) :: t
        real(real64) :: tb
        real(real64) :: tt
        integer      :: i
        integer      :: db(this % ndof / 2)
        integer      :: dt(this % ndof / 2)
        real(real64) :: ub(this % ndof / 2)
        real(real64) :: ut(this % ndof / 2)
        real(real64) :: uzb(this % ndof / 2)
        real(real64) :: uzt(this % ndof / 2)

        ! vars init
        u = 0.0d0
        t = 0.0d0
        tb = 0.0d0
        tt = 0.0d0
        i = 0
        ub = 0.0d0 
        ut = 0.0d0 
        uzb = 0.0d0
        uzt = 0.0d0
        db = 0 
        dt = 0

        t = 0.0d0
        do i = 1, this % nlayers
            t = t + this % layers(i) % t
        end do
        tb = this % layers(ilay) % z / t
        tt = (this % layers(ilay) % z + this % layers(ilay) % t) / t

        db = (/(i, i=1, this % ndof / 2)/)
        dt = (/(i, i=this % ndof / 2 + 1, this % ndof)/)

        ub = this % u(db)
        ut = this % u(dt)
        uzb = ((1 - tb) * ub) + (tb * ut)
        uzt = ((1 - tt) * ub) + (tt * ut)

        u(db) = uzb
        u(dt) = uzt

    end function

    subroutine layeredelm3D_computeStrainStress(this)
        implicit none
    
        class(LayeredElm3D), intent(inout) :: this
    
        integer :: ilay 

        ! vars init 
        ilay = 0

        if (this % isInit()) then 
            do ilay = 1, this % nlayers 
                call this % layers(ilay) % computeStrainStress( &
                        this % computeLayerDisplacement(ilay)   &
                    ) 
            end do
        else
            call sroXIT(" >>>>>> layeredelm3D_computeStrainStress <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine layeredelm3D_computeStiffness(this)
        implicit none
    
        class(LayeredElm3D), intent(inout) :: this
    
        integer :: ilay 

        ! vars init
        ilay        = 0
        this % k    = 0.0d0

        if (this % isInit()) then 
            do ilay = 1, this % nlayers 
                call this % layers(ilay) % computeStiffness
                call assign_and_increment(this % k, this % layers(ilay) % k)
            end do
        else
            call sroXIT(" >>>>>> layeredelm3D_computeStiffness <<<<<< Object not initialised ")
        end if
    
    end subroutine

    function layeredelm3D_computeLayerAdjoint(this, ilay) result(adj)
        implicit none

        class(LayeredElm3D)         :: this
        integer,        intent(in)  :: ilay

        real(real64)                :: adj(this % ndof)

        real(real64) :: t
        real(real64) :: tb
        real(real64) :: tt
        integer      :: i
        integer      :: db(this % ndof / 2)
        integer      :: dt(this % ndof / 2)
        real(real64) :: ub(this % ndof / 2)
        real(real64) :: ut(this % ndof / 2)
        real(real64) :: uzb(this % ndof / 2)
        real(real64) :: uzt(this % ndof / 2)

        ! REVIEW: update for adjoint

        ! vars init
        adj = 0.0d0
        t = 0.0d0
        tb = 0.0d0
        tt = 0.0d0
        i = 0
        ub = 0.0d0 
        ut = 0.0d0 
        uzb = 0.0d0
        uzt = 0.0d0
        db = 0 
        dt = 0

        t = 0.0d0
        do i = 1, this % nlayers
            t = t + this % layers(i) % t
        end do
        tb = this % layers(ilay) % z / t
        tt = (this % layers(ilay) % z + this % layers(ilay) % t) / t

        db = (/(i, i=1, this % ndof / 2)/)
        dt = (/(i, i=this % ndof / 2 + 1, this % ndof)/)

        ub = this % adj(db)
        ut = this % adj(dt)
        uzb = ((1 - tb) * ub) + (tb * ut)
        uzt = ((1 - tt) * ub) + (tt * ut)

        adj(db) = uzb
        adj(dt) = uzt

    end function

    subroutine layeredelm3D_computeStrainStressAdjoint(this)
        implicit none
    
        class(LayeredElm3D), intent(inout) :: this
    
        integer :: ilay 

        ! vars init 
        ilay = 0

        if (this % isInit()) then 
            do ilay = 1, this % nlayers 
                call this % layers(ilay) % computeStrainStressAdjoint( &
                        this % computeLayerAdjoint(ilay)   &
                    ) 
            end do
        else
            call sroXIT(" >>>>>> layeredelm3D_computeStrainStressAdjoint <<<<<< Object not initialised ")
        end if
    
    end subroutine
    
    function layeredelm3D_isInit(this) result(b)
        implicit none
    
        logical :: b
        class(LayeredElm3D) :: this

        ! vars init
        b = .false.
    
        if ( (this % nlayers > 0) .and. &
             allocated(this % localCS) .and. &
             allocated(this % layers) .and. &
             allocated(this % adj) .and. &
             this % areLayersInit() .and. &
             this % Element % isInit()) then 

            b = .true.
        else
            call sroLog(" >> Sro::LayeredElm3D << | isInit | Attempted use before initialisation ")
            b = .false.
        end if
    
    end function

    function layeredelm3D_areLayersInit(this) result(b)
        implicit none 

        logical :: b 
        class(LayeredElm3D) :: this 

        integer :: ilay 

        ! vars init
        ilay    = 0
        b       = .true.

        do ilay = 1, this % nlayers
            if (.not. this % layers(ilay) % isInit()) b = .false.
        end do
    end function

    subroutine layeredelm3D_final(this)
        implicit none
    
        class(LayeredElm3D), intent(inout) :: this

        integer :: ilay
    
        ! vars init
        ilay = 0

        this % nlayers = -1
        this % volume = 0.0d0 
        if (allocated(this % localCS)) deallocate(this % localCS) 
        if (allocated(this % layers)) deallocate(this % layers)
        if (allocated(this % adj)) deallocate(this % adj)
        
        do ilay = 1, this % nlayers 
            call this % layers(ilay) % finalize 
        end do

        call this % Element % finalize
    
    end subroutine

    function layeredelm3D_computeMeanStrain(this) result(e)
        implicit none
    
        class(LayeredElm3D) :: this
        real(real64) :: e(this % nstr)
    
        integer :: ilay, istr, ipt
        real(real64) :: emax

        ! vars init
        ilay    = 0
        istr    = 0
        e       = 0.0d0
        emax    = 0.0d0
        ipt     = 0
        
        if (this % isInit()) then
            
            do ilay = 1, this % nlayers
                e = e + (sum(this % layers(ilay) % e, 1) / this % npts)
            end do
            e = e / this % nlayers

            do istr = 1, this % nstr 
                if (e(istr) /= e(istr)) then 
                    call sroXIT(" >>>>>> layeredelm3D_computeMeanStrain <<<<<< NaN detected")
                end if
            end do

        else 
            call sroXIT(" >>>>>> layeredelm3D_computeMeanStrain <<<<<< Object not initialised ")
        end if

    end function

    function layeredelm3D_computeMeanStress(this) result(s)
        implicit none
    
        class(LayeredElm3D) :: this
        real(real64) :: s(this % nstr)
    
        integer :: ilay, istr

        ! vars init
        ilay    = 0
        istr    = 0
        s       = 0.0d0
        
        if (this % isInit()) then
            
            do ilay = 1, this % nlayers
                s = s + (sum(this % layers(ilay) % s, 1) / this % npts)
            end do
            s = s / this % nlayers
            do istr = 1, this % nstr 
                if (s(istr) /= s(istr)) then
                    call sroXIT(" >>>>>> layeredelm3D_computeMeanStress <<<<<< NaN detected")
                end if
            end do

        else 
            call sroXIT(" >>>>>> layeredelm3D_computeMeanStress <<<<<< Object not initialised ")
        end if

    end function

    function layeredelm3D_computeMeanStrainAdjoint(this) result(e)
        implicit none
    
        class(LayeredElm3D) :: this
        real(real64) :: e(this % nstr)
    
        integer :: ilay, istr, ipt
        real(real64) :: emax

        ! vars init
        ilay    = 0
        istr    = 0
        e       = 0.0d0
        emax    = 0.0d0
        ipt     = 0
        
        if (this % isInit()) then
            
            do ilay = 1, this % nlayers
                e = e + (sum(this % layers(ilay) % eadj, 1) / this % npts)
            end do
            e = e / this % nlayers

            do istr = 1, this % nstr 
                if (e(istr) /= e(istr)) then 
                    call sroXIT(" >>>>>> layeredelm3D_computeMeanStrainAdjoint <<<<<< NaN detected")
                end if
            end do

        else 
            call sroXIT(" >>>>>> layeredelm3D_computeMeanStrainAdjoint <<<<<< Object not initialised ")
        end if

    end function

    function layeredelm3D_computeMeanStressAdjoint(this) result(s)
        implicit none
    
        class(LayeredElm3D) :: this
        real(real64) :: s(this % nstr)
    
        integer :: ilay, istr

        ! vars init
        ilay    = 0
        istr    = 0
        s       = 0.0d0
        
        if (this % isInit()) then
            
            do ilay = 1, this % nlayers
                s = s + (sum(this % layers(ilay) % sadj, 1) / this % npts)
            end do
            s = s / this % nlayers
            do istr = 1, this % nstr 
                if (s(istr) /= s(istr)) then
                    call sroXIT(" >>>>>> layeredelm3D_computeMeanStressAdjoint <<<<<< NaN detected")
                end if
            end do

        else 
            call sroXIT(" >>>>>> layeredelm3D_computeMeanStressAdjoint <<<<<< Object not initialised ")
        end if

    end function


end module