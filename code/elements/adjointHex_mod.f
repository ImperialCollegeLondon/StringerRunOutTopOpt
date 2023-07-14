!DEC$ FREEFORM

module adjointHex_mod
    use iso_fortran_env
    use element_mod
    implicit none
    
    type, extends(Element) :: AdjHex
 
        real(real64), allocatable   :: adj(:)        

        real(real64)                :: G = 0.0d0
        real(real64), allocatable   :: dkda(:,:)

        integer, allocatable        :: dofMap(:)
    
    contains

        procedure :: print => adjhex_print
        
        procedure :: isInit => adjhex_isInit
        procedure :: finalize => adjhex_final

        procedure :: getDisplacement => adjhex_getDisplacement
        procedure :: getAdjoint => adjhex_getAdjoint        

        procedure :: setDof => adjhex_setDof
        procedure :: getDof => adjhex_getDof

        procedure :: computeG => adjhex_computeG

        procedure :: update => adjhex_update

        final :: adjhex_dtor
    
    end type

    interface AdjHex
        procedure :: adjhex_ctor
    end interface
    
contains

    function adjhex_ctor(id, coords, con, k, dkda) result(this)
        implicit none

        type(AdjHex) :: this

        integer, intent(in) :: id 
        real(real64), intent(in) :: coords(:,:)
        integer, intent(in) :: con(:)
        real(real64), intent(in) :: k(:,:)
        real(real64), intent(in) :: dkda(:,:)

        this % name = "AdjHex"

        this % nnodes       = 8
        this % ndof         = 24
        this % ndofnode     = 3
        this % ndim         = 3
        this % nstr         = 0
        this % npts         = 0
        this % nprops       = 1

        call this % Element % initialize(id, coords, con, (/0.0d0/))

        allocate(this % adj(this % ndof), source = 0.0d0)
        allocate(this % dkda(this % ndof, this % ndof), source=0.0d0)

        call assign(this % k, k)
        call assign(this % dkda, dkda)

        allocate(                   &
            this % dofMap(48),  &
            source = (/         &
                 1,  2,  3,     &
                 7,  8,  9,     &
                13, 14, 15,     & 
                19, 20, 21,     &
                25, 26, 27,     &
                31, 32, 33,     &
                37, 38, 39,     &
                43, 44, 45,     &
                 4,  5,  6,     &
                10, 11, 12,     &
                16, 17, 18,     &
                22, 23, 24,     &
                28, 29, 30,     &
                34, 35, 36,     &
                40, 41, 42,     &
                46, 47, 48      &
            /))

    end function

    subroutine adjhex_dtor(this)
        implicit none
    
        type(AdjHex), intent(inout) :: this
    
        call this % finalize 
    end subroutine

    subroutine adjhex_final(this)
        implicit none
    
        class(AdjHex), intent(inout) :: this
    
        if (allocated(this % dkda)) deallocate(this % dkda)
        if (allocated(this % dofMap)) deallocate(this % dofMap)
        if (allocated(this % adj)) deallocate(this % adj)

        call this % Element % finalize
    
    end subroutine
    
    function adjhex_isInit(this) result(b)
        implicit none 

        logical :: b 
        class(AdjHex) :: this

        if (allocated(this % dkda) .and. &
            allocated(this % dofMap) .and. &
            allocated(this % adj) .and. &
            this % Element % isInit() ) then 
            b = .true.
        else
            call sroLog(" >> Sro::AdjHex << | isInit | Attempted use before initialisation ")
            b = .false.
        end if

    end function

    function adjhex_getDisplacement(this) result(u)
        implicit none 

        class(AdjHex) :: this
        real(real64) :: u(this % nnodes, this % ndim)

        integer :: i, ids(this % nnodes)

        i = 0
        ids = 0

        ids = (/ (this % ndim * (i - 1) + 1, i = 1, this % nnodes) /)

        u(:,1) = this % u(ids)
        u(:,2) = this % u(ids + 1)
        u(:,3) = this % u(ids + 2)

    end function


    function adjhex_getAdjoint(this) result(p)
        implicit none 

        class(AdjHex) :: this
        real(real64) :: p(this % nnodes, this % ndim)

        integer :: i, ids(this % nnodes)

        i = 0
        ids = 0

        ids = (/ (this % ndim * (i - 1) + 1, i = 1, this % nnodes) /)

        p(:,1) = this % adj(ids)
        p(:,2) = this % adj(ids + 1)
        p(:,3) = this % adj(ids + 2)

    end function

    subroutine adjhex_setDof(this, q)
        implicit none
    
        class(AdjHex), intent(inout) :: this
        real(real64) :: q(:)
    
        integer :: i, udof(24), adof(24)

        udof = this % dofMap((/ (i, i=1     , 24    ) /))
        adof = this % dofMap((/ (i, i=1+24  , 24+24 ) /))

        call assign(this % u, q(udof))
        call assign(this % adj, q(adof))
    
    end subroutine

    function adjhex_getDof(this) result(q)
        implicit none

        class(AdjHex) :: this
        real(real64) :: q(this % ndof * 2)

        integer :: i, udof(24), adof(24)

        udof = this % dofMap((/ (i, i=1     , 24    ) /))
        adof = this % dofMap((/ (i, i=1+24  , 24+24 ) /))

        q(udof) = this % u 
        q(adof) = this % adj

    end function

    subroutine adjhex_computeG(this)
        implicit none
    
        class(AdjHex), intent(inout) :: this
    
        this % G = -0.5d0 * dot_product(this % u, matmul(this % dkda, this % u))

    end subroutine

    subroutine adjhex_update(this, u, k, rhs, analysisType, unit)
        implicit none
    
        class(AdjHex),      intent(inout)   :: this
        real(real64),       intent(in)      :: u(:)
        real(real64),       intent(inout)   :: k(:,:)
        real(real64),       intent(inout)   :: rhs(:)
        character(len=*),   intent(in)      :: analysisType
        integer, optional,  intent(in)      :: unit

        
        real(real64) :: kel(48, 48), fel(48), fex(48), eye(24, 24)

        integer :: i, j, isub, idx(1)
        integer :: udof(24), adof(24)

        ! vars init
        kel     = 0.0d0
        eye     = 0.0d0
        fel     = 0.0d0
        fex     = 0.0d0
        i       = 0
        j       = 0
        isub    = 0
        udof    = 0
        adof    = 0

        udof = this % dofMap((/ (i, i=1     , 24    ) /))
        adof = this % dofMap((/ (i, i=1+24  , 24+24 ) /))

        do i = 1, 24
            eye(i, i) = 1.0d0
        end do
        kel(adof, adof) = eye

        idx = findloc(sro % adjelmap, this % id)
        sro % adjdof(idx(1), :) = u 
        ! if (idx(1) /= 0) then
        !     sro % adjdof(idx(1), :) = u 
        ! else
        !     call sroXIT(" >>>>>> adjhex_update <<<<<< Unable to match element index with sro.elmap (" // int2str(this % id) // ")")
        ! end if

        ! if (this % id == 1) then
        !     call sroLog(" >> Sro::AdjHex    << | update | " // analysisType)
        !     call sroLog(" >> Sro::AdjHex    << | update | u_max = " // float2str(maxval(u(udof))))
        !     call sroLog(" >> Sro::AdjHex    << | update | u_min = " // float2str(minval(u(udof))))
        !     call sroLog(" >> Sro::AdjHex    << | update | a_max = " // float2str(maxval(u(adof))))
        !     call sroLog(" >> Sro::AdjHex    << | update | a_min = " // float2str(minval(u(adof))))
        ! end if

        if (analysisType == "displacement") then

            sro % adjdof(idx(1), adof) = this % adj
            ! if (idx(1) /= 0) then 
            !     sro % adjdof(idx(1), adof) = this % adj
            ! end if
            
            fex(adof)       = this % adj

        else if (analysisType == "adjoint") then

            sro % adjdof(idx(1), udof) = this % u
            ! if (idx(1) /= 0) then
            !     sro % adjdof(idx(1), udof) = this % u
            ! end if

            kel(adof, adof) = this % k
            fex(adof)       = matmul(this % dkda, this % u)
            
        else

            sro % adjdof(idx(1), udof) = this % u
            sro % adjdof(idx(1), adof) = this % adj
            ! if (idx(1) /= 0) then
            !     sro % adjdof(idx(1), udof) = this % u
            !     sro % adjdof(idx(1), adof) = this % adj
            ! end if
            
            fex(adof)       = this % adj

        end if

        fel     = fex - matmul(kel, u)

        k       = kel
        rhs     = fel
    end subroutine

    subroutine adjhex_print(this, varName, unit)
        implicit none
    
        class(AdjHex),      intent(inout)   :: this
        character(len=*),   intent(in)      :: varName
        integer, optional,  intent(in)      :: unit
    
        ! TODO print 
    
    end subroutine
   
end module