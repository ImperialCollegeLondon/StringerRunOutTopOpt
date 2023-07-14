!DEC$ FREEFORM

module energyReleaseRateHex_mod
    use iso_fortran_env
    use element_mod
    implicit none
    
    type, extends(Element) :: EnrrHex 

        real(real64)              :: G = 0.0d0
        real(real64), allocatable :: dkda(:,:)

    contains 

        procedure :: print      => enrrhex_print 
        
        procedure :: finalize   => enrrhex_final

        procedure :: computeG   => enrrhex_computeG
        procedure :: update     => enrrhex_update

        procedure :: isInit     => enrrhex_isInit 

        final :: enrrhex_dtor

    end type

    interface EnrrHex
        procedure :: enrrhex_ctor
    end interface
    
contains

    function enrrhex_ctor(id, coords, con, dkda) result(this)
        implicit none 

        type(EnrrHex)                       :: this 
        integer,            intent(in)      :: id
        real(real64),       intent(in)      :: coords(:,:)
        integer,            intent(in)      :: con(:)
        real(real64),       intent(in)      :: dkda(:,:)

        this % name = "EnrrHex"

        this % nnodes       = 8
        this % ndof         = 24
        this % ndofnode     = 3
        this % ndim         = 3
        this % nstr         = 6
        this % npts         = 5
        this % nprops       = 1

        call this % initialize(id, coords, con, (/0.0d0/))

        allocate(this % dkda(this % ndof, this % ndof), source = dkda)

    end function 

    subroutine enrrhex_dtor(this)
        implicit none
    
        type(EnrrHex), intent(inout) :: this
    
        call this % finalize
    end subroutine

    subroutine enrrhex_print(this, varName, unit)
        implicit none
        
        class(EnrrHex),     intent(inout)   :: this
        character(len=*),   intent(in)      :: varName
        integer, optional,  intent(in)      :: unit
        
        integer :: uid

        ! vars init
        uid = 0
        
        uid = 6
        if (present(unit)) uid = unit 
        
        if (this % isInit()) then
            call this % Element % print(varName, uid) 
            write(uid, fmt100) ""
            call printVar(this % dkda, "dkda", uid)
            write(uid, fmt100) ""
        else 
            call sroXIT(" >>>>>> enrrhex_print <<<<<< Object not initialised ")
        end if 

    end subroutine

    subroutine enrrhex_final(this)
        implicit none
    
        class(EnrrHex), intent(inout) :: this
    
        if (allocated(this % dkda)) deallocate(this % dkda)

        call this % Element % finalize
    end subroutine

    function enrrhex_isInit(this) result(b)
        implicit none 

        logical :: b 
        class(EnrrHex) :: this 

        ! vars init
        b = .false.

        if (allocated(this % dkda) .and. &
            this % Element % isInit()) then 

            b = .true.
        else
            call sroLog(" >> Sro::EnrrHex << | isInit | Attempted use before initialisation ")
            b = .false.
        end if

    end function

    subroutine enrrhex_computeG(this)
        implicit none
    
        class(EnrrHex), intent(inout) :: this
    
        ! vars init
        this % G = 0.0d0

        if (this % isInit()) then 
            this % G = dot_product(this % u, matmul(this % u, this % dkda))
        else 
            call sroXIT(" >>>>>> enrrhex_computeG <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine enrrhex_update(this, u)
        implicit none
    
        class(EnrrHex),     intent(inout)   :: this
        real(real64),       intent(in)      :: u(:)

        if (this % isInit()) then
            call assign(this % u, u)
            call this % computeG
        else 
            call sroXIT(" >>>>>> enrrhex_update <<<<<< Object not initialised ")
        end if
    
    end subroutine
   
end module