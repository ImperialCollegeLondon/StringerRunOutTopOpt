!DEC$ FREEFORM

module visualizationHex_mod
    use iso_fortran_env
    use element_mod
    implicit none
    
    type, extends(Element) :: VizHex 
    contains 
        procedure :: print      => vizhex_print 
        
        procedure :: finalize   => vizhex_final

        procedure :: update     => vizhex_update

        procedure :: isInit     => vizhex_isInit 

        final :: vizhex_dtor

    end type

    interface VizHex
        procedure :: vizhex_ctor
    end interface
    
contains

    function vizhex_ctor(id, coords, con) result(this)
        implicit none 

        type(VizHex)                        :: this 
        integer,            intent(in)      :: id
        real(real64),       intent(in)      :: coords(:,:)
        integer,            intent(in)      :: con(:)

        this % name = "VizHex"

        this % nnodes       = 8
        this % ndof         = 24
        this % ndofnode     = 3
        this % ndim         = 3
        this % nstr         = 6
        this % npts         = 5
        this % nprops       = 1

        call this % initialize(id, coords, con, (/0.0d0/))

    end function 

    subroutine vizhex_dtor(this)
        implicit none
    
        type(VizHex), intent(inout) :: this
    
        call this % finalize
    end subroutine

    subroutine vizhex_print(this, varName, unit)
        implicit none
        
        class(VizHex), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid

        ! vars init
        uid = 0
        
        uid = 6
        if (present(unit)) uid = unit 
        
        if (this % isInit()) then
            call this % Element % print(varName, uid)
        else 
            call sroXIT(" >>>>>> vizhex_print <<<<<< Object not initialised ")
        end if 

    end subroutine

    subroutine vizhex_final(this)
        implicit none
    
        class(VizHex), intent(inout) :: this

        call this % Element % finalize
    end subroutine

    function vizhex_isInit(this) result(b)
        implicit none 

        logical :: b 
        class(VizHex) :: this 

        b = .false.

        if (this % Element % isInit()) then 
            b = .true.
        else
            call sroLog(" >> Sro::VizHex << | isInit | Attempted use before initialisation ")
            b = .false.
        end if

    end function

    subroutine vizhex_update(this, u, analysisType)
        implicit none
    
        class(VizHex),  intent(inout)   :: this
        real(real64),   intent(in)      :: u(:)
        character(len=*),   intent(in)      :: analysisType

        if (this % isInit()) then 
            select case(analysisType)
            case("displacement")
                call assign(this % u, u)
            end select
        else 
            call sroXIT(" >>>>>> vizhex_update <<<<<< Object not initialised ")
        end if
    
    end subroutine
   
end module