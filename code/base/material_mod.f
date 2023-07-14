!DEC$ FREEFORM

module material_mod
    use iso_fortran_env
    use component_mod 
    implicit none
   
    type, extends(Component), abstract :: Material
        
        integer                     :: nprops = -1
        integer                     :: nstr   = -1

        real(real64), allocatable   :: props(:)
        real(real64), allocatable   :: D(:,:)

    contains

        procedure(mat_computeD), deferred :: computeD
        
        procedure :: print      => mat_print

        procedure :: initialize => mat_init
        procedure :: finalize   => mat_final
        
        procedure :: isInit     => mat_isInit 

    end type
   
    interface
        subroutine mat_computeD(this)
            import Material
            implicit none
        
            class(Material), intent(inout) :: this
        
        end subroutine
    end interface

contains

    subroutine mat_print(this, varName, unit)
        implicit none
        
        class(Material), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid

        ! vars init
        uid = 0
        
        uid = 6
        if (present(unit)) uid = unit 
        
        if (this % isInit()) then 
            write(uid, fmt100) ""
            write(uid, fmt100) varName // " < Material | " // this % name // " >"
            write(uid, fmt100) ""
            write(uid, fmt101) "nprops", this % nprops
            write(uid, fmt101) "nstr", this % nstr
            write(uid, fmt100) ""
            call printVar(this % props, "props", uid)
            call printVar(this % D, "D", uid)
        else 
            call sroXIT(" >>>>>> mat_print <<<<<< Object not initialised ")
        end if

    end subroutine

    subroutine mat_init(this, props)
        implicit none
    
        class(Material),    intent(inout)   :: this
        real(real64),       intent(in)      :: props(:)
    
        if ( (this % nstr /= -1) .and. &
             (this % nprops /= -1) .and. &
             (size(props) == this % nprops)) then 

            allocate(this % props(this % nprops), source = props)
            allocate(this % D(this % nstr, this % nstr), source = 0.0d0)

            call this % computeD
        else
            call sroXIT(" >>>>>> mat_init <<<<<< Bad input data ")
        end if
    
    end subroutine

    subroutine mat_final(this)
        implicit none
    
        class(Material), intent(inout) :: this
    
        this % nstr = -1 
        this % nprops = -1 
        if (allocated(this % props)) deallocate(this % props)
        if (allocated(this % D)) deallocate(this % D)
    
    end subroutine

    function mat_isInit(this) result(b)
        implicit none
        logical :: b
        class(Material) :: this

        ! vars init
        b = .false.

        if ( (this % nstr /= -1) .and. &
             (this % nprops /= -1) .and. &
             allocated(this % props) .and. &
             allocated(this % D)) then 

            b = .true.
        else
            call sroLog(" >> Sro::Material << | isInit | Attempted use before initialisation ")
            b = .false.
        end if
        
    end function
   
end module