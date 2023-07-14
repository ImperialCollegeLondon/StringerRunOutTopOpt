!DEC$ FREEFORM

module integrator_mod
    use iso_fortran_env
    use component_mod
    implicit none
   
    type, extends(Component), abstract :: Integrator
        
        integer                     :: npts = -1
        integer                     :: ndim = -1
        
        real(real64), allocatable   :: pts(:,:) 
        real(real64), allocatable   :: wts(:,:)

    contains

        procedure(igt_setPointsWeights), deferred :: setPointsWeights

        procedure :: print      => igt_print

        procedure :: initialize => igt_init
        procedure :: finalize   => igt_final 

        procedure :: isInit     => igt_isInit 

    end type

    interface

        subroutine igt_setPointsWeights(this)
            import Integrator
            implicit none
            class(Integrator), intent(inout) :: this 
        end subroutine

    end interface
   
contains
   
    subroutine igt_print(this, varName, unit)
        implicit none
        
        class(Integrator),  intent(inout)   :: this
        character(len=*),   intent(in)      :: varName
        integer, optional,  intent(in)      :: unit
        
        integer :: uid
        
        ! vars init
        uid = 0

        uid = 6
        if (present(unit)) uid = unit 
        
        if (this % isInit()) then
            write(uid, fmt100) ""
            write(uid, fmt100) varName // " < Integrator | " // this % name // " >"
            write(uid, fmt100) ""
            write(uid, fmt101) "npts", this % npts 
            write(uid, fmt101) "ndim", this % ndim 
            write(uid, fmt100) ""
            call printVar(this % pts, "pts", uid)
            call printVar(this % wts, "wts", uid)
        else 
            call sroXIT(" >>>>>> igt_print <<<<<< Object not initialised")
        end if
        
    end subroutine

    function igt_isInit(this) result(b)
        implicit none

        logical             :: b
        class(Integrator)   :: this

        ! vars init
        b = .false.
    
        if ( (this % npts /= -1) .and. &
             (this % ndim /= -1) .and. &
             allocated(this % pts) .and. &
             allocated(this % wts) ) then 

            b = .true. 
        else 
            call sroLog(" >> Sro::Integrator << | isInit | Attempted use before initialisation")
            b = .false. 
        end if    
    end function

    subroutine igt_init(this, numPts, numDim)
        implicit none
    
        class(Integrator), intent(inout)    :: this
        integer, intent(in)                 :: numPts, numDim
    
        this % npts = numPts 
        this % ndim = numDim

        allocate(this % pts(this % npts, this % ndim), source = 0.0d0)
        allocate(this % wts(this % npts, this % ndim), source = 0.0d0)

        call this % setPointsWeights
    end subroutine

    subroutine igt_final(this)
        implicit none
    
        class(Integrator), intent(inout) :: this
    
        this % npts = -1
        this % ndim = -1
        if (allocated(this % pts)) deallocate(this % pts)
        if (allocated(this % wts)) deallocate(this % wts)
    end subroutine

end module