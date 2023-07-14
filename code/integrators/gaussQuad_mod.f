!DEC$ FREEFORM

module gaussQuad_mod
    use iso_fortran_env
    use integrator_mod
    implicit none
    
    type, extends(Integrator) :: GaussIntegrator2DQuad
    contains
        procedure :: setPointsWeights => gaussQuad_setPointsWeights
    end type
    
    interface GaussIntegrator2DQuad
        procedure :: gaussQuad_ctor
    end interface
    
contains

    function gaussQuad_ctor(numPts) result(this)
        implicit none 

        type(GaussIntegrator2DQuad)     :: this 
        integer, intent(in), optional   :: numPts

        integer :: npts
        
        ! vars init
        npts = 5

        this % name = "GaussIntegrator2DQuad"

        if (present(numPts)) npts = numPts 

        call this % initialize(npts, 2)
    end function 

    subroutine gaussQuad_setPointsWeights(this)
        implicit none
    
        class(GaussIntegrator2DQuad), intent(inout) :: this
    
        if (this % isInit()) then

            ! vars init
            this % pts = 0.0d0
            this % wts = 0.0d0

            select case (this % npts)
            case(1)
                this % pts(1,:) = 0.0d0
                this % wts(1,:) = 2.0d0
            case(4)
                this % pts(1,:) = (/ -1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0) /)
                this % pts(2,:) = (/  1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0) /)
                this % pts(3,:) = (/  1.0d0/sqrt(3.0d0),  1.0d0/sqrt(3.0d0) /)
                this % pts(4,:) = (/ -1.0d0/sqrt(3.0d0),  1.0d0/sqrt(3.0d0) /)
                this % wts(:,:) = 1.0d0
            case default
                this % pts(1,:) = (/ -1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0) /)
                this % pts(2,:) = (/  1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0) /)
                this % pts(3,:) = (/  1.0d0/sqrt(3.0d0),  1.0d0/sqrt(3.0d0) /)
                this % pts(4,:) = (/ -1.0d0/sqrt(3.0d0),  1.0d0/sqrt(3.0d0) /)
                this % wts(:,:) = 1.0d0
            end select 
        else 
            call sroXIT(" >>>>>> gaussQuad_setPointsWeights <<<<<< Object not initialised ")
        end if
    end subroutine
   
end module