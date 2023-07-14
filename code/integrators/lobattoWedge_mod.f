!DEC$ FREEFORM

module lobattoWedge_mod
    use iso_fortran_env
    use integrator_mod
    implicit none

    type, extends(Integrator) :: LobattoIntegrator3DWedge
    contains
        procedure :: setPointsWeights => lobattoWed_setPointsWeights
    end type
    
    interface LobattoIntegrator3DWedge
        procedure :: lobattoWed_ctor
    end interface
    
contains

    function lobattoWed_ctor(numPts) result(this)
        implicit none 

        type(LobattoIntegrator3DWedge)   :: this
        integer, intent(in), optional    :: numPts

        integer :: npts
        
        ! vars init
        npts = 5

        this % name = "LobattoIntegrator3DWedge"

        if (present(numPts)) npts = numPts

        call this % initialize(npts, 3)

    end function

    subroutine lobattoWed_setPointsWeights(this)
        implicit none
    
        class(LobattoIntegrator3DWedge), intent(inout) :: this
    
        if (this % isInit()) then

            ! vars init
            this % pts = 0.0d0
            this % wts = 0.0d0

            this % pts(:,1:2) = 1.0d0/3.0d0
            ! this % wts(:,1:2) = 1.0d0/2.0d0
            this % wts(:,1:2) = 1.0d0

            select case (this % npts)
            case(1)
                this % pts(1,3) = 0.0d0
                this % wts(1,3) = 2.0d0
            case(3)
                this % pts(1,3) = -1.0d0
                this % pts(2,3) =  0.0d0
                this % pts(3,3) =  1.0d0
                this % wts(1,3) = 1.0d0 / 3.0d0
                this % wts(2,3) = 4.0d0 / 3.0d0
                this % wts(3,3) = 1.0d0 / 3.0d0
            case(4)
                this % pts(1,3) = -1.0d0
                this % pts(2,3) = -(1.0d0/5.0d0) * sqrt(5.0d0)
                this % pts(3,3) =  (1.0d0/5.0d0) * sqrt(5.0d0)
                this % pts(4,3) =  1.0d0
                this % wts(1,3) = 1.0d0 / 6.0d0
                this % wts(2,3) = 5.0d0 / 6.0d0
                this % wts(3,3) = 5.0d0 / 6.0d0
                this % wts(4,3) = 1.0d0 / 6.0d0
            case(5)
                this % pts(1,3) = -1.0d0
                this % pts(2,3) = -(1.0d0/7.0d0) * sqrt(21.0d0)
                this % pts(3,3) =  0.0d0
                this % pts(4,3) =  (1.0d0/7.0d0) * sqrt(21.0d0)
                this % pts(5,3) =  1.0d0
                this % wts(1,3) =  1.0d0 / 10.0d0
                this % wts(2,3) = 49.0d0 / 90.0d0
                this % wts(3,3) = 32.0d0 / 45.0d0
                this % wts(4,3) = 49.0d0 / 90.0d0
                this % wts(5,3) =  1.0d0 / 10.0d0
            case(6)
                this % pts(1,3) = -1.0d0
                this % pts(2,3) = -sqrt((1.0d0/21.0d0) * (7.0d0 + 2.0d0 * sqrt(7.0d0)))
                this % pts(3,3) = -sqrt((1.0d0/21.0d0) * (7.0d0 - 2.0d0 * sqrt(7.0d0)))
                this % pts(4,3) =  sqrt((1.0d0/21.0d0) * (7.0d0 - 2.0d0 * sqrt(7.0d0)))
                this % pts(5,3) =  sqrt((1.0d0/21.0d0) * (7.0d0 + 2.0d0 * sqrt(7.0d0)))
                this % pts(6,3) =  1.0d0
                this % wts(1,3) =  1.0d0 / 15.0d0                        
                this % wts(2,3) = (1.0d0 / 30.0d0) * (14.0d0 - sqrt(7.0d0))
                this % wts(3,3) = (1.0d0 / 30.0d0) * (14.0d0 + sqrt(7.0d0))
                this % wts(4,3) = (1.0d0 / 30.0d0) * (14.0d0 + sqrt(7.0d0))
                this % wts(5,3) = (1.0d0 / 30.0d0) * (14.0d0 - sqrt(7.0d0))
                this % wts(6,3) =  1.0d0 / 15.0d0  
            case default
                this % pts(1,3) = -1.0d0
                this % pts(2,3) = -(1.0d0/7.0d0) * sqrt(21.0d0)
                this % pts(3,3) =  0.0d0
                this % pts(4,3) =  (1.0d0/7.0d0) * sqrt(21.0d0)
                this % pts(5,3) =  1.0d0
                this % wts(1,3) =  1.0d0 / 10.0d0
                this % wts(2,3) = 49.0d0 / 90.0d0
                this % wts(3,3) = 32.0d0 / 45.0d0
                this % wts(4,3) = 49.0d0 / 90.0d0
                this % wts(5,3) =  1.0d0 / 10.0d0             
            end select
        else
            call sroXIT(" >>>>>> lobattoWed_setPointsWeights <<<<<< Object not initialised ")
        end if
    
    end subroutine
   
end module