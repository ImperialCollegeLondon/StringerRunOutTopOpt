!DEC$ FREEFORM

module laminaPSMat_mod
    use iso_fortran_env
    use material_mod
    implicit none
    
    type, extends(Material) :: LaminaPlaneStress

    contains 
        procedure :: computeD => lamPS_computeD

        final :: lamPS_dtor
    end type

    interface LaminaPlaneStress
        procedure :: lamPS_ctor
    end interface

contains

    function lamPS_ctor(properties) result(this)
        implicit none 

        type(LaminaPlaneStress) :: this 
        real(real64) :: properties(:)

        this % name = "LaminaPlaneStress"

        this % nprops = size(properties)
        this % nstr = 6 

        call this % initialize(properties)
    end function

    subroutine lamPS_dtor(this)
        implicit none
    
        type(LaminaPlaneStress), intent(inout) :: this
    
        call this % finalize
    end subroutine

    subroutine lamPS_computeD(this)
        implicit none
    
        class(LaminaPlaneStress), intent(inout) :: this
    
        real(real64) :: dmat(this % nstr, this % nstr)
        real(real64) :: e1, e2, nu12, g12, g13, g23
        real(real64) :: d11, d22, d33, d44, d55, d66, d12
        real(real64) :: aux1

        ! from Abaqus 2019 docs: Linear Elastic Behaviour > Defining orthotropic elasticity in plane stress

        ! vars init
        dmat    = 0.0d0
        e1      = 0.0d0
        e2      = 0.0d0
        nu12    = 0.0d0
        g12     = 0.0d0
        g13     = 0.0d0
        g23     = 0.0d0
        d11     = 0.0d0
        d22     = 0.0d0
        d33     = 0.0d0
        d44     = 0.0d0
        d55     = 0.0d0
        d66     = 0.0d0
        d12     = 0.0d0
        aux1    = 0.0d0

        e1      = this % props(1)
        e2      = this % props(2)
        nu12    = this % props(3)
        g12     = this % props(4)
        g13     = this % props(5)
        g23     = this % props(6)

        aux1 = e1 - (e2 * (nu12 ** 2.0d0))

        d11 = (e1 ** 2.0d0) / aux1
        d22 = (e1 * e2) / aux1
        d33 = e1!min(e1, e2)
        d44 = g12
        d55 = g13
        d66 = g23
        d12 = (e1 * e2 * nu12) / aux1
        
        ! hardcoded values for apparent properties of web layup
        ! dmat(1,1) = 81023.85887d0
        ! dmat(1,2) = 24815.01251d0
        ! dmat(2,1) = 24815.01251d0
        ! dmat(2,2) = 49518.44665d0
        ! dmat(3,3) = 165000.0d0   
        ! dmat(4,4) = 27542.10519d0
        ! dmat(5,5) = 3920.0d0     
        ! dmat(6,6) = 3780.0d0     

        dmat(1,1) = d11
        dmat(1,2) = d12
        dmat(2,1) = d12
        dmat(2,2) = d22
        dmat(3,3) = d33
        dmat(4,4) = d44
        dmat(5,5) = d55
        dmat(6,6) = d66

        call assign(this % D, dmat)
        
    end subroutine
   
end module