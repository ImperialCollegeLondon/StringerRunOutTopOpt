!DEC$ FREEFORM

module isoparametricWedge_mod
    use iso_fortran_env
    use shapeFunctionLayered_mod
    implicit none

    type, extends(ShapeFunctionLayered) :: WedgeIsoSF
    contains 
        procedure :: N_at  => wediso_N_at
        procedure :: dN_at => wediso_dN_at

        final     :: wediso_dtor
    end type

    interface WedgeIsoSF
        procedure :: wediso_ctor
    end interface
    
contains

    function wediso_ctor(x, xl, ptsl) result(this)
        implicit none 

        type(WedgeIsoSF)           :: this
        real(real64), intent(in)   :: x(:,:)
        real(real64), intent(in)   :: xl(:,:), ptsl(:,:)

        this % name = "WedgeIsoSF"

        this % nnodes = 6 
        this % ndim = 3 
        this % ndof = 18
        this % nstr = 6 
        this % npts = size(ptsl,1)
    
        call this % initialize(x, xl, ptsl)
    end function

    subroutine wediso_dtor(this)
        implicit none
    
        type(WedgeIsoSF), intent(inout) :: this
    
        call this % finalize
    end subroutine

    function wediso_N_at(this, xi, eta, zeta, n) result(nmat)
        implicit none
    
        class(WedgeIsoSF), intent(inout) :: this
        integer,         intent(in)      :: n
        real(real64),    intent(in)      :: xi, eta, zeta
        real(real64)                     :: nmat(1, n)

        ! vars init
        nmat = 0.0d0
        
        nmat(1,1) = (1.0d0 / 2.0d0) * (1 - zeta) * (1 - xi - eta)
        nmat(1,2) = (1.0d0 / 2.0d0) * (1 - zeta) * xi
        nmat(1,3) = (1.0d0 / 2.0d0) * (1 - zeta) * eta
        nmat(1,4) = (1.0d0 / 2.0d0) * (1 + zeta) * (1 - xi - eta)
        nmat(1,5) = (1.0d0 / 2.0d0) * (1 + zeta) * xi
        nmat(1,6) = (1.0d0 / 2.0d0) * (1 + zeta) * eta
    end function

    function wediso_dN_at(this, xi, eta, zeta, n, m) result(dnmat)
        implicit none
    
        class(WedgeIsoSF), intent(inout) :: this
        integer,         intent(in)      :: n, m
        real(real64),    intent(in)      :: xi, eta, zeta
        real(real64)                     :: dnmat(n, m)

        ! vars init
        dnmat = 0.0d0
        
        dnmat(1,1) = -(1.0d0 / 2.0d0) * (1 - zeta)
        dnmat(1,2) =  (1.0d0 / 2.0d0) * (1 - zeta)
        dnmat(1,3) =  0.0d0
        dnmat(1,4) = -(1.0d0 / 2.0d0) * (1 + zeta)
        dnmat(1,5) =  (1.0d0 / 2.0d0) * (1 + zeta)
        dnmat(1,6) =  0.0d0
        dnmat(2,1) = -(1.0d0 / 2.0d0) * (1 - zeta)
        dnmat(2,2) =  0.0d0
        dnmat(2,3) =  (1.0d0 / 2.0d0) * (1 - zeta)
        dnmat(2,4) = -(1.0d0 / 2.0d0) * (1 + zeta)
        dnmat(2,5) =  0.0d0
        dnmat(2,6) =  (1.0d0 / 2.0d0) * (1 + zeta)
        dnmat(3,1) = -(1.0d0 / 2.0d0) * (1 - xi - eta)
        dnmat(3,2) = -(1.0d0 / 2.0d0) * xi
        dnmat(3,3) = -(1.0d0 / 2.0d0) * eta
        dnmat(3,4) =  (1.0d0 / 2.0d0) * (1 - xi - eta)
        dnmat(3,5) =  (1.0d0 / 2.0d0) * xi
        dnmat(3,6) =  (1.0d0 / 2.0d0) * eta
    end function
   
end module