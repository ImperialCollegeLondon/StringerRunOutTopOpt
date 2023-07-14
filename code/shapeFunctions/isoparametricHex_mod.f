!DEC$ FREEFORM 

module isoparametricHex_mod
    use iso_fortran_env
    use shapeFunctionLayered_mod
    implicit none
    
    type, extends(ShapeFunctionLayered) :: HexIsoSF
    contains
        procedure :: N_at  => hexiso_N_at
        procedure :: dN_at => hexiso_dN_at

        final     :: hexiso_dtor
    end type

    interface HexIsoSF
        procedure :: hexiso_ctor
    end interface
    
contains

    function hexiso_ctor(x, xl, ptsl) result(this)
        implicit none 

        type(HexIsoSF)           :: this
        real(real64), intent(in) :: x(:,:)
        real(real64), intent(in) :: xl(:,:), ptsl(:,:)

        this % name = "HexIsoSF"

        this % nnodes = 8 
        this % ndim = 3 
        this % ndof = 24
        this % nstr = 6 
        this % npts = size(ptsl,1)

        call this % initialize(x, xl, ptsl)
    end function

    subroutine hexiso_dtor(this)
        implicit none
        
        type(HexIsoSF), intent(inout) :: this
        
        call this % finalize
    end subroutine

    function hexiso_N_at(this, xi, eta, zeta, n) result(nmat)
        implicit none
    
        class(HexIsoSF), intent(inout) :: this
        integer,         intent(in)    :: n
        real(real64),    intent(in)    :: xi, eta, zeta
        real(real64)                   :: nmat(1, n)

        ! vars init
        nmat = 0.0d0
        
        nmat(1,1) = (1.0d0 / 8.0d0) * (1 - xi) * (1 - eta) * (1 - zeta)
        nmat(1,2) = (1.0d0 / 8.0d0) * (1 + xi) * (1 - eta) * (1 - zeta)
        nmat(1,3) = (1.0d0 / 8.0d0) * (1 + xi) * (1 + eta) * (1 - zeta)
        nmat(1,4) = (1.0d0 / 8.0d0) * (1 - xi) * (1 + eta) * (1 - zeta)
        nmat(1,5) = (1.0d0 / 8.0d0) * (1 - xi) * (1 - eta) * (1 + zeta)
        nmat(1,6) = (1.0d0 / 8.0d0) * (1 + xi) * (1 - eta) * (1 + zeta)
        nmat(1,7) = (1.0d0 / 8.0d0) * (1 + xi) * (1 + eta) * (1 + zeta)
        nmat(1,8) = (1.0d0 / 8.0d0) * (1 - xi) * (1 + eta) * (1 + zeta)
    end function

    function hexiso_dN_at(this, xi, eta, zeta, n, m) result(dnmat)
        implicit none
    
        class(HexIsoSF), intent(inout) :: this
        integer,         intent(in)    :: n, m
        real(real64),    intent(in)    :: xi, eta, zeta
        real(real64)                   :: dnmat(n, m)

        ! vars init
        dnmat = 0.0d0
        
        dnmat(1, 1) = -(1.0d0 / 8.0d0) * (1 - eta) * (1 - zeta)
        dnmat(1, 2) =  (1.0d0 / 8.0d0) * (1 - eta) * (1 - zeta)
        dnmat(1, 3) =  (1.0d0 / 8.0d0) * (1 + eta) * (1 - zeta)
        dnmat(1, 4) = -(1.0d0 / 8.0d0) * (1 + eta) * (1 - zeta)
        dnmat(1, 5) = -(1.0d0 / 8.0d0) * (1 - eta) * (1 + zeta)
        dnmat(1, 6) =  (1.0d0 / 8.0d0) * (1 - eta) * (1 + zeta)
        dnmat(1, 7) =  (1.0d0 / 8.0d0) * (1 + eta) * (1 + zeta)
        dnmat(1, 8) = -(1.0d0 / 8.0d0) * (1 + eta) * (1 + zeta)
        dnmat(2, 1) = -(1.0d0 / 8.0d0) * (1 - xi) * (1 - zeta)
        dnmat(2, 2) = -(1.0d0 / 8.0d0) * (1 + xi) * (1 - zeta)
        dnmat(2, 3) =  (1.0d0 / 8.0d0) * (1 + xi) * (1 - zeta)
        dnmat(2, 4) =  (1.0d0 / 8.0d0) * (1 - xi) * (1 - zeta)
        dnmat(2, 5) = -(1.0d0 / 8.0d0) * (1 - xi) * (1 + zeta)
        dnmat(2, 6) = -(1.0d0 / 8.0d0) * (1 + xi) * (1 + zeta)
        dnmat(2, 7) =  (1.0d0 / 8.0d0) * (1 + xi) * (1 + zeta)
        dnmat(2, 8) =  (1.0d0 / 8.0d0) * (1 - xi) * (1 + zeta)
        dnmat(3, 1) = -(1.0d0 / 8.0d0) * (1 - xi) * (1 - eta)
        dnmat(3, 2) = -(1.0d0 / 8.0d0) * (1 + xi) * (1 - eta)
        dnmat(3, 3) = -(1.0d0 / 8.0d0) * (1 + xi) * (1 + eta)
        dnmat(3, 4) = -(1.0d0 / 8.0d0) * (1 - xi) * (1 + eta)
        dnmat(3, 5) =  (1.0d0 / 8.0d0) * (1 - xi) * (1 - eta)
        dnmat(3, 6) =  (1.0d0 / 8.0d0) * (1 + xi) * (1 - eta)
        dnmat(3, 7) =  (1.0d0 / 8.0d0) * (1 + xi) * (1 + eta)
        dnmat(3, 8) =  (1.0d0 / 8.0d0) * (1 - xi) * (1 + eta)
    end function
   
end module