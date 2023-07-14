!DEC$ FREEFORM

module isoparametricQuad_mod
    use iso_fortran_env
    use shapeFunction_mod
    implicit none
    
    type, extends(ShapeFunction) :: QuadIsoSF
    contains
        procedure :: N_at   => quadsf_N_at 
        procedure :: dN_at  => quadsf_dN_at 

        final     :: quadsf_dtor
    end type
    
    interface QuadIsoSF
        procedure :: quadsf_ctor
    end interface
    
contains

    function quadsf_ctor(x, pts, useB) result(this)
        implicit none 

        type(QuadIsoSF) :: this 
        real(real64), intent(in) :: x(:,:), pts(:,:)
        logical, optional, intent(in) :: useB 

        logical :: useBFlag 

        useBFlag = .True.
        if (present(useB)) useBFlag = useB

        this % name = "QuadIsoSF"

        this % nnodes = 4
        this % ndim = 2 
        this % ndof = 4 
        this % nstr = 3 
        this % npts = size(pts,1)

        call this % initialize(x, pts, useBFlag)
    end function

    subroutine quadsf_dtor(this)
        implicit none
    
        type(QuadIsoSF), intent(inout) :: this
    
        call this % finalize
    end subroutine

    function quadsf_N_at(this, xi, eta, zeta, n) result(nmat)
        implicit none
    
        class(QuadIsoSF), intent(inout) :: this
        integer,          intent(in)    :: n
        real(real64),     intent(in)    :: xi, eta, zeta
        real(real64)                    :: nmat(1, n)

        ! vars init
        nmat = 0.0d0
    
        nmat(1,1) = 0.25d0 * (1 - xi) * (1 - eta)
        nmat(1,2) = 0.25d0 * (1 + xi) * (1 - eta)
        nmat(1,3) = 0.25d0 * (1 + xi) * (1 + eta)
        nmat(1,4) = 0.25d0 * (1 - xi) * (1 + eta)
    end function

    function quadsf_dN_at(this, xi, eta, zeta, n, m) result(dnmat)
        implicit none
    
        class(QuadIsoSF), intent(inout) :: this
        integer,          intent(in)    :: n, m
        real(real64),     intent(in)    :: xi, eta, zeta
        real(real64)                    :: dnmat(n, m)

        ! vars init
        dnmat = 0.0d0
    
        dnmat(1,1) = -0.25d0 * (1 - eta)
        dnmat(1,2) = +0.25d0 * (1 - eta)
        dnmat(1,3) = +0.25d0 * (1 + eta)
        dnmat(1,4) = -0.25d0 * (1 + eta)
        dnmat(2,1) = -0.25d0 * (1 - xi)
        dnmat(2,2) = -0.25d0 * (1 + xi)
        dnmat(2,3) = +0.25d0 * (1 + xi)
        dnmat(2,4) = +0.25d0 * (1 - xi)
    end function
   
end module