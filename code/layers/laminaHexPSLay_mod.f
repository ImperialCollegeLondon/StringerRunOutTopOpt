!DEC$ FREEFORM

module laminaHexPSLay_mod
    use iso_fortran_env
    use layer_mod
    use laminaPSMat_mod
    use lobattoHex_mod
    use isoparametricHex_mod
    implicit none
    
    type, extends(Layer) :: S8PSLayer
    contains
        procedure :: computeX => s8pslay_computeX
        procedure :: computeXg => s8pslay_computeXg

        final :: s8pslay_dtor
    end type

    interface S8PSLayer
        procedure :: s8pslay_ctor
    end interface
    
contains

    function s8pslay_ctor(id, z, t, theta, x, props, Tlg, Tgl) result(this)
        implicit none 

        type(S8PSLayer)                             :: this
        integer,                intent(in)          :: id
        real(real64),           intent(in)          :: z, t, theta
        real(real64),           intent(in)          :: x(:,:), props(:)
        class(Transform3D),     intent(in), target  :: Tlg
        class(Transform3D),     intent(in), target  :: Tgl

        this % name = "S8PSLayer"

        this % nnodes = 8
        this % ndof = 24
        this % ndim = 3 
        this % nstr = 6
        this % npts = 5
        this % nprops = size(props)

        call this % initialize(id, z, t, theta, x, Tlg, Tgl)

        this % material = LaminaPlaneStress(props)
        this % integrator = LobattoIntegrator3DHex(this % npts)
        this % shapeFunction = HexIsoSF(x, this % x, this % integrator % pts)
        call assign(this % ipx, this % shapeFunction % computeIPCoordinates(this % x))

    end function

    subroutine s8pslay_dtor(this)
        implicit none
    
        type(S8PSLayer), intent(inout) :: this
    
        call this % finalize
    
    end subroutine

    subroutine s8pslay_computeX(this)
        implicit none
    
        class(S8PSLayer), intent(inout) :: this
    
        real(real64) :: c(this % nnodes, this % ndim)
        real(real64) :: nodalCoords(this % nnodes, this % ndim)
        real(real64) :: z, t, l1, l2, l3, l4
        real(real64) :: e1(this % ndim), e2(this % ndim), e3(this % ndim), e4(this % ndim)
        real(real64) :: v1(this % ndim), v2(this % ndim), v3(this % ndim), v4(this % ndim)

        ! vars init
        c           = 0.0d0
        nodalCoords = 0.0d0
        z           = 0.0d0
        t           = 0.0d0
        l1          = 0.0d0
        l2          = 0.0d0
        l3          = 0.0d0
        l4          = 0.0d0
        e1          = 0.0d0
        e2          = 0.0d0
        e3          = 0.0d0
        e4          = 0.0d0
        v1          = 0.0d0
        v2          = 0.0d0
        v3          = 0.0d0
        v4          = 0.0d0

        nodalCoords = reshape(                                              &
        (/                                                                  &
            -1.0d0,  1.0d0,  1.0d0, -1.0d0, -1.0d0,  1.0d0,  1.0d0, -1.0d0, &
            -1.0d0, -1.0d0,  1.0d0,  1.0d0, -1.0d0, -1.0d0,  1.0d0,  1.0d0, &
            -1.0d0, -1.0d0, -1.0d0, -1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0  &
        /), (/this % nnodes, this % ndim/))
        z = this % z 
        t = this % t

        v1 = nodalCoords(5,:) - nodalCoords(1,:) 
        v2 = nodalCoords(6,:) - nodalCoords(2,:) 
        v3 = nodalCoords(7,:) - nodalCoords(3,:) 
        v4 = nodalCoords(8,:) - nodalCoords(4,:) 
        l1 = vecNorm(v1)
        l2 = vecNorm(v2)
        l3 = vecNorm(v3)
        l4 = vecNorm(v4)
        e1 = v1 / l1
        e2 = v2 / l2
        e3 = v3 / l3
        e4 = v4 / l4

        c(1,:) = nodalCoords(1,:) + (z * l1) * e1
        c(2,:) = nodalCoords(2,:) + (z * l2) * e2
        c(3,:) = nodalCoords(3,:) + (z * l3) * e3
        c(4,:) = nodalCoords(4,:) + (z * l4) * e4
        c(5,:) = c(1,:) + (t * l1) * e1
        c(6,:) = c(2,:) + (t * l2) * e2
        c(7,:) = c(3,:) + (t * l3) * e3
        c(8,:) = c(4,:) + (t * l4) * e4

        call assign(this % x, c)
    
    end subroutine

    subroutine s8pslay_computeXg(this)
        implicit none
    
        class(S8PSLayer), intent(inout) :: this
    
        real(real64) :: c(this % nnodes, this % ndim)
        real(real64) :: nodalCoords(this % nnodes, this % ndim)
        real(real64) :: z, t, l1, l2, l3, l4
        real(real64) :: e1(this % ndim), e2(this % ndim), e3(this % ndim), e4(this % ndim)
        real(real64) :: v1(this % ndim), v2(this % ndim), v3(this % ndim), v4(this % ndim)

        ! vars init
        c           = 0.0d0
        nodalCoords = 0.0d0
        z           = 0.0d0
        t           = 0.0d0
        l1          = 0.0d0
        l2          = 0.0d0
        l3          = 0.0d0
        l4          = 0.0d0
        e1          = 0.0d0
        e2          = 0.0d0
        e3          = 0.0d0
        e4          = 0.0d0
        v1          = 0.0d0
        v2          = 0.0d0
        v3          = 0.0d0
        v4          = 0.0d0

        nodalCoords = this % ndx
        z = this % z 
        t = this % t

        v1 = nodalCoords(5,:) - nodalCoords(1,:) 
        v2 = nodalCoords(6,:) - nodalCoords(2,:) 
        v3 = nodalCoords(7,:) - nodalCoords(3,:) 
        v4 = nodalCoords(8,:) - nodalCoords(4,:) 
        l1 = vecNorm(v1)
        l2 = vecNorm(v2)
        l3 = vecNorm(v3)
        l4 = vecNorm(v4)
        e1 = v1 / l1
        e2 = v2 / l2
        e3 = v3 / l3
        e4 = v4 / l4

        c(1,:) = nodalCoords(1,:) + (z * l1) * e1
        c(2,:) = nodalCoords(2,:) + (z * l2) * e2
        c(3,:) = nodalCoords(3,:) + (z * l3) * e3
        c(4,:) = nodalCoords(4,:) + (z * l4) * e4
        c(5,:) = c(1,:) + (t * l1) * e1
        c(6,:) = c(2,:) + (t * l2) * e2
        c(7,:) = c(3,:) + (t * l3) * e3
        c(8,:) = c(4,:) + (t * l4) * e4

        call assign(this % xg, c)
    
    end subroutine
    
end module