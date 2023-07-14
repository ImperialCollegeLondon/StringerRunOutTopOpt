!DEC$ FREEFORM

module levelsetQuad_mod
    use iso_fortran_env
    use element_mod
    use isoparametricQuad_mod
    use gaussQuad_mod
    implicit none
    
    type, extends(Element) :: Q4LS

        real(real64)                :: h        = -1.0d0
        real(real64)                :: area     = -1.0d0

        real(real64), allocatable   :: fv(:)
        real(real64), allocatable   :: v(:)

        type(Transform3D)           :: T

        type(QuadIsoSF)             :: shapeFunction
        type(GaussIntegrator2DQuad) :: integrator 

    contains

        procedure :: print              => q4ls_print

        procedure :: computeArea        => q4ls_computeArea
        procedure :: computeStiffness   => q4ls_computeStiffness
        procedure :: computeForce       => q4ls_computeForce
        procedure :: computeExtensForce => q4ls_computeExtensForce

        procedure :: initialize         => q4ls_init 
        procedure :: finalize           => q4ls_final

        procedure :: isInit             => q4ls_isInit

        final     :: q4ls_dtor

    end type
    
    interface Q4LS
        procedure :: q4ls_ctor
    end interface
    
contains

    subroutine q4ls_print(this, varName, unit)
        implicit none
        
        class(Q4LS), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid

        ! vars init
        uid = 0
        
        uid = 6
        if (present(unit)) uid = unit 

        if (this % isInit()) then 
            call this % Element % print(varName, uid)
            write(uid, fmt100) ""
            write(uid, fmt102) "h", this % h
            write(uid, fmt102) "area", this % area
            write(uid, fmt100) ""
            call printVar(this % fv, "fv", uid)
            call printVar(this % v, "v", uid)
            write(uid, fmt100) ""
            call this % shapeFunction % print(varName // ".shapeFunction", uid)
            call this % integrator % print(varName // ".integrator", uid)
        else 
            call sroXIT(" >>>>>> q4ls_print <<<<<< Object not initialised ")
        end if

        
    end subroutine
   
    function q4ls_ctor(id, coords, con, props, transform, origin) result(this)
        implicit none
        type(Q4LS)                      :: this
        integer,            intent(in)  :: id, con(:)
        real(real64),       intent(in)  :: coords(:,:), props(:)
        type(Transform3D),  intent(in)  :: transform
        real(real64),       intent(in)  :: origin(:)

        real(real64) :: ct(4,3)

        ! vars init
        ct = 0.0d0
        
        this % name = "Q4LS"

        this % nnodes = 4
        this % ndofnode = 1
        this % ndof = 4
        this % ndim = 2
        this % nstr = 3
        this % npts = 4
        this % nprops = 1

        this % T = transform
        ct = coords
        ct(1,:) = ct(1,:) - origin
        ct(2,:) = ct(2,:) - origin
        ct(3,:) = ct(3,:) - origin
        ct(4,:) = ct(4,:) - origin
        ct = this % T % transformV3x4(ct)

        call this % initialize(id, ct(:,1:2), con, props)
        call this % computeStiffness
    end function

    subroutine q4ls_dtor(this)
        implicit none
    
        type(Q4LS), intent(inout) :: this
    
        call this % finalize
    end subroutine

    subroutine q4ls_init(this, id, coords, con, props)
        implicit none
    
        class(Q4LS),    intent(inout)   :: this
        integer,        intent(in)      :: id
        real(real64),   intent(in)      :: coords(:,:)
        integer,        intent(in)      :: con(:)
        real(real64),   intent(in)      :: props(:)
    
        call this % Element % initialize(id, coords, con, props)

        this % integrator = GaussIntegrator2DQuad(4)
        this % shapeFunction = QuadIsoSF(coords, this % integrator % pts, .False.)

        allocate(this % fv(this % ndof), source = 0.0d0)
        allocate(this % v(this % ndof), source = 0.0d0)
    end subroutine

    subroutine q4ls_final(this)
        implicit none
    
        class(Q4LS), intent(inout) :: this
    
        if (allocated(this % fv)) deallocate(this % fv)
        if (allocated(this % v)) deallocate(this % v)
        call this % shapeFunction % finalize
        call this % integrator % finalize 
        call this % Element % finalize
    
    end subroutine

    function q4ls_isInit(this) result(b)
        implicit none 

        logical :: b 
        class(Q4LS) :: this 

        ! vars init 
        b = .false.

        if ( allocated(this % fv) .and. &
             allocated(this % v) .and. &
             this % shapeFunction % isInit() .and. &
             this % integrator % isInit() .and. &
             this % Element % isInit()) then
            
            b = .true.
        else 
            call sroLog(" >> Sro::Q4LS << | isInit | Attempted use before initialisation ")
            b = .false.
        end if

    end function

    subroutine q4ls_computeArea(this)
        implicit none
    
        class(Q4LS), intent(inout) :: this

        real(real64) :: e1, e2, e3, e4

        ! vars init
        e1          = 0.0d0
        e2          = 0.0d0
        e3          = 0.0d0
        e4          = 0.0d0
        this % area = 0.0d0
    
        this % area = areaQuadrilateral(    &
            this % x(1,:),        &
            this % x(2,:),        &
            this % x(3,:),        &
            this % x(4,:)         &
        )

        e1 = vecNorm(this % x(1,:) - this % x(2,:))
        e2 = vecNorm(this % x(2,:) - this % x(3,:))
        e3 = vecNorm(this % x(3,:) - this % x(4,:))
        e4 = vecNorm(this % x(4,:) - this % x(1,:))

        this % h = min(min(min(e1, e2), e3), e4)
    end subroutine

    subroutine q4ls_computeStiffness(this)
        implicit none
    
        class(Q4LS), intent(inout) :: this
    
        integer :: igp 

        real(real64) :: kaux(this % ndof, this % ndof)
        real(real64) :: nvec(this % nnodes)
        real(real64) :: aux

        ! vars init
        igp         = 0
        kaux        = 0.0d0
        nvec        = 0.0d0
        aux         = 0.0d0
        this % k    = 0.0d0

        do igp = 1, this % npts
            aux = 0.0d0
            aux = this % shapeFunction % detJ(igp)
            aux = aux * vecProd(this % integrator % wts(igp,:))
            nvec = this % shapeFunction % N(1,:,igp)
            kaux = vec4mul(nvec, nvec) * aux
            call assign_and_increment(this % k, matLump(kaux, this % ndof))
        end do
    end subroutine

    subroutine q4ls_computeForce(this, dt)
        implicit none
    
        class(Q4LS),    intent(inout)   :: this
        real(real64),   intent(in)      :: dt
    
        integer         :: igp, i
        real(real64)    :: fvec(this % ndof)
        real(real64)    :: phi(this % ndof)
        real(real64)    :: vel(this % ndof)
        real(real64)    :: invjmat(this % ndim, this % ndim)
        real(real64)    :: nmat(this % nnodes)
        real(real64)    :: dnmat(this % ndim, this % nnodes)
        real(real64)    :: detjmat
        real(real64)    :: dnxy(this % ndim, this % nnodes)
        real(real64)    :: p, dp(this % ndim), normdp
        real(real64)    :: vn, normal(this % ndim), v(this % ndim)
        real(real64)    :: beta
        real(real64)    :: W(this % nnodes)
        real(real64)    :: r1(this % ndof), r2(this % ndof)

        ! vars init
        igp     = 0
        i       = 0
        fvec    = 0.0d0
        phi     = 0.0d0
        vel     = 0.0d0
        invjmat = 0.0d0
        nmat    = 0.0d0
        dnmat   = 0.0d0
        detjmat = 0.0d0
        dnxy    = 0.0d0
        p       = 0.0d0
        dp      = 0.0d0
        normdp  = 0.0d0
        vn      = 0.0d0
        normal  = 0.0d0
        v       = 0.0d0
        beta    = 0.0d0
        W       = 0.0d0
        r1      = 0.0d0
        r2      = 0.0d0
        
        phi = this % u
        vel = this % v

        do igp = 1, this % npts
            invjmat = this % shapeFunction % Jinv(:,:,igp)
            nmat    = this % shapeFunction % N(1,:,igp)
            dnmat   = this % shapeFunction % dN(:,:,igp)
            detjmat = this % shapeFunction % detJ(igp)

            dnxy    = matmul(invjmat, dnmat)
            p       = dot_product(nmat, phi)
            dp      = matmul(dnxy, phi)
            normdp  = vecNorm(dp)
            vn      = dot_product(nmat, vel)
            normal  = -dp / normdp
            v       = vn * normal 

            beta = max( 0.0d0,                                                      &
                1.0d0 / (2.0d0 * sqrt(dt**(-2) + vecNorm(matmul(invjmat, v))**2))   &
            )
            W = nmat + beta * matmul(v, dnxy)

            r1 = W * vn * normdp * detjmat * vecProd(this % integrator % wts(igp,:))
            r2 = nmat * p * detjmat * vecProd(this % integrator % wts(igp,:))
            
            do i = 1, 4
                ! if is NaN (not equal itself) set to 0
                if (r1(i) /= r1(i)) then
                    r1(i) = 0.0d0
                    !call sroLog(" --[WARN]--  >>>>>> q4ls_computeForce <<<<<< NaN detected")
                end if
                if (r2(i) /= r2(i)) then
                    r2(i) = 0.0d0 
                    !call sroLog(" --[WARN]--  >>>>>> q4ls_computeForce <<<<<< NaN detected")
                end if
            end do

            fvec = fvec + dt*r1 + r2
        end do 

        call assign(this % f, fvec)
    end subroutine

    subroutine q4ls_computeExtensForce(this, dt)
        implicit none
    
        class(Q4LS), intent(inout) :: this
        real(real64), intent(in) :: dt
    
        integer         :: igp, i
        real(real64)    :: h
        real(real64)    :: fvec(this % ndof)
        real(real64)    :: phi(this % ndof)
        real(real64)    :: vel(this % ndof)
        real(real64)    :: invjmat(this % ndim, this % ndim)
        real(real64)    :: nmat(this % nnodes)
        real(real64)    :: dnmat(this % ndim, this % nnodes)
        real(real64)    :: detjmat
        real(real64)    :: dnxy(this % ndim, this % nnodes)
        real(real64)    :: p, dp(this % ndim), normdp
        real(real64)    :: vn, normal(this % ndim), v(this % ndim)
        real(real64)    :: dvn(this % ndim)
        real(real64)    :: s
        real(real64)    :: w(this % ndim)
        real(real64)    :: r1(this % ndof), r2(this % ndof)

        ! vars init
        igp     = 0
        i       = 0
        h       = 0.0d0
        fvec    = 0.0d0
        phi     = 0.0d0
        vel     = 0.0d0
        invjmat = 0.0d0
        nmat    = 0.0d0
        dnmat   = 0.0d0
        detjmat = 0.0d0
        dnxy    = 0.0d0
        p       = 0.0d0
        dp      = 0.0d0
        normdp  = 0.0d0
        vn      = 0.0d0
        normal  = 0.0d0
        v       = 0.0d0
        dvn     = 0.0d0
        s       = 0.0d0
        w       = 0.0d0
        r1      = 0.0d0
        r2      = 0.0d0
        
        phi = this % u
        vel = this % v
        h   = this % h

        do igp = 1, this % npts
            invjmat = this % shapeFunction % Jinv(:,:,igp)
            nmat    = this % shapeFunction % N(1,:,igp)
            dnmat   = this % shapeFunction % dN(:,:,igp)
            detjmat = this % shapeFunction % detJ(igp)

            dnxy    = matmul(invjmat, dnmat)
            p       = dot_product(nmat, phi)
            dp      = matmul(dnxy, phi)
            normdp  = vecNorm(dp) ! can be zero if phi is equal in all nodes
            vn      = dot_product(nmat, vel)
            dvn     = matmul(dnxy, vel)
            normal  = -dp / normdp ! can be NaN if normdp is zero
            s       = p / sqrt(p**2 + h**2 * normdp**2)
            w       = s * normal
            v       = -w

            r1 = nmat * vn * detjmat * vecProd(this % integrator % wts(igp,:))
            r2 = nmat * dot_product(v, dvn) * detjmat * vecProd(this % integrator % wts(igp,:))

            do i = 1, 4
                ! if is NaN (not equal itself) set to 0
                if (r1(i) /= r1(i)) r1(i) = 0.0d0
                if (r2(i) /= r2(i)) r2(i) = 0.0d0 
            end do
            
            fvec = fvec + r1 - dt*r2
        end do 

        call assign(this % fv, fvec)
    end subroutine

end module