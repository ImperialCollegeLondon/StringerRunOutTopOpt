!DEC$ FREEFORM

module shapeFunction_mod
    use iso_fortran_env
    use mathUtils_mod
    use component_mod
    implicit none
   
    type, extends(Component), abstract :: ShapeFunction

        integer                     :: nnodes = -1
        integer                     :: ndof   = -1
        integer                     :: ndim   = -1
        integer                     :: nstr   = -1
        integer                     :: npts   = -1
        
        real(real64), allocatable   :: N(:, :, :)
        real(real64), allocatable   :: dN(:, :, :)

        real(real64), allocatable   :: J(:, :, :)
        real(real64), allocatable   :: Jinv(:, :, :)
        real(real64), allocatable   :: detJ(:)

        real(real64), allocatable   :: JinvdN(:, :, :)
        real(real64), allocatable   :: B(:, :, :)

    contains

        procedure(sf_N_at),  deferred :: N_at 
        procedure(sf_dN_at), deferred :: dN_at

        procedure :: print      => sf_print

        procedure :: computeN   => sf_computeN
        procedure :: computeJ   => sf_computeJ
        procedure :: computeB   => sf_computeB

        procedure :: initialize => sf_init
        procedure :: finalize   => sf_final

        procedure :: isInit     => sf_isInit

        procedure :: getPointX  => sf_getPointX

    end type

    interface
        function sf_N_at(this, xi, eta, zeta, n) result(nmat)
            use iso_fortran_env
            import ShapeFunction
            implicit none 

            class(ShapeFunction), intent(inout) :: this
            real(real64), intent(in) :: xi, eta, zeta 
            integer,      intent(in) :: n
            real(real64) :: nmat(1, n)
        end function

        function sf_dN_at(this, xi, eta, zeta, n, m) result(dnmat)
            use iso_fortran_env
            import ShapeFunction
            implicit none 

            class(ShapeFunction), intent(inout) :: this
            real(real64), intent(in) :: xi, eta, zeta 
            integer, intent(in) :: n, m
            real(real64) :: dnmat(n, m)
        end function
    end interface
   
contains

    subroutine sf_print(this, varName, unit)
        implicit none
        
        class(ShapeFunction), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid, ipt

        ! vars init
        uid = 0
        ipt = 0
        
        uid = 6
        if (present(unit)) uid = unit 

        if (this % isInit()) then
            write(uid, fmt100) ""
            write(uid, fmt100) varName // " < ShapeFunction | " // this % name // " >"
            write(uid, fmt100) ""
            write(uid, fmt101) "nnodes", this % nnodes
            write(uid, fmt101) "ndof", this % ndof
            write(uid, fmt101) "ndim", this % ndim
            write(uid, fmt101) "nstr", this % nstr
            write(uid, fmt101) "npts", this % npts
            write(uid, fmt100) ""
            do ipt = 1, this % npts
                write(uid, fmt101) "Integration point", ipt
                call printVar(this % N(:,:,ipt), "N", uid)
                call printVar(this % dN(:,:,ipt), "dN", uid)
                call printVar(this % J(:,:,ipt), "J", uid)
                call printVar(this % Jinv(:,:,ipt), "Jinv", uid)
                write(uid, fmt102) "detJ", this % detJ(ipt)
                call printVar(this % JinvdN(:,:,ipt), "JinvdN", uid)
                call printVar(this % B(:,:,ipt), "B", uid)
            end do
        else 
            call sroXIT(" >>>>>> sf_print <<<<<< Object not initialised ")
        end if
    end subroutine

    subroutine sf_computeN(this, pts)
        implicit none
    
        class(ShapeFunction),   intent(inout)   :: this
        real(real64),           intent(in)      :: pts(:,:)
    
        integer         :: ipt, n
        real(real64)    :: xi, eta, zeta
        real(real64)    :: nmat(1, this % nnodes, this % npts)
        real(real64)    :: dnmat(this % ndim, this % nnodes, this % npts)

        ! vars init
        ipt     = 0
        n       = 0
        xi      = 0.0d0
        eta     = 0.0d0
        zeta    = 0.0d0
        nmat    = 0.0d0
        dnmat   = 0.0d0

        if (this % isInit()) then 

            do ipt = 1, this % npts 
                select case(this % ndim)
                case(2)
                    xi = pts(ipt, 1)
                    eta = pts(ipt, 2)
                    zeta = 0.0d0
                case(3)
                    xi = pts(ipt, 1)
                    eta = pts(ipt, 2)
                    zeta = pts(ipt, 3)
                end select

                nmat(:,:,ipt) = this % N_at(xi, eta, zeta, this % nnodes)
                dnmat(:,:,ipt) = this % dN_at(xi, eta, zeta, this % ndim, this % nnodes)
            end do

            call assign(this % N, nmat)
            call assign(this % dN, dnmat)

        else 
            call sroXIT(" >>>>>> sf_computeN <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine sf_computeJ(this, coords)
        implicit none
    
        class(ShapeFunction),   intent(inout)   :: this
        real(real64),           intent(in)      :: coords(:,:)

        integer         :: ipt
        real(real64)    :: jmat(this % ndim, this % ndim, this % npts)
        real(real64)    :: jinvmat(this % ndim, this % ndim, this % npts)
        real(real64)    :: detjmat(this % npts)

        ! vars init
        ipt         = 0
        jmat        = 0.0d0
        jinvmat     = 0.0d0
        detjmat     = 0.0d0
    
        if (this % isInit()) then 
            
            do ipt = 1, this % npts 
                jmat(:, :, ipt)     = matmul(this % dN(:, :, ipt), coords)
                if (this % ndim == 3) then
                    detjmat(ipt)        = mat3det(jmat(:, :, ipt))
                    jinvmat(:, :, ipt)  = mat3inv(jmat(:, :, ipt), detjmat(ipt))
                else if (this % ndim == 2) then 
                    detjmat(ipt)        = mat2det(jmat(:, :, ipt))
                    jinvmat(:, :, ipt)  = mat2inv(jmat(:, :, ipt))
                else 
                    call sroXIT(" >>>>>> sf_computeJ <<<<<< Invalid shapeFunction dimension ")
                end if
            end do

            call assign(this % J, jmat)
            call assign(this % Jinv, jinvmat)
            call assign(this % detJ, detjmat)

        else 
            call sroXIT(" >>>>>> sf_computeJ <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine sf_computeB(this)
        implicit none
    
        class(ShapeFunction), intent(inout) :: this

        real(real64)    :: dnmat(this % ndim, this % nnodes)
        real(real64)    :: jinvmat(this % ndim, this % ndim)
        real(real64)    :: jinvdnmat(this % ndim, this % nnodes)
        real(real64)    :: bmat(this % nstr, this % ndof, this % npts)
        integer         :: ipt, ind, i1, i2

        ! vars init
        dnmat       = 0.0d0
        jinvmat     = 0.0d0
        jinvdnmat   = 0.0d0
        bmat        = 0.0d0
        ipt         = 0
        ind         = 0
        i1          = 0
        i2          = 0
    
        if (this % isInit()) then 

            do ipt = 1, this % npts
                dnmat = this % dN(:, :, ipt)
                jinvmat = this % Jinv(:, :, ipt)
                jinvdnmat = matmul(jinvmat, dnmat)
                this % JinvdN(:, :, ipt) = jinvdnmat
                if (this % ndim == 3) then
                    do ind = 1, this % nnodes
                        i1 = this % ndim * ind - (this % ndim - 1)
                        i2 = this % ndim * ind
                        bmat(1, i1:i2, ipt) = [jinvdnmat(1, ind),             0.0d0,             0.0d0]
                        bmat(2, i1:i2, ipt) = [            0.0d0, jinvdnmat(2, ind),             0.0d0]
                        bmat(3, i1:i2, ipt) = [            0.0d0,             0.0d0, jinvdnmat(3, ind)]
                        bmat(4, i1:i2, ipt) = [jinvdnmat(2, ind), jinvdnmat(1, ind),             0.0d0]
                        bmat(5, i1:i2, ipt) = [jinvdnmat(3, ind),             0.0d0, jinvdnmat(1, ind)]
                        bmat(6, i1:i2, ipt) = [            0.0d0, jinvdnmat(3, ind), jinvdnmat(2, ind)]
                    end do
                else if (this % ndim == 2) then 
                    do ind = 1, this % nnodes
                        i1 = this % ndim * ind - (this % ndim - 1)
                        i2 = this % ndim * ind
                        bmat(1, i1:i2, ipt) = [jinvdnmat(1, ind),             0.0d0]
                        bmat(2, i1:i2, ipt) = [            0.0d0, jinvdnmat(2, ind)]
                        bmat(3, i1:i2, ipt) = [jinvdnmat(2, ind), jinvdnmat(1, ind)]
                    end do
                else 
                    call sroXIT(" >>>>>> sf_computeB <<<<<< Invalid shapeFunction dimension ")
                end if
            end do

            call assign(this % B, bmat)
        else 
            call sroXIT(" >>>>>> sf_computeB <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine sf_init(this, x, pts, useB)
        implicit none
    
        class(ShapeFunction),   intent(inout)   :: this
    
        real(real64),           intent(in)      :: x(:,:)
        real(real64),           intent(in)      :: pts(:,:)
        logical, optional,      intent(in)      :: useB

        logical :: useBFlag

        useBFlag = .True.
        if (present(useB)) useBFlag = useB


        if ( (this % nnodes /= -1) .and. &
             (this % ndof /= -1) .and. &
             (this % ndim /= -1) .and. &
             (this % nstr /= -1) .and. &
             (this % npts /= -1) .and. &
             all(shape(x) == (/this % nnodes, this % ndim/)) .and. &
             all(shape(pts) == (/this % npts, this % ndim/))) then 

            allocate(this % N(1, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % dN(this % ndim, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % J(this % ndim, this % ndim, this % npts), source = 0.0d0)
            allocate(this % Jinv(this % ndim, this % ndim, this % npts), source = 0.0d0)
            allocate(this % detJ(this % npts), source = 0.0d0)
            allocate(this % JinvdN(this % ndim, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % B(this % nstr, this % ndof, this % npts), source = 0.0d0)

            call this % computeN(pts)
            call this % computeJ(x)
            if (useBFlag) call this % computeB

        else 
            call sroXIT(" >>>>>> sf_init <<<<<< Bad input data ")
        end if
    
    end subroutine

    subroutine sf_final(this)
        implicit none
    
        class(ShapeFunction), intent(inout) :: this
    
        this % nnodes = -1
        this % ndof = -1
        this % ndim = -1
        this % nstr = -1
        this % npts = -1

        if (allocated(this % N)) deallocate(this % N)
        if (allocated(this % dN)) deallocate(this % dN)
        if (allocated(this % J)) deallocate(this % J)
        if (allocated(this % Jinv)) deallocate(this % Jinv)
        if (allocated(this % detJ)) deallocate(this % detJ)
        if (allocated(this % JinvdN)) deallocate(this % JinvdN)
        if (allocated(this % B)) deallocate(this % B)

    end subroutine

    function sf_isInit(this) result(b)
        implicit none
        logical :: b
        class(ShapeFunction) :: this

        ! vars init
        b = .false.

        if ( (this % nnodes /= -1) .and. &
             (this % ndof /= -1) .and. &
             (this % ndim /= -1) .and. &
             (this % nstr /= -1) .and. &
             (this % npts /= -1) .and. &
             allocated(this % N) .and. &
             allocated(this % dN) .and. &
             allocated(this % J) .and. &
             allocated(this % Jinv) .and. &
             allocated(this % detJ) .and. &
             allocated(this % JinvdN) .and. &
             allocated(this % B) ) then 

            b = .true.
        else 
            call sroLog(" >> Sro::ShapeFunction << | isInit | Attempted use before initialisation ")
            b = .false.
        end if
        
    end function

    function sf_getPointX(this, coords, N) result(c)
        implicit none
        
        class(ShapeFunction)        :: this
        real(real64)                :: c(this % ndim)
        real(real64), intent(in)    :: coords(this % nnodes, this % ndim)
        real(real64), intent(in)    :: N(1, this % nnodes)

        integer :: i

        ! vars init
        c = 0.0d0
        i = 0

        do i = 1, this % ndim
            c(i) = dot_product(N(1,:), coords(:,i))
        end do
    end function
   
end module