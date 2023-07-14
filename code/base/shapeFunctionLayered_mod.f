!DEC$ FREEFORM

module shapeFunctionLayered_mod
    use iso_fortran_env
    use mathUtils_mod
    use component_mod
    implicit none
   
    type, extends(Component), abstract :: ShapeFunctionLayered

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

        real(real64), allocatable   :: Nl(:, :, :)
        real(real64), allocatable   :: dNl(:, :, :)

        real(real64), allocatable   :: Jl(:, :, :)
        real(real64), allocatable   :: Jinvl(:, :, :)
        real(real64), allocatable   :: detJl(:)

        real(real64), allocatable   :: JinvdN(:, :, :)
        real(real64), allocatable   :: B(:, :, :)

    contains

        procedure(sfl_N_at),  deferred :: N_at 
        procedure(sfl_dN_at), deferred :: dN_at

        procedure :: print      => sfl_print

        procedure :: computeN   => sfl_computeN
        procedure :: computeJ   => sfl_computeJ
        procedure :: computeB   => sfl_computeB

        procedure :: initialize => sfl_init
        procedure :: finalize   => sfl_final

        procedure :: isInit     => sfl_isInit

        procedure :: getPointX  => sfl_getPointX
        procedure :: computeIPCoordinates => sfl_computeIPCoordinates

    end type

    interface
        function sfl_N_at(this, xi, eta, zeta, n) result(nmat)
            use iso_fortran_env
            import ShapeFunctionLayered
            implicit none 

            class(ShapeFunctionLayered), intent(inout) :: this
            real(real64), intent(in) :: xi, eta, zeta 
            integer,      intent(in) :: n
            real(real64) :: nmat(1, n)
        end function

        function sfl_dN_at(this, xi, eta, zeta, n, m) result(dnmat)
            use iso_fortran_env
            import ShapeFunctionLayered
            implicit none 

            class(ShapeFunctionLayered), intent(inout) :: this
            real(real64), intent(in) :: xi, eta, zeta 
            integer, intent(in) :: n, m
            real(real64) :: dnmat(n, m)
        end function
    end interface
   
contains

    subroutine sfl_print(this, varName, unit)
        implicit none
        
        class(ShapeFunctionLayered), intent(inout) :: this
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
            write(uid, fmt100) varName // " < ShapeFunctionLayered | " // this % name // " >"
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
                call printVar(this % Nl(:,:,ipt), "N", uid)
                call printVar(this % dNl(:,:,ipt), "dN", uid)
                call printVar(this % Jl(:,:,ipt), "J", uid)
                call printVar(this % Jinvl(:,:,ipt), "Jinv", uid)
                write(uid, fmt102) "detJ", this % detJl(ipt)
                call printVar(this % JinvdN(:,:,ipt), "JinvdN", uid)
                call printVar(this % B(:,:,ipt), "B", uid)
            end do
        else 
            call sroXIT(" >>>>>> sfl_print <<<<<< Object not initialised ")
        end if
    end subroutine

    subroutine sfl_computeN(this, pts, nmat, dnmat)
        implicit none
    
        class(ShapeFunctionLayered),   intent(inout)   :: this
        real(real64),           intent(in)      :: pts(:,:)
        real(real64),           intent(inout)   :: nmat(:,:,:)
        real(real64),           intent(inout)   :: dnmat(:,:,:)
    
        integer         :: ipt, n
        real(real64)    :: xi, eta, zeta

        ! vars init
        ipt     = 0
        n       = 0
        xi      = 0.0d0
        eta     = 0.0d0
        zeta    = 0.0d0

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
    
    end subroutine

    subroutine sfl_computeJ(this, coords, dnmat, jmat, jinvmat, detjmat)
        implicit none
    
        class(ShapeFunctionLayered),   intent(inout)   :: this
        real(real64),           intent(in)      :: coords(:,:)
        real(real64),           intent(in)      :: dnmat(:,:,:)
        real(real64),           intent(inout)   :: jmat(:,:,:)
        real(real64),           intent(inout)   :: jinvmat(:,:,:)
        real(real64),           intent(inout)   :: detjmat(:)

        integer         :: ipt

        ! vars init
        ipt         = 0
            
        do ipt = 1, this % npts 
            jmat(:, :, ipt)         = matmul(dnmat(:, :, ipt), coords)
            if (this % ndim == 3) then
                detjmat(ipt)        = mat3det(jmat(:, :, ipt))
                jinvmat(:, :, ipt)  = mat3inv(jmat(:, :, ipt), detjmat(ipt))
            else if (this % ndim == 2) then 
                detjmat(ipt)        = mat2det(jmat(:, :, ipt))
                jinvmat(:, :, ipt)  = mat2inv(jmat(:, :, ipt))
            else 
                call sroXIT(" >>>>>> sfl_computeJ <<<<<< Invalid ShapeFunctionLayered dimension ")
            end if
        end do
    
    end subroutine

    subroutine sfl_computeB(this)
        implicit none
    
        class(ShapeFunctionLayered), intent(inout) :: this

        real(real64)    :: dnmat(this % ndim, this % nnodes)
        real(real64)    :: jinvmat(this % ndim, this % ndim)
        real(real64)    :: jinvmatl(this % ndim, this % ndim)
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
                dnmat = this % dNl(:, :, ipt)
                jinvmat = this % Jinv(:, :, ipt)
                jinvmatl = this % Jinvl(:, :, ipt)
                jinvdnmat = matmul(jinvmat, matmul(jinvmatl, dnmat))
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
                    call sroXIT(" >>>>>> sfl_computeB <<<<<< Invalid ShapeFunctionLayered dimension ")
                end if
            end do

            call assign(this % B, bmat)
        else 
            call sroXIT(" >>>>>> sfl_computeB <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine sfl_init(this, x, xl, ptsl)
        implicit none
    
        class(ShapeFunctionLayered),   intent(inout)   :: this
    
        real(real64),           intent(in)      :: x(:,:)
        real(real64),           intent(in)      :: xl(:,:)
        real(real64),           intent(in)      :: ptsl(:,:)

        real(real64) :: nmat(1, this % nnodes, this % npts)
        real(real64) :: dnmat(this % ndim, this % nnodes, this % npts)
        real(real64) :: jmat(this % ndim, this % ndim, this % npts)
        real(real64) :: jinvmat(this % ndim, this % ndim, this % npts)
        real(real64) :: detjmat(this % npts)
        real(real64) :: pts(this % npts, this % ndim)
            
        ! vars init 
        nmat = 0.0d0
        dnmat = 0.0d0
        jmat = 0.0d0
        jinvmat = 0.0d0
        detjmat = 0.0d0
        pts = 0.0d0

        if ( (this % nnodes /= -1) .and. &
             (this % ndof /= -1) .and. &
             (this % ndim /= -1) .and. &
             (this % nstr /= -1) .and. &
             (this % npts /= -1) .and. &
             all(shape(x) == (/this % nnodes, this % ndim/)) .and. &
             all(shape(xl) == (/this % nnodes, this % ndim/)) .and. &
             all(shape(ptsl) == (/this % npts, this % ndim/))) then 

            allocate(this % N(1, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % dN(this % ndim, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % J(this % ndim, this % ndim, this % npts), source = 0.0d0)
            allocate(this % Jinv(this % ndim, this % ndim, this % npts), source = 0.0d0)
            allocate(this % detJ(this % npts), source = 0.0d0)
            allocate(this % Nl(1, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % dNl(this % ndim, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % Jl(this % ndim, this % ndim, this % npts), source = 0.0d0)
            allocate(this % Jinvl(this % ndim, this % ndim, this % npts), source = 0.0d0)
            allocate(this % detJl(this % npts), source = 0.0d0)
            allocate(this % JinvdN(this % ndim, this % nnodes, this % npts), source = 0.0d0)
            allocate(this % B(this % nstr, this % ndof, this % npts), source = 0.0d0)

            call this % computeN(ptsl, nmat, dnmat)
            call this % computeJ(xl, dnmat, jmat, jinvmat, detjmat)
            call assign(this % Nl, nmat)
            call assign(this % dNl, dnmat)
            call assign(this % Jl, jmat)
            call assign(this % Jinvl, jinvmat)
            call assign(this % detJl, detjmat)
            
            pts = this % computeIPCoordinates(xl)
            nmat = 0.0d0
            dnmat = 0.0d0
            jmat = 0.0d0
            jinvmat = 0.0d0
            detjmat = 0.0d0
            call this % computeN(pts, nmat, dnmat)
            call this % computeJ(x, dnmat, jmat, jinvmat, detjmat)
            call assign(this % N, nmat)
            call assign(this % dN, dnmat)
            call assign(this % J, jmat)
            call assign(this % Jinv, jinvmat)
            call assign(this % detJ, detjmat)

            call this % computeB

        else 
            call sroXIT(" >>>>>> sfl_init <<<<<< Bad input data ")
        end if
    
    end subroutine

    subroutine sfl_final(this)
        implicit none
    
        class(ShapeFunctionLayered), intent(inout) :: this
    
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
        if (allocated(this % Nl)) deallocate(this % Nl)
        if (allocated(this % dNl)) deallocate(this % dNl)
        if (allocated(this % Jl)) deallocate(this % Jl)
        if (allocated(this % Jinvl)) deallocate(this % Jinvl)
        if (allocated(this % detJl)) deallocate(this % detJl)
        if (allocated(this % JinvdN)) deallocate(this % JinvdN)
        if (allocated(this % B)) deallocate(this % B)

    end subroutine

    function sfl_isInit(this) result(b)
        implicit none
        logical :: b
        class(ShapeFunctionLayered) :: this

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
             allocated(this % Nl) .and. &
             allocated(this % dNl) .and. &
             allocated(this % Jl) .and. &
             allocated(this % Jinvl) .and. &
             allocated(this % detJl) .and. &
             allocated(this % JinvdN) .and. &
             allocated(this % B) ) then 

            b = .true.
        else 
            call sroLog(" >> Sro::ShapeFunctionLayered << | isInit | Attempted use before initialisation ")
            b = .false.
        end if
        
    end function

    function sfl_getPointX(this, coords, N) result(c)
        implicit none
        
        class(ShapeFunctionLayered)        :: this
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

    function sfl_computeIPCoordinates(this, x) result(ipxi)
        implicit none
    
        class(ShapeFunctionLayered), intent(inout) :: this
        real(real64), intent(in) :: x(:,:)
        real(real64) :: ipxi(this % npts, this % ndim)

        real(real64) :: nmat(this % npts, this % nnodes)
        integer      :: ipt

        ! vars init
        ipxi    = 0.0d0
        nmat    = 0.0d0
        ipt     = 0

        do ipt = 1, this % npts 
            nmat(ipt, :) = this % Nl(1,:,ipt)
        end do 

        ipxi = matmul(nmat, x) 
    end function
   
end module