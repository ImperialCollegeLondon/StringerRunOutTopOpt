!DEC$ FREEFORM

module element_mod
    use iso_fortran_env
    use component_mod
    implicit none

    type, extends(Component) :: Element

        integer                     :: id       = -1

        integer                     :: nnodes   = -1
        integer                     :: ndof     = -1
        integer                     :: ndofnode = -1
        integer                     :: ndim     = -1
        integer                     :: nstr     = -1
        integer                     :: npts     = -1
        integer                     :: nprops   = -1

        real(real64), allocatable   :: x(:,:)
        real(real64), allocatable   :: k(:,:)
        real(real64), allocatable   :: f(:)
        real(real64), allocatable   :: u(:)
        real(real64), allocatable   :: props(:)
        integer,      allocatable   :: con(:)
        integer,      allocatable   :: dof(:)

    contains

        procedure :: print      => elm_print

        procedure :: initialize => elm_init
        procedure :: finalize   => elm_final

        procedure :: isInit     => elm_isInit

        procedure :: updateDof  => elm_updateDof

    end type
   
contains

    subroutine elm_print(this, varName, unit)
        implicit none
        
        class(Element), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid, ii
        character(3) :: str

        ! vars init
        uid = 0
        ii  = 0
        str = "   "
        
        uid = 6
        if (present(unit)) uid = unit 

        if (this % isInit()) then
            write(uid, fmt100) ""
            write(uid, fmt100) varName // " < Element | " // this % name // " >"
            write(uid, fmt100) ""
            write(uid, fmt101) "id", this % id
            write(uid, fmt100) ""
            write(uid, fmt101) "nnodes", this % nnodes
            write(uid, fmt101) "ndof", this % ndof
            write(uid, fmt101) "ndofnode", this % ndofnode
            write(uid, fmt101) "ndim", this % ndim
            write(uid, fmt101) "nstr", this % nstr
            write(uid, fmt101) "npts", this % npts
            write(uid, fmt101) "nprops", this % nprops
            write(uid, fmt100) ""
            call printVar(this % x, "x", uid)
            call printVar(this % con, "con", uid)
            call printVar(this % dof, "dof", uid)
            call printVar(this % props, "props", uid)
            call printVar(this % k, "k", uid)
            call printVar(this % f, "f", uid)
            call printVar(this % u, "u", uid)
        else
            call sroXIT(" >>>>>> elm_print <<<<<< Object not initialised")
        end if
    end subroutine

    subroutine elm_init(this, id, coords, con, props)
        implicit none
    
        class(Element), intent(inout)   :: this
        integer,        intent(in)      :: id
        real(real64),   intent(in)      :: coords(:,:)
        integer,        intent(in)      :: con(:)
        real(real64),   intent(in)      :: props(:)
    
        if ( (this % nnodes /= -1) .and. &
             (this % ndof /= -1) .and. &
             (this % ndofnode /= -1) .and. &
             (this % ndim /= -1) .and. &
             (this % nstr /= -1) .and. &
             (this % npts /= -1) .and. &
             (this % nprops /= -1) .and. &
             all(shape(coords) == (/this % nnodes, this % ndim/)) .and. &
             (size(con) == this % nnodes) .and. &
             (size(props) == this % nprops)) then

            this % id = id

            allocate(this % x(this % nnodes, this % ndim), source=coords)
            allocate(this % k(this % ndof, this % ndof), source=0.0d0)
            allocate(this % f(this % ndof), source=0.0d0)
            allocate(this % u(this % ndof), source=0.0d0)
            allocate(this % props(this % nprops), source=props)
            allocate(this % con(this % nnodes), source=con)
            allocate(this % dof(this % ndof), source=0)
            call this % updateDof
        else 
            call sroXIT(" >>>>>> elm_init <<<<<< Bad input data")
        end if
    end subroutine

    subroutine elm_final(this)
        implicit none
    
        class(Element), intent(inout) :: this
    
        this % id = -1
        this % nnodes = -1
        this % ndof = -1
        this % ndofnode = -1
        this % ndim = -1
        this % nstr = -1
        this % npts = -1
        this % nprops = -1
        if (allocated(this % x)) deallocate(this % x)
        if (allocated(this % k)) deallocate(this % k)
        if (allocated(this % f)) deallocate(this % f)
        if (allocated(this % u)) deallocate(this % u)
        if (allocated(this % props)) deallocate(this % props)
        if (allocated(this % con)) deallocate(this % con)
        if (allocated(this % dof)) deallocate(this % dof)
    end subroutine

    function elm_isInit(this) result(b)
        implicit none
        logical         :: b
        class(Element)  :: this

        ! vars init
        b = .false.
        
        if ( (this % id         /= -1) .and. &
             (this % nnodes     /= -1) .and. &
             (this % ndof       /= -1) .and. &
             (this % ndofnode   /= -1) .and. &
             (this % ndim       /= -1) .and. &
             (this % nstr       /= -1) .and. &
             (this % npts       /= -1) .and. &
             (this % nprops     /= -1) .and. &
             allocated(this % x) .and. &
             allocated(this % k) .and. &
             allocated(this % f) .and. &
             allocated(this % u) .and. &
             allocated(this % props) .and. &
             allocated(this % con) .and. &
             allocated(this % dof)) then 
        
            b = .true. 
        else 
            call sroLog(" >> Sro::Element << | isInit | Attempted use before initialisation ")
            b = .false.
        end if
    end function

    subroutine elm_updateDof(this)
        implicit none
    
        class(Element), intent(inout) :: this
        
        integer :: ii, jj, kk

        ! vars init
        ii = 0
        jj = 0 
        kk = 0
        
        kk = 1
        do ii = 1, this % nnodes
            do jj = 1, this % ndofnode
                this % dof(kk) = this % ndofnode * (this % con(ii) - 1) + jj
                kk = kk + 1
            end do
        end do
    end subroutine

end module