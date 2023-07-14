!DEC$ FREEFORM

module layer_mod
    use iso_fortran_env
    use component_mod
    use integrator_mod
    use shapeFunctionLayered_mod
    use material_mod
    implicit none
    
    type, extends(Component), abstract :: Layer   

        integer                             :: id       = -1

        integer                             :: nnodes   = -1
        integer                             :: ndof     = -1
        integer                             :: ndim     = -1
        integer                             :: nstr     = -1
        integer                             :: npts     = -1
        integer                             :: nprops   = -1

        real(real64)                        :: z        = 0.0d0
        real(real64)                        :: t        = 0.0d0
        real(real64)                        :: theta    = 0.0d0

        real(real64),           allocatable :: x(:, :)
        real(real64),           allocatable :: xg(:, :)
        real(real64),           allocatable :: k(:, :)
        real(real64),           allocatable :: e(:, :)
        real(real64),           allocatable :: s(:, :)
        real(real64),           allocatable :: ndx(:,:)
        real(real64),           allocatable :: ipx(:,:)

        real(real64),           allocatable :: eadj(:,:)
        real(real64),           allocatable :: sadj(:,:)

        type(Transform3D)                   :: Tlm
        type(Transform3D)                   :: Tml
        type(Transform3D)                   :: Tlg
        type(Transform3D)                   :: Tgl

        class(Integrator),              allocatable :: integrator
        class(ShapeFunctionLayered),    allocatable :: shapeFunction
        class(Material),                allocatable :: material

    contains

        procedure(lay_computeX), deferred   :: computeX
        procedure(lay_computeXg), deferred   :: computeXg

        procedure :: print                  => lay_print

        procedure :: initialize             => lay_init
        procedure :: finalize               => lay_final

        procedure :: isInit                 => lay_isInit

        procedure :: updateLocalCS          => lay_updateLocalCS

        procedure :: computeStrainStress    => lay_computeStrainStress
        procedure :: computeStiffness       => lay_computeStiffness

        procedure :: computeStrainStressAdjoint    => lay_computeStrainStressAdjoint

        procedure :: globalCoordinates      => lay_globalCoordinates

    end type

    interface
        subroutine lay_computeX(this)
            import Layer 
            implicit none 
            class(Layer), intent(inout) :: this 
        end subroutine
        
        subroutine lay_computeXg(this)
            import Layer 
            implicit none 
            class(Layer), intent(inout) :: this 
        end subroutine
    end interface
   
contains

    subroutine lay_print(this, varName, unit)
        implicit none
        
        class(Layer),       intent(inout)   :: this
        character(len=*),   intent(in)      :: varName
        integer, optional,  intent(in)      :: unit
        
        integer :: uid

        ! vars init
        uid = 0
        
        uid = 6
        if (present(unit)) uid = unit 
        
        if (this % isInit()) then 
            write(uid, fmt100) ""
            write(uid, fmt100) varName // " < Layer | " // this % name // " >"
            write(uid, fmt100) ""
            write(uid, fmt101) "nnodes", this % nnodes
            write(uid, fmt101) "ndof", this % ndof
            write(uid, fmt101) "ndim", this % ndim
            write(uid, fmt101) "nstr", this % nstr 
            write(uid, fmt101) "npts", this % npts
            write(uid, fmt101) "nprops", this % nprops
            write(uid, fmt100) ""
            write(uid, fmt102) "z", this % z
            write(uid, fmt102) "t", this % t
            write(uid, fmt102) "theta", this % theta
            write(uid, fmt100) ""
            call printVar(this % x, "x", uid)
            call printVar(this % ndx, "ndx", uid)
            call printVar(this % ipx, "ipx", uid)
            call printVar(this % k, "k", uid)
            call printVar(this % e, "e", uid)
            call printVar(this % s, "s", uid)
            write(uid, fmt100) ""
            call this % Tlm % print(varName // ".Tlm", uid)
            call this % Tml % print(varName // ".Tml", uid)
            call this % Tlg % print(varName // ".Tlg", uid)
            call this % Tgl % print(varName // ".Tgl", uid)
            write(uid, fmt100) ""
            call this % integrator % print(varName // ".integrator", uid)
            call this % shapeFunction % print(varName // ".shapeFunction", uid)
            call this % material % print(varName // ".material", uid)
        else 
            call sroXIT(" >>>>>> lay_print <<<<<< Object not initialised")
        end if

        
    end subroutine

    subroutine lay_init(this, id, z, t, theta, x, Tlg, Tgl)
        implicit none
    
        class(Layer), intent(inout) :: this
    
        integer,                intent(in)         :: id
        real(real64),           intent(in)         :: z, t, theta
        real(real64),           intent(in)         :: x(:,:)
        class(Transform3D),     intent(in)         :: Tlg
        class(Transform3D),     intent(in)         :: Tgl

        if ( (this % nnodes /= -1) .and. &
             (this % ndof /= -1) .and. &
             (this % ndim /= -1) .and. &
             (this % nstr /= -1) .and. &
             (this % npts /= -1) .and. &
             (this % nprops /= -1) .and. &
             (id >= 0) .and. &
             (z >= 0) .and. &
             (t >= 0)) then

            this % id       = id
            this % z        = z 
            this % t        = t 
            this % theta    = theta
            
            this % Tlm = Transform3D(this % theta)
            this % Tml = Transform3D(this % Tlm % T3, Imat)
            this % Tlg = Tlg 
            this % Tgl = Tgl
            
            allocate(this % x(this % nnodes, this % ndim),  source = 0.0d0)
            allocate(this % xg(this % nnodes, this % ndim),  source = 0.0d0)
            allocate(this % k(this % ndof, this % ndof),    source = 0.0d0)
            allocate(this % e(this % npts, this % nstr),    source = 0.0d0)
            allocate(this % s(this % npts, this % nstr),    source = 0.0d0)
            allocate(this % ipx(this % npts, this % ndim),  source = 0.0d0)
            allocate(this % ndx(this % nnodes, this % ndim),source = x)

            allocate(this % eadj(this % npts, this % nstr), source = 0.0d0)
            allocate(this % sadj(this % npts, this % nstr), source = 0.0d0)

            call this % computeX
            call this % computeXg
        else 
            call sroXIT(" >>>>>> lay_init <<<<<< Bad input data ")
        end if
    end subroutine

    subroutine lay_final(this)
        implicit none
    
        class(Layer), intent(inout) :: this
    
        this % id = -1
        this % z = -1.0d0 
        this % t = -1.0d0
        this % theta = -1.0d0
        this % nnodes = -1
        this % ndof = -1
        this % ndim = -1
        this % nstr = -1
        this % npts = -1
        this % nprops = -1
        if (allocated(this % x)) deallocate(this % x)
        if (allocated(this % xg)) deallocate(this % xg)
        if (allocated(this % k)) deallocate(this % k)
        if (allocated(this % e)) deallocate(this % e)
        if (allocated(this % s)) deallocate(this % s)
        if (allocated(this % ipx)) deallocate(this % ipx)
        if (allocated(this % ndx)) deallocate(this % ndx)
        if (allocated(this % eadj)) deallocate(this % eadj)
        if (allocated(this % sadj)) deallocate(this % sadj)
    
    end subroutine

    function lay_isInit(this) result(b)
        implicit none
        logical :: b
        class(Layer) :: this

        ! vars init
        b = .false.

        if ( (this % id         /= -1) .and. &
             (this % z          >= 0.0d0) .and. &
             (this % t          >= 0.0d0) .and. &
             (this % nnodes     /= -1) .and. &
             (this % ndof       /= -1) .and. &
             (this % ndim       /= -1) .and. &
             (this % nstr       /= -1) .and. &
             (this % npts       /= -1) .and. &
             (this % nprops     /= -1) .and. &
             allocated(this % x) .and. &
             allocated(this % xg) .and. &
             allocated(this % k) .and. &
             allocated(this % e) .and. &
             allocated(this % s) .and. &
             allocated(this % ipx) .and. &
             allocated(this % ndx) .and. &
             allocated(this % sadj) .and. &
             allocated(this % eadj) .and. &
             this % shapeFunction % isInit() .and. &
             this % integrator % isInit() .and. &
             this % material % isInit()) then 
            
            b = .true. 
        else
            call sroLog(" >> Sro::Layer << | isInit | Attempted use before initialisation ")
            b = .false. 
        end if
        
    end function

    subroutine lay_updateLocalCS(this, Tlg, Tgl)
        implicit none
    
        class(Layer),       intent(inout)   :: this
        type(Transform3D),  intent(in)      :: Tlg, Tgl
    
        this % Tlg = Tlg 
        this % Tgl = Tgl
    
    end subroutine

    subroutine lay_computeStrainStress(this, u)
        implicit none
    
        class(Layer), intent(inout) :: this
        real(real64), intent(in)    :: u(:)

        integer         :: ipt
        real(real64)    :: dmat(this % nstr, this % nstr)
        real(real64)    :: eps(this % nstr), sig(this % nstr) 
        real(real64)    :: ut(this % ndof)

        ! vars init
        ipt     = 0
        dmat    = 0.0d0
        eps     = 0.0d0
        sig     = 0.0d0
        ut      = 0.0d0
    
        if (this % isInit()) then 

            this % e = 0.0d0
            this % s = 0.0d0
            
            dmat = this % Tlg % transformM6(    &
                    this % Tml % transformM6(   &
                        this % material % D     &
                    )                           &
                   ) 

            do ipt = 1, this % npts 
                eps = matmul(this % shapeFunction % B(:,:,ipt), u)
                sig = matmul(dmat, eps)
                this % e(ipt,:) = this % Tgl % TransformV6(eps)
                this % s(ipt,:) = this % Tgl % TransformV6(sig)
            end do

        else 
            call sroXIT(" >>>>>> lay_computeStrainStress <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine lay_computeStrainStressAdjoint(this, adj)
        implicit none
    
        class(Layer), intent(inout) :: this
        real(real64), intent(in)    :: adj(:)

        integer         :: ipt
        real(real64)    :: dmat(this % nstr, this % nstr)
        real(real64)    :: eps(this % nstr), sig(this % nstr) 
        real(real64)    :: adjt(this % ndof)

        ! vars init
        ipt     = 0
        dmat    = 0.0d0
        eps     = 0.0d0
        sig     = 0.0d0
        adjt      = 0.0d0
    
        if (this % isInit()) then 

            this % eadj = 0.0d0
            this % sadj = 0.0d0
            
            dmat = this % Tlg % transformM6(    &
                    this % Tml % transformM6(   &
                        this % material % D     &
                    )                           &
                   ) 

            do ipt = 1, this % npts 
                eps = matmul(this % shapeFunction % B(:,:,ipt), adj)
                sig = matmul(dmat, eps)
                this % eadj(ipt,:) = this % Tgl % TransformV6(eps)
                this % sadj(ipt,:) = this % Tgl % TransformV6(sig)
            end do

        else 
            call sroXIT(" >>>>>> lay_computeStrainStressAdjoint <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine lay_computeStiffness(this)
        implicit none
    
        class(Layer), intent(inout) :: this
        
        integer         :: ipt
        real(real64)    :: kmat(this % ndof, this % ndof)
        real(real64)    :: kmataux(this % ndof, this % ndof)
        real(real64)    :: aux 
        real(real64)    :: bmat(this % nstr, this % ndof, this % npts)
        real(real64)    :: bmatt(this % ndof, this % nstr)
        real(real64)    :: dmat(this % nstr, this % nstr)
        real(real64)    :: detJl(this % npts), detJe(this % npts)
        real(real64)    :: w(this % npts, this % ndim)

        ! vars init
        ipt     = 0
        kmat    = 0.0d0
        kmataux = 0.0d0
        aux     = 0.0d0 
        bmat    = 0.0d0
        bmatt   = 0.0d0
        dmat    = 0.0d0
        detJl   = 0.0d0
        detJe   = 0.0d0
        w       = 0.0d0
    
        if (this % isInit()) then 
            
            dmat = this % Tlg % transformM6(    &
                    this % Tml % transformM6(   &
                        this % material % D     &
                    )                           &
                   )

            detJl   = this % shapeFunction % detJl 
            detJe   = this % shapeFunction % detJ
            bmat    = this % shapeFunction % B 
            w       = this % integrator % wts

            do ipt = 1, this % npts 
                aux = detJe(ipt) * detJl(ipt) * vecProd(w(ipt,:))
                if (this % nnodes == 6) then
                    aux = aux * 0.5d0
                end if
                bmatt = transpose(bmat(:,:,ipt)) 
                kmataux = matmul(bmatt, matmul(dmat, bmat(:,:,ipt))) * aux

                kmat = kmat + matmul(bmatt, matmul(dmat, bmat(:,:,ipt))) * aux
            end do

            call assign(this % k, kmat)
        else 
            call sroXIT(" >>>>>> lay_computeStiffness <<<<<< Object not initialised ")
        end if
    
    end subroutine

    function lay_globalCoordinates(this) result(xg)
        implicit none 

        class(Layer) :: this
        real(real64) :: xg(this % nnodes, this % ndim)
        
        integer      :: ind
        real(real64) :: nc(this % nnodes, this % ndim)
        real(real64) :: cl(this % nnodes, this % ndim) 
        real(real64) :: nmat(1, this % nnodes)

        ! vars init
        xg      = 0.0d0
        ind     = 0
        nc      = 0.0d0
        cl      = 0.0d0
        nmat    = 0.0d0
    
        if (this % isInit()) then 
            nc = this % ndx 
            cl = this % x
            do ind = 1, this % nnodes 
                nmat = this % shapeFunction % N_at(cl(ind,1), cl(ind,2), cl(ind,3), this % nnodes)
                xg(ind,:) = this % shapeFunction % getPointX(nc, nmat)
            end do
        else 
            call sroXIT(" >>>>>> lay_globalCoordinates <<<<<< Object not initialised ") 
        end if
    end function
   
end module