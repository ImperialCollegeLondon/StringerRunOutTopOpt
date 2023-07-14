#ifndef SRO_DEBUG
#define SRO_DEBUG
#endif
!DEC$ FREEFORM

module solidShellWedge_mod
    use iso_fortran_env
    use layered3D_mod
    use laminaWedgePSLay_mod
    implicit none

    type, extends(LayeredElm3D) :: S6PS
    contains

        procedure :: computeVolume  => s6ps_computeVolume
        procedure :: computeLocalCS => s6ps_computeLocalCS
        procedure :: setupLayers    => s6ps_setupLayers

        final :: s6ps_dtor
        
    end type
    
    interface S6PS
        procedure :: s6ps_ctor
    end interface
    
contains

    function s6ps_ctor(id, coords, con, props, orientations, thicknesses, localCS, Tlg, Tgl) result(this)
        implicit none 

        type(S6PS)                      :: this
        integer,            intent(in)  :: id
        integer,            intent(in)  :: con(:)
        real(real64),       intent(in)  :: coords(:,:)
        real(real64),       intent(in)  :: props(:)
        real(real64),       intent(in)  :: orientations(:)
        real(real64),       intent(in)  :: thicknesses(:)
        type(Transform3D),  intent(in)  :: Tlg, Tgl 
        real(real64),       intent(in)  :: localCS(:,:)

        this % name = "S6PS"

        if (size(orientations) == size(thicknesses)) then

            this % nnodes       = 6
            this % ndof         = 18
            this % ndofnode     = 3
            this % ndim         = 3
            this % nstr         = 6
            this % npts         = 5
            this % nprops       = size(props)

            this % nlayers      = size(orientations)

            call this % initialize(id, coords, con, props) 
            ! call this % computeLocalCS(transform % toCS)
            call assign(this % localCS, localCS)
            this % Tgl = Tgl
            this % Tlg = Tlg
            call this % setupLayers(orientations, thicknesses)
            call this % computeVolume 
        else 
            call sroXIT(" >>>>>> s6ps_ctor <<<<<< Bad input data ")
        end if
    end function

    subroutine s6ps_dtor(this)
        implicit none
    
        type(S6PS), intent(inout) :: this
    
        call this % finalize
    end subroutine

    subroutine s6ps_computeVolume(this)
        implicit none
    
        class(S6PS), intent(inout) :: this

        ! vars init
        this % volume = 0.0d0
    
        this % volume = volWedge(this % x)
    end subroutine

    subroutine s6ps_computeLocalCS(this, csys)
        implicit none
    
        class(S6PS), intent(inout) :: this
        real(real64), intent(in) :: csys(:,:)
    
        real(real64) :: r21(this % ndim), r31(this % ndim), r32(this % ndim)
        real(real64) :: s1(this % ndim), s3(this % ndim)
        real(real64) :: e1(this % ndim), e2(this % ndim), e3(this % ndim)
        real(real64) :: lcs(this % ndim, this % ndim)
        real(real64) :: c(this % nnodes, this % ndim)
        real(real64) :: xdir(this % ndim), xprojdir(this % ndim)

        ! vars init
        r21         = 0.0d0
        r31         = 0.0d0
        s1          = 0.0d0
        s3          = 0.0d0
        e1          = 0.0d0
        e2          = 0.0d0
        e3          = 0.0d0
        lcs         = 0.0d0
        c           = 0.0d0
        xdir        = 0.0d0
        xprojdir    = 0.0d0
        
        if (this % Element % isInit()) then 
            c(:,1) = this%x(:,1) + this%u((/1,4,7,10,13,16/))
            c(:,2) = this%x(:,2) + this%u((/2,5,8,11,14,17/))
            c(:,3) = this%x(:,3) + this%u((/3,6,9,12,15,18/))

            r31 = c(3, :) - c(1, :)
            r32 = c(3, :) - c(2, :)
            
            s3 = vec3cross(r31, r32)
            e3 = s3 / vecNorm(s3)

            xdir = csys(:,1)
            xprojdir = xdir - (dot_product(xdir, e3) * e3)
            e1 = xprojdir / vecNorm(xprojdir)
            
            ! r21 = c(2, :) - c(1, :)
            ! s1 = r21 - (dot_product(r21, e3) * e3)
            ! e1 = s1 / vecNorm(s1)
            
            e2 = vec3cross(e3, e1)
            
            lcs(:,:) = 0.0d0 
            lcs(:,1) = e1
            lcs(:,2) = e2
            lcs(:,3) = e3

            call assign(this % localCS, lcs)
            this % Tgl = Transform3D(this % localCS, Imat)
            this % Tlg = Transform3D(this % Tgl % T3, Imat)
        else 
            call sroXIT(" >>>>>> s6ps_computeLocalCS <<<<<< Object not initialised ")
        end if
    
    end subroutine

    subroutine s6ps_setupLayers(this, orientations, thicknesses)
        implicit none
    
        class(S6PS),    intent(inout)   :: this
        real(real64),   intent(in)      :: orientations(:)
        real(real64),   intent(in)      :: thicknesses(:)
        
        type(S6PSLayer) :: dummy(this % nlayers)
        integer         :: ilay
        real(real64)    :: z
        
        ! vars init
        z       = 0.0d0
        ilay    = 0

#ifdef SRO_DEBUG
        do ilay = 1, this % nlayers 
            dummy(ilay) = S6PSLayer(    &
                ilay,                   &
                z,                      &
                thicknesses(ilay),      &
                orientations(ilay),     &
                this % x,               &
                this % props,           &
                this % Tlg,             &
                this % Tgl              &
            )
            z = z + thicknesses(ilay)
        end do

        allocate(this % layers(this % nlayers), source = dummy)
#else
        allocate(this % layers(this % nlayers), source = dummy)
        do ilay = 1, this % nlayers 
            this % layers(ilay) = S6PSLayer(    &
                ilay,                   &
                z,                      &
                thicknesses(ilay),      &
                orientations(ilay),     &
                this % x,               &
                this % props,           &
                this % Tlg,             &
                this % Tgl              &
            )
            z = z + thicknesses(ilay)
        end do
#endif

    end subroutine
   
end module

#ifdef SRO_DEBUG
#undef SRO_DEBUG
#endif