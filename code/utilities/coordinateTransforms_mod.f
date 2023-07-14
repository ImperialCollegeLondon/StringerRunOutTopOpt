!DEC$ FREEFORM

module coordinateTransforms_mod
    use iso_fortran_env
    use mathUtils_mod
    use printUtils_mod
    implicit none 

    real(real64), parameter :: toler = 1.0d-19

    type :: Transform3D
        real(real64) :: toCS(3,3)   = 0.0d0
        real(real64) :: fromCS(3,3) = 0.0d0
        real(real64) :: T3(3,3)     = 0.0d0
        real(real64) :: T6(6,6)     = 0.0d0
    contains
        procedure :: computeT3
        procedure :: computeT6
        procedure :: transformM3
        procedure :: transformM6
        procedure :: transformV3
        procedure :: transformV6
        procedure :: transformV3x3
        procedure :: transformV3x8
        procedure :: transformV3x4
        procedure :: transformV24

        procedure :: print => transform3D_print
    end type

    interface Transform3D
        procedure :: transform3D_init_default
        procedure :: transform3D_init
        procedure :: transform3D_init_rot
    end interface

contains

    function transform3D_init_default() result(instance)
        implicit none 

        type(Transform3D) :: instance

        instance % toCS(:,:) = 0.0d0
        instance % fromCS(:,:) = 0.0d0
        instance % T3(:,:) = 0.0d0 
        instance % T6(:,:) = 0.0d0
    end function

    function transform3D_init(toCS, fromCS) result(instance)
        implicit none 

        type(Transform3D) :: instance 
        real(real64), intent(in) :: toCS(3,3)
        real(real64), intent(in) :: fromCS(3,3)
        integer :: i, j

        ! vars init
        i = 0
        j = 0

        instance % toCS = toCS
        instance % fromCS = fromCS
        do i = 1, 3
            do j = 1, 3
                if (abs(instance%toCS(i,j)) .lt. toler) then
                    instance%toCS(i,j) = 0.0d0
                end if
                if (abs(instance%fromCS(i,j)) .lt. toler) then
                    instance%fromCS(i,j) = 0.0d0
                end if
            end do
        end do
        call instance % computeT3 
        call instance % computeT6
    end function

    function transform3D_init_rot(theta_deg) result(instance)
        implicit none 

        type(Transform3D) :: instance 
        real(real64), intent(in) :: theta_deg
        real(real64) :: I(3,3), CS(3,3)
        real(real64) :: t

        ! vars init
        I = 0.0d0
        CS = 0.0d0
        t = 0.0d0

        t = theta_deg * pi / 180.0d0

        I = reshape((/1.0d0,0.0d0,0.0d0, &
                      0.0d0,1.0d0,0.0d0, &
                      0.0d0,0.0d0,1.0d0 /), shape(I))

        CS = reshape ((/ cos(t), -sin(t),0.0d0, &
                         sin(t),  cos(t),0.0d0, &
                         0.0d0,  0.0d0,1.0d0 /), shape(CS))

        instance % toCS = CS 
        instance % fromCS = I 
        call instance % computeT3 
        call instance % computeT6
    end function

    subroutine computeT3(this)
        implicit none 

        class(Transform3D), intent(inout) :: this
        real(real64) :: T(3,3), tcs(3,3), fcs(3,3)
        integer :: i, j

        ! vars init
        T   = 0.0d0 
        tcs = 0.0d0
        fcs = 0.0d0
        i   = 0
        j   = 0

        tcs = this % toCS
        fcs = this % fromCS

        T(1,1) = dot_product(tcs(:,1), fcs(:,1))
        T(1,2) = dot_product(tcs(:,1), fcs(:,2))
        T(1,3) = dot_product(tcs(:,1), fcs(:,3))
        T(2,1) = dot_product(tcs(:,2), fcs(:,1))
        T(2,2) = dot_product(tcs(:,2), fcs(:,2))
        T(2,3) = dot_product(tcs(:,2), fcs(:,3))
        T(3,1) = dot_product(tcs(:,3), fcs(:,1))
        T(3,2) = dot_product(tcs(:,3), fcs(:,2))
        T(3,3) = dot_product(tcs(:,3), fcs(:,3))

        do i = 1, 3
            do j = 1, 3
                if (abs(T(i,j)) .lt. toler) then
                    T(i,j) = 0.0d0
                end if
            end do
        end do
        
        this % T3 = T
    end subroutine

    subroutine computeT6(this)
        implicit none 

        class(Transform3D), intent(inout) :: this 
        real(real64) :: R(6,6), T(3,3)
        integer :: i, j

        ! vars init
        R = 0.0d0
        T = 0.0d0
        i = 0
        j = 0

        R(:,:) = 0.0d0 
        T = this%T3

        R(1,1) = T(1,1)**2
        R(1,2) = T(2,1)**2
        R(1,3) = T(3,1)**2
        R(1,4) = T(1,1)*T(2,1)
        R(1,5) = T(1,1)*T(3,1)
        R(1,6) = T(2,1)*T(3,1)

        R(2,1) = T(1,2)**2
        R(2,2) = T(2,2)**2
        R(2,3) = T(3,2)**2
        R(2,4) = T(1,2)*T(2,2)
        R(2,5) = T(1,2)*T(3,2)
        R(2,6) = T(2,2)*T(3,2)

        R(3,1) = T(1,3)**2
        R(3,2) = T(2,3)**2
        R(3,3) = T(3,3)**2
        R(3,4) = T(1,3)*T(2,3)
        R(3,5) = T(1,3)*T(3,3)
        R(3,6) = T(2,3)*T(3,3)

        R(4,1) = 2*T(1,1)*T(1,2)
        R(4,2) = 2*T(2,1)*T(2,2)
        R(4,3) = 2*T(3,1)*T(3,2)
        R(4,4) = T(1,2)*T(2,1)+T(1,1)*T(2,2)
        R(4,5) = T(1,2)*T(3,1)+T(1,1)*T(3,2)
        R(4,6) = T(2,2)*T(3,1)+T(2,1)*T(3,2)

        R(5,1) = 2*T(1,1)*T(1,3)
        R(5,2) = 2*T(2,1)*T(2,3)
        R(5,3) = 2*T(3,1)*T(3,3)
        R(5,4) = T(2,1)*T(1,3)+T(1,1)*T(2,3)
        R(5,5) = T(3,1)*T(1,3)+T(1,1)*T(3,3)
        R(5,6) = T(3,1)*T(2,3)+T(2,1)*T(3,3)

        R(6,1) = 2*T(1,2)*T(1,3)
        R(6,2) = 2*T(2,2)*T(2,3)
        R(6,3) = 2*T(3,2)*T(3,3)
        R(6,4) = T(2,2)*T(1,3)+T(1,2)*T(2,3)
        R(6,5) = T(3,2)*T(1,3)+T(1,2)*T(3,3)
        R(6,6) = T(3,2)*T(2,3)+T(2,2)*T(3,3)

        do i = 1, 6
            do j = 1, 6
                if (abs(R(i,j)) .lt. toler) then
                    R(i,j) = 0.0d0
                end if
            end do
        end do

        this % T6 = transpose(R)
    end subroutine

    function transformM3(this, M3) result(M)
        implicit none 
        real(real64) :: M(3,3)
        class(Transform3D), intent(in) :: this
        real(real64), intent(in) :: M3(3,3)

        ! vars init
        M = 0.0d0

        M = matmul(matmul(this%T3, M3), transpose(this%T3))
    end function

    function transformM6(this, M6) result(M)
        implicit none 
        real(real64) :: M(6,6)
        class(Transform3D), intent(in) :: this
        real(real64), intent(in) :: M6(6,6)
        
        ! vars init
        M = 0.0d0

        M = matmul(matmul(this%T6, M6), transpose(this%T6))
    end function

    function transformV3(this, V3) result(V)
        implicit none 
        real(real64) :: V(3)
        class(Transform3D), intent(in) :: this
        real(real64), intent(in) :: V3(3)
        
        ! vars init
        V = 0.0d0

        V = matmul(this%T3, V3)
    end function

    function transformV6(this, V6) result(V)
        implicit none 
        real(real64) :: V(6)
        class(Transform3D), intent(in) :: this
        real(real64), intent(in) :: V6(6)
        
        ! vars init
        V = 0.0d0

        V = matmul(this%T6, V6)
    end function

    function transformV3x3(this, M3) result(M)
        implicit none 
        real(real64) :: M(3,3)
        class(Transform3D), intent(in) :: this
        real(real64), intent(in) :: M3(3,3)

        ! vars init
        M = 0.0d0
        
        M = matmul(this%T3, M3)
    end function
    
    function transformV3x8(this, M8x3) result(M)
        implicit none 
        real(real64) :: M(8,3)
        class(Transform3D), intent(in) :: this 
        real(real64), intent(in) :: M8x3(8,3)
        real(real64) :: aux1(3,8), aux2(3,8)

        ! vars init
        aux1 = 0.0d0
        aux2 = 0.0d0
        M    = 0.0d0
        
        aux1 = transpose(M8x3)
        aux2 = matmul(this%T3, aux1)
        M = transpose(aux2)
    end function

    function transformV3x4(this, M4x3) result(M)
        implicit none 
        real(real64) :: M(4,3)
        class(Transform3D), intent(in) :: this 
        real(real64), intent(in) :: M4x3(4,3)
        real(real64) :: aux1(3,4), aux2(3,4)

        ! vars init
        M    = 0.0d0
        aux1 = 0.0d0
        aux2 = 0.0d0
        
        aux1 = transpose(M4x3)
        aux2 = matmul(this%T3, aux1)
        M = transpose(aux2)
    end function

    function transformV24(this, V24) result(V)
        implicit none 
        real(real64) :: V(24)
        class(Transform3D), intent(in) :: this 
        real(real64), intent(in) :: V24(24)

        real(real64) :: M38(3, 8), M38transf(3,8)

        ! vars init
        M38         = 0.0d0
        M38transf   = 0.0d0
        V           = 0.0d0

        M38 = reshape(V24, shape(M38))
        M38transf = matmul(this % T3, M38)
        V = reshape(M38, shape(V))
    end function

    subroutine transform3D_print(this, varName, unit)
        implicit none
        class(Transform3D), intent(inout) :: this
        character(len=*),      intent(in) :: varName
        integer, optional,     intent(in) :: unit

        integer :: uid, ii

        ! vars init
        uid = 0
        ii = 0
        
        uid = 6 
        if (present(unit)) uid = unit
 
        write(uid, fmt100) ""
        write(uid, fmt100) sep1
        write(uid, fmt100) trim(varName)//"<Transform3D>"
        
        write(uid, fmt100) ""
        write(uid, fmt100) "toCS"
        do ii = 1, 3
            write(uid, fmt203) this % toCS(ii,:)
        end do 
        write(uid, fmt100) "fromCS"
        do ii = 1, 3
            write(uid, fmt203) this % fromCS(ii,:)
        end do 
        
        write(uid, fmt100) ""
        write(uid, fmt100) "T3"
        do ii = 1, 3
            write(uid, fmt203) this % T3(ii,:)
        end do 
        write(uid, fmt100) "T6"
        do ii = 1, 6
            write(uid, fmt206) this % T6(ii,:)
        end do
    end subroutine

end module