!DEC$ FREEFORM

module allocatable_mod
    use iso_fortran_env
    implicit none
    
    interface assign
        procedure :: assign_r1, assign_i1
        procedure :: assign_r2, assign_i2
        procedure :: assign_r3, assign_i3
        procedure :: assign_r4, assign_i4
    end interface

    interface assign_and_increment
        procedure :: assignincr_r1, assignincr_i1
        procedure :: assignincr_r2, assignincr_i2
        procedure :: assignincr_r3, assignincr_i3
        procedure :: assignincr_r4, assignincr_i4
    end interface
    
contains

    subroutine assign_r1(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:)
        real(real64), intent(in)    :: src(:)

        integer :: i

        ! vals init
        i = 0

        do i = 1, size(dst)
            dst(i) = src(i)
        end do
    end subroutine

    subroutine assign_i1(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:)
        integer, intent(in)    :: src(:)

        integer :: i

        ! vals init
        i = 0

        do i = 1, size(dst)
            dst(i) = src(i)
        end do
    end subroutine

    subroutine assign_r2(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:,:)
        real(real64), intent(in)    :: src(:,:)

        integer :: i, j

        ! vals init
        i = 0
        j = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                dst(i, j) = src(i, j)
            end do
        end do
    end subroutine

    subroutine assign_i2(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:,:)
        integer, intent(in)    :: src(:,:)

        integer :: i, j

        ! vals init
        i = 0 
        j = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                dst(i, j) = src(i, j)
            end do
        end do
    end subroutine

    subroutine assign_r3(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:,:,:)
        real(real64), intent(in)    :: src(:,:,:)

        integer :: i, j, k

        ! vals init
        i = 0 
        j = 0
        k = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    dst(i,j,k) = src(i,j,k)
                end do
            end do
        end do
    end subroutine

    subroutine assign_i3(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:,:,:)
        integer, intent(in)    :: src(:,:,:)

        integer :: i, j, k

        ! vals init
        i = 0 
        j = 0
        k = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    dst(i,j,k) = src(i,j,k)
                end do
            end do
        end do
    end subroutine

    subroutine assign_r4(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:,:,:,:)
        real(real64), intent(in)    :: src(:,:,:,:)

        integer :: i, j, k, l

        ! vals init
        i = 0 
        j = 0
        k = 0
        l = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    do l = 1, size(dst, dim=4)
                        dst(i,j,k,l) = src(i,j,k,l)
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine assign_i4(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:,:,:,:)
        integer, intent(in)    :: src(:,:,:,:)

        integer :: i, j, k, l

        ! vals init
        i = 0 
        j = 0
        k = 0
        l = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    do l = 1, size(dst, dim=4)
                        dst(i,j,k,l) = src(i,j,k,l)
                    end do
                end do
            end do
        end do
    end subroutine




    subroutine assignincr_r1(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:)
        real(real64), intent(in)    :: src(:)

        integer :: i

        ! vals init
        i = 0

        do i = 1, size(dst)
            dst(i) = dst(i) + src(i)
        end do
    end subroutine

    subroutine assignincr_i1(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:)
        integer, intent(in)    :: src(:)

        integer :: i

        ! vals init
        i = 0

        do i = 1, size(dst)
            dst(i) = dst(i) + src(i)
        end do
    end subroutine

    subroutine assignincr_r2(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:,:)
        real(real64), intent(in)    :: src(:,:)

        integer :: i, j

        ! vals init
        i = 0 
        j = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                dst(i, j) = dst(i, j) + src(i, j)
            end do
        end do
    end subroutine

    subroutine assignincr_i2(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:,:)
        integer, intent(in)    :: src(:,:)

        integer :: i, j

        ! vals init
        i = 0 
        j = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                dst(i, j) = dst(i, j) + src(i, j)
            end do
        end do
    end subroutine

    subroutine assignincr_r3(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:,:,:)
        real(real64), intent(in)    :: src(:,:,:)

        integer :: i, j, k

        ! vals init
        i = 0 
        j = 0
        k = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    dst(i,j,k) = dst(i,j,k) + src(i,j,k)
                end do
            end do
        end do
    end subroutine

    subroutine assignincr_i3(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:,:,:)
        integer, intent(in)    :: src(:,:,:)

        integer :: i, j, k

        ! vals init
        i = 0 
        j = 0
        k = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    dst(i,j,k) = dst(i,j,k) + src(i,j,k)
                end do
            end do
        end do
    end subroutine

    subroutine assignincr_r4(dst, src)
        implicit none
    
        real(real64), intent(inout) :: dst(:,:,:,:)
        real(real64), intent(in)    :: src(:,:,:,:)

        integer :: i, j, k, l

        ! vals init
        i = 0 
        j = 0
        k = 0
        l = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    do l = 1, size(dst, dim=4)
                        dst(i,j,k,l) = dst(i,j,k,l) + src(i,j,k,l)
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine assignincr_i4(dst, src)
        implicit none
    
        integer, intent(inout) :: dst(:,:,:,:)
        integer, intent(in)    :: src(:,:,:,:)

        integer :: i, j, k, l

        ! vals init
        i = 0 
        j = 0
        k = 0
        l = 0

        do i = 1, size(dst, dim=1)
            do j = 1, size(dst, dim=2)
                do k = 1, size(dst, dim=3)
                    do l = 1, size(dst, dim=4)
                        dst(i,j,k,l) = dst(i,j,k,l) + src(i,j,k,l)
                    end do
                end do
            end do
        end do
    end subroutine
   
end module