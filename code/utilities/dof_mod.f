!DEC$ FREEFORM

module dof_mod
    use iso_fortran_env
    implicit none
    
    type :: SroGlobalData
        real(real64), allocatable :: dof(:,:)
        real(real64), allocatable :: adjdof(:,:)
        
        integer, allocatable :: elmap(:)
        integer, allocatable :: adjelmap(:)
    contains
        procedure :: allocateDof => sroglobaldata_allocatedof
        procedure :: allocateAdjdof => sroglobaldata_allocateadjdof

        procedure :: allocateElmap => sroglobaldata_allocateelmap
        procedure :: allocateAdjelmap => sroglobaldata_allocateadjelmap
        
        final :: SroGlobalData_dtor
    end type

    interface SroGlobalData 
        procedure :: sroglobaldata_ctor
    end interface
    
    type(SroGlobalData) :: sro

contains

    function sroglobaldata_ctor() result(this)
        type(SroGlobalData) :: this

        call this % allocateDof(1,1)
        call this % allocateAdjdof(1,1)

        call this % allocateElmap(1)
        call this % allocateAdjelmap(1)
        
    end function

    subroutine sroGlobalData_dtor(this)
        type(SroGlobalData) :: this

        if (allocated(this % dof)) deallocate(this % dof)
        if (allocated(this % adjdof)) deallocate(this % adjdof)

        if (allocated(this % elmap)) deallocate(this % elmap)
        if (allocated(this % adjelmap)) deallocate(this % adjelmap)
    end subroutine

    subroutine sroglobaldata_allocatedof(this, n, m)
        class(SroGlobalData), intent(inout) :: this
        integer, intent(in) :: n, m

        if (allocated(this % dof)) deallocate(this % dof)
        allocate(this % dof(n, m), source = 0.0d0)
    end subroutine

    subroutine sroglobaldata_allocateadjdof(this, n, m)
        class(SroGlobalData), intent(inout) :: this
        integer, intent(in) :: n, m

        if (allocated(this % adjdof)) deallocate(this % adjdof)
        allocate(this % adjdof(n, m), source = 0.0d0)
    end subroutine

    subroutine sroglobaldata_allocateelmap(this, n)
        class(SroGlobalData), intent(inout) :: this
        integer, intent(in) :: n

        if (allocated(this % elmap)) deallocate(this % elmap)
        allocate(this % elmap(n), source = 0)
    end subroutine

    subroutine sroglobaldata_allocateadjelmap(this, n)
        class(SroGlobalData), intent(inout) :: this
        integer, intent(in) :: n

        if (allocated(this % adjelmap)) deallocate(this % adjelmap)
        allocate(this % adjelmap(n), source = 0)
    end subroutine
   
end module