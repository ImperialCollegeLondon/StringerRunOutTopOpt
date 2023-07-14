!DEC$ FREEFORM

module partitioning_mod
    use iso_fortran_env
    use mathUtils_mod
    use printUtils_mod
    implicit none 

    type :: PartitionCase
        integer                 :: id                   = -1
        integer                 :: nsubs                = -1
        integer                 :: nnodes               = -1
        integer, allocatable    :: connectivity(:,:)
        integer, allocatable    :: boundaryNds(:)

    contains
        procedure :: print => partition_print

        procedure :: finalize => partition_final
        final :: partition_dtor
    end type

    interface PartitionCase 
        procedure :: partition_ctor
    end interface 

contains 

    function partition_ctor(id, con, bnd) result(this)
        implicit none 

        type(PartitionCase) :: this 

        integer, intent(in) :: id
        integer :: con(:,:), bnd(:)
        integer :: s(2), sbnd

        ! vars init
        s       = 0
        sbnd    = 0

        this % id = id
        s = shape(con)
        this % nsubs = s(1) 
        this % nnodes = s(2) 
        allocate(this % connectivity(this % nsubs, this % nnodes), source=con)
        sbnd = size(bnd)
        allocate(this % boundaryNds(sbnd), source=bnd)
    end function

    subroutine partition_dtor(this)
        implicit none 

        type(PartitionCase), intent(inout) :: this

        call this % finalize
    end subroutine

    subroutine partition_final(this)
        implicit none 

        class(PartitionCase), intent(inout) :: this 

        this % id = -1
        this % nsubs = 0
        this % nnodes = 0
        if (allocated(this % connectivity)) deallocate(this % connectivity)
        if (allocated(this % boundaryNds)) deallocate(this % boundaryNds)
    end subroutine

    subroutine partition_print(this, varName, unit)
        implicit none
        
        class(PartitionCase), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid

        ! vars init
        uid = 0
        
        uid = 6
        if (present(unit)) uid = unit 
        
        write(uid, fmt100) ""
        write(uid, fmt100) varName // " <PartitionCase>"
        write(uid, fmt100) ""

        write(uid, fmt100) ""
        write(uid, fmt101) "id", this % id
        write(uid, fmt101) "nsubs", this % nsubs
        write(uid, fmt101) "nnodes", this % nnodes
        call printVar(this % connectivity, "connectivity", uid)
        call printVar(this % boundaryNds, "boundaryNds", uid)
    end subroutine
end module 