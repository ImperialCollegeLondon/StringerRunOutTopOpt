!DEC$ FREEFORM 

module vtkLegacy_mod
    use iso_fortran_env
    use stringUtils_mod
    use printUtils_mod
    implicit none

    type :: VTKFile

        integer                         :: unit     = -1
        character(len=:), allocatable   :: name     
        character(len=:), allocatable   :: desc     
        logical                         :: isOpen   = .false.

    contains 

        procedure :: open                   => vtkfile_open
        procedure :: close                  => vtkfile_close

        procedure :: writeUnstructuredGrid  => vtkfile_writeug

        procedure :: openPointData          => vtkfile_openpointdata
        procedure :: openCellData           => vtkfile_opencelldata

        procedure :: writeVector          => vtkfile_writevector
        procedure :: writeScalar          => vtkfile_writescalar
        procedure :: writeScalarInt       => vtkfile_writescalarint

    end type

    interface VTKFile
        procedure :: vtkf_ctor
    end interface
    
contains

    function vtkf_ctor(filename, description) result(this)
        implicit none

        type(VTKFile)                   :: this
        character(len=*), intent(in)    :: filename
        character(len=*), intent(in)    :: description

        this % name = filename
        this % desc = description
        this % isOpen = .false.
    end function

    subroutine vtkfile_open(this)
        implicit none
    
        class(VTKFile), intent(inout) :: this
    
        open(newunit=this % unit, file=this % name)
        this % isOpen = .true.
    end subroutine

    subroutine vtkfile_close(this)
        implicit none

        class(VTKFile), intent(inout) :: this

        if (this % isOpen) then
            close(this % unit)
            this % isOpen = .false.
        end if
    end subroutine

    subroutine vtkfile_writeug(this, points, cells, offsets, types)
        implicit none
    
        class(VTKFile), intent(inout) :: this

        real(real64), intent(in)    :: points(:,:)
        integer     , intent(in)    :: cells(:)
        integer     , intent(in)    :: offsets(:)
        integer     , intent(in)    :: types(:)

        integer :: i, j, k
        integer :: npts, ncells, nlencells

        ! vars init
        i           = 0
        j           = 0
        k           = 0
        npts        = 0
        ncells      = 0
        nlencells   = 0

        if (this % isOpen) then
            npts = size(points, dim=1)
            nlencells = size(cells)
            ncells = size(offsets)
            
            write(this % unit, fmt100) "# vtk DataFile Version 2.0"
            write(this % unit, fmt100) this % desc
            write(this % unit, fmt100) "ASCII"
            write(this % unit, fmt100) "DATASET UNSTRUCTURED_GRID"

            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "POINTS " // int2str(npts) // " double"
            do i = 1, npts 
                write(this % unit, fmt100) floatArray2str(points(i, :))
            end do

            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "CELLS " // int2str(ncells) // " " // int2str(nlencells + ncells)
            j = 0
            k = 0
            do i = 1, ncells
                j = k + 1
                k = offsets(i)
                write(this % unit, fmt100) int2str(k-j+1) // ' ' // intArray2str(cells(j:k))
            end do

            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "CELL_TYPES " // int2str(ncells)
            do i = 1, ncells
                write(this % unit, fmt100) int2str(types(i))
            end do
        end if
    end subroutine

    subroutine vtkfile_openpointdata(this, npts)
        implicit none
    
        class(VTKFile), intent(inout) :: this
        
        integer, intent(in) :: npts

        if (this % isOpen) then 
            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "POINT_DATA " // int2str(npts)
        end if
    end subroutine

    subroutine vtkfile_opencelldata(this, ncells)
        implicit none
    
        class(VTKFile), intent(inout) :: this

        integer, intent(in) :: ncells
    
        if (this % isOpen) then 
            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "CELL_DATA " // int2str(ncells)
        end if
    end subroutine

    subroutine vtkfile_writevector(this, v, name)
        implicit none
    
        class(VTKFile), intent(inout)   :: this
        real(real64),   intent(in)      :: v(:,:)
        character(len=*), intent(in)    :: name

        integer :: i

        ! vars init
        i = 0
    
        if (this % isOpen) then 
            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "VECTORS " // name // " double"
            do i = 1, size(v, dim=1)
                write(this % unit, fmt100) floatArray2str(v(i,:))
            end do
        end if
    end subroutine

    subroutine vtkfile_writescalar(this, v, name)
        implicit none
    
        class(VTKFile), intent(inout) :: this
        real(real64),   intent(in)      :: v(:)
        character(len=*), intent(in)    :: name

        integer :: i

        ! vars init
        i = 0
    
        if (this % isOpen) then 
            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "SCALARS " // name // " double"
            write(this % unit, fmt100) "LOOKUP_TABLE default"
            do i = 1, size(v)
                write(this % unit, fmt100) float2str(v(i))
            end do
        end if
    end subroutine
    
    subroutine vtkfile_writescalarint(this, v, name)
        implicit none
    
        class(VTKFile), intent(inout) :: this
        integer,   intent(in)      :: v(:)
        character(len=*), intent(in)    :: name

        integer :: i

        ! vars init
        i = 0
    
        if (this % isOpen) then 
            write(this % unit, fmt100) ""
            write(this % unit, fmt100) "SCALARS " // name // " int"
            write(this % unit, fmt100) "LOOKUP_TABLE default"
            do i = 1, size(v)
                write(this % unit, fmt100) int2str(v(i))
            end do
        end if
    end subroutine

   
end module