!DEC$ FREEFORM

module component_mod
    use iso_fortran_env
    use dof_mod
    use abaqusXIT_mod
    use printUtils_mod
    use stringUtils_mod
    use mathUtils_mod
    use abaqusUtils_mod
    use partitioning_mod
    use coordinateTransforms_mod
    use allocatable_mod
    implicit none
   
    type, abstract :: Component
        character(len=:), allocatable :: name
    contains
        procedure(component_print), deferred :: print
    end type

    interface
        subroutine component_print(this, varName, unit)
            import Component
            implicit none
        
            class(Component), intent(inout) :: this
            character(len=*), intent(in) :: varName
            integer, optional, intent(in) :: unit
        end subroutine
    end interface

end module