!DEC$ FREEFORM

module fnmLsAnalysis_mod
    use iso_fortran_env
    use component_mod
    implicit none
    
    type, extends(Component) :: FnmLsSettings
        logical         :: useEnergyReleaseRate = .false.
        logical         :: useConstantVelocity = .false.
        logical         :: useForcedPartition = .false.

        real(real64)    :: constantVelocity = 0.0d0
        integer         :: partitionMode = 0
        integer         :: velocityIter = 0

    contains
        procedure :: print => fnmlsset_print

    end type

    interface FnmLsSettings
        procedure :: fnmlssett_ctor
    end interface
    
contains

    function fnmlssett_ctor() result(this)
        implicit none

        type(FnmLsSettings) :: this

        this % name = "FnmLsSetting"

        this % useEnergyReleaseRate = .false.
        this % useConstantVelocity = .false.
        this % useForcedPartition = .false.

        this % constantVelocity = 0.0d0
        this % partitionMode = 0
        this % velocityIter = 0
    end function

    subroutine fnmlsset_print(this, varName, unit)
        implicit none
        
        class(FnmLsSettings), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid
        
        uid = 6
        if (present(unit)) uid = unit 
        
        write(uid, fmt100) ""
        write(uid, fmt100) varName // " <FnmLsSettings>"
        write(uid, fmt100) ""
        write(uid, fmt100) "useEnergyReleaseRate" // logical2str(this % useEnergyReleaseRate)
        write(uid, fmt100) "useConstantVelocity" &
            // logical2str(this % useConstantVelocity) &
            // " - " // float2str(this % constantVelocity)
        write(uid, fmt100) "useForcedPartition" &
            // logical2str(this % useForcedPartition) &
            // " - " // int2str(this % partitionMode)
        
    end subroutine
   
end module