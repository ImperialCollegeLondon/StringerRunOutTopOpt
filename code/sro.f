#ifdef ABQ_WIN86_64
#define SLASH "\"
#else
#define SLASH "/"
#endif

!DEC$ FREEFORM

#include "utilities/dof_mod.f"
#include "utilities/timer_mod.f"
#include "utilities/abaqusXIT_mod.f"
#include "utilities/stringUtils_mod.f"
#include "utilities/printUtils_mod.f"
#include "utilities/mathUtils_mod.f"
#include "utilities/abaqusUtils_mod.f"
#include "utilities/partitioning_mod.f"
#include "utilities/coordinateTransforms_mod.f"
#include "utilities/vtkLegacy_mod.f"
#include "utilities/allocatable_mod.f"

#include "base/component_mod.f"
#include "base/material_mod.f"
#include "base/integrator_mod.f"
#include "base/shapeFunction_mod.f"
#include "base/shapeFunctionLayered_mod.f"
#include "base/layer_mod.f"
#include "base/element_mod.f"
#include "base/layered3D_mod.f"

#include "materials/laminaPSMat_mod.f"

#include "integrators/lobattoHex_mod.f"
#include "integrators/lobattoWedge_mod.f"
#include "integrators/gaussQuad_mod.f"

#include "shapeFunctions/isoparametricHex_mod.f"
#include "shapeFunctions/isoparametricWedge_mod.f"
#include "shapeFunctions/isoparametricQuad_mod.f"

#include "layers/laminaHexPSLay_mod.f"
#include "layers/laminaWedgePSLay_mod.f"

#include "elements/levelsetQuad_mod.f"
#include "elements/solidShellHex_mod.f"
#include "elements/solidShellWedge_mod.f"
#include "elements/floatingNodeLevelSetHex_mod.f"
#include "elements/adjointHex_mod.f"

#include "optimizers/levelsetOptimization_mod.f"

#include "meshes/fnmLsAnalysis_mod.f"
#include "meshes/fnmLsMesh_mod.f"

module stringerRunoutOptimisation

    use fnmLsMesh_mod

    ! logical                         :: sroDebug     = .true.
    logical                         :: sroDebug     = .false.
    logical                         :: sroIsInit    = .false.
    integer                         :: sroUnit      = -1
    character(len=:), allocatable   :: sroJobname
    character(len=:), allocatable   :: sroOutdir
    character(len=21)               :: sroStages(7) =   &
        [                                               &
            character(len=21) ::                        &
            "analysisStart",                            &
            "incrementStart",                           &
            "incrementEnd",                             &
            "analysisEnd",                              &
            "restartAnalysisStart",                     &
            "stepStart",                                &
            "stepEnd"                                   &
        ]

    type(FnmLsMesh)                 :: sroMesh

contains

    pure function getStage(lopId) result(stage)
        implicit none 

        character(len=:), allocatable :: stage
        integer, intent(in) :: lopId

        stage = ""
        if (lopId >= 1 .and. lopId <= 7) then
            stage = trim(adjustl(sroStages(lopId)))
        end if
    end function

end module