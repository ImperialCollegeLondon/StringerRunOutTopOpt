#ifndef SRO_VEL_SCALING
#define SRO_VEL_SCALING 5
#endif

#define SRO_CONSTRAINT_STAB

!DEC$ FREEFORM
#include "sro.f"

module wing_nd_mod
    use iso_fortran_env
    implicit none

    integer :: matrixUnit = -1
    integer :: currentStep = -1
    integer :: ribOffsetId = -1054

contains
end module

subroutine UEL( RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,                     &
    &           PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,       &
    &           KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,     &
    &           LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    
    use iso_fortran_env
    use stringerRunoutOptimisation
    use wing_nd_mod
    include 'aba_param.inc'
    
    real(real64) :: RHS(MLVARX, 1), AMATRX(NDOFEL, NDOFEL), PROPS(NPROPS),      &       
        &           DU(MLVARX, *), ENERGY(8), COORDS(MCRD, NNODE), U(NDOFEL),   &      
        &           SVARS(*), V(NDOFEL), A(NDOFEL), PREDEF(2, NPREDF, NNODE),   &
        &           ADLMAG(MDLOAD, *), DDLMAG(MDLOAD, *), PARAMS(*), TIME(2)
    real(real64) :: PNEWDT, DTIME, PERIOD
    integer      :: JPROPS(*), JDLTYP(MDLOAD, *), LFLAGS(*)
    integer      :: MLVARX, NDOFEL, MCRD, NNODE, MDLOAD, NPREDF, NRHS, NSVARS,  &
        &           NPROPS, JTYPE, KSTEP, KINC, JELEM, NDLOAD, NJPROP
    
    integer :: i, j

    if (LFLAGS(3) == 1) then

        call sroMesh % updateElement(JELEM + ribOffsetId, JTYPE, KSTEP, KINC, U, AMATRX, RHS(:,1), ENERGY)
        
    end if

end subroutine UEL
    

subroutine UEXTERNALDB( LOP, LRESTART, TIME, DTIME, KSTEP, KINC )

    use iso_fortran_env
    use stringerRunoutOptimisation
    use wing_nd_mod
    include 'aba_param.inc'

    real(real32) :: TIME(2), DTIME
    integer      :: LOP, LRESTART, KSTEP, KINC

    character*256   :: outdir, jobname
    integer         :: debugVar = 1, lenoutdir, lenjobname
    integer         :: iel
    integer         :: isub
    integer         :: ilay

    if (sroDebug) then 
        do while( debugVar /= 2 )
            debugVar = 1
        end do 
    end if

    if (.not. sroIsInit) then

        call GETOUTDIR(outdir, lenoutdir)
        sroOutdir = trim(outdir)
        
        call GETJOBNAME(jobname, lenjobname)
        sroJobname = trim(jobname)
        
        open(newunit=sroUnit, file=sroOutdir // SLASH // sroJobname // ".sro")
        call sroSetLogFile(sroUnit)
        call sroInitTimer()

        open(newunit=matrixUnit, file=sroOutdir // SLASH // sroJobname // "_MAT_.sro")

        call sroLog(" >> Sro::Wing << | uexternaldb | Initialising")

        call sroLog(" >> Sro::Wing << | uexternaldb | Job details: SRO_VEL_SCALING=" // int2str(SRO_VEL_SCALING))
#ifdef SRO_CONSTRAINT_STAB
        call sroLog(" >> Sro::Wing << | uexternaldb | Job details: SRO_CONSTRAINT_STAB=true")
#else
        call sroLog(" >> Sro::Wing << | uexternaldb | Job details: SRO_CONSTRAINT_STAB=false")
#endif

        call sroLog(" >> Sro::Wing << | uexternaldb | Creating mesh")
        sroMesh = FnmLsMesh(sroJobname, sroOutdir // SLASH, sroUnit)
        call sroMesh % registerSteps(   &
            [                           &
                character(len=50) ::    &
                "displacement",         &
                "adjoint",              &
                "velocity 1",           &
                "velocity 2",           &
                "velocity 3",           &
                "velocity 4",           &
                "velocity 5",           &
                "levelset"              &
            ]                           &         
        )
        call sroLog(" >> Sro::Wing << | uexternaldb | Initialising levelset")
        call sroMesh % initialiseLs()

        call sroMesh % updateAnalysisStatus
        sroIsInit = .true.

    end if

    call sroMesh % update(getStage(LOP + 1), KSTEP, KINC)

end subroutine UEXTERNALDB