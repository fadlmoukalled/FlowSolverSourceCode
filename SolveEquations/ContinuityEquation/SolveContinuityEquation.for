c
c#############################################################################################
c
      SUBROUTINE SolveContinuity
c
c#############################################################################################
c
      use User0, only: Lcompressible,RRFContinuity,ASIterContinuity,
     *                 ASSolverContinuity,LMultigridContinuity,
     *                 LsolveEnergy,Linviscid,LfixPressure,
     *                 LFreeSurfaceFlow
      use Variables1, only: PressureC
      use BoundaryConditions2, only: LcalculateBeta,
     *             LRotationalPeriodicity,LTranslationalPeriodicity
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE SolveEquation(Variable,FiT,rrF,itmax,
     *                                       solver,LMultigrid)
c--------------------------------------------------------------
          character*10 Variable
          character*6 solver
          logical LMultigrid
          integer itmax  
          double precision rrF
          double precision, dimension(:)  :: FiT
c--------------------------------------------------------------
        end SUBROUTINE SolveEquation
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      Variable='pressc'
c
      call InitializeFluxes
      call InitializeCoefficients
      call InitializePressureCorrection
c
c---- Calculate mass flow rates at outflow boundary 
c
      if(.not.Lcompressible) then
c
        call FlowIn
	  call FlowOut
	  call AdjustFlow
c
      endif
c
c---- Calculate mass flow rates at internal and boundary cell faces 
c
      call AssembleMdot
c
      call UpdateBoundaryForTotalConditions
      call CalculateBoundaryMassFlowRates
      call UpdateOutletPressureForResistanceConditions
c
      if(LFreeSurfaceFlow) then
c
        call AssemblePressureCorrectionrField
c
      else
c
        call AssemblePressureCorrection
c
      endif
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
      if(LfixPressure) call FixPressure
c
	call SolveEquation(Variable,PressureC,RRFContinuity,
     *         ASIterContinuity,ASSolverContinuity,LMultigridContinuity)
c
      call UpdateBoundaryPressureCorrection
      call CorrectPressure
      if(LCompressible) call CorrectDensity
      call CorrectVelocity
      call Correctmdot
      call UpdateBoundaryPressure
      call CorrectBoundarymdot
c
      if(LTranslationalPeriodicity.and.LcalculateBeta) 
     *                                     call updateperiodicBeta
c
      If(Lcompressible.and.LsolveEnergy) call UpdateBoundaryTemperature
c
      call CalculateResidualsContinuity
c
c---- Update values along symmetry boundaries
c
      call UpdateSymmetry
c
c---- Update velocity field along slip walls
c
      call UpdateSlipWallBoundary
C
C----   If the flow is inviscid and temperature is not solved perform the following: 
C
      if(Lcompressible.and..not.LsolveEnergy.and.Linviscid) then
c
        call CalculateTemperature  
        call Calculatedrhodp  
        call CalculateDensity  
        call CalculateFaceDensity  
        call CalculateMachNumber  
c
      endif
c
	return
      end