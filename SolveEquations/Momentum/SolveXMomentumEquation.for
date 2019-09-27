c
C#############################################################################################
c
      SUBROUTINE SolveMomentumX
c
C#############################################################################################
c
      use User0
      use Variables1
      use Variables4
      use PhysicalProperties1
      use variables2, only: dphi
      use BoundaryConditions1
      use BoundaryConditions2, only: LTranslationalPeriodicity
      use Coriolis1
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
c********************************************************************************************
c--- Interfaces
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------
          character*10 Variable
          integer :: MethodCalcGradient,nIterGradientPhi
          logical :: LimitGradient
          integer :: LimitGradientMethod
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------
        SUBROUTINE InterpolateGradientToFace
     *   (GradientInterpolationScheme,FiT,dfidxT,dfidyT,dfidzT,
     *                                  dfidxfT,dfidyfT,dfidzfT)
c--------------------------------------------------------------
          character*16 GradientInterpolationScheme
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:) :: dfidxfT
          double precision, dimension(:) :: dfidyfT
          double precision, dimension(:) :: dfidzfT
c--------------------------------------------------------------
        end SUBROUTINE InterpolateGradientToFace
c--------------------------------------------------------------
        SUBROUTINE AssembleDiffusionTerm(Variable,
     *       Gam,BGam,FiT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *                           dfidxfT,dfidyfT,dfidzfT)
c--------------------------------------------------------------
          character*10 Variable
          double precision, dimension(:) :: FiT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
          double precision, dimension(:) :: dfidxfT
          double precision, dimension(:) :: dfidyfT
          double precision, dimension(:) :: dfidzfT
          double precision, dimension(:) :: Gam
          double precision, dimension(:,:) :: Bgam
c--------------------------------------------------------------
        end SUBROUTINE AssembleDiffusionTerm
c--------------------------------------------------------------
        SUBROUTINE AssembleConvectionTerm(Variable,Bleed,
     *        ConvectionScheme,HRFramework,FiT,BFiT,dfidxT,
     *              dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------
          character*10 Variable
          character*20 ConvectionScheme
          character*4 HRFramework
          double precision :: Bleed
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------
        end SUBROUTINE AssembleConvectionTerm
c--------------------------------------------------------------
        SUBROUTINE AssembleSourceTerm(Variable,FiT,Sc,Sb)
          character*10 Variable
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: Sc
          double precision, dimension(:) :: Sb
c--------------------------------------------------------------
        end SUBROUTINE AssembleSourceTerm
c--------------------------------------------------------------
        SUBROUTINE AssemblePointSourceTerm(FiT,ScPointSource,
     *                                                  SbPointSource)
c--------------------------------------------------------------
          double precision, dimension(:)  :: FiT
          double precision, dimension(:)  :: ScPointSource
          double precision, dimension (:) :: SbPointSource
c--------------------------------------------------------------
        end SUBROUTINE AssemblePointSourceTerm
c--------------------------------------------------------------
      SUBROUTINE AssembleTransientTerm(Variable,FiT,FiTold,FiToldold)
c--------------------------------------------------------------
        character*10 Variable
        double precision, dimension(:) :: FiT
        double precision, dimension(:) :: FiTold
        double precision, dimension(:) :: FiToldold
c--------------------------------------------------------------
      end SUBROUTINE AssembleTransientTerm
c--------------------------------------------------------------
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
        SUBROUTINE AssemblePressureGradientTerm(dPdxy)
c--------------------------------------------------------------
          double precision, dimension(:) :: dPdxy
c--------------------------------------------------------------
        end SUBROUTINE AssemblePressureGradientTerm
c--------------------------------------------------------------
        SUBROUTINE AssembleCoriolisSourceTerm(phi1,phi2,phi3,phi4)
c--------------------------------------------------------------
          double precision :: phi2,phi4
          double precision, dimension(:) :: phi1
          double precision, dimension(:) :: phi3
c--------------------------------------------------------------
        end SUBROUTINE AssembleCoriolisSourceTerm
c--------------------------------------------------------------
        SUBROUTINE AssembleExplicitDiffusionTerm(Variable,Gam,BGam)
c--------------------------------------------------------------
          implicit none
          character*10 Variable
          double precision, dimension(:) :: Gam
          double precision, dimension(:,:) :: BGam
c--------------------------------------------------------------
        end SUBROUTINE AssembleExplicitDiffusionTerm
c--------------------------------------------------------------
        SUBROUTINE AssembleBulkViscosityTerm(Variable,Gam,BGam)
c--------------------------------------------------------------
          implicit none
          character*10 Variable
          double precision, dimension(:) :: Gam
          double precision, dimension(:,:) :: BGam
c--------------------------------------------------------------
        end SUBROUTINE AssembleBulkViscosityTerm
c--------------------------------------------------------------
        SUBROUTINE storeDvalues(dphi1,dphi2,dphi3,dphi4)
c--------------------------------------------------------------
          double precision, dimension(:) :: dphi1,dphi2,dphi3,dphi4
c--------------------------------------------------------------
        end SUBROUTINE storeDvalues
c--------------------------------------------------------------
      end interface
C*********************************************************************************************
c
      Variable='velx'
c
      call InitializeFluxes
      call InitializeCoefficients
c
c---- Calculate Gradients of the u,v,w, and p variable
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelocity,uVelGradx,uVelGrady,uVelGradz,BuVelocity,
     *       BuVelGradx,BuVelGrady,BuVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      Variable='vely'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelocity,vVelGradx,vVelGrady,vVelGradz,BvVelocity,
     *        BvVelGradx,BvVelGrady,BvVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      Variable='velz'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelocity,wVelGradx,wVelGrady,wVelGradz,BwVelocity,
     *        BwVelGradx,BwVelGrady,BwVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      Variable='press'
c
      call Gradient(Variable,MethodCalcGradientContinuity,
     *   Pressure,PressGradx,PressGrady,PressGradz,BPressure,
     *      BPressGradx,BPressGrady,BPressGradz,nIterGradientContinuity,
     *           LimitGradientContinuity,LimitGradientContinuityMethod)
c
      Variable='velx'
c
      call InterpolateGradientToFace
     *  (GradientInterpolationSchemeMomentum,uVelocity,uVelGradx,
     *             uVelGrady,uVelGradz,uVelGradfx,uVelGradfy,uVelGradfz)
      call InterpolateGradientToFace
     *  (GradientInterpolationSchemeMomentum,vVelocity,vVelGradx,
     *             vVelGrady,vVelGradz,vVelGradfx,vVelGradfy,vVelGradfz)
      call InterpolateGradientToFace
     *  (GradientInterpolationSchemeMomentum,wVelocity,wVelGradx,
     *             wVelGrady,wVelGradz,wVelGradfx,wVelGradfy,wVelGradfz)
      call InterpolateGradientToFace
     *     (GradientInterpolationSchemeContinuity,
     *          Pressure,PressGradx,PressGrady,PressGradz,
     *                   PressGradfx,PressGradfy,PressGradfz)
c
c---- Assemble terms
c
      if(.not.Linviscid) then
c
        call CalculateEffectiveDiffusionCoefficient(Variable)
        call AssembleDiffusionTerm(Variable,      
     *        eDiffCoefficient,BeDiffCoefficient,uVelocity,BuVelocity,
     *             BuVelGradx,BuVelGrady,BuVelGradz,
     *                     uVelGradfx,uVelGradfy,uVelGradfz)
c
        call AssembleExplicitDiffusionTerm
     *                  (Variable,eDiffCoefficient,BeDiffCoefficient)
c
        call AssembleBulkViscosityTerm
     *                  (Variable,eDiffCoefficient,BeDiffCoefficient)
c
      endif
c
      call AssembleConvectionTerm(Variable,BleedMomentum,
     *    ConvectionSchemeMomentum,HRFrameworkMomentum,
     * uVelocity,BuVelocity,uVelGradx,uVelGrady,uVelGradz,
     *                    BuVelGradx,BuVelGrady,BuVelGradz)
c
      if(LUnsteady) then
c
        call AssembleTransientTerm
     *        (Variable,uVelocity,uVelocityOld,uVelocityOldOld)
c
      endif
c
      if(LBuoyancy.and..not.LsolveEnergy) 
     *                         call reDistributeTheBuoyancyTerm
      if(LBuoyancy) call AssembleBuoyancyTerm(Variable)
c
      if(LTranslationalPeriodicity) 
     *                call AssembleLinearPressureTerm(Variable)
c
      call AssembleSourceTerm(Variable,uVelocity,
     *                            ScMomentumx,SbMomentumx)
      call AssemblePressureGradientTerm(PressGradx)
      if(LTurbulentFlow) call MinusTwoThirdsRhoK(Variable)
      if(LCoriolis) then
c
        call AssembleCentrifugalForceSourceTerm(Variable)
        Call AssembleCoriolisSourceTerm(WVelocity,
     *      AngularVelocityY,vVelocity,AngularVelocityZ)
c
      endif       
c
      if(LFreeSurfaceFlow.and.LSurfaceTension) 
     *           call AssembleSurfaceTensionSourceTerm(Variable)
c
c---- Underrelax using false transient
c
      if(LFalseTransientMomentum) 
     *        call AssembleFalseTransient(Variable,FalsedtMomentum)
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxMomentum) call UnderRelaxEquation(urfMomentum)
c
c---- Store Volume/ac and u*
c
      call storeDvalues(Du1Velocity,Du2Velocity,uVelocity,uVelocityStar)
c
c---- Solve equations
c
	call SolveEquation(Variable,uVelocity,RRFMomentum,
     *           ASIterMomentum,ASSolverMomentum,LMultigridMomentum)
c
      call UpdatePeriodic
c
	return
	end
