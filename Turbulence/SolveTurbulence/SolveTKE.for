c
c#############################################################################################
c
      SUBROUTINE SolveTurbulenceKineticEnergy
c
C#############################################################################################
c
      use User0
      use PhysicalProperties1
      use Variables1
      use Variables4
c********************************************************************************************
c
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
     *        ConvectionScheme,NVF,TVD,FiT,BFiT,dfidxT,
     *              dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------
          character*10 Variable
          character*20 ConvectionScheme
          logical NVF,TVD
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
      SUBROUTINE AssembleTransientTerm(Variable,FiT,FiTold,FiToldold)
c--------------------------------------------------------------
        character*10 Variable
        double precision, dimension(:) :: FiT
        double precision, dimension(:) :: FiTold
        double precision, dimension(:) :: FiToldold
c--------------------------------------------------------------
      end SUBROUTINE AssembleTransientTerm
c--------------------------------------------------------------
        SUBROUTINE updateValuesFromBoundaryConditions
     *               (Variable,Gam,BGam,FiT,BFiT,BdfidxT,
     *                  BdfidyT,BdfidzT,dfidxfT,dfidyfT,dfidzfT)
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
        end SUBROUTINE updateValuesFromBoundaryConditions
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
        SUBROUTINE updateFarFieldTurbulenceValues(Variable,FiT,BFiT)
c--------------------------------------------------------------
          character*10 Variable
          double precision, dimension(:) :: FiT
          double precision, dimension(:,:) :: BFiT
c--------------------------------------------------------------
        end SUBROUTINE updateFarFieldTurbulenceValues
c--------------------------------------------------------------
      end interface
C*********************************************************************************************
c
      Variable='tke'
c
      call InitializeFluxes
      call InitializeCoefficients
c
      call CalculateInletTurbulence
      call updateFarFieldTurbulenceValues
     *                (Variable,TurbulentKE,BTurbulentKE)
c
c---- Calculate Gradients of the phi variable
c
      call Gradient(Variable,MethodCalcGradientTKE,
     *      TurbulentKE,TKEGradx,TKEGrady,TKEGradz,BTurbulentKE,
     *            BTKEGradx,BTKEGrady,BTKEGradz,nIterGradientTKE,
     *             LimitGradientTKE,LimitGradientTKEMethod)
c
      call InterpolateGradientToFace(GradientInterpolationSchemeTKE,
     *            TurbulentKE,TKEGradx,TKEGrady,TKEGradz,
     *                         TKEGradfx,TKEGradfy,TKEGradfz)
c
c---- Assemble terms
c
      if(TurbulenceModel.eq.'komega2006'.or.
     *               TurbulenceModel.eq.'komega2006lrn') then
        call ModifyEffectiveDiffusionCoefficient(Variable)
      else
        call CalculateEffectiveDiffusionCoefficient(Variable)
      endif
      call AssembleTurbulentKESources
      if(LBuoyancy.and.LBuoyancyCorrection) 
     *            call AssembleTurbulentKEBuoyancySources
      if(LCompressible.and.LCompressibilityCorrection) 
     *            call AssembleTurbulentKECompressibilitySources
      call AssembleDiffusionTerm(Variable,      
     *         eDiffCoefficient,BeDiffCoefficient,TurbulentKE,
     *             BTurbulentKE,BTKEGradx,BTKEGrady,BTKEGradz,
     *                           TKEGradfx,TKEGradfy,TKEGradfz)
c
       call AssembleConvectionTerm(Variable,BleedTKE,
     *        ConvectionSchemeTKE,LNVFtke,LTVDtke,TurbulentKE,
     *          BTurbulentKE,TKEGradx,TKEGrady,TKEGradz,
     *                      BTKEGradx,BTKEGrady,BTKEGradz)
c
      if(LUnsteady) then
c
        call AssembleTransientTerm
     *      (Variable,TurbulentKE,TurbulentKEOld,TurbulentKEOldOld)
c
      endif
c
      call AssembleSourceTerm(Variable,TurbulentKE,ScTKE,SbTKE)
c
c---- Underrelax using false transient
c
      if(LFalseTransientTKE) 
     *        call AssembleFalseTransient(Variable,FalsedtTKE)
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxTKE) call UnderRelaxEquation(urfTKE)
c
c---- Solve equations
c
	call SolveEquation(Variable,TurbulentKE,rrfTKE,
     *                  ASIterTKE,ASSolverTKE,LMultigridTKE)
c
      call BoundTurbulentKE
c
      call updateValuesFromBoundaryConditions(Variable,      
     *      eDiffCoefficient,BeDiffCoefficient,TurbulentKE,
     *           BTurbulentKE,BTKEGradx,BTKEGrady,BTKEGradz,
     *                           TKEGradfx,TKEGradfy,TKEGradfz)
c
	return
      end