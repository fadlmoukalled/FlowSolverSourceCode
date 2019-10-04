c
c#############################################################################################
c
      SUBROUTINE SolveTurbulenceV2fEquation
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
      SUBROUTINE UpdateTurbulentV2ValueAtInlet(BFiT)
c--------------------------------------------------------------
        double precision, dimension(:,:) :: BFiT
c--------------------------------------------------------------
        end SUBROUTINE UpdateTurbulentV2ValueAtInlet
c********************************************************************************************
      end interface
C*********************************************************************************************
c
      Variable='tv2'
c
      call InitializeFluxes
      call InitializeCoefficients
      call updateFarFieldTurbulenceValues
     *                (Variable,TurbulentV2,BTurbulentV2)
      call UpdateTurbulentV2ValueAtInlet(BTurbulentV2)
c
c---- Calculate Gradients of the phi variable
c
      call Gradient(Variable,MethodCalcGradientTurbulentV2,
     *   TurbulentV2,TurbulentV2Gradx,TurbulentV2Grady,
     *     TurbulentV2Gradz,BTurbulentV2,BTurbulentV2Gradx,
     *             BTurbulentV2Grady,BTurbulentV2Gradz,
     *         nIterGradientTurbulentV2,LimitGradientTurbulentV2,
     *                              LimitGradientTurbulentV2Method)
c
      call InterpolateGradientToFace
     *     (GradientInterpolationSchemeTurbulentV2,
     *        TurbulentV2,TurbulentV2Gradx,TurbulentV2Grady,
     *         TurbulentV2Gradz,TurbulentV2Gradfx,TurbulentV2Gradfy,
     *                               TurbulentV2Gradfz)
c
c---- Assemble terms
c
      call CalculateEffectiveDiffusionCoefficient(Variable)
      call AssembleTurbulentV2Sources
      call AssembleDiffusionTerm(Variable,eDiffCoefficient,      
     *       BeDiffCoefficient,TurbulentV2,BTurbulentV2,
     *        BTurbulentV2Gradx,BTurbulentV2Grady,BTurbulentV2Gradz,
     *           TurbulentV2Gradfx,TurbulentV2Gradfy,TurbulentV2Gradfz)
c
       call AssembleConvectionTerm(Variable,BleedTurbulentV2,
     *     ConvectionSchemeTurbulentV2,HRFrameworkTurbulentV2,
     *       TurbulentV2,BTurbulentV2,TurbulentV2Gradx,
     *         TurbulentV2Grady,TurbulentV2Gradz,BTurbulentV2Gradx,
     *                           BTurbulentV2Grady,BTurbulentV2Gradz)
c
      if(LUnsteady) then
c
        call AssembleTransientTerm
     *       (Variable,TurbulentV2,TurbulentV2Old,TurbulentV2OldOld)
c
      endif
c
      call AssembleSourceTerm
     *            (Variable,TurbulentV2,ScTurbulentV2,SbTurbulentV2)
c
c---- Underrelax using false transient
c
      if(LFalseTransientTurbulentV2) 
     *        call AssembleFalseTransient(Variable,FalsedtTurbulentV2)
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxTurbulentV2) call UnderRelaxEquation(urfTurbulentV2)
c
c---- Solve equations
c
	call SolveEquation(Variable,TurbulentV2,rrfTurbulentV2,
     *    ASIterTurbulentV2,ASSolverTurbulentV2,LMultigridTurbulentV2)
c
      call BoundTurbulentV2
c
      call updateValuesFromBoundaryConditions(Variable,      
     *      eDiffCoefficient,BeDiffCoefficient,TurbulentV2,
     *        BTurbulentV2,BTurbulentV2Gradx,BTurbulentV2Grady,
     *           BTurbulentV2Gradz,TurbulentV2Gradfx,
     *                      TurbulentV2Gradfy,TurbulentV2Gradfz)
c
	return
	end