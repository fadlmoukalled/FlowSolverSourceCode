c
c#############################################################################################
c
      SUBROUTINE SolveTurbulenceSpecificDissipationRate
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
      Variable='tomega'
c
      call InitializeFluxes
      call InitializeCoefficients
      call updateFarFieldTurbulenceValues
     *                (Variable,TurbulentOmega,BTurbulentOmega)
c
c---- Calculate Gradients of the phi variable
c
      call Gradient(Variable,MethodCalcGradientTOmega,
     *    TurbulentOmega,TOmegaGradx,TOmegaGrady,TOmegaGradz,
     *     BTurbulentOmega,BTOmegaGradx,BTOmegaGrady,BTOmegaGradz,
     *         nIterGradientTOmega,LimitGradientTOmega,
     *                              LimitGradientTOmegaMethod)
c
      call InterpolateGradientToFace(GradientInterpolationSchemeTOmega,
     *        TurbulentOmega,TOmegaGradx,TOmegaGrady,TOmegaGradz,
     *                        TOmegaGradfx,TOmegaGradfy,TOmegaGradfz)
c
c---- Assemble terms
c
      if(TurbulenceModel.eq.'komega2006'.or.
     *           TurbulenceModel.eq.'komega2006lrn') then
        call ModifyEffectiveDiffusionCoefficient(Variable)
      else
        call CalculateEffectiveDiffusionCoefficient(Variable)
      endif
      call AssembleTurbulentOmegaSources
      if(LBuoyancy.and.LBuoyancyCorrection) 
     *          call AssembleTurbulentOmegaBuoyancySources
      if(LCompressible.and.LCompressibilityCorrection) 
     *          call AssembleTOmegaCompressibilitySources
      call AssembleDiffusionTerm(Variable,eDiffCoefficient,      
     *         BeDiffCoefficient,TurbulentOmega,BTurbulentOmega,
     *          BTOmegaGradx,BTOmegaGrady,BTOmegaGradz,
     *                TOmegaGradfx,TOmegaGradfy,TOmegaGradfz)
c
      call AssembleConvectionTerm(Variable,BleedTOmega,
     *     ConvectionSchemeTOmega,LNVFTOmega,LTVDTOmega,TurbulentOmega,
     *       BTurbulentOmega,TOmegaGradx,TOmegaGrady,TOmegaGradz,
     *                      BTOmegaGradx,BTOmegaGrady,BTOmegaGradz)
c
      if(LUnsteady) then
c
        call AssembleTransientTerm
     *  (Variable,TurbulentOmega,TurbulentOmegaOld,TurbulentOmegaOldOld)
c
      endif
c
      call AssembleSourceTerm(Variable,TurbulentOmega,ScTOmega,SbTOmega)
c
c---- Underrelax using false transient
c
      if(LFalseTransientTOmega) 
     *        call AssembleFalseTransient(Variable,FalsedtTOmega)
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxTOmega) call UnderRelaxEquation(urfTOmega)
c
      call ModifyCoefficientsAlongWalls
c
c---- Solve equations
c
	call SolveEquation(Variable,TurbulentOmega,rrfTOmega,
     *                  ASIterTOmega,ASSolverTOmega,LMultigridTOmega)
c
      call BoundTurbulentOmega
c
      call updateValuesFromBoundaryConditions(Variable,      
     *      eDiffCoefficient,BeDiffCoefficient,TurbulentOmega,
     *        BTurbulentOmega,BTOmegaGradx,BTOmegaGrady,BTOmegaGradz,
     *                         TOmegaGradfx,TOmegaGradfy,TOmegaGradfz)
c
	return
	end