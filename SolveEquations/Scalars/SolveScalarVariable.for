c
c#############################################################################################
c
      SUBROUTINE SolveScalar(iscalar)
c
C#############################################################################################
c
      use User0
      use Variables1
      use Scalar1
      use Scalar2
      use TransferVariable1
      use Variables4
      use variables2, only: dphi
      use BoundaryConditions1
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: iScalar
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
        SUBROUTINE AssembleAnisotropicDiffusionTerm
     *       (Variable,FiT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *                            dfidxfT,dfidyfT,dfidzfT)
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
c--------------------------------------------------------------
        end SUBROUTINE AssembleAnisotropicDiffusionTerm
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
      end interface
C*********************************************************************************************
c
      iScalarVariable=iScalar
      Variable=ScalarName(iScalar)
      call TransferVariable(iScalar)
c
      call InitializeFluxes
      call InitializeCoefficients
c
c---- Calculate Gradients of the phi variable
c
      call Gradient(Variable,MethodCalcGradientScalar(iScalar),
     *     ScalarT,ScalarGradxT,ScalarGradyT,ScalarGradzT,BScalarT,
     *     BScalarGradxT,BScalarGradyT,
     *     BScalarGradzT,nIterGradientScalar(iScalar),
     *     LimitGradientScalar(iScalar),
     *     LimitGradientScalarMethod(iScalar))
c
      call InterpolateGradientToFace
     *      (GradientInterpolationSchemeScalar(iScalar),ScalarT,
     *        ScalarGradxT,ScalarGradyT,ScalarGradzT,ScalarGradfxT,
     *        ScalarGradfyT,ScalarGradfzT)
c
c---- Assemble terms
c
      if(LanisotropicDiffusion) then
c
        call AssembleAnisotropicDiffusionTerm(Variable,ScalarT,BScalarT,
     *            BScalarGradxT,BScalarGradyT,BScalarGradzT,
     *                  ScalarGradfxT,ScalarGradfyT,ScalarGradfzT)
c
      else
c
        call CalculateEffectiveDiffusionCoefficient(Variable)
        call AssembleDiffusionTerm(Variable,DiffusionCoefficientT,
     *       BDiffusionCoefficientT,ScalarT,BScalarT,BScalarGradxT,
     *                   BScalarGradyT,BScalarGradzT,ScalarGradfxT,
     *                                  ScalarGradfyT,ScalarGradfzT)
c
      endif
c
      if(LConvectScalar.or.LSolveMomentum)
     *   call AssembleConvectionTerm(Variable,BleedScalar(iScalar),
     *       ConvectionSchemeScalar(iScalar),LNVFScalar(iScalar),
     *       LTVDScalar(iScalar),ScalarT,BScalarT,ScalarGradxT,
     *       ScalarGradyT,ScalarGradzT,BScalarGradxT,
     *                          BScalarGradyT,BScalarGradzT)
c
      if(LUnsteady) then
c
        call AssembleTransientTerm(Variable,ScalarT,
     *                                  ScalarOldT,ScalarOldOldT)
c
      endif
c
      call AssembleSourceTerm(Variable,ScalarT,ScScalarT,SbScalarT)
      call AssemblePointSourceTerm(ScalarT,
     *                   ScPointSourceScalarT,SbPointSourceScalarT)
c
c---- Underrelax using false transient
c
      if(LFalseTransientScalar(iScalar)) 
     *     call AssembleFalseTransient(Variable,FalsedtScalar(iScalar))
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxScalar(iScalar)) 
     *        call UnderRelaxEquation(urfScalar(iScalar))
c
c---- Solve equations
c
	call SolveEquation(Variable,ScalarT,RRFScalar(iScalar),
     *         ASIterScalar(iScalar),ASSolverScalar(iScalar),
     *                                 LMultigridScalar(iScalar))
c
      call updateValuesFromBoundaryConditions(Variable,      
     *     DiffusionCoefficientT,BDiffusionCoefficientT,
     *     ScalarT,BScalarT,BScalarGradxT,BScalarGradyT,BScalarGradzT,
     *     ScalarGradfxT,ScalarGradfyT,ScalarGradfzT)
c
      call TransferVariableBack(iScalar)
c
	return
      end
