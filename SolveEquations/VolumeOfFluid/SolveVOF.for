c
c#############################################################################################
c
      SUBROUTINE SolveVolumeOfFluid(irField)
c
C#############################################################################################
c
      use User0
      use TransferrField1
      use VolumeOfFluid1
      use VolumeOfFluid2, only: irFieldVariable
      use TransferVariable1
      use Variables4
      use variables2, only: dphi
      use BoundaryConditions1
      use PhysicalProperties1, only: Viscosity,Bviscosity
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: irField
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
        SUBROUTINE AssembleConvectionTermrField(Variable,Bleed,
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
        end SUBROUTINE AssembleConvectionTermrField
c--------------------------------------------------------------
      SUBROUTINE AssembleTransientTermrField
     *                       (Variable,FiT,FiTold,FiToldold)
c--------------------------------------------------------------
        character*10 Variable
        double precision, dimension(:) :: FiT
        double precision, dimension(:) :: FiTold
        double precision, dimension(:) :: FiToldold
c--------------------------------------------------------------
      end SUBROUTINE AssembleTransientTermrField
c--------------------------------------------------------------
        SUBROUTINE updateValuesFromBoundaryConditionsrField
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
        end SUBROUTINE updateValuesFromBoundaryConditionsrField
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
      if(.not.LSolveMomentum.and..not.LConvectScalar) then
c
        write(*,*) 'cannot solve the volume of fluid equation' 
        write(*,*) 'because the momentum equation is not solved.'
        write(*,*) 'Program has stopped'
        stop
c
      endif      
c
      irFieldVariable=irField
      Variable=rFieldName(irField)
      call TransferrField(irField)
c
      call InitializeFluxes
      call InitializeCoefficients
c
c---- Calculate Gradients of the phi variable
c
      call Gradient(Variable,MethodCalcGradientrField(irField),
     *     rFieldT,rFieldGradxT,rFieldGradyT,rFieldGradzT,BrFieldT,
     *     BrFieldGradxT,BrFieldGradyT,
     *     BrFieldGradzT,nIterGradientrField(irField),
     *     LimitGradientrField(irField),
     *     LimitGradientrFieldMethod(irField))
c
      call InterpolateGradientToFace
     *      (GradientInterpolationSchemerField(irField),rFieldT,
     *        rFieldGradxT,rFieldGradyT,rFieldGradzT,rFieldGradfxT,
     *        rFieldGradfyT,rFieldGradfzT)
c
      call AssembleConvectionTermrField(Variable,BleedrField(irField),
     *     ConvectionSchemerField(irField),LNVFrField(irField),
     *       LTVDrField(irField),rFieldT,BrFieldT,rFieldGradxT,
     *          rFieldGradyT,rFieldGradzT,BrFieldGradxT,
     *                          BrFieldGradyT,BrFieldGradzT)
c
      if(LUnsteady) call AssembleTransientTermrField
     *                 (Variable,rFieldT,rFieldOldT,rFieldOldOldT)
c
c---- Underrelax using false transient
c
      if(LFalseTransientrField(irField)) 
     *     call AssembleFalseTransientrField
     *                     (Variable,FalsedtrField(irField))
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxrField(irField)) 
     *        call UnderRelaxEquation(urfrField(irField))
c
c---- Solve equations
c
	call SolveEquation(Variable,rFieldT,RRFrField(irField),
     *         ASIterrField(irField),ASSolverrField(irField),
     *                                 LMultigridrField(irField))
c
      call updateValuesFromBoundaryConditionsrField(Variable,      
     *     Viscosity,BViscosity,rFieldT,BrFieldT,BrFieldGradxT,
     *        BrFieldGradyT,BrFieldGradzT,rFieldGradfxT,
     *                          rFieldGradfyT,rFieldGradfzT)
c
      call TransferrFieldBack(irField)
c
	return
      end