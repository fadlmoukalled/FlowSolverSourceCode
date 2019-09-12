c
c#############################################################################################
c
      SUBROUTINE CalculateSurfaceTensionTerm(irField)
c
c#############################################################################################
c
      use User0, only: rFieldName,MethodCalcGradientrField,
     *                 nIterGradientrField,LimitGradientrField,
     *                 LimitGradientrFieldMethod
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFacesMax
      use VolumeOfFluid1, only: Curvature,BCurvature,
     *                          delrFieldMagnitude,BdelrFieldMagnitude,
     *                          delrFieldMagnitudeGradx,
     *                          BdelrFieldMagnitudeGradx,
     *                          delrFieldMagnitudeGrady,
     *                          BdelrFieldMagnitudeGrady,
     *                          delrFieldMagnitudeGradz,
     *                          BdelrFieldMagnitudeGradz,
     *                          rFieldGradxxT,rFieldGradxyT,
     *                          rFieldGradxzT,BrFieldGradxxT,
     *                          BrFieldGradxyT,BrFieldGradxzT,
     *                          rFieldGradyyT,BrFieldGradyyT,
     *                          rFieldGradyzT,BrFieldGradyzT,
     *                          rFieldGradzzT,BrFieldGradzzT
      use TransferrField1, only: rFieldT,BrFieldT,rFieldGradxT,
     *                           rFieldGradyT,rFieldGradzT,
     *                           BrFieldGradxT,BrFieldGradyT,
     *                           BrFieldGradzT
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 :: Variable
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
c--------------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------------
c
      delrFieldMagnitude=0.
      BdelrFieldMagnitude=0.
      delrFieldMagnitudeGradx=0.
      BdelrFieldMagnitudeGradx=0.
      delrFieldMagnitudeGrady=0.
      BdelrFieldMagnitudeGrady=0.
      delrFieldMagnitudeGradz=0.
      BdelrFieldMagnitudeGradz=0.
c
      rFieldGradxxT=0.
      rFieldGradxyT=0.
      rFieldGradxzT=0.
      BrFieldGradxxT=0.
      BrFieldGradxyT=0.
      BrFieldGradxzT=0.
      rFieldGradyyT=0.
      BrFieldGradyyT=0.
      rFieldGradyzT=0.
      BrFieldGradyzT=0.
      rFieldGradzzT=0.
      BrFieldGradzzT=0.
c
      Variable=rFieldName(irField)
c
      call Gradient(Variable,MethodCalcGradientrField(irField),
     *     rFieldT,rFieldGradxT,rFieldGradyT,rFieldGradzT,BrFieldT,
     *     BrFieldGradxT,BrFieldGradyT,
     *     BrFieldGradzT,nIterGradientrField(irField),
     *     LimitGradientrField(irField),
     *     LimitGradientrFieldMethod(irField))
      call Gradient(Variable,MethodCalcGradientrField(irField),
     *     rFieldGradxT,rFieldGradxxT,rFieldGradxyT,rFieldGradxzT,
     *     BrFieldGradxT,BrFieldGradxxT,BrFieldGradxyT,
     *     BrFieldGradxzT,nIterGradientrField(irField),
     *     LimitGradientrField(irField),
     *     LimitGradientrFieldMethod(irField))
      call Gradient(Variable,MethodCalcGradientrField(irField),
     *     rFieldGradyT,rFieldGradxyT,rFieldGradyyT,rFieldGradyzT,
     *     BrFieldGradyT,BrFieldGradxyT,BrFieldGradyyT,
     *     BrFieldGradyzT,nIterGradientrField(irField),
     *     LimitGradientrField(irField),
     *     LimitGradientrFieldMethod(irField))
      call Gradient(Variable,MethodCalcGradientrField(irField),
     *     rFieldGradzT,rFieldGradxzT,rFieldGradyzT,rFieldGradzzT,
     *     BrFieldGradzT,BrFieldGradxzT,BrFieldGradyzT,
     *     BrFieldGradzzT,nIterGradientrField(irField),
     *     LimitGradientrField(irField),
     *     LimitGradientrFieldMethod(irField))
c
      delrFieldMagnitude=dsqrt(rFieldGradxT*rFieldGradxT+
     *     rFieldGradyT*rFieldGradyT+rFieldGradzT*rFieldGradzT)    
      BdelrFieldMagnitude=dsqrt(BrFieldGradxT*BrFieldGradxT+
     *     BrFieldGradyT*BrFieldGradyT+BrFieldGradzT*BrFieldGradzT)    
c
      call Gradient(Variable,MethodCalcGradientrField(irField),
     *     delrFieldMagnitude,delrFieldMagnitudeGradx,
     *     delrFieldMagnitudeGrady,delrFieldMagnitudeGradz,
     *     BdelrFieldMagnitude,BdelrFieldMagnitudeGradx,
     *     BdelrFieldMagnitudeGrady,BdelrFieldMagnitudeGradz,
     *     nIterGradientrField(irField),
     *     LimitGradientrField(irField),
     *     LimitGradientrFieldMethod(irField))
c
      Curvature=((rFieldGradxT*delrFieldMagnitudeGradx+
     *       rFieldGradyT*delrFieldMagnitudeGrady+
     *       rFieldGradzT*delrFieldMagnitudeGradz)/
     *       (delrFieldMagnitude+tiny)-
     *       (rFieldGradxxT+rFieldGradyyT+rFieldGradzzT))/
     *       (delrFieldMagnitude+tiny)
      BCurvature=((BrFieldGradxT*BdelrFieldMagnitudeGradx+
     *       BrFieldGradyT*BdelrFieldMagnitudeGrady+
     *       BrFieldGradzT*BdelrFieldMagnitudeGradz)/
     *       (BdelrFieldMagnitude+tiny)-
     *       (BrFieldGradxxT+BrFieldGradyyT+BrFieldGradzzT))/
     *       (BdelrFieldMagnitude+tiny)
c
      return
      end
c