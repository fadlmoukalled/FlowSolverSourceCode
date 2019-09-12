c
c#############################################################################################
c
      SUBROUTINE TransferVariable(iscalar)
c
C#############################################################################################
c
      use User0, only: LUnsteady
      use Scalar1
      use TransferVariable1
      use BoundaryConditions1, only: NumberOfPointSources
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use PhysicalProperties1, only: DiffusionCoefficient,
     *                               BDiffusionCoefficient
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,iScalar
c********************************************************************************************
c
      do i=1,NumberOfElements
c      
        ScalarT(i)=Scalar(i,iScalar)
        ScalarGradxT(i)=ScalarGradx(i,iScalar)
        ScalarGradyT(i)=ScalarGrady(i,iScalar)
        ScalarGradzT(i)=ScalarGradz(i,iScalar)
        ScScalarT(i)=ScScalar(i,iScalar)
        SbScalarT(i)=SbScalar(i,iScalar)
        DiffusionCoefficientT(i)=DiffusionCoefficient(i,iScalar)
c
      enddo      
c
      if(NumberOfPointSources.gt.0) then
c
        do i=1,NumberOfPointSources
c      
          ScPointSourceScalarT(i)=ScPointSourceScalar(i,iScalar)
          SbPointSourceScalarT(i)=SbPointSourceScalar(i,iScalar)
c
        enddo      
c      
      endif
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BScalarT(i,j)=BScalar(i,j,iScalar)
          BScalarGradxT(i,j)=BScalarGradx(i,j,iScalar)
          BScalarGradyT(i,j)=BScalarGrady(i,j,iScalar)
          BScalarGradzT(i,j)=BScalarGradz(i,j,iScalar)
          BDiffusionCoefficientT(i,j)=BDiffusionCoefficient(i,j,iScalar)
c
        enddo
      enddo      
c      
      do i=1,NIFaces
c
        ScalarGradfxT(i)=ScalarGradfx(i,iScalar)
        ScalarGradfyT(i)=ScalarGradfy(i,iScalar)
        ScalarGradfzT(i)=ScalarGradfz(i,iScalar)
c
      enddo      
c
      if(LUnsteady) then
c
        do i=1,NumberOfElements
c      
        ScalarOldT(i)=ScalarOld(i,iScalar)
        ScalarOldOldT(i)=ScalarOldOld(i,iScalar)
c
        enddo      
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BScalarOldT(i,j)=BScalarOld(i,j,iScalar)
            BScalarOldOldT(i,j)= BScalarOldOld(i,j,iScalar)
c
          enddo
        enddo      
c      
      endif
c
      return
      end
