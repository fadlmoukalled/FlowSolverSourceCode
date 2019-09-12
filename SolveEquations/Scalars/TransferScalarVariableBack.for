c
c#############################################################################################
c
      SUBROUTINE TransferVariableBack(iscalar)
c
C#############################################################################################
c
      use User0, only: LUnsteady
      use Scalar1
      use TransferVariable1
      use BoundaryConditions1, only: NumberOfPointSources
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,iScalar
c********************************************************************************************
c
      do i=1,NumberOfElements
c      
        Scalar(i,iScalar)=ScalarT(i)
        ScalarGradx(i,iScalar)=ScalarGradxT(i)
        ScalarGrady(i,iScalar)=ScalarGradyT(i)
        ScalarGradz(i,iScalar)=ScalarGradzT(i)
c
      enddo      
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BScalar(i,j,iScalar)=BScalarT(i,j)
          BScalarGradx(i,j,iScalar)=BScalarGradxT(i,j) 
          BScalarGrady(i,j,iScalar)=BScalarGradyT(i,j)
          BScalarGradz(i,j,iScalar)=BScalarGradzT(i,j)
c
        enddo
      enddo      
c      
      do i=1,NIFaces
c
        ScalarGradfx(i,iScalar)=ScalarGradfxT(i)
        ScalarGradfy(i,iScalar)=ScalarGradfyT(i)
        ScalarGradfz(i,iScalar)=ScalarGradfzT(i)
c
      enddo      
c
      if(LUnsteady) then
c
        do i=1,NumberOfElements
c      
        ScalarOld(i,iScalar)=ScalarOldT(i)
        ScalarOldOld(i,iScalar)=ScalarOldOldT(i)
c
        enddo      
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BScalarOld(i,j,iScalar)=BScalarOldT(i,j)
            BScalarOldOld(i,j,iScalar)=BScalarOldOldT(i,j) 
c
          enddo
        enddo      
c      
      endif
c
      return
      end