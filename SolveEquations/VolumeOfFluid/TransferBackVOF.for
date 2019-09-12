c
c#############################################################################################
c
      SUBROUTINE TransferrFieldBack(irField)
c
C#############################################################################################
c
      use User0, only: LUnsteady
      use VolumeOfFluid1
      use TransferrField1
      use TransferVariable1, only: DiffusionCoefficientT,
     *                             BDiffusionCoefficientT
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,irField
c********************************************************************************************
c
      do i=1,NumberOfElements
c      
        rField(i,irField)=rFieldT(i)
        rFieldGradx(i,irField)=rFieldGradxT(i)
        rFieldGrady(i,irField)=rFieldGradyT(i)
        rFieldGradz(i,irField)=rFieldGradzT(i)
c
      enddo      
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BrField(i,j,irField)=BrFieldT(i,j)
          BrFieldGradx(i,j,irField)=BrFieldGradxT(i,j) 
          BrFieldGrady(i,j,irField)=BrFieldGradyT(i,j)
          BrFieldGradz(i,j,irField)=BrFieldGradzT(i,j)
c
        enddo
      enddo      
c      
      do i=1,NIFaces
c
        rFieldGradfx(i,irField)=rFieldGradfxT(i)
        rFieldGradfy(i,irField)=rFieldGradfyT(i)
        rFieldGradfz(i,irField)=rFieldGradfzT(i)
c
      enddo      
c
      if(LUnsteady) then
c
        do i=1,NumberOfElements
c      
        rFieldOld(i,irField)=rFieldOldT(i)
        rFieldOldOld(i,irField)=rFieldOldOldT(i)
c
        enddo      
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BrFieldOld(i,j,irField)=BrFieldOldT(i,j)
            BrFieldOldOld(i,j,irField)=BrFieldOldOldT(i,j) 
c
          enddo
        enddo      
c      
      endif
c
      return
      end