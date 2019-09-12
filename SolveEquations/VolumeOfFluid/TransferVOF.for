c
c#############################################################################################
c
      SUBROUTINE TransferrField(irField)
c
C#############################################################################################
c
      use User0, only: LUnsteady
      use VolumeOfFluid1
      use TransferrField1
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use PhysicalProperties1, only: DiffusionCoefficient,
     *                               BDiffusionCoefficient
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,irField
c********************************************************************************************
c
      do i=1,NumberOfElements
c      
        rFieldT(i)=rField(i,irField)
        rFieldGradxT(i)=rFieldGradx(i,irField)
        rFieldGradyT(i)=rFieldGrady(i,irField)
        rFieldGradzT(i)=rFieldGradz(i,irField)
c
      enddo      
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BrFieldT(i,j)=BrField(i,j,irField)
          BrFieldGradxT(i,j)=BrFieldGradx(i,j,irField)
          BrFieldGradyT(i,j)=BrFieldGrady(i,j,irField)
          BrFieldGradzT(i,j)=BrFieldGradz(i,j,irField)
c
        enddo
      enddo      
c      
      do i=1,NIFaces
c
        rFieldGradfxT(i)=rFieldGradfx(i,irField)
        rFieldGradfyT(i)=rFieldGradfy(i,irField)
        rFieldGradfzT(i)=rFieldGradfz(i,irField)
c
      enddo      
c
      if(LUnsteady) then
c
        do i=1,NumberOfElements
c      
        rFieldOldT(i)=rFieldOld(i,irField)
        rFieldOldOldT(i)=rFieldOldOld(i,irField)
c
        enddo      
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BrFieldOldT(i,j)=BrFieldOld(i,j,irField)
            BrFieldOldOldT(i,j)= BrFieldOldOld(i,j,irField)
c
          enddo
        enddo      
c      
      endif
c
      return
      end