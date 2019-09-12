c
c#############################################################################################
c
      SUBROUTINE SetDensity
c
C#############################################################################################
c
      Use User0, only: ConstantDensity
      Use PhysicalProperties1, only: Density,BDensity
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,irField
c********************************************************************************************
c      
      if(nIter.gt.1) return
c
      do i=1,NumberOfElements
c
        Density(i)=ConstantDensity
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BDensity(i,j)=ConstantDensity
c
        enddo 
      enddo 
c      
      return
      end