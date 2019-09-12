c
c#############################################################################################
c
      SUBROUTINE CalculateSpecificHeat
c
C#############################################################################################
c
      Use User0, only: LSolveEnergy,ConstantSpecificHeat
      Use PhysicalProperties1, only: SpecificHeat,BSpecificHeat,
     *                               RGas,GammaGas,EquationOfState
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c      
      if(.not.LSolveEnergy) return
c
      if(nIter.gt.1) return
c
      if(EquationOfState.eq.'constant') then
c      
        do i=1,NumberOfElements
c
          SpecificHeat(i)=ConstantSpecificHeat
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeat(i,j)=ConstantSpecificHeat
c
          enddo 
        enddo 
c
      elseif(EquationOfState.eq.'tait') then
c      
        do i=1,NumberOfElements
c
          SpecificHeat(i)=ConstantSpecificHeat
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeat(i,j)=ConstantSpecificHeat
c
          enddo 
        enddo 
c      
      elseif(EquationOfState.eq.'idealgas') then
c
        do i=1,NumberOfElements
c
          SpecificHeat(i)=GammaGas*Rgas/(GammaGas-1.)
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeat(i,j)=GammaGas*Rgas/(GammaGas-1.)
c
          enddo 
        enddo 
c
      endif
c
      return
      end
