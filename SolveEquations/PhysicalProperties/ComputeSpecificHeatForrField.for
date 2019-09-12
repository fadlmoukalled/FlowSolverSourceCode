c
c#############################################################################################
c
      SUBROUTINE CalculateSpecificHeatrField(irField)
c
C#############################################################################################
c
      Use User0, only: LSolveEnergy,ConstantSpecificHeatrField
      Use PhysicalProperties1, only: SpecificHeatrField,
     *                               BSpecificHeatrField,
     *                               RGasrField,GammaGasrField,
     *                               EquationOfStaterField
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,irField
c********************************************************************************************
c      
      if(.not.LSolveEnergy) return
c
      if(nIter.gt.1) return
c
      if(EquationOfStaterField(irField).eq.'constant') then
c      
        do i=1,NumberOfElements
c
          SpecificHeatrField(i,irField)=
     *                   ConstantSpecificHeatrField(irField)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeatrField(i,j,irField)=
     *                   ConstantSpecificHeatrField(irField)
c
          enddo 
        enddo 
c
      elseif(EquationOfStaterField(irField).eq.'tait') then
c      
        do i=1,NumberOfElements
c
          SpecificHeatrField(i,irField)=
     *                   ConstantSpecificHeatrField(irField)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeatrField(i,j,irField)=
     *                   ConstantSpecificHeatrField(irField)
c
          enddo 
        enddo 
c      
      elseif(EquationOfStaterField(irField).eq.'idealgas') then
c
        do i=1,NumberOfElements
c
          SpecificHeatrField(i,irField)=GammaGasrField(irField)*
     *                  RGasrField(irField)/(GammaGasrField(irField)-1.)
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeatrField(i,j,irField)=GammaGasrField(irField)*
     *                  RGasrField(irField)/(GammaGasrField(irField)-1.)
c
          enddo 
        enddo 
c
      endif
c
      return
      end
