c
c#############################################################################################
c
      SUBROUTINE CalculateDensityrField(irField)
c
C#############################################################################################
c
      Use User0, only: ConstantDensityrField,BulkModulus
      Use PhysicalProperties1, only: DensityrField,BDensityrField,
     *                               RGasrField,EquationOfStaterField,
     *                               ReferencePressure,ReferenceDensity,
     *                               ThetaEOS
      use Variables1, only: Pressure,BPressure,Temperature,BTemperature
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,irField
c********************************************************************************************
c      
c
      if(EquationOfStaterField(irField).eq.'constant') then
c      
        if(nIter.gt.1) return
c
        do i=1,NumberOfElements
c
          DensityrField(i,irField)=ConstantDensityrField(irField)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BDensityrField(i,j,irField)=ConstantDensityrField(irField)
c
          enddo 
        enddo 
c      
      elseif(EquationOfStaterField(irField).eq.'tait') then
c
        do i=1,NumberOfElements
c
          DensityrField(i,irField)=
     *     ((1.+(Pressure(i)-ReferencePressure)*(ThetaEOS/
     *        BulkModulus))*(ReferenceDensity**ThetaEOS))**(1./ThetaEOS)
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BDensityrField(i,j,irField)=
     *     ((1.+(BPressure(i,j)-ReferencePressure)*(ThetaEOS/
     *        BulkModulus))*(ReferenceDensity**ThetaEOS))**(1./ThetaEOS)
c
          enddo 
        enddo 
c
      elseif(EquationOfStaterField(irField).eq.'idealgas') then
c
        do i=1,NumberOfElements
c
          DensityrField(i,irField)=Pressure(i)/
     *              (RGasrField(irField)*dmax1(Temperature(i),tiny))
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BDensityrField(i,j,irField)=BPressure(i,j)/
     *             (RGasrField(irField)*dmax1(BTemperature(i,j),tiny))
c
          enddo 
        enddo 
c
      endif
c
      return
      end