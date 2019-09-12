c
c#############################################################################################
c
      SUBROUTINE CalculateDensity
c
c#############################################################################################
c
      use User0, only: ConstantDensity,BulkModulus
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: Pressure,BPressure
      use PhysicalProperties1, only: drhodP,BdrhodP,Density,BDensity,
     *                               ReferencePressure,ReferenceDensity,
     *                               ThetaEOS,EquationOfState
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,j
c
c********************************************************************************************
c
      if(EquationOfState.eq.'constant') then
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
      elseif(EquationOfState.eq.'tait') then
c
        do i=1,NumberOfElements
c
          Density(i)=
     *     ((1.+(Pressure(i)-ReferencePressure)*(ThetaEOS/
     *        BulkModulus))*(ReferenceDensity**ThetaEOS))**(1./ThetaEOS)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BDensity(i,j)=
     *       ((1.+(BPressure(i,j)-ReferencePressure)*(ThetaEOS/
     *        BulkModulus))*(ReferenceDensity**ThetaEOS))**(1./ThetaEOS)
c
          enddo
        enddo
c
      elseif(EquationOfState.eq.'idealgas') then
c
        do i=1,NumberOfElements
c
          Density(i)=Pressure(i)*drhodP(i)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BDensity(i,j)=BPressure(i,j)*BdrhodP(i,j)
c
          enddo
        enddo
c
      endif
c
      return
      end