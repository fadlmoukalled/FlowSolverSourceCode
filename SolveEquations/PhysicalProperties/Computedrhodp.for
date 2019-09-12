c
c#############################################################################################
c
      SUBROUTINE Calculatedrhodp
c
c#############################################################################################
c
      use User0, only: BulkModulus
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: Temperature,BTemperature
      use PhysicalProperties1, only: drhodP,BdrhodP,RGas,
     *                               EquationOfState,ThetaEOS,
     *                               ReferenceDensity,ReferencePressure,
     *                               Density,BDensity
      use Constants1, only: tiny
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,j
c
c********************************************************************************************
c
      if(EquationOfState.eq.'idealgas') then
c
        do i=1,NumberOfElements
c
          drhodP(i)=1./(RGas*Temperature(i)+tiny)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BdrhodP(i,j)=1./(RGas*BTemperature(i,j)+tiny)
c
          enddo
        enddo
c
      elseif(EquationOfState.eq.'tait') then
c
        do i=1,NumberOfElements
c
          drhodP(i)=(ReferenceDensity**ThetaEOS)/
     *             (BulkModulus*(Density(i)**(ThetaEOS-1)))
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BdrhodP(i,j)=(ReferenceDensity**ThetaEOS)/
     *             (BulkModulus*(BDensity(i,j)**(ThetaEOS-1)))
c
          enddo
        enddo
c
      elseif(EquationOfState.eq.'constant') then
c
        do i=1,NumberOfElements
c
          drhodP(i)=0.
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BdrhodP(i,j)=0.
c
          enddo
        enddo
c
      endif
c
      return
      end