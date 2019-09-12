c
c#############################################################################################
c
      SUBROUTINE CalculateMachNumber
c
c#############################################################################################
c
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      Temperature,BTemperature,MachNumber,
     *                      BMachNumber
      use PhysicalProperties1, only: SpecificHeat,BSpecificHeat,RGas
      use Constants1, only: tiny
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: velocity2,cpr
c
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        velocity2=uVelocity(i)**2+vVelocity(i)**2+wVelocity(i)**2
        cpr=SpecificHeat(i)*RGas/(SpecificHeat(i)-RGas)
        MachNumber(i)=dsqrt(velocity2/(cpr*Temperature(i)+tiny))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          velocity2=BuVelocity(i,j)**2+
     *                    BvVelocity(i,j)**2+BwVelocity(i,j)**2
          cpr=BSpecificHeat(i,j)*RGas/(BSpecificHeat(i,j)-RGas)
          BMachNumber(i,j)=dsqrt(velocity2/(cpr*BTemperature(i,j)+tiny))
c
        enddo
      enddo
c
      return
      end
