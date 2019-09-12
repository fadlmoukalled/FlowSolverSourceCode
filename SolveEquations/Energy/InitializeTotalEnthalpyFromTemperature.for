c
c#############################################################################################
c
      SUBROUTINE InitializeHTotalFromTemperature
c
c#############################################################################################
c
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      Temperature,BTemperature,HTotal,
     *                      BHTotal
      use PhysicalProperties1, only: SpecificHeat,BSpecificHeat,
     *                               RGas,ReferenceTemperature
      use Constants1, only: tiny
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: velocity2,Tref
c
c********************************************************************************************
c
      Tref=ReferenceTemperature
c
      do i=1,NumberOfElements
c
        Htotal(i)=SpecificHeat(i)*(Temperature(i)-Tref)+
     *        0.5*(uVelocity(i)**2+vVelocity(i)**2+wVelocity(i)**2)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BHtotal(i,j)=BSpecificHeat(i,j)*(BTemperature(i,j)-Tref)+
     *              0.5*(BuVelocity(i,j)**2+
     *                    BvVelocity(i,j)**2+BwVelocity(i,j)**2)
c
        enddo
      enddo
c
      return
      end
