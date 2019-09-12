c
C#############################################################################################
c
      SUBROUTINE updateTemperatureFromHtotal
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Variables1, only: Temperature,Htotal,
     *                      uVelocity,vVelocity,wVelocity
      use PhysicalProperties1, only: SpecificHeat,ReferenceTemperature
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: Tref
c********************************************************************************************
c
      Tref=ReferenceTemperature
c
c--- Update static temperature at centroids of internal elements
c
      do i=1,NumberOfElements
c
        Temperature(i)=Tref+(Htotal(i)-0.5*(uVelocity(i)**2+
     *            vVelocity(i)**2+wVelocity(i)**2))/SpecificHeat(i)
c       
      enddo
c
      return
      end