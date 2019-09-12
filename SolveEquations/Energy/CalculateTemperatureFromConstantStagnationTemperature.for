c
c#############################################################################################
c
      SUBROUTINE CalculateTemperature
c
c#############################################################################################
c
      use User0, only: constantStagnationTemperature,urfEnergy
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: uVelocity,vVelocity,wVelocity,Temperature,
     *                    BuVelocity,BvVelocity,BwVelocity,BTemperature
      use PhysicalProperties1, only: SpecificHeat,BSpecificHeat
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: velocity2,urfT
c********************************************************************************************
c
      urfT=1.-urfEnergy
      do i=1,NumberOfElements
c
        velocity2=uVelocity(i)**2+vVelocity(i)**2+wVelocity(i)**2
        velocity2=velocity2/(2.*SpecificHeat(i))
        Temperature(i)=urfT*Temperature(i)+
     *          urfEnergy*(constantStagnationTemperature-velocity2)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          velocity2=BuVelocity(i,j)**2+
     *               BvVelocity(i,j)**2+BwVelocity(i,j)**2
          velocity2=velocity2/(2.*BSpecificHeat(i,j))
          BTemperature(i,j)=urfT*BTemperature(i,j)+
     *          urfEnergy*(constantStagnationTemperature-velocity2)
c
        enddo
      enddo
c
      return
      end
