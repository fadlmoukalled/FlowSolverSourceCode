c
c#############################################################################################
c
      SUBROUTINE UpdateStaticTemperatureFromTotal
c
C#############################################################################################
c
      use BoundaryConditions2, only: 
     *                             IinletSpecifiedStagnationTemperature,
     *                        IinletSpecifiedStagnationTemperatureOwner,
     *               IinletSpecifiedStagnationTemperatureNumberOfBCSets,
     *                      IinletSpecifiedStagnationTemperatureNBFaces
      use Variables1, only: BuVelocity,BvVelocity,BwVelocity,
     *                      BTemperature,BStagnationTemperature
      use Geometry3, only: NIFaces,NBFaces
      use PhysicalProperties1, only: BSpecificHeat
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j
      double precision :: velocitySquared
c
c********************************************************************************************
c
        do i=1,IinletSpecifiedStagnationTemperature
c
          i1=IinletSpecifiedStagnationTemperatureOwner(i)
          i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          velocitySquared=BuVelocity(i2,i3)**2+
     *                       BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2
          BTemperature(i2,i3)=BStagnationTemperature(i2,i3)
     *                       -0.5*velocitySquared/BSpecificHeat(i2,i3)
c
        enddo
c
      return
      end
