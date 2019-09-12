c
c#############################################################################################
c
      SUBROUTINE CalculateBoundaryHtotal
c
C#############################################################################################
c
      use BoundaryConditions2
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      BTemperature,BStagnationTemperature,
     *                      BHtotal,Htotal
      use PhysicalProperties1, only: BSpecificHeat,ReferenceTemperature
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: Tref      
c********************************************************************************************
c      
      Tref= ReferenceTemperature
c      
      do i=1,IWallDirichlet
c
        i1=IWallDirichletOwner(i)
        i2=IWallDirichletNumberOfBCSets(i)
        i3=IWallDirichletNBFaces(i)
c
        BHtotal(i2,i3)=BSpecificHeat(i2,3)*(BTemperature(i2,i3)-Tref)+
     *             0.5*(BuVelocity(i2,i3)**2+
     *                  BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2)
c
      enddo
c----------------------------------------------------------------------
      do i=1,IWallVonNeumann
c
        i1=IWallVonNeumannOwner(i)
        i2=IWallVonNeumannNumberOfBCSets(i)
        i3=IWallVonNeumannNBFaces(i)
c
        BHtotal(i2,i3)=BSpecificHeat(i2,3)*(BTemperature(i2,i3)-Tref)+
     *               0.5*(BuVelocity(i2,i3)**2+
     *                    BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2) 
c
      enddo
c----------------------------------------------------------------------
      do i=1,IWallRobin
c
        i2=IWallRobinNumberOfBCSets(i)
        i3=IWallRobinNBFaces(i)
c
        BHtotal(i2,i3)=BSpecificHeat(i2,3)*(BTemperature(i2,i3)-Tref)+
     *               0.5*(BuVelocity(i2,i3)**2+
     *                    BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2) 
c  
      enddo
!c----------------------------------------------------------------------
      do i=1,IinletSupersonic
c
        i2=IinletSupersonicNumberOfBCSets(i)
        i3=IinletSupersonicNBFaces(i)
c
        BHtotal(i2,i3)=BSpecificHeat(i2,3)*(BTemperature(i2,i3)-Tref)+
     *               0.5*(BuVelocity(i2,i3)**2+
     *                    BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2)
c
      enddo
c----------------------------------------------------------------------
      do i=1,IinletSpecifiedStaticTemperature
c
        i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
        i3=IinletSpecifiedStaticTemperatureNBFaces(i)
c
        BHtotal(i2,i3)=BSpecificHeat(i2,3)*(BTemperature(i2,i3)-Tref)+
     *               0.5*(BuVelocity(i2,i3)**2+
     *                    BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2)
c
      enddo
c----------------------------------------------------------------------
      do i=1,IinletSpecifiedStagnationTemperature
c
        i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
        i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
c
        BHtotal(i2,i3)=BSpecificHeat(i2,3)*
     *                       (BStagnationTemperature(i2,i3)-Tref)
        BTemperature(i2,i3)=BStagnationTemperature(i2,i3)-
     *             0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                  BwVelocity(i2,i3)**2)/BSpecificHeat(i2,i3)
c
      enddo
c
      return
      end
