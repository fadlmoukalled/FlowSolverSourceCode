c
C#############################################################################################
c
      SUBROUTINE AssembleKineticEnergyWorkTerm
c
C#############################################################################################
c
      use Variables3, only: FluxTE
      use variables1, only: uVelocity,vVelocity,wVelocity,
     *                      uVelGradx,vVelGrady,wVelGradz
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use PhysicalProperties1, only: Density
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i
      double precision FluxCElocal,divVel,kineticEnergy
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        divVel=uVelGradx(i)+vVelGrady(i)+wVelGradz(i)
        kineticEnergy=0.5*Density(i)*(uVelocity(i)*uVelocity(i)+
     *             vVelocity(i)*vVelocity(i)+ wVelocity(i)*wVelocity(i))
        FluxCElocal=-kineticEnergy*divVel*Volume(i)
        FluxTE(i)=FluxTE(i)+FluxCElocal
c        
      enddo

      return
      end