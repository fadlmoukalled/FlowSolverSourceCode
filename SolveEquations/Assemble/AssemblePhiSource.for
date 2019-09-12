c
C#############################################################################################
      SUBROUTINE AssemblePhiSourceTerm
C#############################################################################################
      use User0, only: LTurbulentFlow
      use Variables1, only: uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz
      use Variables3, only: FluxTE
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use PhysicalProperties1, only:Viscosity,TurbulentViscosity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: FluxTElocal,term1,visc
c********************************************************************************************
c
      if(LTurbulentFlow) then
c
        do i=1,NumberOfElements
c
          visc=Viscosity(i)+TurbulentViscosity(i)
          term1=2.*(uVelGradx(i)**2+
     *               vVelGrady(i)**2+wVelGradz(i)**2)+
     *                   (uVelGrady(i)+vVelGradx(i))**2+
     *                      (uVelGradz(i)+wVelGradx(i))**2+
     *                          (vVelGradz(i)+wVelGrady(i))**2
          FluxTElocal=-visc*Volume(i)*term1
     *                          
c
          FluxTE(i)=FluxTE(i)+FluxTElocal
c
        enddo
c
      else
c
        do i=1,NumberOfElements
c
          visc=Viscosity(i)
          term1=2.*(uVelGradx(i)**2+
     *               vVelGrady(i)**2+wVelGradz(i)**2)+
     *                   (uVelGrady(i)+vVelGradx(i))**2+
     *                      (uVelGradz(i)+wVelGradx(i))**2+
     *                          (vVelGradz(i)+wVelGrady(i))**2
          FluxTElocal=-visc*Volume(i)*term1
     *                          
c
          FluxTE(i)=FluxTE(i)+FluxTElocal
c
        enddo
c
      endif
c
      return
      end