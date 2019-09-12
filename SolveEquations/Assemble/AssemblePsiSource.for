c
C#############################################################################################
      SUBROUTINE AssemblePsiSourceTerm
C#############################################################################################
      use User0, only: LSolveTurbulenceKineticEnergy,LTurbulentFlow
      use Variables1, only: uVelGradx,vVelGrady,wVelGradz,TurbulentKE
      use Variables3, only: FluxTE
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use PhysicalProperties1, only:Viscosity,Density,TurbulentViscosity
      use Constants1, only: twothird
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: FluxTElocal,tke1,visc
c********************************************************************************************
c
      if(LTurbulentFlow) then
c
        do i=1,NumberOfElements
c
          visc=Viscosity(i)+TurbulentViscosity(i)
          FluxTElocal=-(-twothird)*visc*Volume(i)*
     *                   ((uVelGradx(i)+vVelGrady(i)+wVelGradz(i))**2)
c
          FluxTE(i)=FluxTE(i)+FluxTElocal
c
        enddo
c
        if(LSolveTurbulenceKineticEnergy) then
c
          do i=1,NumberOfElements    
c
            tke1=dmax1(TurbulentKE(i),0.)
c
            FluxTElocal=-(-twothird)*Density(i)*tke1*Volume(i)*
     *                      (uVelGradx(i)+vVelGrady(i)+wVelGradz(i))
            FluxTE(i)=FluxTE(i)+FluxTElocal
c          
          enddo
c
        endif
c
      else
c
        do i=1,NumberOfElements
c
          visc=Viscosity(i)
          FluxTElocal=-(-twothird)*visc*Volume(i)*
     *                   ((uVelGradx(i)+vVelGrady(i)+wVelGradz(i))**2)
c
          FluxTE(i)=FluxTE(i)+FluxTElocal
c
        enddo
c
      endif
c
      return
      end