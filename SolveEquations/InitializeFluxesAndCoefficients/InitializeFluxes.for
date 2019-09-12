c
C#############################################################################################
c
      SUBROUTINE InitializeFluxes
c
C#############################################################################################
c
      use User0, only: LUnsteady
      use Variables2
      use Variables3
c********************************************************************************************
      implicit none
c********************************************************************************************
c
c---  Initialize face fluxes
c
      FluxCf=0.
      FluxFf=0.
      FluxVf=0.
      FluxTf=0.
c
c---  Initialize element fluxes
c
      FluxCE=0.
      FluxVE=0.
      FluxTE=0.
c
      if(LUnsteady) then
        FluxCEold=0.
        FluxCEoldold=0.
      endif
c
      return
      end
