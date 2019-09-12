c
c#############################################################################################
c
      SUBROUTINE calculateShearStress
c
C#############################################################################################
c
      use User0, only: Lcompressible,LSolveTurbulenceKineticEnergy,
     *                 TurbulenceModel
      use Variables1, only: TurbulentKE,TurbulentKL,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz
      use PhysicalProperties1, only: TurbulentViscosity,
     *                               Density,Viscosity
      use Turbulence1, only: S11,S12,S13,S22,S23,S33,
     *                       Tau11,Tau12,Tau13,Tau22,Tau23,Tau33
      use Constants1, only: tiny,twothird
      use Geometry1, only: NumberOfElements
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1,tke1,effectiveViscosity
c********************************************************************************************
c
      if(TurbulenceModel.eq.'kepsilon'.or.
     *       TurbulenceModel.eq.'realizable'.or.
     *           TurbulenceModel.eq.'kepsilonrng') then
c      
        do i=1,NumberOfElements    
c
          effectiveViscosity=TurbulentViscosity(i)+Viscosity(i)
c
          Tau11(i)=effectiveViscosity*2.*S11(i)
          Tau12(i)=effectiveViscosity*2.*S12(i)
          Tau13(i)=effectiveViscosity*2.*S13(i)
          Tau22(i)=effectiveViscosity*2.*S22(i)
          Tau23(i)=effectiveViscosity*2.*S23(i)
          Tau33(i)=effectiveViscosity*2.*S33(i)
c          
        enddo
c
        do i=1,NumberOfElements    
c
          tke1=dmax1(TurbulentKE(i),0.)
c
          Tau11(i)=Tau11(i)-twothird*Density(i)*tke1
          Tau22(i)=Tau22(i)-twothird*Density(i)*tke1
          Tau33(i)=Tau33(i)-twothird*Density(i)*tke1
c          
        enddo
c
        if(Lcompressible) then
c
          do i=1,NumberOfElements    
c
            term1=uVelGradx(i)+vVelGrady(i)+wVelGradz(i)
            effectiveViscosity=TurbulentViscosity(i)+Viscosity(i)
c
            Tau11(i)=Tau11(i)-twothird*effectiveViscosity*term1
            Tau22(i)=Tau22(i)-twothird*effectiveViscosity*term1
            Tau33(i)=Tau33(i)-twothird*effectiveViscosity*term1
c          
          enddo
c
        endif
c
      else
c
        do i=1,NumberOfElements    
c
          Tau11(i)=TurbulentViscosity(i)*2.*S11(i)
          Tau12(i)=TurbulentViscosity(i)*2.*S12(i)
          Tau13(i)=TurbulentViscosity(i)*2.*S13(i)
          Tau22(i)=TurbulentViscosity(i)*2.*S22(i)
          Tau23(i)=TurbulentViscosity(i)*2.*S23(i)
          Tau33(i)=TurbulentViscosity(i)*2.*S33(i)
c          
        enddo
c
        if(LSolveTurbulenceKineticEnergy) then
c
          if(TurbulenceModel.eq.'kklomega') then
c
            do i=1,NumberOfElements    
c
              tke1=dmax1(TurbulentKE(i)+TurbulentKL(i),0.)
c
              Tau11(i)=Tau11(i)-twothird*Density(i)*tke1
              Tau22(i)=Tau22(i)-twothird*Density(i)*tke1
              Tau33(i)=Tau33(i)-twothird*Density(i)*tke1
c          
            enddo
c
          else
c
            do i=1,NumberOfElements    
c
              tke1=dmax1(TurbulentKE(i),0.)
c
              Tau11(i)=Tau11(i)-twothird*Density(i)*tke1
              Tau22(i)=Tau22(i)-twothird*Density(i)*tke1
              Tau33(i)=Tau33(i)-twothird*Density(i)*tke1
c          
            enddo
c
          endif
c
        endif
c
        if(Lcompressible) then
c
          do i=1,NumberOfElements    
c
            term1=uVelGradx(i)+vVelGrady(i)+wVelGradz(i)
c
            Tau11(i)=Tau11(i)-twothird*TurbulentViscosity(i)*term1
            Tau22(i)=Tau22(i)-twothird*TurbulentViscosity(i)*term1
            Tau33(i)=Tau33(i)-twothird*TurbulentViscosity(i)*term1
c          
          enddo
c
        endif
c
      endif
c
      return
      end
