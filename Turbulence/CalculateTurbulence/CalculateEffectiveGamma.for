c
c#############################################################################################
c
      SUBROUTINE CalculateTGammaEff
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Variables1, only: TurbulentKE,TurbulentOmega,TGamma,TReTheta,
     *                      uVelocity,vVelocity,wVelocity,TGammaEff
      use PhysicalProperties1, only: Density,Viscosity
      use Turbulence1, only: s1,ce2,Vorticity,StrainRate
      use Constants1, only: tiny
      use WallDistance1, only: WallDistance
c********************************************************************************************
      Implicit none
c********************************************************************************************
      integer :: i
      double precision :: ReV,RT,ReW,Fwake,Uvel2,delta,Fthetat,
     *                    FReattach,GammaSep
c********************************************************************************************
      Interface
c********************************************************************************************
        FUNCTION ReThetaC(x)
c*********************************************************************************************
          double precision :: x
          double precision :: ReThetaC
c*********************************************************************************************
        end FUNCTION ReThetaC
c*********************************************************************************************
      end interface
c*********************************************************************************************
c
      do i=1,NumberOfElements
c      
        ReV=Density(i)*StrainRate(i)*WallDistance(i)*
     *                WallDistance(i)/Viscosity(i)
        RT=Density(i)*dmax1(TurbulentKE(i),0.)/
     *            (Viscosity(i)*dmax1(TurbulentOmega(i),tiny))
        ReW=Density(i)*dmax1(TurbulentOmega(i),0.)*
     *         WallDistance(i)*WallDistance(i)/Viscosity(i)
        Fwake=dexp(-ReW*ReW/1.d10)
        Uvel2=uVelocity(i)*uVelocity(i)+vVelocity(i)*
     *                    vVelocity(i)+wVelocity(i)*wVelocity(i)
        Uvel2=dmax1(Uvel2,tiny)
        delta=375.*Vorticity(i)*Viscosity(i)*TReTheta(i)*
     *        WallDistance(i)/(Density(i)*Uvel2)
        if(delta.eq.0.) delta=tiny
        Fthetat=dmin1(dmax1(
     *        Fwake*dexp(-(wallDistance(i)/delta)**4),
     *                         1.-((ce2*TGamma(i)-1.)/(ce2-1.))**2),1.)
        FReattach=dexp(-RT*RT*RT*RT/160000.)
        GammaSep=dmin1(s1*dmax1(0.,ReV/(3.235*
     *     dmax1(ReThetaC(TReTheta(i)),tiny))-1.)*FReattach,2.)*Fthetat
        TGammaEff(i)=dmax1(TGamma(i),GammaSep)     
c       
      enddo
c
      return
      end