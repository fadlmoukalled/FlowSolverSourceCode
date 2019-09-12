c
c#############################################################################################
c
      SUBROUTINE AssembleYapCorrection
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use WallDistance1, only: WallDistance
      use Variables1, only: TurbulentKE,TurbulentED
      use Variables3, only: FluxTE
      use PhysicalProperties1, only: Density
      use Turbulence1, only: cappa,cmu75,LTKE,ModelNumber
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: dNorm,Le,tke,ted,term1,Source,tedt,tke1
c********************************************************************************************
c
c-----------------------------------------------------------------------------
      CalculateYapCorrection: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(1) CalculateYapCorrection           !kepsilon
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(2) CalculateYapCorrection           !kepsilonchien
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i)+LTKE(i)*tke/Density(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(3) CalculateYapCorrection           !kepsilonsharma
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i)+LTKE(i)/Density(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(4) CalculateYapCorrection           !kepsilonchc
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(5) CalculateYapCorrection           !kepsilonkasagi
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(6) CalculateYapCorrection           !kepsilontagawa
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(7) CalculateYapCorrection           !kepsilonhishida
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i)+LTKE(i)/Density(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(8) CalculateYapCorrection           !kelambremhorst
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(9) CalculateYapCorrection           !kelambremhorstm
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            dNorm=WallDistance(i)
c
            Le=cappa*dnorm/cmu75
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i),tiny)
            term1=(tke**1.5)/(ted*Le)
c
            if(term1.gt.1.) then
c
              tedt=dmax1(TurbulentED(i),0.)
              tke1=dmax1(TurbulentKE(i),tiny)
              Source=0.83*Density(i)*((tedt**2)/tke1)*
     *                         (term1-1.)*(term1**2)*Volume(i)
              FluxTE(i)=FluxTE(i)-Source
c
            endif
c
          enddo 
c
c-----------------------------------------------------------------------------
        case(10) CalculateYapCorrection           !realizable
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(22) CalculateYapCorrection           !kepsilonrt
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(25) CalculateYapCorrection           !kepsilonrng
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(26) CalculateYapCorrection           !kepsilonv2f
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(27) CalculateYapCorrection           !kepsilonzetaf
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
      end select CalculateYapCorrection 
c-----------------------------------------------------------------------------
c
      return
      end