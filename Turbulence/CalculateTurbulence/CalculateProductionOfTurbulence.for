c
c#############################################################################################
c
      SUBROUTINE CalculateTurbulenceProduction
c
C#############################################################################################
c
      use User0, only: TurbulenceModel,
     *                 LimitTurbulenceProduction,
     *                 LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LSolveTurbulentKL
      use Geometry1, only: NumberOfElements
      use Variables1, only: TurbulentKE,TurbulentED,
     *                      TurbulentKL,TurbulentOmega,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      TurbulenceProduction
      use PhysicalProperties1, only: Density
      use Turbulence1, only: cmu,cmu75,c1limiter,ModelNumber,
     *                       Tau11,Tau12,Tau13,
     *                       Tau22,Tau23,Tau33
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term11,term12,term13,
     *                    term21,term22,term23,
     *                    term31,term32,term33,
     *                    tke1,ProductionTemp
c********************************************************************************************
c
      CalculatePk: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(1) CalculatePk           !kepsilon
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(2) CalculatePk           !kepsilonchien
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(3) CalculatePk           !kepsilonsharma
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(4) CalculatePk           !kepsilonchc
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(5) CalculatePk           !kepsilonkasagi
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(6) CalculatePk           !kepsilontagawa
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(7) CalculatePk           !kepsilonhishida
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(8) CalculatePk           !kelambremhorst
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(9) CalculatePk           !kelambremhorstm
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(10) CalculatePk           !realizable
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(11) CalculatePk           !komega
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*cmu*Density(i)*
     *                    dmax1(TurbulentKE(i)*TurbulentOmega(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(12) CalculatePk           !komegaepsilon
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*cmu*Density(i)*
     *                    dmax1(TurbulentKE(i)*TurbulentOmega(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(13) CalculatePk           !komegabsl
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*cmu*Density(i)*
     *                    dmax1(TurbulentKE(i)*TurbulentOmega(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(14) CalculatePk           !komegasst
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*cmu*Density(i)*
     *                    dmax1(TurbulentKE(i)*TurbulentOmega(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(15) CalculatePk           !sstgamaretheta
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*cmu*Density(i)*
     *                    dmax1(TurbulentKE(i)*TurbulentOmega(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(16) CalculatePk           !komega2006
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*cmu*Density(i)*
     *                    dmax1(TurbulentKE(i)*TurbulentOmega(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(17) CalculatePk           !komega2006lrn
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*cmu*Density(i)*
     *                    dmax1(TurbulentKE(i)*TurbulentOmega(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(18) CalculatePk           !kklmodel
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              tke1=dmax1(TurbulentKE(i),0.)
              ProductionTemp=c1limiter*cmu75*Density(i)*(tke1**2.5)/
     *                                       dmax1(TurbulentKL(i),tiny)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(21) CalculatePk           !kklomega
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(22) CalculatePk           !kepsilonrt
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(23) CalculatePk           !sstgama
c-----------------------------------------------------------------------------
c
          return
c-----------------------------------------------------------------------------
        case(25) CalculatePk           !kepsilonrng
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(26) CalculatePk           !kepsilonv2f
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(27) CalculatePk           !kepsilonzetaf
c-----------------------------------------------------------------------------
c
          call CalculatePkTerm
c
          if(LimitTurbulenceProduction) then
c
            do i=1,NumberOfElements    
c
              ProductionTemp=c1limiter*Density(i)*
     *                            dmax1(TurbulentED(i),0.)
c
              TurbulenceProduction(i)=
     *            dmin1(TurbulenceProduction(i),ProductionTemp)
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------------------------      
      end select CalculatePk 
c-----------------------------------------------------------------------------------------------      
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculatePkTerm
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Variables1, only: uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      TurbulenceProduction
      use PhysicalProperties1, only: Density
      use Turbulence1, only: cmu,cmu75,c1limiter,
     *                       Tau11,Tau12,Tau13,
     *                       Tau22,Tau23,Tau33
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term11,term12,term13,
     *                    term21,term22,term23,
     *                    term31,term32,term33
c********************************************************************************************
c
      call CalculateShearStress
c
      do i=1,NumberOfElements    
c
        term11=tau11(i)*uVelGradx(i)
        term12=tau12(i)*uVelGrady(i)
        term13=tau13(i)*uVelGradz(i)
        term21=tau12(i)*vVelGradx(i)
        term22=tau22(i)*vVelGrady(i)
        term23=tau23(i)*vVelGradz(i)
        term31=tau13(i)*wVelGradx(i)
        term32=tau23(i)*wVelGrady(i)
        term33=tau33(i)*wVelGradz(i)
c
        TurbulenceProduction(i)=dmax1(term11+term12+term13+
     *            term21+term22+term23+term31+term32+term33,0.)
c
      enddo
c
      return
      end
