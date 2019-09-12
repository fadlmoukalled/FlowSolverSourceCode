c
c#############################################################################################
c
      SUBROUTINE ModifiyTViscosityByTemperatureCorrection
c
C#############################################################################################
c
      use User0, only: urfTViscosity,LCompressibilityCorrection,
     *                 MaximumViscosityRatio
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: Temperature,TempGradx,TempGrady,TempGradz,
     *                      uVelocity,vVelocity,wVelocity,
     *                      uVelGradx,uVelGrady,uVelGradz, 
     *                      vVelGradx,vVelGrady,vVelGradz, 
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BTemperature,BTempGradx,BTempGrady,
     *                      BTempGradz,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      BuVelGradx,BuVelGrady,BuVelGradz, 
     *                      BvVelGradx,BvVelGrady,BvVelGradz, 
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
     *                      TurbulentKE,TurbulentED,TurbulentOmega,
     *                      BTurbulentKE,BTurbulentED,BTurbulentOmega,
     *                      ModifiedED,BModifiedED
      use PhysicalProperties1, only: SpecificHeat,Density,
     *                               BSpecificHeat,BDensity,RGas,
     *                               TurbulentViscosity,Viscosity,
     *                               BTurbulentViscosity,BViscosity
      use Turbulence1, only: sigKW,sigKE,f1WA,fmuWA,Bf1WA,BfmuWA,cmu,
     *                       StrainRate,BStrainRate,ModelNumber
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: Ttotal,termx,termy,termz,delTtotal,tke,ted,Tg,
     *                    ratio,a2,Mt2,FMt,term,TVisOld,sigR
      double precision, save :: Mt02=0.01
c********************************************************************************************
c
      CalculateCompressibleTemperatureCorrection: 
     *                                    select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(1) CalculateCompressibleTemperatureCorrection    !kepsilon
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            Ttotal=Temperature(i)+(uVelocity(i)**2+vVelocity(i)**2+
     *                          wVelocity(i)**2)/(2.*SpecificHeat(i))
            termx=TempGradx(i)+
     *         (uVelocity(i)*uVelGradx(i)+vVelocity(i)*vVelGradx(i)+
     *                       wVelocity(i)*wVelGradx(i))/SpecificHeat(i)
            termy=TempGrady(i)+
     *         (uVelocity(i)*uVelGrady(i)+vVelocity(i)*vVelGrady(i)+
     *                       wVelocity(i)*wVelGrady(i))/SpecificHeat(i)
            termz=TempGradz(i)+
     *         (uVelocity(i)*uVelGradz(i)+vVelocity(i)*vVelGradz(i)+
     *                       wVelocity(i)*wVelGradz(i))/SpecificHeat(i)
            delTtotal=dsqrt(termx**2+termy**2+termz**2)
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentED(i),tiny)
            Tg=(tke**1.5)*delTtotal/(ted*Ttotal)
c        
            ratio=SpecificHeat(i)/(SpecificHeat(i)-RGas)
            a2=ratio*RGas*Temperature(i)
            Mt2=2.*tke/a2    
            FMt=0.
            if(Mt2.gt.Mt02) FMt=Mt2-Mt02
            term=(1.+Tg**3/(0.041+FMt))
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         urfTViscosity*cmu*term*Density(i)*tke**2/ted
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Ttotal=BTemperature(i,j)+
     *               (BuVelocity(i,j)**2+BvVelocity(i,j)**2+
     *                    BwVelocity(i,j)**2)/(2.*BSpecificHeat(i,j))
              termx=BTempGradx(i,j)+
     *              (BuVelocity(i,j)*BuVelGradx(i,j)+
     *                  BvVelocity(i,j)*BvVelGradx(i,j)+
     *                    BwVelocity(i,j)*BwVelGradx(i,j))/
     *                                          BSpecificHeat(i,j)
              termy=BTempGrady(i,j)+
     *              (BuVelocity(i,j)*BuVelGrady(i,j)+
     *                  BvVelocity(i,j)*BvVelGrady(i,j)+
     *                    BwVelocity(i,j)*BwVelGrady(i,j))/
     *                                          BSpecificHeat(i,j)
              termz=BTempGradz(i,j)+
     *              (BuVelocity(i,j)*BuVelGradz(i,j)+
     *                  BvVelocity(i,j)*BvVelGradz(i,j)+
     *                    BwVelocity(i,j)*BwVelGradz(i,j))/
     *                                          BSpecificHeat(i,j)
              delTtotal=dsqrt(termx**2+termy**2+termz**2)
c
              tke=dmax1(BTurbulentKE(i,j),0.)
              ted=dmax1(BTurbulentED(i,j),tiny)
              Tg=(tke**1.5)*delTtotal/(ted*Ttotal)
c        
              ratio=BSpecificHeat(i,j)/(BSpecificHeat(i,j)-RGas)
              a2=ratio*RGas*BTemperature(i,j)
              Mt2=2.*tke/a2    
              FMt=0.
              if(Mt2.gt.Mt02) FMt=Mt2-Mt02
              term=(1.+Tg**3/(0.041+FMt))
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *         urfTViscosity*cmu*term*BDensity(i,j)*tke**2/ted
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------
        case(11:14,16:17) CalculateCompressibleTemperatureCorrection !k-omega based
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            Ttotal=Temperature(i)+(uVelocity(i)**2+vVelocity(i)**2+
     *                            wVelocity(i)**2)/(2.*SpecificHeat(i))
            termx=TempGradx(i)+
     *         (uVelocity(i)*uVelGradx(i)+vVelocity(i)*vVelGradx(i)+
     *                       wVelocity(i)*wVelGradx(i))/SpecificHeat(i)
            termy=TempGrady(i)+
     *         (uVelocity(i)*uVelGrady(i)+vVelocity(i)*vVelGrady(i)+
     *                       wVelocity(i)*wVelGrady(i))/SpecificHeat(i)
            termz=TempGradz(i)+
     *         (uVelocity(i)*uVelGradz(i)+vVelocity(i)*vVelGradz(i)+
     *                       wVelocity(i)*wVelGradz(i))/SpecificHeat(i)
            delTtotal=dsqrt(termx**2+termy**2+termz**2)
c
            tke=dmax1(TurbulentKE(i),0.)
            ted=dmax1(TurbulentOmega(i),tiny)
            Tg=delTtotal*dsqrt(tke)/(ted*Ttotal)
c        
            ratio=SpecificHeat(i)/(SpecificHeat(i)-RGas)
            a2=ratio*RGas*Temperature(i)
            Mt2=2.*tke/a2    
            FMt=0.
            if(Mt2.gt.Mt02) FMt=Mt2-Mt02
            if(LCompressibilityCorrection) then
              term=cmu*(1.+340.*Tg**3)
            else
              term=cmu*(1.+13.5*Tg**3/(0.039+FMt))
            endif
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *                            urfTViscosity*Density(i)*term*tke/ted
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Ttotal=BTemperature(i,j)+
     *               (BuVelocity(i,j)**2+BvVelocity(i,j)**2+
     *                    BwVelocity(i,j)**2)/(2.*BSpecificHeat(i,j))
              termx=BTempGradx(i,j)+
     *              (BuVelocity(i,j)*BuVelGradx(i,j)+
     *                  BvVelocity(i,j)*BvVelGradx(i,j)+
     *                    BwVelocity(i,j)*BwVelGradx(i,j))/
     *                                          BSpecificHeat(i,j)
              termy=BTempGrady(i,j)+
     *              (BuVelocity(i,j)*BuVelGrady(i,j)+
     *                  BvVelocity(i,j)*BvVelGrady(i,j)+
     *                    BwVelocity(i,j)*BwVelGrady(i,j))/
     *                                          BSpecificHeat(i,j)
              termz=BTempGradz(i,j)+
     *              (BuVelocity(i,j)*BuVelGradz(i,j)+
     *                  BvVelocity(i,j)*BvVelGradz(i,j)+
     *                    BwVelocity(i,j)*BwVelGradz(i,j))/
     *                                          BSpecificHeat(i,j)
              delTtotal=dsqrt(termx**2+termy**2+termz**2)
c
              tke=dmax1(BTurbulentKE(i,j),0.)
              ted=dmax1(BTurbulentOmega(i,j),tiny)
              Tg=delTtotal*dsqrt(tke)/(ted*Ttotal)
c        
              ratio=BSpecificHeat(i,j)/(BSpecificHeat(i,j)-RGas)
              a2=ratio*RGas*BTemperature(i,j)
              Mt2=2.*tke/a2    
              FMt=0.
              if(Mt2.gt.Mt02) FMt=Mt2-Mt02
              if(LCompressibilityCorrection) then
                term=cmu*(1.+340.*Tg**3)
              else
                term=cmu*(1.+13.5*Tg**3/(0.039+FMt))
              endif
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *                            urfTViscosity*Density(i)*term*tke/ted
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------
        case(20) CalculateCompressibleTemperatureCorrection !Wray-Agarwal
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            Ttotal=Temperature(i)+(uVelocity(i)**2+vVelocity(i)**2+
     *                          wVelocity(i)**2)/(2.*SpecificHeat(i))
            termx=TempGradx(i)+
     *        (uVelocity(i)*uVelGradx(i)+vVelocity(i)*vVelGradx(i)+
     *                      wVelocity(i)*wVelGradx(i))/SpecificHeat(i)
            termy=TempGrady(i)+
     *         (uVelocity(i)*uVelGrady(i)+vVelocity(i)*vVelGrady(i)+
     *                      wVelocity(i)*wVelGrady(i))/SpecificHeat(i)
            termz=TempGradz(i)+
     *         (uVelocity(i)*uVelGradz(i)+vVelocity(i)*vVelGradz(i)+
     *                      wVelocity(i)*wVelGradz(i))/SpecificHeat(i)
            delTtotal=dsqrt(termx**2+termy**2+termz**2)
c
            sigR=f1WA(i)*(sigKW-sigKE)+sigKE
            Tg=dsqrt(sigR*ModifiedED(i)/StrainRate(i))*delTtotal/Ttotal
            term=1.+18.*Tg**3
c        
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *             urfTViscosity*Density(i)*term*fmuWA(i)*ModifiedED(i)
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Ttotal=BTemperature(i,j)+
     *               (BuVelocity(i,j)**2+BvVelocity(i,j)**2+
     *                    BwVelocity(i,j)**2)/(2.*BSpecificHeat(i,j))
              termx=BTempGradx(i,j)+
     *              (BuVelocity(i,j)*BuVelGradx(i,j)+
     *                  BvVelocity(i,j)*BvVelGradx(i,j)+
     *                    BwVelocity(i,j)*BwVelGradx(i,j))/
     *                                          BSpecificHeat(i,j)
              termy=BTempGrady(i,j)+
     *              (BuVelocity(i,j)*BuVelGrady(i,j)+
     *                  BvVelocity(i,j)*BvVelGrady(i,j)+
     *                    BwVelocity(i,j)*BwVelGrady(i,j))/
     *                                          BSpecificHeat(i,j)
              termz=BTempGradz(i,j)+
     *              (BuVelocity(i,j)*BuVelGradz(i,j)+
     *                  BvVelocity(i,j)*BvVelGradz(i,j)+
     *                    BwVelocity(i,j)*BwVelGradz(i,j))/
     *                                          BSpecificHeat(i,j)
              delTtotal=dsqrt(termx**2+termy**2+termz**2)
c
              sigR=Bf1WA(i,j)*(sigKW-sigKE)+sigKE
              Tg=dsqrt(sigR*BModifiedED(i,j)/
     *                              BStrainRate(i,j))*delTtotal/Ttotal
              term=1.+18.*Tg**3
c        
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *                             urfTViscosity*BDensity(i,j)*term*
     *                                     BfmuWA(i,j)*BModifiedED(i,j)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateCompressibleTemperatureCorrection 
c-----------------------------------------------------------------------------------------------      
c
c      call ApplyBConTurbulentViscosity ! to apply it I should not calculate mut at the boundary
c
c---- Update turbulent viscosity at inlet
c
      call UpdateInletTurbulentViscosity
c
c--- Limit the turbulent viscosity
c
      do i=1,NumberOfElements
c
        TurbulentViscosity(i)=dmin1(TurbulentViscosity(i),
     *                       MaximumViscosityRatio*Viscosity(i))
c
      enddo 
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentViscosity(i,j)=dmin1(BTurbulentViscosity(i,j),
     *         MaximumViscosityRatio*BViscosity(i,j))
c
        enddo 
      enddo 
c
      return
      end