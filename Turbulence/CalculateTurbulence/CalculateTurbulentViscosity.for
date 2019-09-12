c
c#############################################################################################
c
      SUBROUTINE CalculateTurbulentViscosity
c
C#############################################################################################
      use User0, only: TurbulenceModel,urfTViscosity,WallTreatment,
     *                 MaximumViscosityRatio,LNegativeSpalartAllmaras
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use PhysicalProperties1, only: TurbulentViscosity,
     *                               BTurbulentViscosity,
     *                               Density,BDensity,
     *                               Viscosity,BViscosity
      use Variables1, only: TurbulentKE,BTurbulentKE,
     *                      TurbulentED,BTurbulentED,
     *                      TurbulentOmega,BTurbulentOmega,
     *                      uVelGradx,vVelGradx,wVelGradx,
     *                      BuVelGradx,BvVelGradx,BwVelGradx,
     *                      uVelGrady,vVelGrady,wVelGrady,
     *                      BuVelGrady,BvVelGrady,BwVelGrady,
     *                      uVelGradz,vVelGradz,wVelGradz,
     *                      BuVelGradz,BvVelGradz,BwVelGradz,
     *                      ModifiedED,BModifiedED,
     *                      TurbulentKL,BTurbulentKL,
     *                      TurbulentV2,BTurbulentV2,
     *                      TurbulentZeta,BTurbulentZeta
      use Turbulence1, only: cmu,a1sst,F2factor,BF2factor,cmu25,
     *                       F3factor,BF3factor,Clim,
     *                       S11,S12,S13,S22,S23,S33,
     *                       BS11,BS12,BS13,BS22,BS23,BS33,
     *                       StrainRate,BStrainRate,
     *                       alfaStar,BalfaStar,fmuCoefficient,
     *                       BfmuCoefficient,fv1Coefficient,
     *                       Bfv1Coefficient,fmuWA,BfmuWA,
     *                       cmuR,BcmuR,ModelNumber,
     *                       TurbulentViscosityTs,TurbulentViscosityTl,
     *                       BTurbulentViscosityTs,
     *                       BTurbulentViscosityTl,
     *                       TurbulentViscosityKE,
     *                       BTurbulentViscosityKE,
     *                       TurbulentViscosityRT,
     *                       BTurbulentViscosityRT,
     *                       TScale,BTScale
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: TVisOld,term1,term2,term3,term4,term5,sum,
     *                    Smagnitude,OmegaTilda,
     *                    S11bar,S12bar,S13bar,S21bar,S22bar,S23bar,
     *                    S31bar,S32bar,S33bar,tke1,tkl1
c********************************************************************************************
c
c----  Select case number 
c          
!      if(TurbulenceModel.eq.'kepsilon') then
!c
!        ModelNumber=1
!c
!      elseif(TurbulenceModel.eq.'kepsilonchien') then
!c
!        ModelNumber=2
!c
!      elseif(TurbulenceModel.eq.'kepsilonsharma') then
!c
!        ModelNumber=3
!c
!      elseif(TurbulenceModel.eq.'kepsilonchc') then
!c
!        ModelNumber=4
!c
!      elseif(TurbulenceModel.eq.'kepsilonkasagi') then
!c
!        ModelNumber=5
!c
!      elseif(TurbulenceModel.eq.'kepsilontagawa') then
!c
!        ModelNumber=6
!c
!      elseif(TurbulenceModel.eq.'kepsilonhishida') then
!c
!        ModelNumber=7
!c
!      elseif(TurbulenceModel.eq.'kelambremhorst') then
!c
!        ModelNumber=8
!c
!      elseif(TurbulenceModel.eq.'kelambremhorstm') then
!c
!        ModelNumber=9
!c
!      elseif(TurbulenceModel.eq.'realizable') then
!c
!        ModelNumber=10
!c
!      elseif(TurbulenceModel.eq.'komega') then
!c
!        ModelNumber=11
!c
!      elseif(TurbulenceModel.eq.'komegaepsilon') then
!c
!        ModelNumber=12
!c
!      elseif(TurbulenceModel.eq.'komegabsl') then
!c
!        ModelNumber=13
!c
!      elseif(TurbulenceModel.eq.'komegasst') then
!c
!        ModelNumber=14
!c
!      elseif(TurbulenceModel.eq.'sstgamaretheta') then
!c
!        ModelNumber=15
!c
!      elseif(TurbulenceModel.eq.'komega2006') then
!c
!        ModelNumber=16
!c
!      elseif(TurbulenceModel.eq.'komega2006lrn') then
!c
!        ModelNumber=17
!c
!      elseif(TurbulenceModel.eq.'kklmodel') then
!c
!        ModelNumber=18
!c
!      elseif(TurbulenceModel.eq.'spalartallmaras') then
!c
!        ModelNumber=19
!c
!      elseif(TurbulenceModel.eq.'wrayagarwal') then
!c
!        ModelNumber=20
!c
!      elseif(TurbulenceModel.eq.'kklomega') then
!c
!        ModelNumber=21
!c
!      elseif(TurbulenceModel.eq.'kepsilonrt') then
!c
!        ModelNumber=22
!c
!      elseif(TurbulenceModel.eq.'sstgama') then
!c
!        ModelNumber=23
!c
!      elseif(TurbulenceModel.eq.'nut92') then
!c
!        ModelNumber=24
!c
!      elseif(TurbulenceModel.eq.'kepsilonrng') then
!c
!        ModelNumber=25
!c
!      elseif(TurbulenceModel.eq.'kepsilonv2f') then
!c
!        ModelNumber=26
!c
!      elseif(TurbulenceModel.eq.'kepsilonzetaf') then
!c
!        ModelNumber=27
!c
!      endif      
c
c-----------------------------------------------------------------------------------------------      
      CalculateMuT:select case (ModelNumber)
c-----------------------------------------------------------------------------------------------      
        case(1) CalculateMuT           !kepsilon
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*urfTViscosity*Density(i)*TurbulentKE(i)**2/
     *                               dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *          cmu*urfTViscosity*BDensity(i,j)*BTurbulentKE(i,j)**2/
     *                                 dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(2) CalculateMuT           !kepsilonchien
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(3) CalculateMuT           !kepsilonsharma
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(4) CalculateMuT           !kepsilonchc
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(5) CalculateMuT           !kepsilonkasagi
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(6) CalculateMuT           !kepsilontagawa
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(7) CalculateMuT           !kepsilonhishida
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(8) CalculateMuT           !kelambremhorst
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(9) CalculateMuT           !kelambremhorstm
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*fmuCoefficient(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            cmu*BfmuCoefficient(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(10) CalculateMuT           !realizable
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmuR(i)*urfTViscosity*Density(i)*
     *                    TurbulentKE(i)**2/dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *            BcmuR(i,j)*urfTViscosity*BDensity(i,j)*
     *                BTurbulentKE(i,j)**2/dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(11) CalculateMuT           !komega
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *            urfTViscosity*Density(i)*dmax1(TurbulentKE(i),0.)/
     *                                    dmax1(TurbulentOmega(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *          urfTViscosity*BDensity(i,j)*dmax1(BTurbulentKE(i,j),0.)/
     *                                 dmax1(BTurbulentOmega(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(12) CalculateMuT           !komegaepsilon
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *            urfTViscosity*Density(i)*dmax1(TurbulentKE(i),0.)/
     *                                    dmax1(TurbulentOmega(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *          urfTViscosity*BDensity(i,j)*dmax1(BTurbulentKE(i,j),0.)/
     *                                  dmax1(BTurbulentOmega(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(13) CalculateMuT           !komegabsl
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *            urfTViscosity*Density(i)*dmax1(TurbulentKE(i),0.)/
     *                                    dmax1(TurbulentOmega(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*BDensity(i,j)*dmax1(BTurbulentKE(i,j),0.)/
     *                                 dmax1(BTurbulentOmega(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(14) CalculateMuT           !komegasst
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            Smagnitude=StrainRate(i)
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*Density(i)*a1sst*dmax1(TurbulentKE(i),0.)/
     *           dmax1(a1sst*TurbulentOmega(i),
     *                 Smagnitude*F2factor(i)*F3Factor(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Smagnitude=BStrainRate(i,j)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *           urfTViscosity*BDensity(i,j)*a1sst*
     *              dmax1(BTurbulentKE(i,j),0.)/dmax1(a1sst*
     *               BTurbulentOmega(i,j),Smagnitude*BF2factor(i,j)*
     *                                            BF3factor(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(15) CalculateMuT           !sstgamaretheta
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            Smagnitude=StrainRate(i)
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*Density(i)*a1sst*dmax1(TurbulentKE(i),0.)/
     *           dmax1(a1sst*TurbulentOmega(i),
     *                 Smagnitude*F2factor(i)*F3Factor(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Smagnitude=BStrainRate(i,j)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *           urfTViscosity*BDensity(i,j)*a1sst*
     *              dmax1(BTurbulentKE(i,j),0.)/dmax1(a1sst*
     *               BTurbulentOmega(i,j),Smagnitude*BF2factor(i,j)*
     *                                            BF3factor(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(16) CalculateMuT           !komega2006
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            S11bar=S11(i)-(uvelGradx(i)+vVelGrady(i)+wVelGradz(i))/3.
            S12bar=S12(i)
            S13bar=S13(i)
            S21bar=S12(i)
            S22bar=S22(i)-(uvelGradx(i)+vVelGrady(i)+wVelGradz(i))/3.
            S23bar=S23(i)
            S31bar=S13(i)
            S32bar=S23(i)
            S33bar=S33(i)-(uvelGradx(i)+vVelGrady(i)+wVelGradz(i))/3.
            sum=2.*(S11bar*S11bar+S12bar*S12bar+S13bar*S13bar+
     *            S21bar*S21bar+S22bar*S22bar+S23bar*S23bar+
     *            S31bar*S31bar+S32bar*S32bar+S33bar*S33bar)
            Smagnitude=dsqrt(sum/cmu)
            OmegaTilda=dmax1(TurbulentOmega(i),Clim*Smagnitude,tiny)
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *                                urfTViscosity*Density(i)*
     *                          dmax1(TurbulentKE(i),0.)/OmegaTilda
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              S11bar=BS11(i,j)-(BuvelGradx(i,j)+BvVelGrady(i,j)+
     *                                          BwVelGradz(i,j))/3.
              S12bar=BS12(i,j)
              S13bar=BS13(i,j)
              S21bar=BS12(i,j)
              S22bar=BS22(i,j)-(BuvelGradx(i,j)+BvVelGrady(i,j)+
     *                                          BwVelGradz(i,j))/3.
              S23bar=BS23(i,j)
              S31bar=BS13(i,j)
              S32bar=BS23(i,j)
              S33bar=BS33(i,j)-(BuvelGradx(i,j)+BvVelGrady(i,j)+
     *                                          BwVelGradz(i,j))/3.
              sum=2.*(S11bar*S11bar+S12bar*S12bar+S13bar*S13bar+
     *              S21bar*S21bar+S22bar*S22bar+S23bar*S23bar+
     *              S31bar*S31bar+S32bar*S32bar+S33bar*S33bar)
              Smagnitude=dsqrt(sum/cmu)
              OmegaTilda=
     *            dmax1(BTurbulentOmega(i,j),Clim*Smagnitude,tiny)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *             urfTViscosity*BDensity(i,j)*
     *                dmax1(BTurbulentKE(i,j),0.)/OmegaTilda
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(17) CalculateMuT           !komega2006lrn
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            S11bar=S11(i)-(uvelGradx(i)+vVelGrady(i)+wVelGradz(i))/3.
            S12bar=S12(i)
            S13bar=S13(i)
            S21bar=S12(i)
            S22bar=S22(i)-(uvelGradx(i)+vVelGrady(i)+wVelGradz(i))/3.
            S23bar=S23(i)
            S31bar=S13(i)
            S32bar=S23(i)
            S33bar=S33(i)-(uvelGradx(i)+vVelGrady(i)+wVelGradz(i))/3.
            sum=2.*(S11bar*S11bar+S12bar*S12bar+S13bar*S13bar+
     *            S21bar*S21bar+S22bar*S22bar+S23bar*S23bar+
     *            S31bar*S31bar+S32bar*S32bar+S33bar*S33bar)
            Smagnitude=dsqrt(sum/(cmu/alfaStar(i)))
            OmegaTilda=dmax1(TurbulentOmega(i),Clim*Smagnitude,tiny)
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *                     urfTViscosity*alfaStar(i)*Density(i)*
     *                          dmax1(TurbulentKE(i),0.)/OmegaTilda
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              S11bar=BS11(i,j)-(BuvelGradx(i,j)+BvVelGrady(i,j)+
     *                                          BwVelGradz(i,j))/3.
              S12bar=BS12(i,j)
              S13bar=BS13(i,j)
              S21bar=BS12(i,j)
              S22bar=BS22(i,j)-(BuvelGradx(i,j)+BvVelGrady(i,j)+
     *                                          BwVelGradz(i,j))/3.
              S23bar=BS23(i,j)
              S31bar=BS13(i,j)
              S32bar=BS23(i,j)
              S33bar=BS33(i,j)-(BuvelGradx(i,j)+BvVelGrady(i,j)+
     *                                          BwVelGradz(i,j))/3.
              sum=2.*(S11bar*S11bar+S12bar*S12bar+S13bar*S13bar+
     *              S21bar*S21bar+S22bar*S22bar+S23bar*S23bar+
     *              S31bar*S31bar+S32bar*S32bar+S33bar*S33bar)
              Smagnitude=dsqrt(sum/(cmu/BalfaStar(i,j)))
              OmegaTilda=
     *            dmax1(BTurbulentOmega(i,j),Clim*Smagnitude,tiny)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *             urfTViscosity*BDensity(i,j)*BalfaStar(i,j)*
     *                  dmax1(BTurbulentKE(i,j),0.)/OmegaTilda
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(18) CalculateMuT           !kklmodel
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            tke1=dmax1(turbulentKE(i),tiny)
            tkl1=dmax1(turbulentKL(i),0.)
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*cmu25*Density(i)*tkl1/dsqrt(tke1)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              tke1=dmax1(BturbulentKE(i,j),tiny)
              tkl1=dmax1(BTurbulentKL(i,j),0.)
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *               urfTViscosity*cmu25*BDensity(i,j)*tkl1/dsqrt(tke1)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(19) CalculateMuT           !spalartallmaras
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*fv1Coefficient(i)*Density(i)*ModifiedED(i)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *                   urfTViscosity*Bfv1Coefficient(i,j)*
     *                                 BDensity(i,j)*BModifiedED(i,j)
c
            enddo 
          enddo 
c
          if(LNegativeSpalartAllmaras) then
c
            do i=1,NumberOfElements
c
              if(ModifiedED(i).lt.0.)  TurbulentViscosity(i)=0.
c
            enddo 
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                if(BModifiedED(i,j).lt.0.)  BTurbulentViscosity(i,j)=0.
c
              enddo 
            enddo 
c
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(20) CalculateMuT           !wrayagarwal
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*fmuWA(i)*Density(i)*ModifiedED(i)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *                   urfTViscosity*BfmuWA(i,j)*
     *                                 BDensity(i,j)*BModifiedED(i,j)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(21) CalculateMuT           !Transitional k-kl-w
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=
     *         (1.-urfTViscosity)*TVisOld+Density(i)*urfTViscosity*
     *                (TurbulentViscosityTs(i)+TurbulentViscosityTl(i))
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=
     *          (1.-urfTViscosity)*TVisOld+BDensity(i,j)*urfTViscosity*
     *           (BTurbulentViscosityTs(i,j)+BTurbulentViscosityTl(i,j))
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(22) CalculateMuT           !k epsilon Rt
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *            urfTViscosity*dmax1(TurbulentViscosityKE(i),
     *                                  TurbulentViscosityRT(i))
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *              urfTViscosity*dmax1(BTurbulentViscosityKE(i,j),
     *                                     BTurbulentViscosityRT(i,j))
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(23) CalculateMuT           !sstgama
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            Smagnitude=StrainRate(i)
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*Density(i)*a1sst*dmax1(TurbulentKE(i),0.)/
     *           dmax1(a1sst*TurbulentOmega(i),
     *                 Smagnitude*F2factor(i)*F3Factor(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Smagnitude=BStrainRate(i,j)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *           urfTViscosity*BDensity(i,j)*a1sst*
     *              dmax1(BTurbulentKE(i,j),0.)/dmax1(a1sst*
     *               BTurbulentOmega(i,j),Smagnitude*BF2factor(i,j)*
     *                                            BF3factor(i,j),tiny)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(24) CalculateMuT           !nut92
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *                           urfTViscosity*Density(i)*ModifiedED(i)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *                     urfTViscosity*BDensity(i,j)*BModifiedED(i,j)
c
            enddo 
          enddo
c
c-----------------------------------------------------------------------------------------------      
        case(25) CalculateMuT           !kepsilonrng
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         cmu*urfTViscosity*Density(i)*TurbulentKE(i)**2/
     *                               dmax1(TurbulentED(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *          cmu*urfTViscosity*BDensity(i,j)*BTurbulentKE(i,j)**2/
     *                                 dmax1(BTurbulentED(i,j),tiny)
c
            enddo 
          enddo 
c
          call calculatesigRNG
c
c-----------------------------------------------------------------------------------------------      
        case(26) CalculateMuT           !kepsilonv2f
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *         urfTViscosity*cmu*Density(i)*TurbulentV2(i)*TScale(i)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *              urfTViscosity*cmu*BDensity(i,j)*
     *                               BTurbulentV2(i,j)*BTScale(i,j)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
        case(27) CalculateMuT           !kepsilonzetaf
c-----------------------------------------------------------------------------------------------      
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity(i)
            TurbulentViscosity(i)=(1.-urfTViscosity)*TVisOld+
     *             urfTViscosity*cmu*Density(i)*TurbulentZeta(i)*
     *                                      TurbulentKE(i)*TScale(i)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *              urfTViscosity*cmu*BDensity(i,j)*
     *              BTurbulentZeta(i,j)*BTurbulentKE(i,j)*BTScale(i,j)
c
            enddo 
          enddo 
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateMuT 
c
c--- Apply boundary conditions on Turbulent Viscosity
c
c      call ApplyBConTurbulentViscosity  ! to apply it I should not calculate mut at the boundary
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
c
c#############################################################################################
c
      SUBROUTINE ApplyBConTurbulentViscosity  !currently not used
c
c#############################################################################################
c
      use User0, only: WallTreatment
      use BoundaryConditions2, only: LRotationalPeriodicity,
     *                               LTranslationalPeriodicity
      use BoundaryConditionsTurbulence2, only: IoutletTurbulence,
     *                                         IoutletTurbulenceOwner,
     *                                 IoutletTurbulenceNumberOfBCSets,
     *                                 IoutletTurbulenceNBFaces
      use PhysicalProperties1, only: TurbulentViscosity,
     *                               BTurbulentViscosity
      use BoundaryConditions2, only: Isymmetry,IsymmetryOwner,
     *                               IsymmetryNumberOfBCSets,
     *                               IsymmetryNBFaces,Iperiodic,
     *                               IperiodicOwner,
     *                               IperiodicNumberOfBCSets,
     *                               IperiodicNBFaces,PeriodicPair,
     *                               Icorrespondingface,a1r,b1r,c1r,
     *                               a2r,b2r,c2r,a3r,b3r,c3r,
     *                               xTranslation,yTranslation,
     *                               zTranslation,Iaxis,IaxisOwner,
     *                               IaxisNumberOfBCSets,IaxisNBFaces,
     *                               IpressureFarField,
     *                               IpressureFarFieldOwner,
     *                               IpressureFarFieldNumberOfBCSets,
     *                               IpressureFarFieldNBFaces,IWallSlip,
     *                               IWallSlipOwner,IWallSlipNBFaces,
     *                               IWallSlipNumberOfBCSets,
     *                               IWallnoSlip,IWallnoSlipOwner,
     *                               IWallnoSlipNumberOfBCSets,
     *                               IWallnoSlipNBFaces
      use Turbulence1, only: ModelNumber
      use Variables1, only: Bmdot
      use Geometry3, only: NBFaceOwner                   
      use Geometry4, only: xc,yc,zc,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,j1,j2,j3
      double precision :: xF1,yF1,zF1,distance1,distance2,GFactCF
c********************************************************************************************
c
      do i=1,IoutletTurbulence
c
        i1=IoutletTurbulenceOwner(i)
        i2=IoutletTurbulenceNumberOfBCSets(i)
        i3=IoutletTurbulenceNBFaces(i)
c
	  BTurbulentViscosity(i2,i3)=TurbulentViscosity(i1)
c
      enddo
c
      do i=1,Isymmetry
c
        i1=IsymmetryOwner(i)
        i2=IsymmetryNumberOfBCSets(i)
        i3=IsymmetryNBFaces(i)
c
        BTurbulentViscosity(i2,i3)=TurbulentViscosity(i1)
c
      enddo
c
      do i=1,Iperiodic
c
        i1=IperiodicOwner(i)
        i2=IperiodicNumberOfBCSets(i)
        i3=IperiodicNBFaces(i)
c
        j2=PeriodicPair(i2)         
        j3=Icorrespondingface(i2,i3)
        j1=NBFaceOwner(j2,j3)
c
        if(LRotationalPeriodicity) then
c
          xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
          yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
          zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
        elseif(LTranslationalPeriodicity) then
c
          xF1=xc(j1)+xTranslation(j2)
          yF1=yc(j1)+yTranslation(j2)
          zF1=zc(j1)+zTranslation(j2)
c
        endif
c
        distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
        distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
        GFactCF=distance2/(distance1+distance2)
c        
        BTurbulentViscosity(i2,i3)=
     *            GFactCF*TurbulentViscosity(i1)+
     *                    (1.-GFactCF)*TurbulentViscosity(j1)
c
      enddo
c
      do i=1,Iaxis
c
        i1=IaxisOwner(i)
        i2=IaxisNumberOfBCSets(i)
        i3=IaxisNBFaces(i)
c
        BTurbulentViscosity(i2,i3)=TurbulentViscosity(i1)
c
      enddo
c
      do i=1,IpressureFarField
c
        i1=IpressureFarFieldOwner(i)
        i2=IpressureFarFieldNumberOfBCSets(i)
        i3=IpressureFarFieldNBFaces(i)
c
        if(Bmdot(i2,i3).gt.0.) then
c
          BTurbulentViscosity(i2,i3)=TurbulentViscosity(i1)
c
        endif
c
      enddo
c
      if(ModelNumber.eq.19.or.ModelNumber.eq.20.or.
     *                            ModelNumber.eq.24) then
c
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
c
          BTurbulentViscosity(i2,i3)=0.
c
        enddo
c      
      else
c
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
c
          BTurbulentViscosity(i2,i3)=TurbulentViscosity(i1)
c
        enddo
c      
c----------------------------------------------------------------------
c
      endif
c
      if(WallTreatment.eq.'lowreynoldsnumber') then
c
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
c
          BTurbulentViscosity(i2,i3)=0.
c
        enddo
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE UpdateInletTurbulentViscosity
c
C#############################################################################################
      use BoundaryConditionsTurbulence1, only: inletTypeT
      use Geometry1, only: NumberOfBCSets     
      use Geometry3, only: NBFaces,NBFaceOwner
      use PhysicalProperties1, only: BTurbulentViscosity,
     *                               TurbulentViscosity  
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
c********************************************************************************************
c      
      do i=1,NumberOfBCSets
c
        if(inletTypeT(i).eq.'specifiedkew') then
c
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
            BTurbulentViscosity(i,j)=TurbulentViscosity(k)
c
          enddo
c
        endif
c
      enddo
c      
      return
      end
