c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentEDSources
c
C#############################################################################################
c
      use User0, only: LstagnationCorrectionKE,
     *                 MomentumWallFunctionType,
     *                 WallTreatment,urfWallEpsilon
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     BgDiff,BFaceTx,BFaceTy,BFaceTz
      use WallDistance1, only: WallDistance,BWallDistance
      use Variables1, only: TurbulentKE,TurbulentED,
     *                      TurbulenceProduction,
     *                      TKEGradx,TKEGrady,TKEGradz,
     *                      TEDGradx,TEDGrady,TEDGradz,
     *                      BTurbulentED,
     *                      BTEDGradx,BTEDGrady,BTEDGradz
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density,BDensity,
     *                               TurbulentViscosity,
     *                               BViscosity,Viscosity,
     *                               BeDiffCoefficient
      use Turbulence1, only: cmu,cmu25,cmu75,cappa,
     *                       ystar,ctrans,ce1,ce2,ce3,
     *                       alpha,betta,LTED,StrainRate,ReT,
     *                       fmuCoefficient,f1Coefficient,f2Coefficient,
     *                       cmuR,c1R,BcmuR,Bc1R,ModelNumber,ustar,
     *                       C2eRNG,TScale,Ce1Coefficient
      use Constants1, only: tiny
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j !,j1,k1
      double precision :: term0,term1,term2,term3,term4,term5,term6,
     *                    term7,term8,term9,term10,sum,TED1,cmu25k,
     *                    dNorm,Ts,gamf,FluxCflocal,FluxFflocal,
     *                    FluxVflocal,dfidxTf,dfidyTf,dfidzTf,
     *                    ted,tke,tvis,nu
c********************************************************************************************
c
      CalculateTurbulentEDSources: select case(ModelNumber)
c-------------------------------------------------------------------------------------------
        case(1) CalculateTurbulentEDSources           !kepsilon
c-----------------------------------------------------------------------------
c
          if(.not.LstagnationCorrectionKE) then
c
            do i=1,NumberOfElements    
c
              term1=Density(i)*cmu*dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
              FluxCE(i)=FluxCE(i)+ce2*term1*Density(i)
              FluxTE(i)=FluxTE(i)-
     *             dmax1(ce1*term1*TurbulenceProduction(i),0.)+
     *                   ce2*term1*Density(i)*dmax1(TurbulentED(i),0.)
c        
            enddo
c
          else
c
            do i=1,NumberOfElements    
c
              sum=StrainRate(i)
c
              Ts=dmin1(dmax1(TurbulentKE(i),tiny)/
     *         dmax1(TurbulentED(i),tiny),alpha/(dsqrt(betta)*cmu*sum))
c
              term1=Volume(i)/dmax1(Ts,tiny)
c
              FluxCE(i)=FluxCE(i)+dmax1(ce2*term1*Density(i),0.)
              FluxTE(i)=FluxTE(i)-
     *             dmax1(ce1*term1*TurbulenceProduction(i),0.)+
     *             dmax1(ce2*term1*Density(i),0.)*TurbulentED(i)
c        
            enddo
c
          endif
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=Walldistance(i1)
c
            if(MomentumWallFunctionType.eq.'standard') then
c
              if(ystar(i).lt.ctrans) then
c
                TurbulentED(i1)=2.*BViscosity(i2,i3)*
     *           dmax1(TurbulentKE(i1),0.)/(Bdensity(i2,i3)*dNorm**2)
c                TurbulentED(i1)=BDensity(i2,i3)*cmu*
c     *                   dmax1(TurbulentKE(i1),0.)/BViscosity(i2,i3)
c
              else
c
                TurbulentED(i1)=cmu75*
     *           ((dmax1(TurbulentKE(i1),0.))**1.5)/(cappa*dNorm)
c
              endif
c
            elseif(MomentumWallFunctionType.eq.'scalable') then
c
              dNorm=ystar(i)*BViscosity(i2,i3)/
     *                           (BDensity(i2,i3)*ustar(i))
              TurbulentED(i1)=cmu75*
     *         ((dmax1(TurbulentKE(i1),0.))**1.5)/(cappa*dNorm)
c
            endif
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(2) CalculateTurbulentEDSources           !kepsilonchien
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)
            FluxCE(i)=FluxCE(i)+LTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)+
     *           LTED(i)*Volume(i)*dmax1(TurbulentED(i),0.)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            BTurbulentED(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(3) CalculateTurbulentEDSources           !kepsilonsharma
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)
            FluxTE(i)=FluxTE(i)-LTED(i)*Volume(i)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            BTurbulentED(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(4) CalculateTurbulentEDSources           !kepsilonchc
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            dNorm=Walldistance(i1)
c
            term1=4.*Viscosity(i1)*TurbulentKE(i1)/
     *                        (Density(i1)*dNorm**2)
c            BTurbulentED(i2,i3)=0.5*BTurbulentED(i2,i3)+
c     *               0.5*dmax1(4.*Viscosity(i1)*TurbulentKE(i1)/
c     *                        (Density(i1)*dNorm**2)-TurbulentED(i1),0.)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
c            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxCf(i4)=FluxCf(i4)+2.*FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
c            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
c     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+2.*FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*term1+FluxVflocal
c
            BTurbulentED(i2,i3)=(1.-urfWallEpsilon)*BTurbulentED(i2,i3)+
     *                    urfWallEpsilon*dmax1(term1-TurbulentED(i1),0.)
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(5) CalculateTurbulentEDSources           !kepsilonkasagi
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            dNorm=Walldistance(i1)
c
            term1=4.*Viscosity(i1)*TurbulentKE(i1)/
     *                        (Density(i1)*dNorm**2)
c            BTurbulentED(i2,i3)=0.5*BTurbulentED(i2,i3)+
c     *               0.5*dmax1(4.*Viscosity(i1)*TurbulentKE(i1)/
c     *                        (Density(i1)*dNorm**2)-TurbulentED(i1),0.)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
c            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxCf(i4)=FluxCf(i4)+2.*FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
c            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
c     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+2.*FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*term1+FluxVflocal
c
            BTurbulentED(i2,i3)=(1.-urfWallEpsilon)*BTurbulentED(i2,i3)+
     *                    urfWallEpsilon*dmax1(term1-TurbulentED(i1),0.)
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(6) CalculateTurbulentEDSources           !kepsilontagawa
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            dNorm=Walldistance(i1)
c
            term1=4.*Viscosity(i1)*TurbulentKE(i1)/
     *                        (Density(i1)*dNorm**2)
c            BTurbulentED(i2,i3)=0.5*BTurbulentED(i2,i3)+
c     *               0.5*dmax1(4.*Viscosity(i1)*TurbulentKE(i1)/
c     *                        (Density(i1)*dNorm**2)-TurbulentED(i1),0.)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *            dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
c            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxCf(i4)=FluxCf(i4)+2.*FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
c            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
c     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+2.*FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*term1+FluxVflocal
c
            BTurbulentED(i2,i3)=(1.-urfWallEpsilon)*BTurbulentED(i2,i3)+
     *                    urfWallEpsilon*dmax1(term1-TurbulentED(i1),0.)
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(7) CalculateTurbulentEDSources           !kepsilonhishida
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)

            FluxCE(i)=FluxCE(i)+fmuCoefficient(i)*LTED(i)*
     *                     Volume(i)/dmax1(TurbulentED(i),tiny)

            FluxTE(i)=FluxTE(i)-(1.-fmuCoefficient(i))*LTED(i)*Volume(i)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            BTurbulentED(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(8) CalculateTurbulentEDSources           !kelambremhorst
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            FluxCflocal=0.
            FluxFflocal=0.
            FluxVflocal=0.
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxVflocal
c
            BTurbulentED(i2,i3)=TurbulentED(i1)
c
          enddo





!c
!          do i=1,IwallTurbulence
!c
!            i1=IWallTurbulenceOwner(i)
!            i2=IWallTurbulenceNumberOfBCSets(i)
!            i3=IWallTurbulenceNBFaces(i)
!            i4=NIFaces
!c
!            do j=1,i2-1
!c
!              i4=i4+NBFaces(j)
!c
!            enddo
!c
!            i4=i4+i3
!c
cc            nx=BFaceAreanx(i2,i3)
cc            ny=BFaceAreany(i2,i3)
cc            nz=BFaceAreanz(i2,i3)
c            dNorm=Walldistance(i1)
cc            dNorm=BDistanceCFx(i2,i3)*nx+
cc     *            BDistanceCFy(i2,i3)*ny+BDistanceCFz(i2,i3)*nz
!            term1=4.*Viscosity(i1)*TurbulentKE(i1)/
!     *                        (Density(i1)*dNorm**2)
!c            BTurbulentED(i2,i3)=0.5*BTurbulentED(i2,i3)+
!c     *               0.5*dmax1(4.*Viscosity(i1)*TurbulentKE(i1)/
!c     *                        (Density(i1)*dNorm**2)-TurbulentED(i1),0.)
!c
!            gamf=BeDiffCoefficient(i2,i3)
!            dfidxTf=BTEDGradx(i2,i3)
!            dfidyTf=BTEDGrady(i2,i3)
!            dfidzTf=BTEDGradz(i2,i3)
!c
!            FluxCflocal= gamf*BgDiff(i2,i3)
!            FluxFflocal=-gamf*BgDiff(i2,i3)
!            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
!     *            dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
!c
!c            FluxCf(i4)=FluxCf(i4)+FluxCflocal
!            FluxCf(i4)=FluxCf(i4)+2.*FluxCflocal
!            FluxFf(i4)=FluxFf(i4)+FluxFflocal
!            FluxVf(i4)=FluxVf(i4)+FluxVflocal
!c            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
!c     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
!            FluxTf(i4)=FluxTf(i4)+2.*FluxCflocal*TurbulentED(i1)+
!     *                      FluxFflocal*term1+FluxVflocal
!c
!            BTurbulentED(i2,i3)=(1.-urfWallEpsilon)*BTurbulentED(i2,i3)+
!     *                    urfWallEpsilon*dmax1(term1-TurbulentED(i1),0.)
!c
!          enddo










c
c-------------------------------------------------------------------------------------------
        case(9) CalculateTurbulentEDSources           !kelambremhorstm
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*fmuCoefficient(i)*
     *               dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+ce2*f2Coefficient(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-ce1*f1Coefficient(i)*term1*
     *         dmax1(TurbulenceProduction(i),0.)+ce2*f2Coefficient(i)*
     *                         term1*Density(i)*dmax1(TurbulentED(i),0.)
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            FluxCflocal=0.
            FluxFflocal=0.
            FluxVflocal=0.
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxVflocal
c
            BTurbulentED(i2,i3)=TurbulentED(i1)
c
          enddo







!c
!          do i=1,IwallTurbulence
!c
!            i1=IWallTurbulenceOwner(i)
!            i2=IWallTurbulenceNumberOfBCSets(i)
!            i3=IWallTurbulenceNBFaces(i)
!            i4=NIFaces
!c
!            do j=1,i2-1
!c
!              i4=i4+NBFaces(j)
!c
!            enddo
!c
!            i4=i4+i3
!c
cc            nx=BFaceAreanx(i2,i3)
cc            ny=BFaceAreany(i2,i3)
cc            nz=BFaceAreanz(i2,i3)
c            dNorm=Walldistance(i1)
cc            dNorm=BDistanceCFx(i2,i3)*nx+
cc     *            BDistanceCFy(i2,i3)*ny+BDistanceCFz(i2,i3)*nz
!            term1=4.*Viscosity(i1)*TurbulentKE(i1)/
!     *                        (Density(i1)*dNorm**2)
!c            BTurbulentED(i2,i3)=0.5*BTurbulentED(i2,i3)+
!c     *               0.5*dmax1(4.*Viscosity(i1)*TurbulentKE(i1)/
!c     *                        (Density(i1)*dNorm**2)-TurbulentED(i1),0.)
!c
!            gamf=BeDiffCoefficient(i2,i3)
!            dfidxTf=BTEDGradx(i2,i3)
!            dfidyTf=BTEDGrady(i2,i3)
!            dfidzTf=BTEDGradz(i2,i3)
!c
!            FluxCflocal= gamf*BgDiff(i2,i3)
!            FluxFflocal=-gamf*BgDiff(i2,i3)
!            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
!     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
!c
!c            FluxCf(i4)=FluxCf(i4)+FluxCflocal
!            FluxCf(i4)=FluxCf(i4)+2.*FluxCflocal
!            FluxFf(i4)=FluxFf(i4)+FluxFflocal
!            FluxVf(i4)=FluxVf(i4)+FluxVflocal
!c            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
!c     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
!            FluxTf(i4)=FluxTf(i4)+2.*FluxCflocal*TurbulentED(i1)+
!     *                      FluxFflocal*term1+FluxVflocal
!c
!            BTurbulentED(i2,i3)=(1.-urfWallEpsilon)*BTurbulentED(i2,i3)+
!     *                    urfWallEpsilon*dmax1(term1-TurbulentED(i1),0.)
!c
!          enddo















c
c-------------------------------------------------------------------------------------------
        case(10) CalculateTurbulentEDSources           !realizable
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            ted=dmax1(TurbulentED(i),tiny)
            tke=dmax1(TurbulentKE(i),tiny)
            tvis=dmax1(Viscosity(i),tiny)/Density(i)
c
            term1=Density(i)*c1R(i)*StrainRate(i)*ted*Volume(i)
            term2=Density(i)*ce2*ted*Volume(i)/(tke+dsqrt(tvis*ted))
c
            FluxCE(i)=FluxCE(i)+term2
            FluxTE(i)=FluxTE(i)-term1+term2*ted
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=Walldistance(i1)
c
            if(MomentumWallFunctionType.eq.'standard') then
c
              if(ystar(i).lt.ctrans) then
c
              TurbulentED(i1)=2.*BViscosity(i2,i3)*
     *         dmax1(TurbulentKE(i1),0.)/(Bdensity(i2,i3)*dNorm**2)
c                TurbulentED(i1)=BDensity(i2,i3)*BcmuR(i2,i3)*
c     *                   dmax1(TurbulentKE(i1),0.)/BViscosity(i2,i3)
c
              else
c
c                TurbulentED(i1)=(BcmuR(i2,i3)**0.75)*
                TurbulentED(i1)=(cmu**0.75)*
     *           ((dmax1(TurbulentKE(i1),0.))**1.5)/(cappa*dNorm)
c
              endif
c
            elseif(MomentumWallFunctionType.eq.'scalable') then
c
              dNorm=ystar(i)*BViscosity(i2,i3)/
     *                           (BDensity(i2,i3)*ustar(i))
c              TurbulentED(i1)=(BcmuR(i2,i3)**0.75)*
              TurbulentED(i1)=(cmu**0.75)*
     *         ((dmax1(TurbulentKE(i1),0.))**1.5)/(cappa*dNorm)
c
            endif
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(22) CalculateTurbulentEDSources           !K-Epsilon-Rt
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            ReT(i)=Density(i)*TurbulentKE(i)*TurbulentKE(i)/
     *                (Viscosity(i)*dmax1(TurbulentED(i),tiny))
            term2=dmax1(1.,dsqrt(2./ReT(i)))
            term0=Volume(i)/term2

            term1=dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+term0*ce2*term1*Density(i)
            FluxTE(i)=FluxTE(i)-
     *           term0*dmax1(ce1*term1*TurbulenceProduction(i),0.)+
     *              term0*ce2*term1*Density(i)*dmax1(TurbulentED(i),0.)
c        
            term3=Density(i)*dsqrt(dmax1(TurbulentKE(i),0.)*term2)
            nu=Viscosity(i)/Density(i)
            term4=(nu*dmax1(TurbulentED(i),0.))**0.25
            term5=dmax1(dsqrt(dmax1(TurbulentKE(i),0.)),term4)
            term6=TKEGradx(i)**2+TKEGrady(i)**2+TKEGradz(i)**2
            term7=TKEGradx(i)*TEDGradx(i)+
     *               TKEGrady(i)*TEDGrady(i)+TKEGradz(i)*TEDGradz(i)
            term8=term6/dmax1(TurbulentED(i),tiny)-term7*
     *            dmax1(TurbulentKE(i),0.)/dmax1(TurbulentED(i)**2,tiny)
            term9=dmax1(term8,0.)
            term10=term3*term5*term9
            FluxTE(i)=FluxTE(i)-term0*ce3*term10*term1
c
          enddo
c
c--- Apply boundary conditions along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            dNorm=Walldistance(i1)
c
            BTurbulentED(i2,i3)=2.*BViscosity(i2,i3)*
     *             dmax1(TurbulentKE(i1),0.)/(BDensity(i2,i3)*dNorm**2)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*dmax1(TurbulentED(i1),0.)+
     *                       FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(25) CalculateTurbulentEDSources           !kepsilonrng
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*cmu*dmax1(TurbulentKE(i),0.)*Volume(i)/
     *                                dmax1(TurbulentViscosity(i),tiny)
            FluxCE(i)=FluxCE(i)+C2eRNG(i)*term1*Density(i)
            FluxTE(i)=FluxTE(i)-
     *          dmax1(ce1*term1*TurbulenceProduction(i),0.)+
     *               C2eRNG(i)*term1*Density(i)*dmax1(TurbulentED(i),0.)
c       
          enddo
c
c--- Modify eddy diffusivity along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              dNorm=Walldistance(i1)
c
              if(MomentumWallFunctionType.eq.'standard') then
c
                if(ystar(i).lt.ctrans) then
c
                  TurbulentED(i1)=2.*Viscosity(i1)*
     *             dmax1(TurbulentKE(i1),0.)/(density(i1)*dNorm**2)
c                  TurbulentED(i1)=BDensity(i2,i3)*cmu*
c     *                   dmax1(TurbulentKE(i1),0.)/BViscosity(i2,i3)
c
                else
c
                  TurbulentED(i1)=cmu75*
     *             ((dmax1(TurbulentKE(i1),0.))**1.5)/(cappa*dNorm)
c
                endif
c
              elseif(MomentumWallFunctionType.eq.'scalable') then
c
                dNorm=ystar(i)*BViscosity(i2,i3)/
     *                           (BDensity(i2,i3)*ustar(i))
                TurbulentED(i1)=cmu75*
     *            ((dmax1(TurbulentKE(i1),0.))**1.5)/(cappa*dNorm)
c
              endif
c
            enddo
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
          endif
c
c-------------------------------------------------------------------------------------------
        case(26) CalculateTurbulentEDSources           !kepsilonv2f
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=ce2*Density(i)*Volume(i)/dmax1(TScale(i),tiny)
            term2=Ce1Coefficient(i)*dmax1(TurbulenceProduction(i),0.)*
     *           Volume(i)/dmax1(TScale(i),tiny)
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentED(i),0.)-term2
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            dNorm=Walldistance(i1)
c
            BTurbulentED(i2,i3)=2.*Viscosity(i1)*
     *             dmax1(TurbulentKE(i1),0.)/(Density(i1)*dNorm**2)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
c
          enddo
c
c-------------------------------------------------------------------------------------------
        case(27) CalculateTurbulentEDSources           !kepsilonzeta
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=ce2*Density(i)*Volume(i)/dmax1(TScale(i),tiny)
            term2=Ce1Coefficient(i)*dmax1(TurbulenceProduction(i),0.)*
     *           Volume(i)/dmax1(TScale(i),tiny)
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentED(i),0.)-term2
c        
          enddo
c
c--- Modify eddy diffusivity along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            dNorm=Walldistance(i1)
c
            BTurbulentED(i2,i3)=2.*Viscosity(i1)*
     *             dmax1(TurbulentKE(i1),0.)/(Density(i1)*dNorm**2)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTEDGradx(i2,i3)
            dfidyTf=BTEDGrady(i2,i3)
            dfidzTf=BTEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentED(i1)+
     *                      FluxFflocal*BTurbulentED(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTurbulentEDSources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end