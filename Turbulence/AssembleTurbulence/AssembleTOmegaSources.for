c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentOmegaSources
c
C#############################################################################################
c
      use User0, only: MomentumWallFunctionType,Lcompressible,
     *                 WallTreatment,LRough
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BgDiff,BFaceTx,BFaceTy,BFaceTz
      use WallDistance1, only: WallDistance,BWallDistance
      use Variables1, only: TurbulentKE,
     *                      uVelocity,vVelocity,wVelocity,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      TempGradx,TempGrady,TempGradz,
     *                      TurbulenceProduction,
     *                      BTurbulenceProduction,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      TurbulentOmega,BTurbulentKE,
     *                      TKEGradx,TKEGrady,TKEGradz,
     *                      TEDGradx,TEDGrady,TEDGradz,
     *                      TOmegaGradx,TOmegaGrady,TOmegaGradz,
     *                      BTurbulentED,BTEDGradx,BTEDGrady,BTEDGradz,
     *                      BTurbulentOmega,BTOmegaGradx,BTOmegaGrady,
     *                      BTOmegaGradz,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
     *                      uvVelGradxy,BuvVelGradxy,
     *                      uwVelGradxz,BuwVelGradxz,
     *                      uVelGrad2x,BuVelGrad2x,uVelGrad2y,
     *                      BuVelGrad2y,uVelGrad2z,BuVelGrad2z,
     *                      vVelGrad2x,BvVelGrad2x,vVelGrad2y,
     *                      BvVelGrad2y,vVelGrad2z,BvVelGrad2z,
     *                      wVelGrad2x,BwVelGrad2x,wVelGrad2y,
     *                      BwVelGrad2y,wVelGrad2z,BwVelGrad2z,
     *                      TurbulentKL,BTurbulentKL,
     *                      BTurbulentKLGradx,BTurbulentKLGrady,
     *                      BTurbulentKLGradz
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: GravityX,GravityY,GravityZ,Density,
     *                               TurbulentViscosity,BDensity,
     *                               CoefficientOfThermalExpansion,
     *                               BViscosity,Viscosity,
     *                               TurbulentViscosity1,
     *                               BeDiffCoefficient
      use Turbulence1, only: SigT,cmu,rhoTED,cmu25,cmu75,cappa,ustar,
     *                       ystar,yplus,BrhoTED,ctrans,ce1,ce2,ce3,cc,
     *                       a1sst,alpha,betta,sigTED,alpha1,betta1,
     *                       sigTED1,alpha2,betta2,sigTED2,
     *                       F1factor,sigdo,KsPlus,
     *                       betta0,StrainRate,F4factor,fr1Coefficient,
     *                       alfa,tau11,tau12,tau13,tau22,tau23,tau33,
     *                       fmuCoefficient,f1Coefficient,f2Coefficient,
     *                       LTED,c11,c12,xi1,xi2,xi3,cd1,ModelNumber,
     *                       W11,W12,W13,W22,W23,W33,
     *                       S11,S12,S13,S22,S23,S33,
     *                       cmuR,c1R,BcmuR,Bc1R,
     *                       Vorticity,BVorticity,
     *                       StrainRate,BStrainRate,
     *                       Ctscrit,Ats,TurbulentKTs,LambdaEff,Clambda,
     *                       TurbulentViscosityTl,C11,C12,Ctl,
     *                       ProductionKL,ReT,Cw1,Cw2,Cw3,Cwr,
     *                       ProductionKT,BLambdaEff,SourceRbp,
     *                       SourceRnat,LambdaT,CoefficientFW,AlfaT
      use Constants1, only: tiny
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j !,j1,k1
      double precision :: term1,term2,term3,dNorm,OmegaViscous,OmegaLog,
     *                    gamf,FluxCflocal,FluxFflocal,FluxVflocal,
     *                    dfidxTf,dfidyTf,dfidzTf,TED1,alphabsl,
     *                    bettabsl,bettaUse,Xiw,fbetta,
     *                    CoefficientFws,SourceDL,ctrans1,
     *                    SourcePkt,Source1,Source2,Source3
c
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c********************************************************************************************
          character*10 Variable
          integer :: MethodCalcGradient,nIterGradientPhi
          logical :: LimitGradient
          integer :: LimitGradientMethod
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c********************************************************************************************
        end SUBROUTINE Gradient
c********************************************************************************************
        FUNCTION TangentialVelocity(i1)
c********************************************************************************************
          integer :: i1
          double precision :: TangentialVelocity
c********************************************************************************************
        end FUNCTION TangentialVelocity
c********************************************************************************************
      end interface
c********************************************************************************************
c
      CalculateTurbulentOmegaSources: select case(ModelNumber)
c-----------------------------------------------------------------------------
        case(11) CalculateTurbulentOmegaSources           !komega
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=dmax1(TurbulentOmega(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxTE(i)=FluxTE(i)-Volume(i)*
     *              dmax1(alpha*term1*TurbulenceProduction(i),0.)
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+2.*dmax1(betta*term1,0.)
            FluxTE(i)=FluxTE(i)+dmax1(betta*term1*TurbulentOmega(i),0.)
c        
          enddo
c
c--- Modify specific eddy diffusivity along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            if(.not.Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                if(MomentumWallFunctionType.eq.'automatic') then
c
                  OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                  OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
                elseif(MomentumWallFunctionType.eq.'standard') then
c
                  ctrans1=ctrans
                  if(LRough) ctrans1=KsPlus(i)/2.
c                  
                  if(ystar(i).lt.dmin1(ctrans,ctrans1)) then
c
                    TurbulentOmega(i1)=6.*Viscosity(i1)/
     *                            (Density(i1)*betta*dNorm**2)
c
                  else
c
                    TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  endif
c
                elseif(MomentumWallFunctionType.eq.'scalable') then
c
                  TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                endif
c
              enddo
c
            elseif(Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
              enddo
c
            endif
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
            call SetTOmegaWallBC
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
c              dNorm=Walldistance(i1)
c
c              BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*betta*dNorm**2)
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTOmegaGradx(i2,i3)
              dfidyTf=BTOmegaGrady(i2,i3)
              dfidzTf=BTOmegaGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(12) CalculateTurbulentOmegaSources           !komegaepsilon
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=dmax1(TurbulentOmega(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxTE(i)=FluxTE(i)-Volume(i)*
     *           dmax1(alpha*term1*TurbulenceProduction(i),0.)
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+2.*dmax1(betta*term1,0.)
            FluxTE(i)=FluxTE(i)+dmax1(betta*term1*TurbulentOmega(i),0.)
c
            term1=2.*sigTED*Density(i)/dmax1(TurbulentOmega(i),tiny)
            term2=TKEGradx(i)*TOmegaGradx(i)+
     *          TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
            term1=term1*term2*Volume(i)
c          
            if(term1.gt.0.) then         
c          
              FluxTE(i)=FluxTE(i)-term1
c          
            else
c
              TED1=dmax1(TurbulentOmega(i),tiny) 
              FluxCE(i)=FluxCE(i)-term1/TED1
              FluxTE(i)=FluxTE(i)-term1
c
            endif         
c
          enddo
c
c--- Modify specific eddy diffusivity along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            if(.not.Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                if(MomentumWallFunctionType.eq.'automatic') then
c
                  OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                  OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
                elseif(MomentumWallFunctionType.eq.'standard') then
c
                  ctrans1=ctrans
                  if(LRough) ctrans1=KsPlus(i)/2.
c                  
                  if(ystar(i).lt.dmin1(ctrans,ctrans1)) then
c
                    TurbulentOmega(i1)=6.*Viscosity(i1)/
     *                            (Density(i1)*betta*dNorm**2)
c
                  else
c
                    TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  endif
c
                elseif(MomentumWallFunctionType.eq.'scalable') then
c
                  TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                endif
c
              enddo
c
            elseif(Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
              enddo
c
            endif
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
            call SetTOmegaWallBC
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
c              dNorm=Walldistance(i1)
c
c              BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*0.075*dNorm**2)
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTOmegaGradx(i2,i3)
              dfidyTf=BTOmegaGrady(i2,i3)
              dfidzTf=BTOmegaGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(13) CalculateTurbulentOmegaSources           !komegabsl
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceInterpolationFactors
c
          do i=1,NumberOfElements    
c
            term1=dmax1(TurbulentOmega(i),0.)/dmax1(TurbulentKE(i),tiny)
            alpha1=betta1/cmu-sigTED1*cappa**2/dsqrt(cmu)
            alpha2=betta2/cmu-sigTED2*cappa**2/dsqrt(cmu)
            alphabsl=F1factor(i)*alpha1+(1.-F1factor(i))*alpha2
            bettabsl=F1factor(i)*betta1+(1.-F1factor(i))*betta2
c
            FluxTE(i)=FluxTE(i)-Volume(i)*
     *          dmax1(alphabsl*term1*TurbulenceProduction(i),0.)
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+2.*dmax1(bettabsl*term1,0.)
            FluxTE(i)=FluxTE(i)+
     *                 dmax1(bettabsl*term1*TurbulentOmega(i),0.)
c
            term1=2.*Density(i)*sigTED2/dmax1(TurbulentOmega(i),tiny)
            term2=TKEGradx(i)*TOmegaGradx(i)+
     *          TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
            term1=term1*term2*Volume(i)
c          
            if(term1.gt.0.) then         
c          
              FluxTE(i)=FluxTE(i)-(1.-F1factor(i))*term1
c          
            else
c
              TED1=dmax1(TurbulentOmega(i),tiny) 
              FluxCE(i)=FluxCE(i)-(1.-F1factor(i))*term1/TED1
              FluxTE(i)=FluxTE(i)-(1.-F1factor(i))*term1
c
            endif         
c        
          enddo
c
c--- Modify specific eddy diffusivity along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            if(.not.Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                if(MomentumWallFunctionType.eq.'automatic') then
c
                  OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                  OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
                elseif(MomentumWallFunctionType.eq.'standard') then
c
                  ctrans1=ctrans
                  if(LRough) ctrans1=KsPlus(i)/2.
c                  
                  if(ystar(i).lt.dmin1(ctrans,ctrans1)) then
c
                    TurbulentOmega(i1)=6.*Viscosity(i1)/
     *                            (Density(i1)*betta*dNorm**2)
c
                  else
c
                    TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  endif
c
                elseif(MomentumWallFunctionType.eq.'scalable') then
c
                  TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                endif
c
              enddo
c
            elseif(Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
              enddo
c
            endif
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
            call SetTOmegaWallBC
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
c              dNorm=Walldistance(i1)
c
c              BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*betta*dNorm**2)
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTOmegaGradx(i2,i3)
              dfidyTf=BTOmegaGrady(i2,i3)
              dfidzTf=BTOmegaGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(14) CalculateTurbulentOmegaSources           !komegasst
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceInterpolationFactors
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*TurbulenceProduction(i)*fr1Coefficient(i)/
     *                               dmax1(TurbulentViscosity(i),1.d-10)
            alphabsl=F1factor(i)*alpha1+(1.-F1factor(i))*alpha2
            FluxTE(i)=FluxTE(i)-Volume(i)*dmax1(alphabsl*term1,0.)
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
            bettabsl=F1factor(i)*betta1+(1.-F1factor(i))*betta2
            bettabsl=F4factor(i)*bettabsl
            FluxCE(i)=FluxCE(i)+2.*dmax1(bettabsl*term1,0.)
            FluxTE(i)=FluxTE(i)+
     *               dmax1(bettabsl*term1*TurbulentOmega(i),0.)
c
            term1=2.*Density(i)*sigTED2/dmax1(TurbulentOmega(i),tiny)
            term2=TKEGradx(i)*TOmegaGradx(i)+
     *          TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
            term1=(1.-F1factor(i))*term1*term2*Volume(i)
c          
            if(term1.gt.0.) then         
c          
              FluxTE(i)=FluxTE(i)-term1
c          
            else
c
              TED1=dmax1(TurbulentOmega(i),tiny) 
              FluxCE(i)=FluxCE(i)-term1/TED1
              FluxTE(i)=FluxTE(i)-term1
c
            endif         
c        
          enddo
c
c--- Modify specific eddy diffusivity along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            if(.not.Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                if(MomentumWallFunctionType.eq.'automatic') then
c
                  OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                  OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
                elseif(MomentumWallFunctionType.eq.'standard') then
c
                  ctrans1=ctrans
                  if(LRough) ctrans1=KsPlus(i)/2.
c                  
                  if(ystar(i).lt.dmin1(ctrans,ctrans1)) then
c
                    TurbulentOmega(i1)=6.*Viscosity(i1)/
     *                            (Density(i1)*betta*dNorm**2)
c
                  else
c
                    TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  endif
c
                elseif(MomentumWallFunctionType.eq.'scalable') then
c
                  TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                endif
c
              enddo
c
            elseif(Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
              enddo
c
            endif
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
c--- Apply boundary conditions along walls
c
c
            call SetTOmegaWallBC
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
c              dNorm=Walldistance(i1)
c
c              BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*betta*dNorm**2)
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTOmegaGradx(i2,i3)
              dfidyTf=BTOmegaGrady(i2,i3)
              dfidzTf=BTOmegaGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(15) CalculateTurbulentOmegaSources           !sstgamaretheta
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceInterpolationFactors
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*TurbulenceProduction(i)*fr1Coefficient(i)/
     *                               dmax1(TurbulentViscosity(i),1.d-10)
            alphabsl=F1factor(i)*alpha1+(1.-F1factor(i))*alpha2
            FluxTE(i)=FluxTE(i)-Volume(i)*dmax1(alphabsl*term1,0.)
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
            bettabsl=F1factor(i)*betta1+(1.-F1factor(i))*betta2
            bettabsl=F4factor(i)*bettabsl
            FluxCE(i)=FluxCE(i)+2.*dmax1(bettabsl*term1,0.)
            FluxTE(i)=FluxTE(i)+
     *               dmax1(bettabsl*term1*TurbulentOmega(i),0.)
c
            term1=2.*Density(i)*sigTED2/dmax1(TurbulentOmega(i),tiny)
            term2=TKEGradx(i)*TOmegaGradx(i)+
     *          TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
            term1=term1*term2*Volume(i)
c          
            if(term1.gt.0.) then         
c          
              FluxTE(i)=FluxTE(i)-(1.-F1factor(i))*term1
c          
            else
c
              TED1=dmax1(TurbulentOmega(i),tiny) 
              FluxCE(i)=FluxCE(i)-(1.-F1factor(i))*term1/TED1
              FluxTE(i)=FluxTE(i)-(1.-F1factor(i))*term1
c
            endif         
c        
          enddo
c
c--- Modify specific eddy diffusivity along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            Print*,'Transitional model is a low-Reynolds number model'
            Print*,'Wall treatment cannot be wall functions'
            Print*,'Program is terminated'
            stop
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
c--- Apply boundary conditions along walls
c
c
            call SetTOmegaWallBC
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
c              dNorm=Walldistance(i1)
c
c              BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*betta*dNorm**2)
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTOmegaGradx(i2,i3)
              dfidyTf=BTOmegaGrady(i2,i3)
              dfidzTf=BTOmegaGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(16) CalculateTurbulentOmegaSources           !komega2006
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*TurbulenceProduction(i)/
     *                              dmax1(TurbulentViscosity1(i),tiny)
            FluxTE(i)=FluxTE(i)-Volume(i)*dmax1(alpha*term1,0.)
            bettaUse=betta0
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
c          
            term3=0.5*(uVelGradx(i)+vVelGrady(i)+wVelGradz(i))
c
c---- The full term
c
!            Xiw=W11(i)*W11(i)*(S11(i)-term3)+W11(i)*W12(i)*S12(i)+
!     *       W11(i)*W13(i)*S13(i)+W12(i)*(-W12(i))*(S11(i)-term3)+
!     *       W12(i)*W22(i)*S12(i)+W12(i)*W23(i)*S13(i)+
!     *       W13(i)*(-W13(i))*(S11(i)-term3)+W13(i)*(-W23(i))*S12(i)+
!     *       W13(i)*W33(i)*S13(i)+(-W12(i))*W11(i)*S12(i)+
!     *       (-W12(i))*W12(i)*(S22(i)-term3)+(-W12(i))*W13(i)*S23(i)+
!     *       W22(i)*(-W12(i))*S12(i)+W22(i)*W22(i)*(S22(i)-term3)+
!     *       W22(i)*W23(i)*S23(i)+W23(i)*(-W13(i))*S12(i)+
!     *       W23(i)*(-W23(i))*(S22(i)-term3)+W23(i)*W33(i)*S23(i)+
!     *      (-W13(i))*W11(i)*S13(i)+(-W13(i))*W12(i)*S23(i)+
!     *      (-W13(i))*W13(i)*(S33(i)-term3)+(-W23(i))*(-W12(i))*S13(i)+
!     *      (-W23(i))*W22(i)*S23(i)+(-W23(i))*W23(i)*(S33(i)-term3)+
!     *       W33(i)*(-W13(i))*S13(i)+W33(i)*(-W23(i))*S23(i)+
!     *       W33(i)*W33(i)*(S33(i)-term3)
c
c---- The nozero terms
c
            Xiw=W12(i)*(-W12(i))*(S11(i)-term3)+W12(i)*W23(i)*S13(i)+
     *       W13(i)*(-W13(i))*(S11(i)-term3)+W13(i)*(-W23(i))*S12(i)+
     *       (-W12(i))*W12(i)*(S22(i)-term3)+(-W12(i))*W13(i)*S23(i)+
     *       W23(i)*(-W13(i))*S12(i)+W23(i)*(-W23(i))*(S22(i)-term3)+
     *      (-W13(i))*W12(i)*S23(i)+(-W13(i))*W13(i)*(S33(i)-term3)+
     *      (-W23(i))*(-W12(i))*S13(i)+(-W23(i))*W23(i)*(S33(i)-term3)
c
            ted1=dmax1(TurbulentOmega(i),tiny)
            Xiw=dabs(Xiw/((cmu*ted1)**3))
c          
            fbetta=(1.+85.*Xiw)/(1.+100.*Xiw)
            bettaUse=betta0*fbetta          
c          
            FluxCE(i)=FluxCE(i)+2.*dmax1(bettaUse*term1,0.)
            FluxTE(i)=FluxTE(i)+
     *                dmax1(bettaUse*term1*TurbulentOmega(i),0.)
c
            term2=TKEGradx(i)*TOmegaGradx(i)+
     *           TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
c
            if(term2.gt.0.) then
c
              term1=Density(i)*sigdo/dmax1(TurbulentOmega(i),tiny)
c          
              FluxTE(i)=FluxTE(i)-term1*term2*Volume(i)
c          
            endif         
c        
          enddo
c
          if(WallTreatment.eq.'wallfunctions') then
c
c--- Modify specific eddy diffusivity along walls
c
            if(.not.Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                if(MomentumWallFunctionType.eq.'automatic') then
c
                  OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                  OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
                elseif(MomentumWallFunctionType.eq.'standard') then
c
                  ctrans1=ctrans
                  if(LRough) ctrans1=KsPlus(i)/2.
c                  
                  if(ystar(i).lt.dmin1(ctrans,ctrans1)) then
c
                    TurbulentOmega(i1)=6.*Viscosity(i1)/
     *                            (Density(i1)*betta*dNorm**2)
c
                  else
c
                    TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                  endif
c
                elseif(MomentumWallFunctionType.eq.'scalable') then
c
                  TurbulentOmega(i1)=BDensity(i2,i3)*ustar(i)*
     *                ustar(i)/(a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                endif
c
              enddo
c
            elseif(Lcompressible) then
c
              do i=1,IwallTurbulence
c
                i1=IWallTurbulenceOwner(i)
                i2=IWallTurbulenceNumberOfBCSets(i)
                i3=IWallTurbulenceNBFaces(i)
c
                dNorm=Walldistance(i1)
c
                OmegaViscous=
     *                 6.*Viscosity(i1)/(betta*Density(i1)*dNorm**2)
                OmegaLog=BDensity(i2,i3)*ustar(i)*ustar(i)/
     *                      (a1sst*cappa*BViscosity(i2,i3)*ystar(i))
c
                TurbulentOmega(i1)=dsqrt(OmegaViscous**2+OmegaLog**2)
c
              enddo
c
            endif
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
            call SetTOmegaWallBC
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
c              dNorm=Walldistance(i1)
c
c              BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*0.075*dNorm**2)
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTOmegaGradx(i2,i3)
              dfidyTf=BTOmegaGrady(i2,i3)
              dfidzTf=BTOmegaGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(17) CalculateTurbulentOmegaSources           !komega2006lrn
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            term1=alfa(i)*dmax1(TurbulentOmega(i),0.)/
     *                           dmax1(TurbulentKE(i),tiny) 
            term1=term1*TurbulenceProduction(i)
            FluxTE(i)=FluxTE(i)-Volume(i)*dmax1(term1,0.)
            bettaUse=betta0
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
c          
            term3=0.5*(uVelGradx(i)+vVelGrady(i)+wVelGradz(i))
c
c---- The full term
c
!            Xiw=W11(i)*W11(i)*(S11(i)-term3)+W11(i)*W12(i)*S12(i)+
!     *       W11(i)*W13(i)*S13(i)+W12(i)*(-W12(i))*(S11(i)-term3)+
!     *       W12(i)*W22(i)*S12(i)+W12(i)*W23(i)*S13(i)+
!     *       W13(i)*(-W13(i))*(S11(i)-term3)+W13(i)*(-W23(i))*S12(i)+
!     *       W13(i)*W33(i)*S13(i)+(-W12(i))*W11(i)*S12(i)+
!     *       (-W12(i))*W12(i)*(S22(i)-term3)+(-W12(i))*W13(i)*S23(i)+
!     *       W22(i)*(-W12(i))*S12(i)+W22(i)*W22(i)*(S22(i)-term3)+
!     *       W22(i)*W23(i)*S23(i)+W23(i)*(-W13(i))*S12(i)+
!     *       W23(i)*(-W23(i))*(S22(i)-term3)+W23(i)*W33(i)*S23(i)+
!     *      (-W13(i))*W11(i)*S13(i)+(-W13(i))*W12(i)*S23(i)+
!     *      (-W13(i))*W13(i)*(S33(i)-term3)+(-W23(i))*(-W12(i))*S13(i)+
!     *      (-W23(i))*W22(i)*S23(i)+(-W23(i))*W23(i)*(S33(i)-term3)+
!     *       W33(i)*(-W13(i))*S13(i)+W33(i)*(-W23(i))*S23(i)+
!     *       W33(i)*W33(i)*(S33(i)-term3)
c
c---- The nozero terms
c
            Xiw=W12(i)*(-W12(i))*(S11(i)-term3)+W12(i)*W23(i)*S13(i)+
     *       W13(i)*(-W13(i))*(S11(i)-term3)+W13(i)*(-W23(i))*S12(i)+
     *       (-W12(i))*W12(i)*(S22(i)-term3)+(-W12(i))*W13(i)*S23(i)+
     *       W23(i)*(-W13(i))*S12(i)+W23(i)*(-W23(i))*(S22(i)-term3)+
     *      (-W13(i))*W12(i)*S23(i)+(-W13(i))*W13(i)*(S33(i)-term3)+
     *      (-W23(i))*(-W12(i))*S13(i)+(-W23(i))*W23(i)*(S33(i)-term3)
c
            ted1=dmax1(TurbulentOmega(i),tiny)
            Xiw=dabs(Xiw/((cmu*ted1)**3))
c          
            fbetta=(1.+85.*Xiw)/(1.+100.*Xiw)
            bettaUse=betta0*fbetta          
c          
            FluxCE(i)=FluxCE(i)+2.*dmax1(bettaUse*term1,0.)
            FluxTE(i)=FluxTE(i)+
     *                dmax1(bettaUse*term1*TurbulentOmega(i),0.)
c
            term2=TKEGradx(i)*TOmegaGradx(i)+
     *            TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
c
            if(term2.gt.0.) then
c
              term1=Density(i)*sigdo/dmax1(TurbulentOmega(i),tiny)
c          
              FluxTE(i)=FluxTE(i)-term1*term2*Volume(i)
c          
            endif         
c        
          enddo
c
          call SetTOmegaWallBC
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
c            dNorm=Walldistance(i1)
c
c            BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*0.075*dNorm**2)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTOmegaGradx(i2,i3)
            dfidyTf=BTOmegaGrady(i2,i3)
            dfidzTf=BTOmegaGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(21) CalculateTurbulentOmegaSources           !transitional k-kl-w model
c-----------------------------------------------------------------------------
c
            do i=1,NumberOfElements    
c
              CoefficientFws=1.-
     *           dexp(-0.41*((LambdaEff(i)/dmax1(LambdaT(i),tiny))**4))
              SourceDL=Cw3*CoefficientFws*AlfaT(i)*CoefficientFW(i)*
     *                CoefficientFW(i)*dsqrt(dmax1(TurbulentKE(i),0.))/
     *                                   (WallDistance(i)**3)
              SourcePkt=Cw1*dmax1(TurbulentOmega(i),0.)*ProductionKT(i)/
     *                dmax1(TurbulentKE(i),tiny)
              Source1=Cwr*dmax1(TurbulentOmega(i),0.)*
     *         (SourceRbp(i)+SourceRnat(i))*dmax1(TurbulentKL(i),0.)/
     *                       dmax1(TurbulentKE(i)*CoefficientFW(i),tiny)
              Source2=(SourceRbp(i)+SourceRnat(i))*
     *               dmax1(TurbulentKL(i),0.)/dmax1(TurbulentKE(i),tiny)
              Source3=Cw2*dmax1(TurbulentOmega(i),0.)*
     *                            CoefficientFW(i)*CoefficientFW(i)
c
              FluxCE(i)=FluxCE(i)+
     *                    Density(i)*Volume(i)*(Source2+Source3)
              FluxTE(i)=FluxTE(i)-
     *              Density(i)*Volume(i)*(SourcePkt+Source1+SourceDL-
     *                   (Source2+Source3)*dmax1(TurbulentOmega(i),0.))
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
              BTurbulentOmega(i2,i3)=TurbulentOmega(i1)
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
            enddo
c
c-----------------------------------------------------------------------------
        case(23) CalculateTurbulentOmegaSources           !sstgama
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceInterpolationFactors
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*TurbulenceProduction(i)*fr1Coefficient(i)/
     *                               dmax1(TurbulentViscosity(i),1.d-10)
            alphabsl=F1factor(i)*alpha1+(1.-F1factor(i))*alpha2
            FluxTE(i)=FluxTE(i)-Volume(i)*dmax1(alphabsl*term1,0.)
c
            term1=Density(i)*dmax1(TurbulentOmega(i),0.)*Volume(i)
            bettabsl=F1factor(i)*betta1+(1.-F1factor(i))*betta2
            bettabsl=F4factor(i)*bettabsl
            FluxCE(i)=FluxCE(i)+2.*dmax1(bettabsl*term1,0.)
            FluxTE(i)=FluxTE(i)+
     *               dmax1(bettabsl*term1*TurbulentOmega(i),0.)
c
            term1=2.*Density(i)*sigTED2/dmax1(TurbulentOmega(i),tiny)
            term2=TKEGradx(i)*TOmegaGradx(i)+
     *          TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
            term1=(1.-F1factor(i))*term1*term2*Volume(i)
c          
            if(term1.gt.0.) then         
c          
              FluxTE(i)=FluxTE(i)-term1
c          
            else
c
              TED1=dmax1(TurbulentOmega(i),tiny) 
              FluxCE(i)=FluxCE(i)-term1/TED1
              FluxTE(i)=FluxTE(i)-term1
c
            endif         
c        
          enddo
c
c--- Apply boundary conditions along walls
c
          call SetTOmegaWallBC
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
c            dNorm=Walldistance(i1)
c
c            BTurbulentOmega(i2,i3)=10.*6.*Viscosity(i1)/
c     *                          (Density(i1)*betta*dNorm**2)
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTOmegaGradx(i2,i3)
            dfidyTf=BTOmegaGrady(i2,i3)
            dfidzTf=BTOmegaGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentOmega(i1)+
     *                 FluxFflocal*BTurbulentOmega(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTurbulentOmegaSources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end