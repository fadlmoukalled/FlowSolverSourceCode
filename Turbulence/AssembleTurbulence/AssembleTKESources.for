c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentKESources
c
C#############################################################################################
c
      use User0, only: WallTreatment,LimitTurbulenceProduction,
     *                 LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LSolveTurbulentKL,MomentumWallFunctionType
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     BgDiff,BFaceTx,BFaceTy,BFaceTz
      use Variables1, only: TurbulentKE,TurbulentED,
     *                      uVelocity,vVelocity,wVelocity,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      TempGradx,TempGrady,TempGradz,
     *                      TurbulenceProduction,
     *                      BTurbulenceProduction,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      BTurbulentKE,TurbulentOmega,
     *                      BTKEGradx,BTKEGrady,BTKEGradz,
     *                      BTurbulentOmega,BTurbulentED,Temperature,
     *                      TurbulentKL,BTurbulentKL,TGammaEff,TGamma
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density,BViscosity,Viscosity,
     *                               TurbulentViscosity,BDensity,
     *                               BeDiffCoefficient,
     *                               BTurbulentViscosity
      use Turbulence1, only: SigT,cmu,rhoTED,cmu25,cmu75,cappa,
     *                       yplus,BrhoTED,ctrans,c1limiter,ustar,uplus,
     *                       duplusdyplus,betta,Tau11,Tau12,Tau13,ystar,
     *                       Tau22,Tau23,Tau33,fr1Coefficient,
     *                       bettaStar,LTKE,ModelNumber,ystar,
     *                       A0,As,Av,Abp,Anat,Ats,Cbpcrit,Cnc,
     *                       Cnatcrit,Cint,Ctscrit,Crnat,Cr,Ca0,Css,
     *                       Ctl,Cwr,Clambda,Ck,Csep,Reoclim,
     *                       coefficientFW,
     *                       TurbulentViscosityTs,
     *                       ProductionKT,SourceRbp,SourceRnat,
     *                       sqrtTurbulentKE,
     *                       sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz,
     *                       BsqrtTurbulentKE,
     *                       BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,
     *                       Vorticity,StrainRate,
     *                       BVorticity,BStrainRate
      use Constants1, only: tiny,twothird
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use WallDistance1, only: WallDistance,BWallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j
      character*10 Variable
      double precision :: term1,term2,term11,term12,term13,
     *                    term21,term22,term23,term31,term32,term33,sum,
     *                    uWall1,vWall1,wWall1,uWall,vWall,wWall,
     *                    WallVelocity,WallVelocity2,
     *                    cmu25k,dfidxTf,dfidyTf,dfidzTf,
     *                    dNorm,nx,ny,nz,ProductionTemp,gamf,
     *                    FluxCflocal,FluxFflocal,FluxVflocal,tke1,uTau,
     *                    dotproduct,phibp,betabp,phinat,betanat,
     *                    SourceDT,fNATcrit,ReOmega,Rev,Fonlim,PkLim,
     *                    Glog,GVisc,factor
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------
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
      CalculateTKESources: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(1) CalculateTKESources           !kepsilon
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
c
            WallVelocity=TangentialVelocity(i1)
            WallVelocity2=WallVelocity*WallVelocity
c
            BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *        duplusdyplus(i)*((ustar(i)/uplus(i))**2)*WallVelocity2/
     *            (Bviscosity(i2,i3))
c
            FluxTE(i1)=FluxTE(i1)+TurbulenceProduction(i1)*Volume(i1)
c
            if(LimitTurbulenceProduction) then
c
              ProductionTemp=c1limiter*BDensity(i2,i3)*
     *                              dmax1(BTurbulentED(i2,i3),0.)
              BTurbulenceProduction(i2,i3)=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
            endif
c            
            TurbulenceProduction(i1)=BTurbulenceProduction(i2,i3)
c
            FluxTE(i1)=FluxTE(i1)-TurbulenceProduction(i1)*Volume(i1)
c
            BTurbulentKE(i2,i3)=TurbulentKE(i1)

          enddo
c
c-----------------------------------------------------------------------------
        case(2) CalculateTKESources           !kepsilonchien
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+LTKE(i)*Volume(i)
            FluxTE(i)=FluxTE(i)+
     *           LTKE(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Modify turbulence production along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(3) CalculateTKESources           !kepsilonsharma
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c
            FluxCE(i)=FluxCE(i)+
     *            LTKE(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)+LTKE(i)*Volume(i)
c        
          enddo
c
c--- Apply boundary condition along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *             dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(4) CalculateTKESources           !kepsilonchc
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Apply boundary condition along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(5) CalculateTKESources           !kepsilonkasagi
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Apply boundary condition along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *              dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(6) CalculateTKESources           !kepsilontagawa
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Apply boundary condition along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(7) CalculateTKESources           !kepsilonhishida
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c
            FluxCE(i)=FluxCE(i)+
     *            LTKE(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)+LTKE(i)*Volume(i)
c        
          enddo
c
c--- Apply boundary condition along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(8) CalculateTKESources           !kelambremhorst
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c
          enddo
c
c--- Apply Boundary Conditions along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(9) CalculateTKESources           !kelambremhorstm
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c
          enddo
c
c--- Apply Boundary Conditions along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(10) CalculateTKESources           !realizable
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
c
            WallVelocity=TangentialVelocity(i1)
            WallVelocity2=WallVelocity*WallVelocity
c
            BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *        duplusdyplus(i)*((ustar(i)/uplus(i))**2)*WallVelocity2/
     *            (Bviscosity(i2,i3))
c
            FluxTE(i1)=FluxTE(i1)+TurbulenceProduction(i1)*Volume(i1)
c
            if(LimitTurbulenceProduction) then
c
              ProductionTemp=c1limiter*BDensity(i2,i3)*
     *                              dmax1(BTurbulentED(i2,i3),0.)
              BTurbulenceProduction(i2,i3)=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
            endif
c            
            TurbulenceProduction(i1)=BTurbulenceProduction(i2,i3)
c
            FluxTE(i1)=FluxTE(i1)-TurbulenceProduction(i1)*Volume(i1)
c
            BTurbulentKE(i2,i3)=TurbulentKE(i1)

          enddo
c
c-----------------------------------------------------------------------------
        case(11) CalculateTKESources           !komega
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=cmu*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *           rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              dNorm=WallDistance(i1)
c
              WallVelocity=TangentialVelocity(i1)
              WallVelocity2=WallVelocity*WallVelocity
c
              FluxTE(i1)=FluxTE(i1)+TurbulenceProduction(i1)*Volume(i1)
c
              if(MomentumWallFunctionType.eq.'automatic') then
c
                factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
                Glog=(BViscosity(i2,i3)+BTurbulentViscosity(i2,i3))*
     *                       ustar(i)*WallVelocity/(cappa*dNorm*dNorm)
                GVisc=BTurbulentViscosity(i2,i3)*
     *                           WallVelocity2/(dNorm*dNorm)
                BTurbulenceProduction(i2,i3)= GVisc*dexp(-factor)+
     *                                            Glog*dexp(-1./factor)
c
              else
c
                BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *                     duplusdyplus(i)*((ustar(i)/uplus(i))**2)*
     *                               WallVelocity2/(Bviscosity(i2,i3))
c
              endif
c
              if(LimitTurbulenceProduction) then
c
                ProductionTemp=c1limiter*cmu*BDensity(i2,i3)*
     *           dmax1(BTurbulentKE(i2,i3)*BTurbulentOmega(i2,i3),0.)
                BTurbulenceProduction(i2,i3)=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
              endif
c            
              TurbulenceProduction(i1)=BTurbulenceProduction(i2,i3)
c
              FluxTE(i1)=FluxTE(i1)-TurbulenceProduction(i1)*Volume(i1)
c
              BTurbulentKE(i2,i3)=TurbulentKE(i1)

            enddo
c
          else
c
            call SetTKEWallBC
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
c              BTurbulentKE(i2,i3)=0.
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTKEGradx(i2,i3)
              dfidyTf=BTKEGrady(i2,i3)
              dfidzTf=BTKEGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(12) CalculateTKESources           !komegaepsilon
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=cmu*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *           rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              dNorm=WallDistance(i1)
c
              WallVelocity=TangentialVelocity(i1)
              WallVelocity2=WallVelocity*WallVelocity
c
              FluxTE(i1)=FluxTE(i1)+TurbulenceProduction(i1)*Volume(i1)
c
              if(MomentumWallFunctionType.eq.'automatic') then
c
                factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
                Glog=(BViscosity(i2,i3)+BTurbulentViscosity(i2,i3))*
     *                       ustar(i)*WallVelocity/(cappa*dNorm*dNorm)
                GVisc=BTurbulentViscosity(i2,i3)*
     *                           WallVelocity2/(dNorm*dNorm)
                BTurbulenceProduction(i2,i3)= GVisc*dexp(-factor)+
     *                                            Glog*dexp(-1./factor)
c
              else
c
                BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *                     duplusdyplus(i)*((ustar(i)/uplus(i))**2)*
     *                               WallVelocity2/(Bviscosity(i2,i3))
c
              endif
c
              if(LimitTurbulenceProduction) then
c
                ProductionTemp=c1limiter*cmu*BDensity(i2,i3)*
     *           dmax1(BTurbulentKE(i2,i3)*BTurbulentOmega(i2,i3),0.)
                BTurbulenceProduction(i2,i3)=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
              endif
c            
              TurbulenceProduction(i1)=BTurbulenceProduction(i2,i3)
c
              FluxTE(i1)=FluxTE(i1)-TurbulenceProduction(i1)*Volume(i1)
c
              BTurbulentKE(i2,i3)=TurbulentKE(i1)

            enddo
c
          else
c
            call SetTKEWallBC
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
c              BTurbulentKE(i2,i3)=0.
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTKEGradx(i2,i3)
              dfidyTf=BTKEGrady(i2,i3)
              dfidzTf=BTKEGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(13) CalculateTKESources           !komegabsl
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=cmu*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *                   rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              dNorm=WallDistance(i1)
c
              WallVelocity=TangentialVelocity(i1)
              WallVelocity2=WallVelocity*WallVelocity
c
              FluxTE(i1)=FluxTE(i1)+TurbulenceProduction(i1)*Volume(i1)
c
              if(MomentumWallFunctionType.eq.'automatic') then
c
                factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
                Glog=(BViscosity(i2,i3)+BTurbulentViscosity(i2,i3))*
     *                       ustar(i)*WallVelocity/(cappa*dNorm*dNorm)
                GVisc=BTurbulentViscosity(i2,i3)*
     *                           WallVelocity2/(dNorm*dNorm)
                BTurbulenceProduction(i2,i3)= GVisc*dexp(-factor)+
     *                                            Glog*dexp(-1./factor)
c
              else
c
                BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *                     duplusdyplus(i)*((ustar(i)/uplus(i))**2)*
     *                               WallVelocity2/(Bviscosity(i2,i3))
c
              endif
c
              if(LimitTurbulenceProduction) then
c
                ProductionTemp=c1limiter*cmu*BDensity(i2,i3)*
     *           dmax1(BTurbulentKE(i2,i3)*BTurbulentOmega(i2,i3),0.)
                BTurbulenceProduction(i2,i3)=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
              endif
c            
              TurbulenceProduction(i1)=BTurbulenceProduction(i2,i3)
c
              FluxTE(i1)=FluxTE(i1)-TurbulenceProduction(i1)*Volume(i1)
c
              BTurbulentKE(i2,i3)=TurbulentKE(i1)

            enddo
c
          else
c
            call SetTKEWallBC
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
c              BTurbulentKE(i2,i3)=0.
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTKEGradx(i2,i3)
              dfidyTf=BTKEGrady(i2,i3)
              dfidzTf=BTKEGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(14) CalculateTKESources           !komegasst
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=cmu*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            ProductionTemp=fr1Coefficient(i)*TurbulenceProduction(i)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-ProductionTemp*Volume(i)+
     *                    rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              dNorm=WallDistance(i1)
c
              WallVelocity=TangentialVelocity(i1)
              WallVelocity2=WallVelocity*WallVelocity
c
              ProductionTemp=fr1Coefficient(i1)*TurbulenceProduction(i1)

              FluxTE(i1)=FluxTE(i1)+ProductionTemp*Volume(i1)
c
              if(MomentumWallFunctionType.eq.'automatic') then
c
                factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
                Glog=(BViscosity(i2,i3)+BTurbulentViscosity(i2,i3))*
     *                       ustar(i)*WallVelocity/(cappa*dNorm*dNorm)
                GVisc=BTurbulentViscosity(i2,i3)*
     *                           WallVelocity2/(dNorm*dNorm)
                BTurbulenceProduction(i2,i3)= GVisc*dexp(-factor)+
     *                                            Glog*dexp(-1./factor)
c
              else
c
                BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *                     duplusdyplus(i)*((ustar(i)/uplus(i))**2)*
     *                               WallVelocity2/(Bviscosity(i2,i3))
c
              endif
c
              ProductionTemp=
     *          c1limiter*cmu*BDensity(i2,i3)*fr1Coefficient(i1)*
     *           dmax1(BTurbulentKE(i2,i3)*BTurbulentOmega(i2,i3),0.)
              ProductionTemp=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
              TurbulenceProduction(i1)=ProductionTemp
c
              FluxTE(i1)=FluxTE(i1)-ProductionTemp*Volume(i1)
c
              BTurbulentKE(i2,i3)=TurbulentKE(i1)

            enddo
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
c--- Apply boundary conditions along walls
c
            call SetTKEWallBC
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
c              BTurbulentKE(i2,i3)=0.
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTKEGradx(i2,i3)
              dfidyTf=BTKEGrady(i2,i3)
              dfidzTf=BTKEGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(15) CalculateTKESources           !sstgamaretheta
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=dmin1(dmax1(TGammaEff(i),0.1),1.)*
     *                    cmu*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            ProductionTemp=fr1Coefficient(i)*TGammaEff(i)*
     *                                     TurbulenceProduction(i)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-ProductionTemp*Volume(i)+
     *                    rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Modify turbulence production along walls
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
            call SetTKEWallBC
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
c              BTurbulentKE(i2,i3)=0.
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTKEGradx(i2,i3)
              dfidyTf=BTKEGrady(i2,i3)
              dfidzTf=BTKEGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(16) CalculateTKESources           !komega2006
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=cmu*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *           rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              dNorm=WallDistance(i1)
c
              WallVelocity=TangentialVelocity(i1)
              WallVelocity2=WallVelocity*WallVelocity
c
              FluxTE(i1)=FluxTE(i1)+TurbulenceProduction(i1)*Volume(i1)
c
              if(MomentumWallFunctionType.eq.'automatic') then
c
                factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
                Glog=(BViscosity(i2,i3)+BTurbulentViscosity(i2,i3))*
     *                       ustar(i)*WallVelocity/(cappa*dNorm*dNorm)
                GVisc=BTurbulentViscosity(i2,i3)*
     *                           WallVelocity2/(dNorm*dNorm)
                BTurbulenceProduction(i2,i3)= GVisc*dexp(-factor)+
     *                                            Glog*dexp(-1./factor)
c
              else
c
                BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *                     duplusdyplus(i)*((ustar(i)/uplus(i))**2)*
     *                               WallVelocity2/(Bviscosity(i2,i3))
c
              endif
c
              if(LimitTurbulenceProduction) then
c
                ProductionTemp=c1limiter*cmu*BDensity(i2,i3)*
     *           dmax1(BTurbulentKE(i2,i3)*BTurbulentOmega(i2,i3),0.)
                BTurbulenceProduction(i2,i3)=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
              endif
c            
              TurbulenceProduction(i1)=BTurbulenceProduction(i2,i3)
c
              FluxTE(i1)=FluxTE(i1)-TurbulenceProduction(i1)*Volume(i1)
c
              BTurbulentKE(i2,i3)=TurbulentKE(i1)
c
            enddo
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
c
c--- Apply boundary conditions
c
            call SetTKEWallBC
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
c              BTurbulentKE(i2,i3)=0.
c
              gamf=BeDiffCoefficient(i2,i3)
              dfidxTf=BTKEGradx(i2,i3)
              dfidyTf=BTKEGrady(i2,i3)
              dfidzTf=BTKEGradz(i2,i3)
c
              FluxCflocal= gamf*BgDiff(i2,i3)
              FluxFflocal=-gamf*BgDiff(i2,i3)
              FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(17) CalculateTKESources           !komega2006lrn
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=
     *          bettaStar(i)*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *           rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),0.)
c        
          enddo
c
c--- Apply boundary conditions
c
          call SetTKEWallBC
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
c            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *             dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(18) CalculateTKESources           !kklmodel
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            tke1=dmax1(TurbulentKE(i),0.)
            term1=cmu75*Density(i)*(tke1**1.5)/
     *                     dmax1(TurbulentKL(i),tiny)+
     *                        2.*Viscosity(i)/(WallDistance(i)**2)
c
            FluxCE(i)=FluxCE(i)+term1*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *                             term1*Volume(i)*tke1
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(21) CalculateTKESources           !kklomega
c-----------------------------------------------------------------------------
c
c--- Calculate turbulence production
c
          do i=1,NumberOfElements    
c
            ProductionKT(i)=TurbulentViscosityTs(i)*
     *                                 StrainRate(i)*StrainRate(i)
            sqrtTurbulentKE(i)=dsqrt(dmax1(TurbulentKE(i),0.))
c     
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BsqrtTurbulentKE(i,j)=dsqrt(dmax1(BTurbulentKE(i,j),0.))
C     
            enddo
          enddo
c
          Variable='TKE05'
c
          call Gradient(Variable,2,sqrtTurbulentKE,sqrtTKEGradx,
     *              sqrtTKEGrady,sqrtTKEGradz,BsqrtTurbulentKE,
     *         BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,2,.false.,1)
c        
          do i=1,NumberOfElements    
c
            phibp=dmax1(Density(i)*dmax1(TurbulentKE(i),0.)/
     *          (Viscosity(i)*dmax1(Vorticity(i),tiny))-Cbpcrit,0.)
            betabp=1.-dexp(-phibp/Abp)
            SourceRbp(i)=Cr*betabp*
     *        dmax1(TurbulentOmega(i),0.)/dmax1(CoefficientFW(i),tiny)
c
            ReOmega=Density(i)*WallDistance(i)*
     *                    WallDistance(i)*Vorticity(i)/Viscosity(i)
            fNATcrit=1.-dexp(-Cnc*dsqrt(dmax1(TurbulentKL(i),0.))*
     *                        WallDistance(i)*Density(i)/Viscosity(i))
            phinat=dmax1(ReOmega-Cnatcrit/fNATcrit,0.)
            betanat=1.-dexp(-phinat/Anat)
            SourceRnat(i)=Crnat*betanat*Vorticity(i)
c
          enddo
c        
          do i=1,NumberOfElements    
c
            SourceDT=2.*(Viscosity(i)/Density(i))*
     *       (sqrtTKEGradx(i)**2+sqrtTKEGrady(i)**2+sqrtTKEGradz(i)**2)
c
            FluxCE(i)=FluxCE(i)+
     *          Density(i)*(dmax1(TurbulentOmega(i),0.)+
     *               SourceDT/dmax1(TurbulentKE(i),tiny))*Volume(i)
            FluxTE(i)=FluxTE(i)-Density(i)*(ProductionKT(i)+
     *                 dmax1(TurbulentKL(i),0.)*(SourceRbp(i)+
     *                               SourceRnat(i)))*Volume(i)+
     *                  Density(i)*(dmax1(TurbulentOmega(i),0.)*
     *                   dmax1(TurbulentKE(i),0.)+SourceDT)*Volume(i)
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *           dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c 
c-----------------------------------------------------------------------------
        case(22) CalculateTKESources           !kepsilonrt
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *                  Density(i)*dmax1(TurbulentED(i),0.)*Volume(i)
c        
          enddo
c
c--- Modify turbulence production along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(23) CalculateTKESources           !sstgama
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            TurbulenceProduction(i)=dmax1(TurbulentViscosity(i)*
     *                               StrainRate(i)*Vorticity(i),0.)
            rhoTED(i)=cmu*Density(i)*dmax1(TurbulentOmega(i),0.)
c
            ProductionTemp=fr1Coefficient(i)*TurbulenceProduction(i)*
     *                 TGamma(i)
c
            Rev=Density(i)*WallDistance(i)*WallDistance(i)*
     *                               StrainRate(i)/Viscosity(i)
            Fonlim=dmin1(dmax1(Rev/(2.2*Reoclim)-1.d0,0.d0),3.d0)
            PkLim=5.*Ck*dmax1(TGamma(i)-0.2,0.)*(1.-Tgamma(i))*Fonlim*
     *            dmax1(3.*Csep*Viscosity(i)-TurbulentViscosity(i),0.)*
     *                      StrainRate(i)*Vorticity(i)
            ProductionTemp=ProductionTemp+PkLim
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)*dmax1(TGamma(i),0.1)
            FluxTE(i)=FluxTE(i)-ProductionTemp*Volume(i)+rhoTED(i)*
     *           Volume(i)*dmax1(TurbulentKE(i),0.)*dmax1(TGamma(i),0.1)
c        
          enddo
c
c--- Apply boundary conditions along walls
c
          call SetTKEWallBC
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
c            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(25) CalculateTKESources           !kepsilonrng
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Modify turbulence production along walls
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              dNorm=WallDistance(i1)
c
              WallVelocity=TangentialVelocity(i1)
              WallVelocity2=WallVelocity*WallVelocity
c
              BTurbulenceProduction(i2,i3)=(BDensity(i2,i3)**2)*
     *          duplusdyplus(i)*((ustar(i)/uplus(i))**2)*WallVelocity2/
     *            (Bviscosity(i2,i3))
c
              FluxTE(i1)=FluxTE(i1)+TurbulenceProduction(i1)*Volume(i1)
c
              if(LimitTurbulenceProduction) then
c
                ProductionTemp=c1limiter*BDensity(i2,i3)*
     *                              dmax1(BTurbulentED(i2,i3),0.)
                BTurbulenceProduction(i2,i3)=
     *           dmin1(BTurbulenceProduction(i2,i3),ProductionTemp)
c        
              endif
c            
              TurbulenceProduction(i1)=BTurbulenceProduction(i2,i3)
c
              FluxTE(i1)=FluxTE(i1)-TurbulenceProduction(i1)*Volume(i1)
c
              BTurbulentKE(i2,i3)=TurbulentKE(i1)

            enddo
c
          elseif(WallTreatment.eq.'lowreynoldsnumber') then
          endif
c
c-----------------------------------------------------------------------------
        case(26) CalculateTKESources           !kepsilonv2f
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Modify turbulence production along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(27) CalculateTKESources           !kepsilonzetaf
c-----------------------------------------------------------------------------
c
          call CalculateTurbulenceProduction
c
          do i=1,NumberOfElements    
c
            rhoTED(i)=Density(i)*dmax1(TurbulentED(i),0.)/
     *                                   dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+rhoTED(i)*Volume(i)
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)+
     *         rhoTED(i)*Volume(i)*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c--- Modify turbulence production along walls
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
            BTurbulentKE(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTKEGradx(i2,i3)
            dfidyTf=BTKEGrady(i2,i3)
            dfidzTf=BTKEGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKE(i1)+
     *                      FluxFflocal*BTurbulentKE(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTKESources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end