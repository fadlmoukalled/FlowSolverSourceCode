c
c#############################################################################################
c
      SUBROUTINE AssembleSpalartAllmarasSources
c
C#############################################################################################
c
      use User0, only: Lcompressible,WallTreatment,LRough,GrainSize,
     *                 LNegativeSpalartAllmaras,LTransitionalSA,
     *                 MethodCalcGradientDensity,nIterGradientDensity,
     *                 LimitGradientDensity,LimitGradientDensityMethod,
     *                 LSpalartCompressibilityCorrection
      use Geometry1, only: NumberOfElements
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BFaceAreanx,BFaceAreany,
     *                     BDistanceCFx,BDistanceCFy,BgDiff,
     *                     BFaceTx,BFaceTy,BFaceTz
      use Variables1, only: uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      TurbulenceProduction,BTurbulenceProduction,
     *                      Temperature,ModifiedED,ModifiedEDGradx,
     *                      ModifiedEDGrady,ModifiedEDGradz,BModifiedED,
     *                      BModifiedEDGradx,BModifiedEDGrady,
     *                      BModifiedEDGradz,TGamma
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density,SpecificHeat,RGas,
     *                               BeDiffCoefficient,Viscosity,
     *                               BViscosity,BturbulentViscosity,
     *                               BDensity,eDiffCoefficient,
     *                               DensGradx,DensGrady,DensGradz,
     *                               BDensGradx,BDensGrady,BDensGradz
      use Turbulence1, only: cb1,cb2,ct3,cw1,sigMED,c5,cappa,Stelda,
     *                       fwCoefficient,ft2Coefficient,uplus,yplus,
     *                       Vorticity,BVorticity,fr1Coefficient
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use WallDistance1, only: WallDistance,iTau
      use constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,i1,i2,i3,i4,j
      double precision :: term1,term2,Turbulencedestruction,gamf,
     *                    dfidxTf,dfidyTf,dfidzTf,FluxCflocal,
     *                    FluxFflocal,FluxVflocal,dNorm
c********************************************************************************************
c--------------------------------------------------------------
c--- Interfaces
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------------------
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
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------------------------
c
      if(.not.LNegativeSpalartAllmaras) then
c
c--- Calculate turbulence production
c
        if(LTransitionalSA) then
c
          do i=1,NumberOfElements     
c
            TurbulenceProduction(i)=cb1*Density(i)*fr1Coefficient(i)*
     *          (1.-ft2Coefficient(i))*Stelda(i)*ModifiedED(i)*TGamma(i)
            TurbulenceProduction(i)=dmax1(TurbulenceProduction(i),0.)
c
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)
c
          enddo
c
        else
c
          do i=1,NumberOfElements     
c
            TurbulenceProduction(i)=cb1*Density(i)*fr1Coefficient(i)*
     *                    (1.-ft2Coefficient(i))*Stelda(i)*ModifiedED(i)
            TurbulenceProduction(i)=dmax1(TurbulenceProduction(i),0.)
c
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)
c
          enddo
c
        endif

        do i=1,NumberOfElements    
c
          dNorm=WallDistance(i)
          if(LRough.and.WallTreatment.eq.'lowreynoldsnumber') 
     *                              dNorm=dNorm+0.03*GrainSize(iTau(i))
          term1=cw1*fwCoefficient(i)-cb1*ft2Coefficient(i)/(cappa*cappa)
          Turbulencedestruction=Density(i)*term1*
     *                   dmax1(ModifiedED(i),0.)/(dNorm**2)
        
c
          FluxCE(i)=FluxCE(i)+Turbulencedestruction*Volume(i)
          FluxTE(i)=FluxTE(i)+Turbulencedestruction*Volume(i)*
     *                                    dmax1(ModifiedED(i),0.)
c
        enddo
c
        do i=1,NumberOfElements    
c
          term1=(ModifiedEDGradx(i)**2+
     *            ModifiedEDGrady(i)**2+ModifiedEDGradz(i)**2)
          term1=term1*Density(i)*Volume(i)*cb2/sigMED 
c
          FluxTE(i)=FluxTE(i)-term1
c
        enddo
c
        if(Lcompressible.and..not.LSpalartCompressibilityCorrection)then
c
c---- Calculate Density Gradients
c
          Variable='density'
          call Gradient(Variable,MethodCalcGradientDensity,
     *      Density,DensGradx,DensGrady,DensGradz,
     *        BDensity,BDensGradx,BDensGrady,
     *          BDensGradz,nIterGradientDensity,
     *            LimitGradientDensity,LimitGradientDensityMethod)
c
c
          do i=1,NumberOfElements    
c
            term1=(DensGradx(i)*ModifiedEDGradx(i)+
     *              DensGrady(i)*ModifiedEDGrady(i)+
     *                 DensGradz(i)*ModifiedEDGradz(i))
            term1=term1*Volume(i)*eDiffCoefficient(i)/Density(i)
c
            FluxCE(i)=FluxCE(i)+
     *                   dmax1(term1,0.)/dmax1(ModifiedED(i),tiny)
            FluxTE(i)=FluxTE(i)+term1
c
          enddo
c
        elseif(Lcompressible.and.LSpalartCompressibilityCorrection) then
c
          do i=1,NumberOfElements    
c
            term1=(uVelGradx(i)**2+uVelGrady(i)**2+uVelGradz(i)**2+
     *       vVelGradx(i)**2+vVelGrady(i)**2+vVelGradz(i)**2+
     *       wVelGradx(i)**2+wVelGrady(i)**2+wVelGradz(i)**2)*Density(i)
            term2=SpecificHeat(i)/(SpecificHeat(i)-RGas)
            term2=term2*RGas*Temperature(i)
            term1=c5*term1*dmax1(ModifiedED(i),0.)*Volume(i)/term2
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(ModifiedED(i),0.)
c        
          enddo
c 
        endif
c      
      elseif(LNegativeSpalartAllmaras) then
c
c--- Calculate turbulence production
c
        do i=1,NumberOfElements     
c
          if(ModifiedED(i).ge.0.) then
c
            TurbulenceProduction(i)=cb1*Density(i)*
     *            (1.-ft2Coefficient(i))*Stelda(i)*ModifiedED(i)
            TurbulenceProduction(i)=dmax1(TurbulenceProduction(i),0.)
c
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)
c
          else
c
            TurbulenceProduction(i)=cb1*Density(i)*
     *            (1.-ct3)*Vorticity(i)*ModifiedED(i)
            TurbulenceProduction(i)=dmax1(TurbulenceProduction(i),0.)
c
            FluxTE(i)=FluxTE(i)-TurbulenceProduction(i)*Volume(i)
c
          endif
c
        enddo
c
        do i=1,NumberOfElements    
c
          if(ModifiedED(i).ge.0.) then
c
            term1=cw1*fwCoefficient(i)-
     *                cb1*ft2Coefficient(i)/(cappa*cappa)
            Turbulencedestruction=Density(i)*term1*
     *                   dmax1(ModifiedED(i),0.)/(WallDistance(i)**2)
        
c
            FluxCE(i)=FluxCE(i)+Turbulencedestruction*Volume(i)
            FluxTE(i)=FluxTE(i)+Turbulencedestruction*Volume(i)*
     *                                    dmax1(ModifiedED(i),0.)
c
          else 
c
            Turbulencedestruction=Density(i)*cw1*
     *                   ((ModifiedED(i)/WallDistance(i))**2)
c
            FluxTE(i)=FluxTE(i)-Turbulencedestruction*Volume(i)
     *                                    
c
          endif
c
        enddo
c
        do i=1,NumberOfElements    
c
          term1=(ModifiedEDGradx(i)**2+
     *            ModifiedEDGrady(i)**2+ModifiedEDGradz(i)**2)
          term1=term1*Density(i)*Volume(i)*cb2/sigMED 
c
          FluxTE(i)=FluxTE(i)-term1
c
        enddo
c
        if(Lcompressible.and..not.LSpalartCompressibilityCorrection)then
c
c
c---- Calculate Density Gradients
c
          call Gradient(Variable,MethodCalcGradientDensity,
     *      Density,DensGradx,DensGrady,DensGradz,
     *        BDensity,BDensGradx,BDensGrady,
     *          BDensGradz,nIterGradientDensity,
     *            LimitGradientDensity,LimitGradientDensityMethod)
c
c
          do i=1,NumberOfElements    
c
            if(ModifiedED(i).ge.0.) then
c
              term1=(DensGradx(i)*ModifiedEDGradx(i)+
     *              DensGrady(i)*ModifiedEDGrady(i)+
     *                 DensGradz(i)*ModifiedEDGradz(i))
              term1=term1*Volume(i)*eDiffCoefficient(i)/Density(i)
c
              FluxCE(i)=FluxCE(i)+
     *                   dmax1(term1,0.)/dmax1(ModifiedED(i),tiny)
              FluxTE(i)=FluxTE(i)+term1
c
            else
c
              term1=-(DensGradx(i)*ModifiedEDGradx(i)+
     *              DensGrady(i)*ModifiedEDGrady(i)+
     *                 DensGradz(i)*ModifiedEDGradz(i))
              term1=term1*Volume(i)*Viscosity(i)/Density(i)
c
              FluxTE(i)=FluxTE(i)-term1
c
            endif
          enddo
c
        elseif(Lcompressible.and.LSpalartCompressibilityCorrection) then
c
          do i=1,NumberOfElements    
c
            term1=(uVelGradx(i)**2+uVelGrady(i)**2+uVelGradz(i)**2+
     *       vVelGradx(i)**2+vVelGrady(i)**2+vVelGradz(i)**2+
     *       wVelGradx(i)**2+wVelGrady(i)**2+wVelGradz(i)**2)*Density(i)
            term2=SpecificHeat(i)/(SpecificHeat(i)-RGas)
            term2=term2*RGas*Temperature(i)
            term1=c5*term1*dmax1(ModifiedED(i),0.)*Volume(i)/term2
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(ModifiedED(i),0.)
c        
          enddo
c 
        endif
c 
      endif
c
      if(WallTreatment.eq.'lowreynoldsnumber') then
c      
        if(LRough) then
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
            term1=0.03*GrainSize(i2)/
     *                (WallDistance(i1)+0.06*GrainSize(i2))
c            term1=0.03*GrainSize(i2)/
c     *                (WallDistance(i1)+0.03*GrainSize(i2))

            BModifiedED(i2,i3)=ModifiedED(i1)*term1
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BModifiedEDGradx(i2,i3)
            dfidyTf=BModifiedEDGrady(i2,i3)
            dfidzTf=BModifiedEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*ModifiedED(i1)+
     *                      FluxFflocal*BModifiedED(i2,i3)+FluxVflocal
c
          enddo
c
        else
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
            BModifiedED(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BModifiedEDGradx(i2,i3)
            dfidyTf=BModifiedEDGrady(i2,i3)
            dfidzTf=BModifiedEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*ModifiedED(i1)+
     *                      FluxFflocal*BModifiedED(i2,i3)+FluxVflocal
c
          enddo
c
        endif
c
      elseif(WallTreatment.eq.'wallfunctions') then
c      
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          BModifiedED(i2,i3)=ModifiedED(i1)
c
        enddo            
c
      endif
c
      return
      end