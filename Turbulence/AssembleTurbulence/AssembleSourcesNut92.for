c
c#############################################################################################
c
      SUBROUTINE AssembleNut92Sources
c
C#############################################################################################
c
      use User0, only: WallTreatment,LRough,GrainSize
      use Geometry1, only: NumberOfElements
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BgDiff,BFaceTx,BFaceTy,BFaceTz
      use Variables1, only: ModifiedED,BModifiedED,ModifiedEDGradx,
     *                      ModifiedEDGrady,ModifiedEDGradz,
     *                      BModifiedEDGradx,BModifiedEDGrady,
     *                      BModifiedEDGradz,
     *                      uVelGradx,vVelGrady,wVelGradz
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BeDiffCoefficient
      use Turbulence1, only: F1Nut92,F2Nut92,G1Nut92,G2Nut92,N1Nut92,
     *                       N2Nut92,ModifiedMutGradx,ModifiedMutGrady,
     *                       ModifiedMutGradz,ModifiedEDGrad2x,
     *                       ModifiedEDGrad2y,ModifiedEDGrad2z,Anut1,
     *                       Anut2,Cnut2,Cnut3,Cnut4,Cnut5,Cnut6,Cnut7,
     *                       ModifiedDensGrad2x,ModifiedDensGrad2y,
     *                       ModifiedDensGrad2z,uTau
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use WallDistance1, only: WallDistance,iTau,jTau
      use Constants1, only: tiny,twothird,onethird
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j
      double precision :: term1,term2,term3,term4,term5,term6,
     *                    nutLocal,nuLocal,Production,aVel,dNorm,
     *                    nutWall,destruction,gamf,
     *                    dfidxTf,dfidyTf,dfidzTf,FluxCflocal,
     *                    FluxFflocal,FluxVflocal
c********************************************************************************************
      interface
c********************************************************************************************
        FUNCTION SpeedOfSound(i)
c********************************************************************************************
          integer :: i
          double precision :: SpeedOfSound
c********************************************************************************************
        end FUNCTION SpeedOfSound
c********************************************************************************************
      end interface
c-------------------------------------------------------------------------------------------
c
      call CalculateNut92Coefficients
c
c---- Calculate production
c
      do i=1,NumberOfElements     
c
        nutLocal=ModifiedED(i)
        nuLocal=Viscosity(i)/Density(i)
c
        term1=Density(i)*Cnut2*F2Nut92(i)       
        term2=nutLocal*G1Nut92(i)
        term3=Anut1*(nutLocal**onethird)*
     *                  (G2Nut92(i)**twothird)
        term4=Density(i)*Cnut2*F2Nut92(i)*Anut2*N1Nut92(i)*
     *    dsqrt((nuLocal+nutLocal)*G1Nut92(i))
c
        Production=(term1*term2+term4)*Volume(i)
        FluxTE(i)=FluxTE(i)-Production
c
        Production=term1*term3*Volume(i)
        FLuxCE(i)=FLuxCE(i)+dmax1(-Production,0.)
        FLuxTE(i)=FLuxTE(i)-Production*nutLocal
c
        term5=ModifiedEDGrad2x(i)+
     *             ModifiedEDGrad2y(i)+ModifiedEDGrad2z(i)
        term6=Density(i)*Cnut3*(term5+N2Nut92(i))*Volume(i)
c
        FLuxCE(i)=FLuxCE(i)+dmax1(-term6,0.)
        FLuxTE(i)=FLuxTE(i)-term6*nutLocal
c
      enddo
c
c---- Calculate destruction
c
      do i=1,NumberOfElements     
c
        nutLocal=ModifiedED(i)
        nuLocal=Viscosity(i)/Density(i)
c
        aVel=SpeedOfSound(i)
c
        dNorm=WallDistance(i)
        nutWall=0.
c
        if(LRough) then
c
          dNorm=dNorm+0.01*GrainSize(iTau(i))
c
          i2=iTau(i)
          i3=jTau(i)
c          i1=NBFaceOwner(i2,i3)
c
c          WallVelocity=TangentialVelocity(i1)
c          TauWall=(BTurbulentViscosity(i2,i3)+
c     *                 BViscosity(i2,i3))*WallVelocity/WallDistance(i1)
c
c          uTauWall=dsqrt(TauWall/BDensity(i2,i3))
c          nutWall=0.02*GrainSize(i2)*uTauWall
          nutWall=BModifiedED(i2,i3)
c
        endif
c
        term1=Cnut5*nutLocal*((G1Nut92(i)/aVel)**2)
        term2=uVelGradx(i)+vVelGrady(i)+wVelGradz(i)+
     *        ModifiedDensGrad2x(i)+ModifiedDensGrad2y(i)+
     *                                 ModifiedDensGrad2z(i)
        term3=Cnut4*(term2+dabs(term2))
        term4=(Cnut6*(N1Nut92(i)*WallDistance(i)+nutWall)+
     *              Cnut7*F1Nut92(i)*nuLocal)/(dNorm*dNorm)
        
        destruction=Density(i)*(term1+term3+term4)*Volume(i)
c
        FluxCE(i)=FluxCE(i)+destruction
        FluxTE(i)=FluxTE(i)+destruction*nutLocal
c
      enddo
c
c---- Additional source term
c
      do i=1,NumberOfElements     
c
        nutLocal=dmax1(ModifiedED(i),tiny)
        term1=(ModifiedMutGradx(i)*ModifiedEDGradx(i)+  
     *        ModifiedMutGrady(i)*ModifiedEDGrady(i)+
     *        ModifiedMutGradz(i)*ModifiedEDGradz(i))*Volume(i)   
        FLuxCE(i)=FLuxCE(i)+dmax1(-term1,0.)/nutLocal
        FLuxTE(i)=FLuxTE(i)-term1
c
      enddo
c
      if(WallTreatment.eq.'lowreynoldsnumber') then
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
          if(LRough) BModifiedED(i2,i3)=0.02*GrainSize(i2)*uTau(i)
c
          gamf=BeDiffCoefficient(i2,i3)
          dfidxTf=BModifiedEDGradx(i2,i3)
          dfidyTf=BModifiedEDGrady(i2,i3)
          dfidzTf=BModifiedEDGradz(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                    dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxCflocal*ModifiedED(i1)+
     *                      FluxFflocal*BModifiedED(i2,i3)+FluxVflocal
c
        enddo
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