c
c#############################################################################################
c
      SUBROUTINE CalculateNut92Coefficients
c
c#############################################################################################
      use User0, only: MethodCalcGradientModifiedED,LRough,
     *                 nIterGradientModifiedED,GrainSize,
     *                 LimitGradientModifiedED,
     *                 LimitGradientModifiedEDMethod,
     *                 MethodCalcGradientDensity,nIterGradientDensity,
     *                 LimitGradientDensity,LimitGradientDensityMethod,
     *                 MethodCalcGradientMomentum,nIterGradientMomentum,
     *                 LimitGradientMomentum,LimitGradientMomentumMethod
      use Variables1, only: uVelGradx,uVelGrady,uVelGradz,
     *                      uVelGrad2x,uVelGrad2y,uVelGrad2z,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BuVelGrad2x,BuVelGrad2y,BuVelGrad2z,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      vVelGrad2x,vVelGrad2y,vVelGrad2z,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BvVelGrad2x,BvVelGrad2y,BvVelGrad2z,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      wVelGrad2x,wVelGrad2y,wVelGrad2z,
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
     *                      BwVelGrad2x,BwVelGrad2y,BwVelGrad2z,
     *                      uvVelGradxy,uwVelGradxz,BuvVelGradxy,
     *                      BuwVelGradxz,ModifiedEDGradx,
     *                      ModifiedEDGrady,ModifiedEDGradz,
     *                      BModifiedEDGradx,BModifiedEDGrady,
     *                      BModifiedEDGradz,ModifiedED,BModifiedED
      use Turbulence1, only: ModifiedEDGrad2x,ModifiedEDGrad2y,
     *                       ModifiedEDGrad2z,BModifiedEDGrad2x,
     *                       BModifiedEDGrad2y,BModifiedEDGrad2z,
     *                       ModifiedMut,BModifiedMut,ModifiedMutGradx,
     *                       ModifiedMutGrady,ModifiedMutGradz,
     *                       BModifiedMutGradx,BModifiedMutGrady,
     *                       BModifiedMutGradz,G1Nut92,G2Nut92,N1Nut92,
     *                       BN1Nut92,N2Nut92,F1Nut92,F2Nut92,
     *                       N1Nut92Gradx,N1Nut92Grady,N1Nut92Gradz,
     *                       BN1Nut92Gradx,BN1Nut92Grady,BN1Nut92Gradz,
     *                       Cnut0,Cnut1,Cnut8,sigT,ModifiedDensGradx,
     *                       ModifiedDensGrady,ModifiedDensGradz,
     *                       BModifiedDensGradx,BModifiedDensGrady,
     *                       BModifiedDensGradz,ModifiedDensGrad2x,
     *                       ModifiedDensGrad2y,ModifiedDensGrad2z,
     *                       BModifiedDensGrad2x,BModifiedDensGrad2y,
     *                       BModifiedDensGrad2z
      use Geometry1, only: NumberOfElements,NumberOfBCSets   
      use Geometry3, only: NBFaces,NBFaceOwner
      use PhysicalProperties1, only: Viscosity,TurbulentViscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               Density,BDensity,
     *                               DensGradx,DensGrady,DensGradz,
     *                               BDensGradx,BDensGrady,BDensGradz
      use WallDistance1, only: WallDistance,iTau,jTau
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,i1,i2,i3,j
      double precision :: term1,term2,term3,term4,term5,dNorm,nutWall,
     *                    Xi,WallVelocity,TauWall,uTauWall
c********************************************************************************************
c--- Interfaces
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
c--------------------------------------------------------------------------------
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
c--- Calculate the Laplacian of u, v, and w velocities
c
      Variable='velx'
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGradx,uVelGrad2x,uvVelGradxy,uwVelGradxz,BuVelGradx,
     *       BuVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGrady,uvVelGradxy,uVelGrad2y,uwVelGradxz,BuVelGrady,
     *       BuvVelGradxy,BuVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGradz,uvVelGradxy,uwVelGradxz,uVelGrad2z,BuVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BuVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGradx,vVelGrad2x,uvVelGradxy,uwVelGradxz,BvVelGradx,
     *       BvVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGrady,uvVelGradxy,vVelGrad2y,uwVelGradxz,BvVelGrady,
     *       BuvVelGradxy,BvVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGradz,uvVelGradxy,uwVelGradxz,vVelGrad2z,BvVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BvVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGradx,wVelGrad2x,uvVelGradxy,uwVelGradxz,BwVelGradx,
     *       BwVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                             LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGrady,uvVelGradxy,wVelGrad2y,uwVelGradxz,BwVelGrady,
     *       BuvVelGradxy,BwVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                             LimitGradientMomentumMethod)
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGradz,uvVelGradxy,uwVelGradxz,wVelGrad2z,BwVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BwVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                            LimitGradientMomentumMethod)
c
c--- Calculate the Laplacian of nut
c
      Variable='med'
      call Gradient(Variable,MethodCalcGradientModifiedED,
     *     ModifiedEDGradx,ModifiedEDGrad2x,uvVelGradxy,uwVelGradxz,
     *     BModifiedEDGradx,BModifiedEDGrad2x,BuvVelGradxy,BuwVelGradxz,
     *     nIterGradientModifiedED,LimitGradientModifiedED,
     *                             LimitGradientModifiedEDMethod)
      call Gradient(Variable,MethodCalcGradientModifiedED,
     *     ModifiedEDGrady,uvVelGradxy,ModifiedEDGrad2y,uwVelGradxz,
     *     BModifiedEDGrady,BuvVelGradxy,BModifiedEDGrad2y,BuwVelGradxz,
     *     nIterGradientModifiedED,LimitGradientModifiedED,
     *                             LimitGradientModifiedEDMethod)
      call Gradient(Variable,MethodCalcGradientModifiedED,
     *     ModifiedEDGradz,uvVelGradxy,uwVelGradxz,ModifiedEDGrad2z,
     *     BModifiedEDGradz,BuvVelGradxy,BuwVelGradxz,BModifiedEDGrad2z,
     *     nIterGradientModifiedED,LimitGradientModifiedED,
     *                             LimitGradientModifiedEDMethod)
c
c--- Calculate the gradient of (-mu+(C1-C0)mut)
c
      do i=1,NumberOfElements
c
        ModifiedMut(i)=-Viscosity(i)+(Cnut1-Cnut0)*TurbulentViscosity(i)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BModifiedMut(i,j)=-BViscosity(i,j)+
     *                   (Cnut1-Cnut0)*BTurbulentViscosity(i,j)
c
        enddo
      enddo
c
      Variable='med'
      call Gradient(Variable,MethodCalcGradientModifiedED,
     *      ModifiedMut,ModifiedMutGradx,ModifiedMutGrady,
     *      ModifiedMutGradz,BModifiedMut,
     *       BModifiedMutGradx,BModifiedMutGrady,BModifiedMutGradz,
     *        nIterGradientModifiedED,LimitGradientModifiedED,
     *                              LimitGradientModifiedEDMethod)
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
      do i=1,NumberOfElements  
c
        ModifiedDensGradx(i)=DensGradx(i)*
     *                  ModifiedED(i)/(sigT*Density(i))
        ModifiedDensGrady(i)=DensGrady(i)*
     *                  ModifiedED(i)/(sigT*Density(i))
        ModifiedDensGradz(i)=DensGradz(i)*
     *                  ModifiedED(i)/(sigT*Density(i))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BModifiedDensGradx(i,j)=BDensGradx(i,j)*
     *                  BModifiedED(i,j)/(sigT*BDensity(i,j))
          BModifiedDensGrady(i,j)=BDensGrady(i,j)*
     *                  BModifiedED(i,j)/(sigT*BDensity(i,j))
          BModifiedDensGradz(i,j)=BDensGradz(i,j)*
     *                  BModifiedED(i,j)/(sigT*BDensity(i,j))
c
        enddo
      enddo
c
      call Gradient(Variable,MethodCalcGradientDensity,
     *    ModifiedDensGradx,ModifiedDensGrad2x,uvVelGradxy,uwVelGradxz,
     *        BModifiedDensGradx,BModifiedDensGrad2x,BuvVelGradxy,
     *          BuwVelGradxz,nIterGradientDensity,
     *            LimitGradientDensity,LimitGradientDensityMethod)
      call Gradient(Variable,MethodCalcGradientDensity,
     *    ModifiedDensGrady,uvVelGradxy,ModifiedDensGrad2y,uwVelGradxz,
     *        BModifiedDensGrady,BuvVelGradxy,BModifiedDensGrad2y,
     *          BuwVelGradxz,nIterGradientDensity,
     *            LimitGradientDensity,LimitGradientDensityMethod)
      call Gradient(Variable,MethodCalcGradientDensity,
     *    ModifiedDensGradz,uvVelGradxy,uwVelGradxz,ModifiedDensGrad2z,
     *        BModifiedDensGradz,BuvVelGradxy,BuwVelGradxz,
     *          BModifiedDensGrad2z,nIterGradientDensity,
     *            LimitGradientDensity,LimitGradientDensityMethod)
c
c---- Calculate coefficients
c
      do i=1,NumberOfElements
c
        term1=uVelGradx(i)*uVelGradx(i)+
     *        vVelGrady(i)*vVelGrady(i)+wVelGradz(i)*wVelGradz(i)
        term2=uVelGradx(i)*uVelGradx(i)+
     *        uVelGrady(i)*uVelGrady(i)+uVelGradz(i)*uVelGradz(i)
        term3=vVelGradx(i)*vVelGradx(i)+
     *        vVelGrady(i)*vVelGrady(i)+vVelGradz(i)*vVelGradz(i)
        term4=wVelGradx(i)*wVelGradx(i)+
     *        wVelGrady(i)*wVelGrady(i)+wVelGradz(i)*wVelGradz(i)
        term5=2.*(uVelGrady(i)*vVelGradx(i)+
     *        uVelGradz(i)*wVelgradx(i)+vVelGradz(i)*wVelGrady(i))
        G1Nut92(i)=dsqrt(dmax1(term1+term2+term3+term4+term5,0.))
c
        term1=uVelGrad2x(i)+uVelGrad2y(i)+uVelGrad2z(i)
        term2=vVelGrad2x(i)+vVelGrad2y(i)+vVelGrad2z(i)
        term3=wVelGrad2x(i)+wVelGrad2y(i)+wVelGrad2z(i)
        G2Nut92(i)=dsqrt(term1*term1+term2*term2+term3*term3)
c
        N1Nut92(i)=dsqrt(ModifiedEDGradx(i)*ModifiedEDGradx(i)+
     *                      ModifiedEDGrady(i)*ModifiedEDGrady(i)+
     *                         ModifiedEDGradz(i)*ModifiedEDGradz(i))
c
        nutWall=0.
c
        if(LRough) then
c
         i2=iTau(i)
         i3=jTau(i)
c         i1=NBFaceOwner(i2,i3)
c
c         dNorm=WallDistance(i1)
c
c         WallVelocity=TangentialVelocity(i1)
c         TauWall=(BTurbulentViscosity(i2,i3)+
c     *                 BViscosity(i2,i3))*WallVelocity/dNorm
c
c          uTauWall=dsqrt(TauWall/BDensity(i2,i3))
c          nutWall=0.02*GrainSize(i2)*uTauWall
          nutWall=BModifiedED(i2,i3)
c
        endif
c
        Xi=TurbulentViscosity(i)/(7.*Viscosity(i))
        F2Nut92(i)=(Xi*Xi+1.3*Xi+0.2)/dmax1((Xi*Xi-1.3*Xi+1.),tiny)
        term1=N1Nut92(i)*WallDistance(i)+
     *                  0.4*Cnut8*Viscosity(i)/Density(i)
        term2=ModifiedED(i)+Cnut8*Viscosity(i)/Density(i)+nutWall
        F1Nut92(i)=term1/dmax1(term2,tiny)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
        BN1Nut92(i,j)=
     *       dsqrt(BModifiedEDGradx(i,j)*BModifiedEDGradx(i,j)+
     *                 BModifiedEDGrady(i,j)*BModifiedEDGrady(i,j)+
     *                     BModifiedEDGradz(i,j)*BModifiedEDGradz(i,j))
c
        enddo
      enddo
c
c---- Calculate gradient of N1Nut92
c
      Variable='med'
      call Gradient(Variable,MethodCalcGradientModifiedED,
     *      N1Nut92,N1Nut92Gradx,N1Nut92Grady,N1Nut92Gradz,BN1Nut92,
     *       BN1Nut92Gradx,BN1Nut92Grady,BN1Nut92Gradz,
     *        nIterGradientModifiedED,LimitGradientModifiedED,
     *                              LimitGradientModifiedEDMethod)
c
c---- Calculate N2
c
      do i=1,NumberOfElements
c
        N2Nut92(i)=dsqrt(N1Nut92Gradx(i)*N1Nut92Gradx(i)+
     *                      N1Nut92Grady(i)*N1Nut92Grady(i)+
     *                            N1Nut92Gradz(i)*N1Nut92Gradz(i))
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculatePseudoEddyViscosityCoefficients
c
c#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Turbulence1, only: Vorticity,BVorticity,fr1Coefficient,
     *                       StrainRate,BStrainRate,
     *                       cmu,cmu25,AmuRT,alphaRT,BettaRT,phiRT,
     *                       cr1,cr2,cr3,
     *                       S11,S12,S13,S22,S23,S33,
     *                       BS11,BS12,BS13,BS22,BS23,BS33,
     *                       S11Old,S12Old,S13Old,S22Old,S23Old,S33Old,
     *                       S11OldOld,S12OldOld,S13OldOld,S22OldOld,
     *                       S23OldOld,S33OldOld,DS11Dt,DS12Dt,DS13Dt,
     *                       DS22Dt,DS23Dt,DS33Dt,
     *                       W11,W12,W13,W22,W23,W33,
     *                       ReT,BReT,f2RT,Bf2RT,fmuKE,BfmuKE,
     *                       fmuRT,BfmuRT,TurbulentViscosityKE,
     *                       BTurbulentViscosityKE,TurbulentViscosityRT,
     *                       BTurbulentViscosityRT
      use Variables1, only: ModifiedED,BModifiedED,MaterialDerivative,
     *                      TurbulentKE,BTurbulentKE,
     *                      TurbulentED,BTurbulentED
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use Constants1, only: tiny,sqrt2
      use User0, only: LKEpsilonRtRotationCurvatureCorrection,
     *                 LRotation,RotationCurvatureMethod
      use Coriolis1, only: AngularVelocityX,AngularVelocityY,
     *                     AngularVelocityZ
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,k
      double precision :: term1,term2,term3,XiLocal,
     *                    rstar,rtelda,D4,
     *                    a1,a2,a3,a4,a5,a6,a7,a8,a9
c
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE CalculateMaterialDerivative(Variable,FiTemp,
     *                                 BFiTemp,FiTempOld,FiTempOldOld)
c--------------------------------------------------------------------------------
          character*10 Variable
          double precision FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
          double precision, dimension(:) :: FiTemp
          double precision, dimension(:) :: FiTempOld
          double precision, dimension(:) :: FiTempOldOld
          double precision, dimension(:,:) :: BFiTemp
c--------------------------------------------------------------------------------
        end SUBROUTINE CalculateMaterialDerivative
c--------------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------------
c
      do i=1,NumberOfElements
c
        ReT(i)=Density(i)*TurbulentKE(i)*TurbulentKE(i)/
     *                (Viscosity(i)*dmax1(TurbulentED(i),tiny))
        XiLocal=Density(i)*dmax1(ModifiedED(i),0.)/(cmu*Viscosity(i))
c
        f2RT(i)=dtanh(dsqrt(XiLocal))/
     *                      dmax1(dtanh(sqrt2*cmu25*XiLocal),tiny)
c
        term1=1.-dexp(-AmuRT*ReT(i))
        term2=1.-dexp(-dsqrt(ReT(i)))
        term3=dmax1(1.,dsqrt(2./dmax1(ReT(i),tiny)))
        fmuKE(i)=term1*term3/dmax1(term2,tiny)
c
        term1=alphaRT*XiLocal*XiLocal   
        term2=BettaRT*XiLocal*XiLocal   
        fmuRT(i)=dtanh(term1)/dmax1(dtanh(term2),tiny)
c
        term1=cmu*fmuKE(i)*Density(i)*TurbulentKE(i)*
     *                TurbulentKE(i)/dmax1(TurbulentED(i),tiny)
        term2=phiRT*Density(i)*
     *          dmax1(TurbulentKE(i),0.)/dmax1(StrainRate(i),tiny)
        TurbulentViscosityKE(i)=dmin1(term1,term2)
c
        TurbulentViscosityRT(i)=
     *           fmuRT(i)*Density(i)*dmax1(ModifiedED(i),0.)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)*BTurbulentKE(i,j)/
     *                  (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
          XiLocal=BDensity(i,j)*
     *                 dmax1(BModifiedED(i,j),0.)/(cmu*BViscosity(i,j))
c
          Bf2RT(i,j)=dtanh(dsqrt(XiLocal))/
     *                      dmax1(dtanh(sqrt2*cmu25*XiLocal),tiny)
c
          term1=1.-dexp(-AmuRT*BReT(i,j))
          term2=1.-dexp(-dsqrt(BReT(i,j)))
          term3=dmax1(1.,dsqrt(2./dmax1(BReT(i,j),tiny)))
          BfmuKE(i,j)=term1*term3/dmax1(term2,tiny)
c
          term1=alphaRT*XiLocal*XiLocal   
          term2=BettaRT*XiLocal*XiLocal   
          BfmuRT(i,j)=dtanh(term1)/dmax1(dtanh(term2),tiny)
c
          term1=cmu*BfmuKE(i,j)*BDensity(i,j)*BTurbulentKE(i,j)*
     *                BTurbulentKE(i,j)/dmax1(BTurbulentED(i,j),tiny)
          term2=phiRT*BDensity(i,j)*
     *          dmax1(BTurbulentKE(i,j),0.)/dmax1(BStrainRate(i,j),tiny)
          BTurbulentViscosityKE(i,j)=dmin1(term1,term2)
c
          BTurbulentViscosityRT(i,j)=
     *           BfmuRT(i,j)*BDensity(i,j)*dmax1(BModifiedED(i,j),0.)
c
        enddo
      enddo
c
c--- Calculate the rotation/curvature correction term
c
      if(LKEpsilonRtRotationCurvatureCorrection) then
c
        call CalculateVorticitytensor
c
        if(RotationCurvatureMethod.eq.'spalartshur') then
c
          allocate(MaterialDerivative(NumberOfElements))
          Variable='S11'
          call CalculateMaterialDerivative(Variable,S11,
     *                               BS11,S11Old,S11OldOld)
          DS11Dt=MaterialDerivative
          Variable='S12'
          call CalculateMaterialDerivative(Variable,S12,
     *                               BS12,S12Old,S12OldOld)
          DS12Dt=MaterialDerivative
          Variable='S13'
          call CalculateMaterialDerivative(Variable,S13,
     *                               BS13,S13Old,S13OldOld)
          DS13Dt=MaterialDerivative
          Variable='S22'
          call CalculateMaterialDerivative(Variable,S22,
     *                               BS22,S22Old,S22OldOld)
          DS22Dt=MaterialDerivative
          Variable='S23'
          call CalculateMaterialDerivative(Variable,S23,
     *                               BS23,S23Old,S23OldOld)
          DS23Dt=MaterialDerivative
          Variable='S33'
          call CalculateMaterialDerivative(Variable,S33,
     *                               BS33,S33Old,S33OldOld)
          DS33Dt=MaterialDerivative
          deallocate(MaterialDerivative)
c           
          term1=0.d0
          term2=0.d0
          term3=0.d0
c
          do i=1,NumberOfElements      
c
            rstar=StrainRate(i)/dmax1(Vorticity(i),tiny)         
c
            a1=W12(i)*S12(i)+W13(i)*S13(i)
            a2=W12(i)*S22(i)+W13(i)*S23(i)
            a3=W12(i)*S23(i)+W13(i)*S33(i)
            a4=-W12(i)*S11(i)+W23(i)*S13(i)
            a5=-W12(i)*S12(i)+W23(i)*S23(i)
            a6=-W12(i)*S13(i)+W23(i)*S33(i)
            a7=-W13(i)*S11(i)-W23(i)*S12(i)
            a8=-W13(i)*S12(i)-W23(i)*S22(i)
            a9=-W13(i)*S13(i)-W23(i)*S23(i)
c
            term1=a1*DS11Dt(i)+a2*DS12Dt(i)+a3*DS13Dt(i)+
     *          a4*DS12Dt(i)+a5*DS22Dt(i)+a6*DS23Dt(i)+
     *          a7*DS13Dt(i)+a8*DS23Dt(i)+a9*DS33Dt(i)
c
            if(LRotation) then
c
              term2=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a3*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a4*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a7*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a8*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
              term3=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a3*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a4*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a7*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a8*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
            endif
c
            D4=(0.5*(StrainRate(i)**2+Vorticity(i)**2))**2
c
            rtelda=2.*(term1+term2+term3)/dmax1(D4,tiny)
c
            fr1Coefficient(i)=(1.+cr1)*(2.*rstar/(1+rstar))*
     *                         (1.-cr3*datan(cr2*rtelda))-cr1
c
          enddo
c
        elseif(RotationCurvatureMethod.eq.'zhangyang') then
c
          do i=1,NumberOfElements      
c
            rstar=StrainRate(i)/dmax1(Vorticity(i),tiny)         
c           
            rtelda=Vorticity(i)/dmax1(StrainRate(i),tiny)
            rtelda=rtelda*(rtelda-1.)
c
            fr1Coefficient(i)=(1.+cr1)*(2.*rstar/(1+rstar))*
     *                         (1.-cr3*datan(cr2*rtelda))-cr1
c
          enddo
c
        endif
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKKLWAlfaTheta
c
c#############################################################################################
      use User0, only: urfTViscosity
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Variables1, only: TurbulentKE,BTurbulentKE,
     *                      TurbulentKL,BTurbulentKL
      use Turbulence1, only: SigT,Ca0,
     *                       LambdaEff,coefficientFW,
     *                       TurbulentViscosityTs,
     *                       BLambdaEff,BcoefficientFW,
     *                       BTurbulentViscosityTs,
     *                       AlfaTheta,BAlfaTheta
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision ::  sum,factor1,factor2,AlfaThetaold
c********************************************************************************************
c
c---- Calculate the small scale eddy viscosity and eddy diffusivity
c        
      do i=1,NumberOfElements    
c
        AlfaThetaold=AlfaTheta(i)
c      
        sum=dmax1(TurbulentKE(i)+TurbulentKL(i),tiny)
        factor1=coefficientFW(i)*dmax1(TurbulentKE(i),0.)*
     *                      TurbulentViscosityTs(i)/(SigT*sum)
        factor2=(1.-coefficientFW(i))*Ca0*
     *             dsqrt(dmax1(TurbulentKE(i),0.))*LambdaEff(i)
c
        AlfaTheta(i)=(1.-urfTViscosity)*AlfaThetaold+
     *                         urfTViscosity*(factor1+factor2)
c     
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          AlfaThetaold=BAlfaTheta(i,j)
c      
          sum=dmax1(BTurbulentKE(i,j)+BTurbulentKL(i,j),tiny)
          factor1=BcoefficientFW(i,j)*dmax1(BTurbulentKE(i,j),0.)*
     *                          BTurbulentViscosityTs(i,j)/(SigT*sum)
          factor2=(1.-BcoefficientFW(i,j))*Ca0*
     *             dsqrt(dmax1(BTurbulentKE(i,j),0.))*BLambdaEff(i,j)
c
          BAlfaTheta(i,j)=(1.-urfTViscosity)*AlfaThetaold+
     *                              urfTViscosity*(factor1+factor2)
c     
         enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKKLWCoefficients
c
c#############################################################################################
      use User0, only: urfTViscosity
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Variables1, only: TurbulentKE,BTurbulentKE,
     *                      TurbulentOmega,BTurbulentOmega,
     *                      TurbulentKL,BTurbulentKL
      use PhysicalProperties1, only: Density,BDensity,
     *                               Viscosity,BViscosity
      use Turbulence1, only: SigT,cmu,
     *                       A0,As,Av,Abp,Anat,Ats,Cbpcrit,Cnc,
     *                       Cnatcrit,Cint,Ctscrit,Crnat,Cr,Ca0,Css,
     *                       Ctl,Cw1,Cw2,Cw3,Cwr,Clambda,C11,C12,Ctl,
     *                       LambdaT,LambdaEff,coefficientFW,
     *                       TurbulentKTs,TurbulentViscosityTs,
     *                       BLambdaT,BLambdaEff,BcoefficientFW,
     *                       BTurbulentKTs,BTurbulentViscosityTs,
     *                       AlfaT,BAlfaT,Vorticity,StrainRate,
     *                       BVorticity,BStrainRate,
     *                       TurbulentViscosityTl,BTurbulentViscosityTl
      use Constants1, only: tiny,twothird
      use WallDistance1, only: WallDistance,BWallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision ::  ReTlocal,coefficientFv,coefficientFSS,
     *                     coefficientCmu,coefficientFint,
     *                     BReTlocal,BcoefficientFv,BcoefficientFSS,
     *                     BcoefficientCmu,BcoefficientFint,ReOmega,
     *                     betaTS,Ktl,denom,fTawl,dEff,AlfaTold
c********************************************************************************************
c
c---- Calculate the small scale eddy viscosity and eddy diffusivity
c        
      do i=1,NumberOfElements    
c
        LambdaT(i)=dsqrt(dmax1(TurbulentKE(i),0.))/
     *                          dmax1(TurbulentOmega(i),tiny)
        LambdaEff(i)=dmin1(Clambda*WallDistance(i),LambdaT(i))
        coefficientFW(i)=
     *             (LambdaEff(i)/dmax1(LambdaT(i),tiny))**twothird
        ReTlocal=coefficientFW(i)*coefficientFW(i)*
     *             dmax1(TurbulentKE(i),0.)* Density(i)/
     *                 (Viscosity(i)*dmax1(TurbulentOmega(i),tiny))
        coefficientFv=1.-dexp(-dsqrt(ReTlocal)/Av)
        coefficientFSS=dexp(-(Css*Viscosity(i)*
     *        Vorticity(i)/(Density(i)*dmax1(TurbulentKE(i),tiny)))**2)
        coefficientCmu=1./(A0+As*StrainRate(i)/
     *                                 dmax1(TurbulentOmega(i),tiny)) 
        coefficientFint=dmin1(dmax1(TurbulentKE(i),0.)/
     *           (Cint*dmax1(TurbulentKE(i)+TurbulentKL(i),tiny)),1.)
        TurbulentKTs(i)=coefficientFSS*
     *                      coefficientFW(i)*dmax1(TurbulentKE(i),0.)
        TurbulentViscosityTs(i)=
     *          coefficientFv*coefficientFint*
     *            coefficientCmu*dsqrt(TurbulentKTs(i))*LambdaEff(i)

        AlfaTold=AlfaT(i)
        AlfaT(i)=(1.-urfTViscosity)*AlfaTold+
     *            urfTViscosity*coefficientFv*cmu*
     *                   LambdaEff(i)*dsqrt(TurbulentKTs(i))
c     
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BLambdaT(i,j)=dsqrt(dmax1(BTurbulentKE(i,j),0.))/
     *                              dmax1(BTurbulentOmega(i,j),tiny)
          BLambdaEff(i,j)=
     *                dmin1(Clambda*BWallDistance(i,j),BLambdaT(i,j))
          BcoefficientFW(i,j)=(BLambdaEff(i,j)/
     *                           dmax1(BLambdaT(i,j),tiny))**twothird
          BReTlocal=BcoefficientFW(i,j)*BcoefficientFW(i,j)*
     *                dmax1(BTurbulentKE(i,j),0.)*BDensity(i,j)/
     *                (BViscosity(i,j)*dmax1(BTurbulentOmega(i,j),tiny))
          BcoefficientFv=1.-dexp(-dsqrt(BReTlocal)/Av)
          BcoefficientFSS=
     *        dexp(-(Css*BViscosity(i,j)*BVorticity(i,j)/
     *               (BDensity(i,j)*dmax1(BTurbulentKE(i,j),tiny)))**2)
          BcoefficientCmu=1./
     *         (A0+As*BStrainRate(i,j)/dmax1(BTurbulentOmega(i,j),tiny))
          BcoefficientFint=dmin1(dmax1(BTurbulentKE(i,j),0.)/
     *        (Cint*dmax1(BTurbulentKE(i,j)+BTurbulentKL(i,j),tiny)),1.)
          BTurbulentKTs(i,j)=BcoefficientFSS*
     *                 BcoefficientFW(i,j)*dmax1(BTurbulentKE(i,j),0.)
          BTurbulentViscosityTs(i,j)=BcoefficientFv*
     *               BcoefficientFint*BcoefficientCmu*
     *                        dsqrt(BTurbulentKTs(i,j))*BLambdaEff(i,j)
          AlfaTold=BAlfaT(i,j)
          BAlfaT(i,j)=(1.-urfTViscosity)*AlfaTold+
     *               urfTViscosity*BcoefficientFv*cmu*
     *                       BLambdaEff(i,j)*dsqrt(BTurbulentKTs(i,j))
C     
         enddo
      enddo
c
c---- Calculate the large scale eddy viscosity
c        
      do i=1,NumberOfElements    
c
        ReOmega=Density(i)*WallDistance(i)*
     *                    WallDistance(i)*Vorticity(i)/Viscosity(i)
        betaTS=1.-dexp(-(dmax1(ReOmega-Ctscrit,0.))**2/Ats)
        Ktl=dmax1(TurbulentKE(i)-TurbulentKTs(i),0.)
        denom=dmax1((LambdaEff(i)*Vorticity(i))**2,tiny)
        fTawl=1.-dexp(-Ctl*Ktl/denom)
        dEff=LambdaEff(i)/Clambda
c          
        TurbulentViscosityTl(i)=
     *    dmin1((Vorticity(i)*Density(i)/Viscosity(i))*
     *         (fTawl*c11*(LambdaEff(i)**3)*dsqrt(Ktl)+
     *                betaTS*C12*(dEff**4)*Vorticity(i)),
     *                        0.5*dmax1(TurbulentKL(i)+Ktl,0.)/
     *                                dmax1(StrainRate(i),tiny))
c
      enddo          
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          ReOmega=BDensity(i,j)*BWallDistance(i,j)*
     *            BWallDistance(i,j)*BVorticity(i,j)/BViscosity(i,j)
          betaTS=1.-dexp(-(dmax1(ReOmega-Ctscrit,0.))**2/Ats)
          Ktl=dmax1(BTurbulentKE(i,j)-BTurbulentKTs(i,j),0.)
          denom=dmax1((BLambdaEff(i,j)*BVorticity(i,j))**2,tiny)
          fTawl=1.-dexp(-Ctl*Ktl/denom)
          dEff=BLambdaEff(i,j)/Clambda
c          
          BTurbulentViscosityTl(i,j)=
     *       dmin1((BVorticity(i,j)*BDensity(i,j)/BViscosity(i,j))*
     *              (fTawl*C11*(BLambdaEff(i,j)**3)*dsqrt(Ktl)+
     *                  betaTS*C12*(dEff**4)*BVorticity(i,j)),
     *                     0.5*dmax1(BTurbulentKL(i,j)+Ktl,0.)/
     *                                dmax1(BStrainRate(i,j),tiny))
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateWrayAgarwalCoefficients
c
c#############################################################################################
      use User0, only: LWrayAgarwalRotationCurvatureCorrection,
     *                 RotationCurvatureMethod,LRotation,
     *                 LRough,GrainSize,LWallDistanceFreeWA
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Turbulence1, only: Vorticity,BVorticity,
     *                       StrainRate,BStrainRate,Cr1WA,
     *                       f1WA,fmuWA,Bf1WA,BfmuWA,cw,cmu,cmu50,
     *                       cr1,cr2,cr3,fr1Coefficient,
     *                       S11,S12,S13,S22,S23,S33,
     *                       BS11,BS12,BS13,BS22,BS23,BS33,
     *                       S11Old,S12Old,S13Old,S22Old,S23Old,S33Old,
     *                       S11OldOld,S12OldOld,S13OldOld,S22OldOld,
     *                       S23OldOld,S33OldOld,DS11Dt,DS12Dt,DS13Dt,
     *                       DS22Dt,DS23Dt,DS33Dt,
     *                       W11,W12,W13,W22,W23,W33
      use Variables1, only: ModifiedED,BModifiedED,MaterialDerivative,
     *                      uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               TurbulentViscosity,BDensity,
     *                               BViscosity,BTurbulentViscosity
      use Constants1, only: tiny
      use Coriolis1, only: AngularVelocityX,AngularVelocityY,
     *                     AngularVelocityZ
      use WallDistance1, only: WallDistance,BWallDistance,iTau,BiTau
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j
      double precision :: term1,term2,term3,term4,sum,Smagnitude,
     *                    Wmagnitude,kWA,wWA,etaWA,arg1WA,rstar,rtelda,
     *                    D4,a1,a2,a3,a4,a5,a6,a7,a8,a9,
     *                    WallVelocity,dNorm,sqrtRS
c
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE CalculateMaterialDerivative(Variable,FiTemp,
     *                                 BFiTemp,FiTempOld,FiTempOldOld)
c********************************************************************************************
          character*10 Variable
          double precision FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
          double precision, dimension(:) :: FiTemp
          double precision, dimension(:) :: FiTempOld
          double precision, dimension(:) :: FiTempOldOld
          double precision, dimension(:,:) :: BFiTemp
c********************************************************************************************
        end SUBROUTINE CalculateMaterialDerivative
c********************************************************************************************
        FUNCTION TangentialVelocity(i1)
c********************************************************************************************
          integer :: i1
          double precision :: TangentialVelocity
c********************************************************************************************
        end FUNCTION TangentialVelocity
c********************************************************************************************
        FUNCTION BTangentialVelocity(i1,i2)
c********************************************************************************************
          integer :: i1,i2
          double precision :: BTangentialVelocity
c********************************************************************************************
        end FUNCTION BTangentialVelocity
c********************************************************************************************
      end interface
c********************************************************************************************
c
c----- Wray Agarwal wall distance free formulatio
c
      if(LWallDistanceFreeWA) then
c
        do i=1,NumberOfElements
c
          term1=Density(i)*ModifiedED(i)/Viscosity(i)
          fmuWA(i)=term1**3/(term1**3+cw**3)
c        
          Wmagnitude=Vorticity(i)
          Smagnitude=dmax1(StrainRate(i),1.d-16)
c
          kWA=TurbulentViscosity(i)*Smagnitude/(density(i)*cmu50)
c
          if(LRough) then
c        
c            WallVelocity=TangentialVelocity(i)
            WallVelocity=dsqrt(uVelocity(i)**2+
     *                      vVelocity(i)**2+WVelocity(i)**2)
            Cr1WA=1.+Density(i)*WallVelocity*
     *                        GrainSize(iTau(i))/Viscosity(i)
            Cr1WA=1./Cr1WA
            kWA=Cr1WA*kWA
c        
          endif
c
          wWA=Smagnitude/cmu50
          etaWA=Smagnitude*dmax1(1.,Wmagnitude/Smagnitude)
          arg1WA=(Viscosity(i)/Density(i)+ModifiedED(i))/2.
          arg1WA=arg1WA*(etaWA*etaWA)/dmax1(cmu*kWA*wWA,tiny)
c
          f1WA(i)=dtanh(arg1WA*arg1WA*arg1WA*arg1WA)   
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            term1=BDensity(i,j)*BModifiedED(i,j)/BViscosity(i,j)
            BfmuWA(i,j)=term1**3/(term1**3+cw**3)
c        
            Wmagnitude=BVorticity(i,j)
            Smagnitude=dmax1(BStrainRate(i,j),1.d-16)
c
            kWA=BTurbulentViscosity(i,j)*Smagnitude/
     *                              (Bdensity(i,j)*cmu50)
c
            if(LRough) then
c        
c              WallVelocity=BTangentialVelocity(i,j) !this requires distance to the wall
              WallVelocity=dsqrt(BuVelocity(i,j)**2+
     *                      BvVelocity(i,j)**2+BWVelocity(i,j)**2)
              Cr1WA=1.+BDensity(i,j)*WallVelocity*
     *                        GrainSize(BiTau(i,j))/BViscosity(i,j)
              Cr1WA=1./Cr1WA
              kWA=Cr1WA*kWA
c        
            endif
c
            wWA=Smagnitude/cmu50
            etaWA=Smagnitude*dmax1(1.,Wmagnitude/Smagnitude)
            arg1WA=(BViscosity(i,j)/BDensity(i,j)+BModifiedED(i,j))/2.
            arg1WA=arg1WA*(etaWA*etaWA)/dmax1(cmu*kWA*wWA,tiny)
c
            Bf1WA(i,j)=dtanh(arg1WA*arg1WA*arg1WA*arg1WA)   
c
          enddo
        enddo
c
      else
c
c----- Wray Agarwal 2017 formulation (with wall distance)
c
        do i=1,NumberOfElements
c
          dNorm=WallDistance(i)
          term1=Density(i)*ModifiedED(i)/Viscosity(i)
c
          if(LRough) then
c
            dNorm=dNorm+0.03*GrainSize(iTau(i))
            term1=term1+Cr1WA*GrainSize(iTau(i))/dNorm
c
          endif
c
          fmuWA(i)=term1**3/(term1**3+cw**3)
c        
          Smagnitude=dmax1(StrainRate(i),1.d-16)
c
          sqrtRS=dsqrt(ModifiedED(i)*Smagnitude)
          etaWA=dNorm*sqrtRS*Density(i)/(20.*Viscosity(i))
          arg1WA=(1.+20.*etaWA)/(1.+(dNorm*Density(i)*
     *         dmax1(sqrtRS,1.5*ModifiedED(i))/(20.*Viscosity(i)))**2)
          f1WA(i)=dmin1(dtanh(arg1WA*arg1WA*arg1WA*arg1WA),0.9)   
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            dNorm=BWallDistance(i,j)
            term1=BDensity(i,j)*BModifiedED(i,j)/BViscosity(i,j)
c
            if(LRough) then
c
              dNorm=dNorm+0.03*GrainSize(BiTau(i,j))
              term1=term1+Cr1WA*GrainSize(BiTau(i,j))/dNorm
c
            endif
c
            BfmuWA(i,j)=term1**3/(term1**3+cw**3)
c        
            Smagnitude=dmax1(BStrainRate(i,j),1.d-16)
c
            sqrtRS=dsqrt(BModifiedED(i,j)*Smagnitude)
            etaWA=dNorm*sqrtRS*BDensity(i,j)/(20.*BViscosity(i,j))
            arg1WA=(1.+20.*etaWA)/(1.+(dNorm*BDensity(i,j)*
     *                     dmax1(sqrtRS,1.5*BModifiedED(i,j))/
     *                                  (20.*BViscosity(i,j)))**2)
            Bf1WA(i,j)=dmin1(dtanh(arg1WA*arg1WA*arg1WA*arg1WA),0.9)   
c
          enddo
        enddo
c
      endif
c
c--- Calculate the rotation/curvature correction term
c
      if(LWrayAgarwalRotationCurvatureCorrection) then
c
        if(RotationCurvatureMethod.eq.'spalartshur') then
c
          allocate(MaterialDerivative(NumberOfElements))
          Variable='S11'
          call CalculateMaterialDerivative(Variable,S11,
     *                               BS11,S11Old,S11OldOld)
          DS11Dt=MaterialDerivative
          Variable='S12'
          call CalculateMaterialDerivative(Variable,S12,
     *                               BS12,S12Old,S12OldOld)
          DS12Dt=MaterialDerivative
          Variable='S13'
          call CalculateMaterialDerivative(Variable,S13,
     *                               BS13,S13Old,S13OldOld)
          DS13Dt=MaterialDerivative
          Variable='S22'
          call CalculateMaterialDerivative(Variable,S22,
     *                               BS22,S22Old,S22OldOld)
          DS22Dt=MaterialDerivative
          Variable='S23'
          call CalculateMaterialDerivative(Variable,S23,
     *                               BS23,S23Old,S23OldOld)
          DS23Dt=MaterialDerivative
          Variable='S33'
          call CalculateMaterialDerivative(Variable,S33,
     *                               BS33,S33Old,S33OldOld)
          DS33Dt=MaterialDerivative
          deallocate(MaterialDerivative)
c           
          term1=0.d0
          term2=0.d0
          term3=0.d0
c
          do i=1,NumberOfElements      
c
            rstar=StrainRate(i)/dmax1(Vorticity(i),tiny)         
c
            a1=W12(i)*S12(i)+W13(i)*S13(i)
            a2=W12(i)*S22(i)+W13(i)*S23(i)
            a3=W12(i)*S23(i)+W13(i)*S33(i)
            a4=-W12(i)*S11(i)+W23(i)*S13(i)
            a5=-W12(i)*S12(i)+W23(i)*S23(i)
            a6=-W12(i)*S13(i)+W23(i)*S33(i)
            a7=-W13(i)*S11(i)-W23(i)*S12(i)
            a8=-W13(i)*S12(i)-W23(i)*S22(i)
            a9=-W13(i)*S13(i)-W23(i)*S23(i)
c
            term1=a1*DS11Dt(i)+a2*DS12Dt(i)+a3*DS13Dt(i)+
     *          a4*DS12Dt(i)+a5*DS22Dt(i)+a6*DS23Dt(i)+
     *          a7*DS13Dt(i)+a8*DS23Dt(i)+a9*DS33Dt(i)
c
            if(LRotation) then
c
              term2=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a3*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a4*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a7*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a8*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
              term3=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a3*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a4*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a7*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a8*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
            endif
c
            D4=(0.5*(StrainRate(i)**2+Vorticity(i)**2))**2
c
            rtelda=2.*(term1+term2+term3)/dmax1(D4,tiny)
c
            fr1Coefficient(i)=(1.+cr1)*(2.*rstar/(1.+rstar))*
     *                         (1.-cr3*datan(cr2*rtelda))-cr1
c
          enddo
c
        elseif(RotationCurvatureMethod.eq.'zhangyang') then
c
          do i=1,NumberOfElements      
c
            rstar=StrainRate(i)/dmax1(Vorticity(i),tiny)         
c           
            rtelda=Vorticity(i)/dmax1(StrainRate(i),tiny)
            rtelda=rtelda*(rtelda-1.)
c
            fr1Coefficient(i)=(1.+cr1)*(2.*rstar/(1+rstar))*
     *                         (1.-cr3*datan(cr2*rtelda))-cr1
c
          enddo
c
        endif
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateSpalartAllmarasCoefficients
c
c#############################################################################################
      use User0, only: LNegativeSpalartAllmaras,LRotation,
     *                 LSpalartAllmarasRotationCorrection,
     *                 LSpalartAllmarasRotationCurvatureCorrection,
     *                 RotationCurvatureMethod,LTransitionalSA,
     *                 WallTreatment,LRough,GrainSize
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use WallDistance1, only: WallDistance,BWallDistance,iTau,BiTau
      use Turbulence1, only: fv1Coefficient,Stelda,fwCoefficient,
     *                       ft2Coefficient,Bfv1Coefficient,
     *                       fnCoefficient,BfnCoefficient,
     *                       Vorticity,BVorticity,fr1Coefficient,
     *                       StrainRate,BStrainRate,cw2,cw3,cv1,
     *                       cappa,c2,c3,ct3,ct4,cn1,crot,
     *                       cr1,cr2,cr3,Cr1SA,
     *                       S11,S12,S13,S22,S23,S33,
     *                       BS11,BS12,BS13,BS22,BS23,BS33,
     *                       S11Old,S12Old,S13Old,S22Old,S23Old,S33Old,
     *                       S11OldOld,S12OldOld,S13OldOld,S22OldOld,
     *                       S23OldOld,S33OldOld,DS11Dt,DS12Dt,DS13Dt,
     *                       DS22Dt,DS23Dt,DS33Dt,
     *                       W11,W12,W13,W22,W23,W33
      use Variables1, only: ModifiedED,BModifiedED,MaterialDerivative,
     *                      TGamma,BTGamma,uVelocity,vVelocity,wVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity,
     *                               TurbulentViscosity
      use Constants1, only: tiny
      use Coriolis1, only: AngularVelocityX,AngularVelocityY,
     *                     AngularVelocityZ
      use ReferenceValues1
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,k
      double precision :: term1,term2,term3,Smagnitude,fv2,Wmagnitude,
     *                    rfactor,gfactor,Sprime,Sbar,rstar,rtelda,D4,
     *                    a1,a2,a3,a4,a5,a6,a7,a8,a9,Rev,ReO,ReOc,Xi1,
     *                    Xi2,Re,Vel,nuBC,dNorm
c
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE CalculateMaterialDerivative(Variable,FiTemp,
     *                                 BFiTemp,FiTempOld,FiTempOldOld)
c--------------------------------------------------------------------------------
          character*10 Variable
          double precision FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
          double precision, dimension(:) :: FiTemp
          double precision, dimension(:) :: FiTempOld
          double precision, dimension(:) :: FiTempOldOld
          double precision, dimension(:,:) :: BFiTemp
c--------------------------------------------------------------------------------
        end SUBROUTINE CalculateMaterialDerivative
c--------------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------------
c
      do i=1,NumberOfElements
c
        if(LRough.and..not.LNegativeSpalartAllmaras.and.
     *                     WallTreatment.eq.'lowreynoldsnumber') then
c
          dNorm=WallDistance(i)+0.03*GrainSize(iTau(i))
          term1=Density(i)*ModifiedED(i)/Viscosity(i)
          term2=Cr1SA*GrainSize(iTau(i))/dNorm
          term1=term1+term2
          fv1Coefficient(i)=term1**3/(term1**3+cv1**3)
          ft2Coefficient(i)=ct3*dexp(-ct4*term1*term1)
          fv2=1.-ModifiedED(i)/(Viscosity(i)/Density(i)+
     *                     fv1Coefficient(i)*ModifiedED(i))
c
        else
c
          dNorm=WallDistance(i)
          term1=Density(i)*ModifiedED(i)/Viscosity(i)
          fv1Coefficient(i)=term1**3/(term1**3+cv1**3)
          ft2Coefficient(i)=ct3*dexp(-ct4*term1*term1)
          fv2=1.-term1/(1.+term1*fv1Coefficient(i))
c
        endif
c
        Sbar=ModifiedED(i)*fv2/((cappa*dNorm)**2)
c
        Wmagnitude=Vorticity(i)
        if(.not.LNegativeSpalartAllmaras.and.
     *                 LSpalartAllmarasRotationCorrection)
     *    Wmagnitude=Wmagnitude+crot*dmin1(0.,StrainRate(i)-Wmagnitude)
c
        if(Sbar.ge.-c2*Wmagnitude) then
c
          Stelda(i)=Wmagnitude+Sbar
c
        elseif(Sbar.lt.-c2*Wmagnitude) then
c
          term2=c2**2*Wmagnitude+c3*Sbar
          term3=(c3-2.*c2)*Wmagnitude-Sbar
          Stelda(i)=Wmagnitude*(1.+term2/term3)
c
        endif
c
c-------------------------------------------------------------------------
c
        if(Stelda(i).eq.0.) then
c        
          rfactor=10.
c
        else        
c
          rfactor=dmin1(ModifiedED(i)/
     *        (Stelda(i)*(cappa*dNorm)**2),10.)
c
        endif
c
        gfactor=rfactor+cw2*(rfactor**6-rfactor)
        term1=cw3**6
        term1=(1.+term1)/(gfactor**6+term1)
        fwCoefficient(i)=gfactor*(term1**(1./6))
c
      enddo
c
      if(LRough.and..not.LNegativeSpalartAllmaras.and.
     *                   WallTreatment.eq.'lowreynoldsnumber') then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            dNorm=BWallDistance(i,j)+0.03*GrainSize(BiTau(i,j))
            term1=BDensity(i,j)*BModifiedED(i,j)/BViscosity(i,j)
            term2=Cr1SA*GrainSize(BiTau(i,j))/dNorm
            term1=term1+term2
            Bfv1Coefficient(i,j)=term1**3/(term1**3+cv1**3)
c
          enddo
        enddo
c
      else
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            term1=BDensity(i,j)*BModifiedED(i,j)/BViscosity(i,j)
            Bfv1Coefficient(i,j)=term1**3/(term1**3+cv1**3)
c
          enddo
        enddo
c
      endif
c
      if(LNegativeSpalartAllmaras) then
c
        do i=1,NumberOfElements
c
          term1=Density(i)*ModifiedED(i)/Viscosity(i) 
          term2=term1**3
          term3=cn1-term2
          if(term3.eq.0.) term3=tiny          
          fnCoefficient(i)=(cn1+term2)/term3
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            term1=BDensity(i,j)*BModifiedED(i,j)/BViscosity(i,j)
            term2=term1**3
            term3=cn1-term2
            if(term3.eq.0.) term3=tiny          
            BfnCoefficient(i,j)=(cn1+term2)/term3
c
          enddo
        enddo
c
      endif
c
c--- Calculate the rotation/curvature correction term
c
      if(LSpalartAllmarasRotationCurvatureCorrection
     *              .and..not.LNegativeSpalartAllmaras.and..not.
     *                    LSpalartAllmarasRotationCorrection) then
c
        if(RotationCurvatureMethod.eq.'spalartshur') then
c
          allocate(MaterialDerivative(NumberOfElements))
          Variable='S11'
          call CalculateMaterialDerivative(Variable,S11,
     *                               BS11,S11Old,S11OldOld)
          DS11Dt=MaterialDerivative
          Variable='S12'
          call CalculateMaterialDerivative(Variable,S12,
     *                               BS12,S12Old,S12OldOld)
          DS12Dt=MaterialDerivative
          Variable='S13'
          call CalculateMaterialDerivative(Variable,S13,
     *                               BS13,S13Old,S13OldOld)
          DS13Dt=MaterialDerivative
          Variable='S22'
          call CalculateMaterialDerivative(Variable,S22,
     *                               BS22,S22Old,S22OldOld)
          DS22Dt=MaterialDerivative
          Variable='S23'
          call CalculateMaterialDerivative(Variable,S23,
     *                               BS23,S23Old,S23OldOld)
          DS23Dt=MaterialDerivative
          Variable='S33'
          call CalculateMaterialDerivative(Variable,S33,
     *                               BS33,S33Old,S33OldOld)
          DS33Dt=MaterialDerivative
          deallocate(MaterialDerivative)
c           
          term1=0.d0
          term2=0.d0
          term3=0.d0
c
          do i=1,NumberOfElements      
c
            rstar=StrainRate(i)/dmax1(Vorticity(i),tiny)         
c
            a1=W12(i)*S12(i)+W13(i)*S13(i)
            a2=W12(i)*S22(i)+W13(i)*S23(i)
            a3=W12(i)*S23(i)+W13(i)*S33(i)
            a4=-W12(i)*S11(i)+W23(i)*S13(i)
            a5=-W12(i)*S12(i)+W23(i)*S23(i)
            a6=-W12(i)*S13(i)+W23(i)*S33(i)
            a7=-W13(i)*S11(i)-W23(i)*S12(i)
            a8=-W13(i)*S12(i)-W23(i)*S22(i)
            a9=-W13(i)*S13(i)-W23(i)*S23(i)
c
            term1=a1*DS11Dt(i)+a2*DS12Dt(i)+a3*DS13Dt(i)+
     *          a4*DS12Dt(i)+a5*DS22Dt(i)+a6*DS23Dt(i)+
     *          a7*DS13Dt(i)+a8*DS23Dt(i)+a9*DS33Dt(i)
c
            if(LRotation) then
c
              term2=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a3*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a4*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a7*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a8*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
              term3=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a3*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a4*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a7*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a8*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
            endif
c
            D4=(0.5*(StrainRate(i)**2+Vorticity(i)**2))**2
c
            rtelda=2.*(term1+term2+term3)/dmax1(D4,tiny)
c
            fr1Coefficient(i)=(1.+cr1)*(2.*rstar/(1+rstar))*
     *                         (1.-cr3*datan(cr2*rtelda))-cr1
c
          enddo
c
        elseif(RotationCurvatureMethod.eq.'zhangyang') then
c
          do i=1,NumberOfElements      
c
            rstar=StrainRate(i)/dmax1(Vorticity(i),tiny)         
c           
            rtelda=Vorticity(i)/dmax1(StrainRate(i),tiny)
            rtelda=rtelda*(rtelda-1.)
c
            fr1Coefficient(i)=(1.+cr1)*(2.*rstar/(1+rstar))*
     *                         (1.-cr3*datan(cr2*rtelda))-cr1
c
          enddo
c
        endif
c
      endif
c
c---- Calculate intermittency coefficient
c
      if(LTransitionalSA) then
c
        do i=1,NumberOfElements
c
          Rev=Density(i)*WallDistance(i)*
     *          WallDistance(i)*Vorticity(i)/Viscosity(i)
          ReO=Rev/2.193
          ReOc=803.73*(TuInfinity+0.6067)**(-1.027)
c
          Xi1=0.002
          Re=RhoInfinity*UInfinity*Lreference/MuInfinity
          Xi2=5./Re
          term1=dmax1(ReO-ReOc,0.)/(Xi1*ReOc)
          Vel=dsqrt(uVelocity(i)**2+vVelocity(i)**2+wVelocity(i)**2)
          nuBC=TurbulentViscosity(i)/(Density(i)*Vel*WallDistance(i))
          term2=dmax1(nuBC-Xi2,0.)/Xi2
c
          TGamma(i)=1.-dexp(-dsqrt(term1)-dsqrt(term2))
c
        enddo  
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
            BTGamma(i,j)=TGamma(k)
c
          enddo
        enddo
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateRNGCoefficients
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Turbulence1, only: StrainRate,C2eRNG,eta0RNG,betaRNG,ce2,cmu
      use Variables1, only: TurbulentKE,TurbulentED
      use Constants1, only: tiny

c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: eta,eta3
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        eta=StrainRate(i)*dmax1(TurbulentKE(i),0.)/
     *                      dmax1(TurbulentED(i),tiny)
c
        eta3=eta*eta*eta
        C2eRNG(i)=ce2+cmu*eta3*(1.-eta/eta0RNG)/(1.+betaRNG*eta3)
        C2eRNG(i)=dmax1(C2eRNG(i),0.)
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateRealizableCoefficients
c
c#############################################################################################
c
      use User0, only: LRotation     
      use Turbulence1, only: 
     *                       S11,S12,S13,S22,S23,S33,      
     *                       BS11,BS12,BS13,BS22,BS23,BS33,      
     *                       W11,W12,W13,W22,W23,W33,      
     *                       BW11,BW12,BW13,BW22,BW23,BW33,
     *                       StrainRate,Vorticity,BStrainRate,
     *                       BVorticity,c1R,Bc1R,cmuR,BcmuR,cmu,A0
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Coriolis1, only: AngularVelocityX,AngularVelocityY,
     *                     AngularVelocityZ
      use Variables1, only: TurbulentKE,TurbulentED,
     *                      BTurbulentKE,BTurbulentED
      use Constants1, only: tiny,sqrt6,sqrt2
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: SijSjkSki,Stelda,W,W12Telda,W13Telda,W23Telda,
     *                    WTelda2,WTelda,Ustar,phi,As,Tke,Ted,eta
c********************************************************************************************
c      
      do i=1,NumberOfElements
c
!     Expanded form of SijSjkSki      
!        SijSjkSki=
!     *   S11(i)*S11(i)*S11(i)+S11(i)*S12(i)*S12(i)+S11(i)*S13(i)*S13(i)+
!     *   S12(i)*S12(i)*S11(i)+S12(i)*S22(i)*S12(i)+S12(i)*S23(i)*S13(i)+
!     *   S13(i)*S13(i)*S11(i)+S13(i)*S23(i)*S12(i)+S13(i)*S33(i)*S13(i)+
!     *   S12(i)*S11(i)*S12(i)+S12(i)*S12(i)*S22(i)+S12(i)*S13(i)*S23(i)+
!     *   S22(i)*S12(i)*S12(i)+S22(i)*S22(i)*S22(i)+S22(i)*S23(i)*S23(i)+
!     *   S23(i)*S13(i)*S12(i)+S23(i)*S23(i)*S22(i)+S23(i)*S33(i)*S23(i)+
!     *   S13(i)*S11(i)*S13(i)+S13(i)*S12(i)*S23(i)+S13(i)*S13(i)*S33(i)+
!     *   S23(i)*S12(i)*S13(i)+S23(i)*S22(i)*S23(i)+S23(i)*S23(i)*S33(i)+
!     *   S33(i)*S13(i)*S13(i)+S33(i)*S23(i)*S23(i)+S33(i)*S33(i)*S33(i)
        SijSjkSki=
     *   S11(i)*S11(i)*S11(i)+S22(i)*S22(i)*S22(i)+S33(i)*S33(i)*S33(i)+
     *   3.*(S11(i)*S12(i)*S12(i)+S11(i)*S13(i)*S13(i)+
     *   S12(i)*S22(i)*S12(i)+2.*S12(i)*S23(i)*S13(i)+
     *   S13(i)*S33(i)*S13(i)+S22(i)*S23(i)*S23(i)+S23(i)*S33(i)*S23(i))
c
        Stelda=StrainRate(i)/sqrt2
        W=SijSjkSki/dmax1(Stelda*Stelda*Stelda,tiny)
c
        if(LRotation) then
c
          W12Telda=W12(i)-2.*AngularVelocityZ
          W13Telda=W13(i)+2.*AngularVelocityY
          W23Telda=W23(i)-2.*AngularVelocityX
c
          WTelda2=  W12Telda*W12Telda+W13Telda*W13Telda+
     *          (-W12Telda)*(-W12Telda)+W23Telda*W23Telda+
     *          (-W13Telda)*(-W13Telda)+(-W23Telda)*(-W23Telda)
c
          Ustar=dsqrt(Stelda*Stelda+WTelda2)       
c        
        else
c
         WTelda=Vorticity(i)/sqrt2
         Ustar=dsqrt(Stelda*Stelda+WTelda*WTelda) 
c
        endif
c
        phi=dacos(dmax1(-1.d0,dmin1(sqrt6*W,1.d0)))/3.
        As=sqrt6*dcos(phi)
        Tke=dmax1(TurbulentKE(i),0.)
        Ted=dmax1(TurbulentED(i),tiny)
c
        cmuR(i)=1./(A0+As*Tke*Ustar/Ted)
        cmuR(i)=dmin1(cmu,cmuR(i))
c
        eta=StrainRate(i)*Tke/Ted
        c1R(i)=dmax1(0.43,eta/(eta+5.))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          SijSjkSki=
     *     BS11(i,j)*BS11(i,j)*BS11(i,j)+BS22(i,j)*BS22(i,j)*BS22(i,j)+
     *     BS33(i,j)*BS33(i,j)*BS33(i,j)+3.*(
     *     BS11(i,j)*BS12(i,j)*BS12(i,j)+BS11(i,j)*BS13(i,j)*BS13(i,j)+
     *     BS12(i,j)*BS22(i,j)*BS12(i,j)+
     *     2.*BS12(i,j)*BS23(i,j)*BS13(i,j)+
     *     BS13(i,j)*BS33(i,j)*BS13(i,j)+BS22(i,j)*BS23(i,j)*BS23(i,j)+
     *     BS23(i,j)*BS33(i,j)*BS23(i,j))
c
          Stelda=BStrainRate(i,j)/sqrt2
          W=SijSjkSki/dmax1(Stelda*Stelda*Stelda,tiny)
c
          if(LRotation) then
c
            W12Telda=BW12(i,j)-2.*AngularVelocityZ
            W13Telda=BW13(i,j)+2.*AngularVelocityY
            W23Telda=BW23(i,j)-2.*AngularVelocityX
c
            WTelda2=  W12Telda*W12Telda+W13Telda*W13Telda+
     *          (-W12Telda)*(-W12Telda)+W23Telda*W23Telda+
     *          (-W13Telda)*(-W13Telda)+(-W23Telda)*(-W23Telda)
c
            Ustar=dsqrt(Stelda*Stelda+WTelda2)       
c        
          else
c
           WTelda=BVorticity(i,j)/sqrt2
           Ustar=dsqrt(Stelda*Stelda+WTelda*WTelda) 
c
          endif
c
          phi=dacos(dmax1(-1.d0,dmin1(sqrt6*W,1.d0)))/3.
          As=sqrt6*dcos(phi)
          Tke=dmax1(BTurbulentKE(i,j),0.)
          Ted=dmax1(BTurbulentED(i,j),tiny)
c
          BcmuR(i,j)=1./(A0+As*Tke*Ustar/Ted)
          BcmuR(i,j)=dmin1(cmu,BcmuR(i,j))
c
          eta=BStrainRate(i,j)*Tke/Ted
          Bc1R(i,j)=dmax1(0.43,eta/(eta+5.))
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculatekelambremhorstCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Turbulence1, only: fmuCoefficient,f1Coefficient,f2Coefficient,
     *                       BfmuCoefficient,ReT,BReT
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use WallDistance1, only: WallDistance,BWallDistance
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: tke1,tde1,term1,term2,f2Old,fmuOld,ReN
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        tke1=dmax1(TurbulentKE(i),0.)
        tde1=dmax1(TurbulentED(i),tiny)
        ReT(i)=Density(i)*tke1**2/(Viscosity(i)*tde1)
        ReN=Density(i)*WallDistance(i)*dsqrt(tke1)/Viscosity(i)
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
c
        f2Coefficient(i)=1.-dexp(-ReT(i)*ReT(i))
c        
        term1=(1.-dexp(-0.0165*ReN))**2
        term2=1.+20.5/dmax1(ReT(i),tiny)
        fmuCoefficient(i)=term1*term2
c
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
c
        f1Coefficient(i)=1.+1.25d-4/dmax1(fmuCoefficient(i)**3,tiny)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          tke1=dmax1(BTurbulentKE(i,j),0.)
          tde1=dmax1(BTurbulentED(i,j),tiny)
          BReT(i,j)=BDensity(i,j)*tke1**2/(BViscosity(i,j)*tde1)
          ReN=BDensity(i,j)*BWallDistance(i,j)*
     *                     dsqrt(tke1)/BViscosity(i,j)
c        
          fmuOld=BfmuCoefficient(i,j)
c
          term1=(1.-dexp(-0.0165*ReN))**2
          term2=1.+20.5/dmax1(BReT(i,j),tiny)
          BfmuCoefficient(i,j)=term1*term2
c
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculatekelambremhorstmCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Turbulence1, only: fmuCoefficient,f1Coefficient,f2Coefficient,
     *                       BfmuCoefficient,ReT,BReT
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use WallDistance1, only: WallDistance,BWallDistance
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: tke1,tde1,term1,term2,f2Old,fmuOld,ReN
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        tke1=dmax1(TurbulentKE(i),0.)
        tde1=dmax1(TurbulentED(i),tiny)
        ReT(i)=Density(i)*tke1**2/(Viscosity(i)*tde1)
        ReN=Density(i)*WallDistance(i)*dsqrt(tke1)/Viscosity(i)
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
c
        f2Coefficient(i)=(1.-0.27*dexp(-ReT(i)*ReT(i)))*(1.-dexp(-ReN))
c        
        fmuCoefficient(i)=dexp(-3.4/((1.+ReT(i)/50.)*(1.+ReT(i)/50.)))
c
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
c
        f1Coefficient(i)=1.+2.744d-3/dmax1(fmuCoefficient(i)**3,tiny)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          tke1=dmax1(BTurbulentKE(i,j),0.)
          tde1=dmax1(BTurbulentED(i,j),tiny)
          BReT(i,j)=BDensity(i,j)*tke1**2/(BViscosity(i,j)*tde1)
c        
          fmuOld=BfmuCoefficient(i,j)
c
          BfmuCoefficient(i,j)=dexp(-3.4/((1.+BReT(i,j)/50.)*
     *                                       (1.+BReT(i,j)/50.)))
c
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKEhishidaCoefficientsOld
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,ZC,BFaceAreanx,BFaceAreany,BFaceAreanZ
      use WallDistance1, only: iTau,jTau,BiTau,BjTau,WallDistance,
     *                         BWallDistance
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,sqrtTurbulentKE,
     *                       BsqrtTurbulentKE,StrainRate,BStrainRate,
     *                       sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz,
     *                       BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,
     *                       SRateGradx,SRateGrady,SRateGradz,
     *                       BSRateGradx,BSRateGrady,BSRateGradz,
     *                       ReT,BReT
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity,
     *                               TurbulentViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,i1,i2,i3
      double precision :: term1,term2,term3,term4,sum,f2Old,fmuOld,
     *                    nx,ny,nz,dNorm,uWall1,vWall1,wWall1,
     *                    WallVelocity,TauWall,uTau,dplus,dotproduct
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
      do i=1,NumberOfElements
c
        sqrtTurbulentKE(i)=dsqrt(dmax1(TurbulentKE(i),0.))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BsqrtTurbulentKE(i,j)=dsqrt(dmax1(BTurbulentKE(i,j),0.))
c
        enddo
      enddo
c
      Variable='TKE05'
      call Gradient(Variable,2,sqrtTurbulentKE,sqrtTKEGradx,
     *              sqrtTKEGrady,sqrtTKEGradz,BsqrtTurbulentKE,
     *         BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,2,.false.,1)
c
      Variable='StrainRate'
      call Gradient(Variable,2,StrainRate,SRateGradx,SRateGrady,
     *              SRateGradz,BStrainRate,BSRateGradx,BSRateGrady,
     *              BSRateGradz,2,.false.,1)
c
      do i=1,NumberOfElements
c
        i2=iTau(i)
        i3=jTau(i)
        i1=NBFaceOwner(i2,i3)
c
        dNorm=WallDistance(i1)
c
        WallVelocity=TangentialVelocity(i1)
        TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
        uTau=dsqrt(TauWall/BDensity(i2,i3))
        dplus=WallDistance(i)*Density(i)*uTau/Viscosity(i)
c        
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=1.-0.3*dexp(-ReT(i)*ReT(i))
        fmuCoefficient(i)=(1.-dexp(-dplus/26.5))*(1.-dexp(-dplus/26.5))
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
c
        LTKE(i)=2.*Viscosity(i)*((sqrtTKEGradx(i)*nx+
     *                     sqrtTKEGrady(i)*ny+sqrtTKEGradz(i)*nz)**2)
        LTED(i)=Viscosity(i)*TurbulentViscosity(i)/Density(i)
        LTED(i)=LTED(i)*((SRateGradx(i)*nx+
     *                     SRateGrady(i)*ny+SRateGradz(i)*nz)**2)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          i2=BiTau(i,j)
          i3=BjTau(i,j)
          i1=NBFaceOwner(i2,i3)
c
          dNorm=WallDistance(i1)
c
          WallVelocity=TangentialVelocity(i1)
          TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
          uTau=dsqrt(TauWall/BDensity(i2,i3))
          dplus=BWallDistance(i,j)*BDensity(i,j)*uTau/BViscosity(i,j)
c        
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=
     *             (1.-dexp(-dplus/26.5))*(1.-dexp(-dplus/26.5))
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKEhishidaCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient,
     *                 MethodCalcGradientMomentum,nIterGradientMomentum,
     *                 LimitGradientMomentum,LimitGradientMomentumMethod
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,ZC,BFaceAreanx,BFaceAreany,BFaceAreanZ
      use WallDistance1, only: iTau,jTau,BiTau,BjTau,WallDistance,
     *                         BWallDistance
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,sqrtTurbulentKE,
     *                       BsqrtTurbulentKE, ReT,BReT,cmu25,
     *                       sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz,
     *                       BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      uVelGrad2x,uVelGrad2y,uVelGrad2z,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BuVelGrad2x,BuVelGrad2y,BuVelGrad2z,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      vVelGrad2x,vVelGrad2y,vVelGrad2z,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BvVelGrad2x,BvVelGrad2y,BvVelGrad2z,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      wVelGrad2x,wVelGrad2y,wVelGrad2z,
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
     *                      BwVelGrad2x,BwVelGrad2y,BwVelGrad2z,
     *                      uvVelGradxy,uwVelGradxz,BuvVelGradxy,
     *                      BuwVelGradxz
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity,
     *                               TurbulentViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,i1,i2,i3
      double precision :: term1,term2,term3,term4,sum,f2Old,fmuOld,
     *                    nx,ny,nz,dNorm,uWall1,vWall1,wWall1,
     *                    WallVelocity,TauWall,uTau,dplus,dotproduct
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
      do i=1,NumberOfElements
c
        sqrtTurbulentKE(i)=dsqrt(dmax1(TurbulentKE(i),0.))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BsqrtTurbulentKE(i,j)=dsqrt(dmax1(BTurbulentKE(i,j),0.))
c
        enddo
      enddo
c
      Variable='TKE05'
      call Gradient(Variable,2,sqrtTurbulentKE,sqrtTKEGradx,
     *              sqrtTKEGrady,sqrtTKEGradz,BsqrtTurbulentKE,
     *         BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,2,.false.,1)
c
      do i=1,NumberOfElements
c
        i2=iTau(i)
        i3=jTau(i)
        i1=NBFaceOwner(i2,i3)
c
        dNorm=WallDistance(i1)
c
        WallVelocity=TangentialVelocity(i1)
        TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
        uTau=dsqrt(TauWall/BDensity(i2,i3))
        uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
        dplus=WallDistance(i)*Density(i)*uTau/Viscosity(i)
c        
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=1.-0.3*dexp(-ReT(i)*ReT(i))
        fmuCoefficient(i)=(1.-dexp(-dplus/26.5))*(1.-dexp(-dplus/26.5))
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
c
        LTKE(i)=2.*Viscosity(i)*(sqrtTKEGradx(i)**2+
     *                     sqrtTKEGrady(i)**2+sqrtTKEGradz(i)**2)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          i2=BiTau(i,j)
          i3=BjTau(i,j)
          i1=NBFaceOwner(i2,i3)
c
          dNorm=WallDistance(i1)
c
          WallVelocity=TangentialVelocity(i1)
          TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
          uTau=dsqrt(TauWall/BDensity(i2,i3))
          uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
          dplus=BWallDistance(i,j)*BDensity(i,j)*uTau/BViscosity(i,j)
c        
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=
     *             (1.-dexp(-dplus/26.5))*(1.-dexp(-dplus/26.5))
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      Variable='velx'
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGradx,uVelGrad2x,uvVelGradxy,uwVelGradxz,BuVelGradx,
     *       BuVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      LTED=uVelGrad2x+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGrady,uvVelGradxy,uVelGrad2y,uwVelGradxz,BuVelGrady,
     *       BuvVelGradxy,BuVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      LTED=LTED+uVelGrad2y+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGradz,uvVelGradxy,uwVelGradxz,uVelGrad2z,BuVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BuVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      LTED=LTED+uVelGrad2z+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGradx,vVelGrad2x,uvVelGradxy,uwVelGradxz,BvVelGradx,
     *       BvVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      LTED=LTED+vVelGrad2x+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGrady,uvVelGradxy,vVelGrad2y,uwVelGradxz,BvVelGrady,
     *       BuvVelGradxy,BvVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      LTED=LTED+vVelGrad2y+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGradz,uvVelGradxy,uwVelGradxz,vVelGrad2z,BvVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BvVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
      LTED=LTED+vVelGrad2z+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGradx,wVelGrad2x,uvVelGradxy,uwVelGradxz,BwVelGradx,
     *       BwVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                             LimitGradientMomentumMethod)
      LTED=LTED+wVelGrad2x+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGrady,uvVelGradxy,wVelGrad2y,uwVelGradxz,BwVelGrady,
     *       BuvVelGradxy,BwVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                             LimitGradientMomentumMethod)
      LTED=LTED+wVelGrad2y+uvVelGradxy+uwVelGradxz
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGradz,uvVelGradxy,uwVelGradxz,wVelGrad2z,BwVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BwVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                            LimitGradientMomentumMethod)
      LTED=LTED+wVelGrad2z+uvVelGradxy+uwVelGradxz
c
      LTED=LTED*LTED
c
      LTED=LTED*Viscosity*TurbulentViscosity/Density
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKEtagawaCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,zc,BFaceAreanx,BFaceAreany,BFaceAreanz
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,ReT,BReT,cmu25
      use WallDistance1, only: iTau,jTau,BiTau,BjTau,WallDistance,
     *                         BWallDistance
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,WVelocity,
     *                      BuVelocity,BvVelocity,BWVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,i1,i2,i3
      double precision :: nx,ny,nz,dNorm,uWall1,vWall1,wWall1,
     *                    WallVelocity,TauWall,uTau,dplus,f2Old,
     *                    fmuOld,dotproduct
c********************************************************************************************
      interface
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
      do i=1,NumberOfElements
c
        i2=iTau(i)
        i3=jTau(i)
        i1=NBFaceOwner(i2,i3)
c
        dNorm=WallDistance(i1)
c
        WallVelocity=TangentialVelocity(i1)
        TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
        uTau=dsqrt(TauWall/BDensity(i2,i3))
        uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
        dplus=WallDistance(i)*Density(i)*uTau/Viscosity(i)
c        
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=(1.-0.3*dexp(-ReT(i)*ReT(i)/42.25))*
     *                         ((1.-dexp(-dplus/6.))**2)
        fmuCoefficient(i)=((1.-dexp(-dplus/26.))**2)*
     *                       (1.+4.1/dmax1(ReT(i)**0.75,tiny))
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          i2=BiTau(i,j)
          i3=BjTau(i,j)
          i1=NBFaceOwner(i2,i3)
c
          dNorm=WallDistance(i1)
c
          WallVelocity=TangentialVelocity(i1)
          TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
          uTau=dsqrt(TauWall/BDensity(i2,i3))
          uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
          dplus=BWallDistance(i,j)*BDensity(i,j)*uTau/BViscosity(i,j)
c        
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=((1.-dexp(-dplus/26.))**2)*
     *                       (1.+4.1/dmax1(BReT(i,j)**0.75,tiny))
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKEkasagiCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,zc,BFaceAreanx,BFaceAreany,BFaceAreanz
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,ReT,BReT,cmu25
      use WallDistance1, only: iTau,jTau,BiTau,BjTau,WallDistance,
     *                         BWallDistance
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,i1,i2,i3
      double precision :: nx,ny,nz,dNorm,uWall1,vWall1,wWall1,
     *                    WallVelocity,TauWall,uTau,dplus,f2Old,fmuOld,
     *                    dotproduct
c********************************************************************************************
      interface
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
      do i=1,NumberOfElements
c
        i2=iTau(i)
        i3=jTau(i)
        i1=NBFaceOwner(i2,i3)
c
        dNorm=WallDistance(i1)
c
        WallVelocity=TangentialVelocity(i1)
        TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
        uTau=dsqrt(TauWall/BDensity(i2,i3))
        uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
        dplus=WallDistance(i)*Density(i)*uTau/Viscosity(i)
c        
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=(1.-(2./9.)*dexp(-ReT(i)*ReT(i)/36.))*
     *                         ((1.-dexp(-dplus/5.))**2)
        fmuCoefficient(i)=(1.-dexp(-dplus/70.))*
     *                       (1.+3.45/dmax1(ReT(i)**0.5,tiny))
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          i2=BiTau(i,j)
          i3=BjTau(i,j)
          i1=NBFaceOwner(i2,i3)
c
          dNorm=WallDistance(i1)
c
          WallVelocity=TangentialVelocity(i1)
          TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
          uTau=dsqrt(TauWall/BDensity(i2,i3))
          uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
          dplus=BWallDistance(i,j)*BDensity(i,j)*uTau/BViscosity(i,j)
c        
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=(1.-dexp(-dplus/70.))*
     *                   (1.+3.45/dmax1(BReT(i,j)**0.5,tiny))
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKEchcCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Turbulence1, only: fmuCoefficient,f2Coefficient,
     *                       BfmuCoefficient,ReT,BReT
      use WallDistance1, only: WallDistance,BWallDistance
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: dNorm,ReK,f2Old,fmuOld
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        dNorm=WallDistance(i)
c        
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
        ReK=dNorm*dsqrt(dmax1(TurbulentKE(i),0.))*
     *                               Density(i)/Viscosity(i)
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=(1.-0.01*dexp(-ReT(i)*ReT(i)))*
     *                                 (1.-dexp(-0.0631*ReK))
        fmuCoefficient(i)=((1.-dexp(-0.0215*ReK))**2)*
     *                       (1.+31.66/dmax1(ReT(i)**1.25,tiny))
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          dNorm=BWallDistance(i,j)
c        
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
          ReK=dNorm*dsqrt(dmax1(BTurbulentKE(i,j),0.))*
     *                         BDensity(i,j)/BViscosity(i,j)
c
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=((1.-dexp(-0.0215*ReK))**2)*
     *                       (1.+31.66/dmax1(BReT(i,j)**1.25,tiny))
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKELaunderSharmaCoefficientsOld
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,zc,BFaceAreanx,BFaceAreany,BFaceAreanz
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,sqrtTurbulentKE,
     *                       BsqrtTurbulentKE,StrainRate,BStrainRate,
     *                       sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz,
     *                       BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,
     *                       SRateGradx,SRateGrady,SRateGradz,
     *                       BSRateGradx,BSRateGrady,BSRateGradz,
     *                       ReT,BReT
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity,
     *                               TurbulentViscosity
      use Constants1, only: tiny
      use WallDistance1, only: iTau,jTau
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,i1,i2
      double precision :: term1,term2,term3,term4,sum,f2Old,fmuOld,
     *                    nx,ny,nz
c********************************************************************************************
c
      interface
c--------------------------------------------------------------------------------
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
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      do i=1,NumberOfElements
c
        sqrtTurbulentKE(i)=dsqrt(dmax1(TurbulentKE(i),0.))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BsqrtTurbulentKE(i,j)=dsqrt(dmax1(BTurbulentKE(i,j),0.))
c
        enddo
      enddo
c
      Variable='tke05'
      call Gradient(Variable,2,sqrtTurbulentKE,sqrtTKEGradx,
     *              sqrtTKEGrady,sqrtTKEGradz,BsqrtTurbulentKE,
     *         BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,2,.false.,1)
c
      Variable='StrainRate'
      call Gradient(Variable,2,StrainRate,SRateGradx,SRateGrady,
     *              SRateGradz,BStrainRate,BSRateGradx,BSRateGrady,
     *              BSRateGradz,2,.false.,1)
c
      do i=1,NumberOfElements
c
        i1=iTau(i)
        i2=jTau(i)
c
        nx=BFaceAreanx(i1,i2)
        ny=BFaceAreany(i1,i2)
        nz=BFaceAreanz(i1,i2)
c
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=1.-0.3*dexp(-ReT(i)*ReT(i))
        fmuCoefficient(i)=dexp(-3.4/((1.+ReT(i)/50.)*(1.+ReT(i)/50.)))
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
        LTKE(i)=2.*Viscosity(i)*((sqrtTKEGradx(i)*nx+
     *                     sqrtTKEGrady(i)*ny+sqrtTKEGradz(i)*nz)**2)
c        LTKE(i)=2.*Viscosity(i)*((TKEGradx(i)*nx+TKEGrady(i)*ny+
c     *              TKEGradz(i)*nz)**2)/(4.*dmax1(TurbulentKE(i),tiny))
c
        LTED(i)=2.*Viscosity(i)*TurbulentViscosity(i)/Density(i)
        LTED(i)=LTED(i)*((SRateGradx(i)*nx+
     *                     SRateGrady(i)*ny+SRateGradz(i)*nz)**2)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=
     *           dexp(-3.4/((1.+BReT(i,j)/50.)*(1.+BReT(i,j)/50.)))
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKELaunderSharmaCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient,
     *                 MethodCalcGradientMomentum,nIterGradientMomentum,
     *                 LimitGradientMomentum,LimitGradientMomentumMethod
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,zc,BFaceAreanx,BFaceAreany,BFaceAreanz
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,sqrtTurbulentKE,
     *                       BsqrtTurbulentKE,StrainRate,BStrainRate,
     *                       sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz,
     *                       BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,
     *                       SRateGradx,SRateGrady,SRateGradz,
     *                       BSRateGradx,BSRateGrady,BSRateGradz,
     *                       ReT,BReT
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      uVelGrad2x,uVelGrad2y,uVelGrad2z,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BuVelGrad2x,BuVelGrad2y,BuVelGrad2z,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      vVelGrad2x,vVelGrad2y,vVelGrad2z,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BvVelGrad2x,BvVelGrad2y,BvVelGrad2z,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      wVelGrad2x,wVelGrad2y,wVelGrad2z,
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
     *                      BwVelGrad2x,BwVelGrad2y,BwVelGrad2z,
     *                      uvVelGradxy,uwVelGradxz,BuvVelGradxy,
     *                      BuwVelGradxz
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity,
     *                               TurbulentViscosity
      use Constants1, only: tiny
      use WallDistance1, only: iTau,jTau
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,i1,i2
      double precision :: term1,term2,term3,term4,sum,f2Old,fmuOld,
     *                    nx,ny,nz
c********************************************************************************************
c
      interface
c--------------------------------------------------------------------------------
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
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      do i=1,NumberOfElements
c
        sqrtTurbulentKE(i)=dsqrt(dmax1(TurbulentKE(i),0.))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BsqrtTurbulentKE(i,j)=dsqrt(dmax1(BTurbulentKE(i,j),0.))
c
        enddo
      enddo
c
      Variable='tke05'
      call Gradient(Variable,2,sqrtTurbulentKE,sqrtTKEGradx,
     *              sqrtTKEGrady,sqrtTKEGradz,BsqrtTurbulentKE,
     *         BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,2,.false.,1)
c
      do i=1,NumberOfElements
c
        i1=iTau(i)
        i2=jTau(i)
c
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=1.-0.3*dexp(-ReT(i)*ReT(i))
        fmuCoefficient(i)=dexp(-3.4/((1.+ReT(i)/50.)*(1.+ReT(i)/50.)))
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
        LTKE(i)=2.*Viscosity(i)*(sqrtTKEGradx(i)**2+
     *                     sqrtTKEGrady(i)**2+sqrtTKEGradz(i)**2)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=
     *           dexp(-3.4/((1.+BReT(i,j)/50.)*(1.+BReT(i,j)/50.)))
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c
        enddo
      enddo
c
      Variable='velx'
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGradx,uVelGrad2x,uvVelGradxy,uwVelGradxz,BuVelGradx,
     *       BuVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
c
      LTED=uVelGrad2x+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGrady,uvVelGradxy,uVelGrad2y,uwVelGradxz,BuVelGrady,
     *       BuvVelGradxy,BuVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
c
      LTED=LTED+uVelGrad2y+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelGradz,uvVelGradxy,uwVelGradxz,uVelGrad2z,BuVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BuVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
c
      LTED=LTED+uVelGrad2z+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGradx,vVelGrad2x,uvVelGradxy,uwVelGradxz,BvVelGradx,
     *       BvVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
c
      LTED=LTED+vVelGrad2x+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGrady,uvVelGradxy,vVelGrad2y,uwVelGradxz,BvVelGrady,
     *       BuvVelGradxy,BvVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
c
      LTED=LTED+vVelGrad2y+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelGradz,uvVelGradxy,uwVelGradxz,vVelGrad2z,BvVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BvVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
c
      LTED=LTED+vVelGrad2z+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGradx,wVelGrad2x,uvVelGradxy,uwVelGradxz,BwVelGradx,
     *       BwVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                             LimitGradientMomentumMethod)
c
      LTED=LTED+wVelGrad2x+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGrady,uvVelGradxy,wVelGrad2y,uwVelGradxz,BwVelGrady,
     *       BuvVelGradxy,BwVelGrad2y,BuwVelGradxz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                             LimitGradientMomentumMethod)
c
      LTED=LTED+wVelGrad2y+uvVelGradxy+uwVelGradxz
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelGradz,uvVelGradxy,uwVelGradxz,wVelGrad2z,BwVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BwVelGrad2z,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                            LimitGradientMomentumMethod)
c
      LTED=LTED+wVelGrad2z+uvVelGradxy+uwVelGradxz
c
      LTED=LTED*LTED
c
      LTED=LTED*2.*Viscosity*TurbulentViscosity/Density
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKEChienCoefficients
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,zc,BFaceAreanx,BFaceAreany,BFaceAreanz
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,ReT,BReT,cmu25
      use WallDistance1, only: iTau,jTau,BiTau,BjTau,WallDistance,
     *                         BWallDistance
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,i1,i2,i3
      double precision :: nx,ny,nz,dNorm,uWall1,vWall1,wWall1,
     *                    WallVelocity,TauWall,uTau,dplus,f2Old,
     *                    fmuOld,dotproduct
c********************************************************************************************
      interface
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
      do i=1,NumberOfElements
c
        i2=iTau(i)
        i3=jTau(i)
        i1=NBFaceOwner(i2,i3)
c
        dNorm=WallDistance(i1)
c
        WallVelocity=TangentialVelocity(i1)
        TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
        uTau=dsqrt(TauWall/BDensity(i2,i3))
        uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
        dplus=WallDistance(i)*Density(i)*uTau/Viscosity(i)
c        
        ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        f2Old=f2Coefficient(i)
        fmuOld=fmuCoefficient(i)
        f2Coefficient(i)=1.-(0.4/1.8)*dexp(-ReT(i)*ReT(i)/36.)
        fmuCoefficient(i)=1.-dexp(-0.0115*dplus)
        f2Coefficient(i)=urff2Coefficient*f2Coefficient(i)+
     *                               (1.-urff2Coefficient)*f2Old
        fmuCoefficient(i)=urffmuCoefficient*fmuCoefficient(i)+
     *                               (1.-urffmuCoefficient)*fmuOld
        LTKE(i)=2.*Viscosity(i)/(WallDistance(i)**2)
        LTED(i)=LTKE(i)*dexp(-dplus/2.)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          i2=BiTau(i,j)
          i3=BjTau(i,j)
          i1=NBFaceOwner(i2,i3)
c
          dNorm=WallDistance(i1)
c
          WallVelocity=TangentialVelocity(i1)
          TauWall=BViscosity(i2,i3)*WallVelocity/dNorm
c
          uTau=dsqrt(TauWall/BDensity(i2,i3))
          uTau=dmax1(cmu25*dsqrt(dmax1(TurbulentKE(i1),0.)),uTau)
          dplus=BWallDistance(i,j)*BDensity(i,j)*uTau/BViscosity(i,j)
c        
          fmuOld=BfmuCoefficient(i,j)
          BfmuCoefficient(i,j)=1.-dexp(-0.0115*dplus)
          BfmuCoefficient(i,j)=urffmuCoefficient*BfmuCoefficient(i,j)+
     *                                     (1.-urffmuCoefficient)*fmuOld
c        
          BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKOmega2006Coefficients
c
c#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Turbulence1, only: alphaStar0,alpha0,ReW,ReB,ReK,cmu,betta0,
     *                      alfaStar,alfa,bettaStar,ReT,BReT,
     *                      BalfaStar,Balfa,BbettaStar
      use Variables1, only: TurbulentKE,TurbulentOmega,
     *                      BTurbulentKE,BTurbulentOmega
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        ReT(i)=Density(i)*dmax1(TurbulentKE(i),0.)
        ReT(i)=ReT(i)/(Viscosity(i)*dmax1(TurbulentOmega(i),tiny))
        alfaStar(i)=(alphaStar0+ReT(i)/ReK)/(1.+ReT(i)/ReK)

        alfa(i)=(13./25.)*(alpha0+ReT(i)/ReW)/(1.+ReT(i)/ReW)
        alfa(i)=alfa(i)/alfaStar(i)
c
        bettaStar(i)=
     *        cmu*(100.*betta0/27.+(ReT(i)/ReB)**4)/(1.+(ReT(i)/ReB)**4)
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BReT(i,j)=BDensity(i,j)*dmax1(BTurbulentKE(i,j),0.)
          BReT(i,j)=BReT(i,j)/
     *           (BViscosity(i,j)*dmax1(BTurbulentOmega(i,j),tiny))
          BalfaStar(i,j)=(alphaStar0+BReT(i,j)/ReK)/(1.+BReT(i,j)/ReK)

          Balfa(i,j)=(13./25.)*(alpha0+BReT(i,j)/ReW)/(1.+BReT(i,j)/ReW)
          Balfa(i,j)=Balfa(i,j)/BalfaStar(i,j)
c
          BbettaStar(i,j)=cmu*(100.*betta0/27.+
     *              (BReT(i,j)/ReB)**4)/(1.+(BReT(i,j)/ReB)**4)
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateV2fCoefficients
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TurbulentED,TurbulentV2,
     *                      BTurbulentKE,BTurbulentED,BTurbulentV2
      use PhysicalProperties1, only: Density,BDensity,
     *                               Viscosity,BViscosity
      use Turbulence1, only: cmu,ce1,alpha,CetaV2f,CLV2f,
     *                       Ce1Coefficient,StrainRate,BStrainRate,
     *                       TScale,BTScale,LScale,BLScale
      use Constants1, only: sqrt6,tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: tke,ted,tv2,StRate,term1,term2,term3
c********************************************************************************************
c
c---- Calculate Ce1 coefficient
c
      do i=1,NumberOfElements
c
        tke=dsqrt(dmax1(TurbulentKE(i),0.))
        tv2=dmax1(dsqrt(dmax1(TurbulentV2(i),0.)),tiny)
        Ce1Coefficient(i)=ce1*(1.+0.05*tke/tv2)
c
      enddo
c
c---- Calculate time scale
c
      do i=1,NumberOfElements
c
        tke=dmax1(TurbulentKE(i),0.)
        ted=dmax1(TurbulentED(i),tiny)
        tv2=dmax1(TurbulentV2(i),tiny)
        StRate=dmax1(StrainRate(i),tiny)
c
        term1=tke/ted
        term2=6.*dsqrt(Viscosity(i)/(Density(i)*ted))
        term3=(alpha/(sqrt6*cmu))*tke/(tv2*StRate)
        TScale(i)=dmax1(dmin1(term1,term3),term2)
c        TScale(i)=dmax1(TScale(i),tiny)
        TScale(i)=dmax1(TScale(i),1.d-10)
c
      enddo
c
      do i=1,NumberOfBCSets      
        do j=1,NBFaces(i)     
c      
          tke=dmax1(BTurbulentKE(i,j),0.)
          ted=dmax1(BTurbulentED(i,j),tiny)
          tv2=dmax1(BTurbulentV2(i,j),tiny)
          StRate=dmax1(BStrainRate(i,j),tiny)
c
          term1=tke/ted
          term2=6.*dsqrt(BViscosity(i,j)/(BDensity(i,j)*ted))
          term3=(alpha/(sqrt6*cmu))*tke/(tv2*StRate)
          BTScale(i,j)=dmax1(dmin1(term1,term3),term2)
c          BTScale(i,j)=dmax1(BTScale(i,j),tiny)
          BTScale(i,j)=dmax1(BTScale(i,j),1.d-10)
c
        enddo
      enddo
c
c---- Calculate length scale
c
      do i=1,NumberOfElements
c
        tke=dmax1(TurbulentKE(i),0.)
        ted=dmax1(TurbulentED(i),tiny)
        tv2=dmax1(TurbulentV2(i),tiny)
        StRate=dmax1(StrainRate(i),tiny)
c
        term1=tke**1.5/ted
        term2=(1./(sqrt6*cmu))*(tke**1.5)/(tv2*StRate)
        term3=CetaV2f*((Viscosity(i)/Density(i))**3/ted)**0.25
        LScale(i)=CLV2f*dmax1(dmin1(term1,term2),term3)
c        LScale(i)=dmax1(LScale(i),tiny)
        LScale(i)=dmax1(LScale(i),1.d-10)
c
      enddo
c
      do i=1,NumberOfBCSets      
        do j=1,NBFaces(i)     
c
          tke=dmax1(BTurbulentKE(i,j),0.)
          ted=dmax1(BTurbulentED(i,j),tiny)
          tv2=dmax1(BTurbulentV2(i,j),tiny)
          StRate=dmax1(BStrainRate(i,j),tiny)
c
          term1=tke**1.5/ted
          term2=(1./(sqrt6*cmu))*(tke**1.5)/(tv2*StRate)
          term3=CetaV2f*((BViscosity(i,j)/BDensity(i,j))**3/ted)**0.25
          BLScale(i,j)=CLV2f*dmax1(dmin1(term1,term2),term3)
c          BLScale(i,j)=dmax1(BLScale(i,j),tiny)
          BLScale(i,j)=dmax1(BLScale(i,j),1.d-10)
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateZetafCoefficients
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TurbulentED,TurbulentZeta,
     *                      BTurbulentKE,BTurbulentED,BTurbulentZeta
      use PhysicalProperties1, only: Density,BDensity,
     *                               Viscosity,BViscosity
      use Turbulence1, only: cmu,ce1,alpha,CetaZeta,CLZeta,
     *                       Ce1Coefficient,StrainRate,BStrainRate,
     *                       TScale,BTScale,LScale,BLScale
      use Constants1, only: sqrt6,tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: tke,ted,tv2,StRate,term1,term2,term3
c********************************************************************************************
c
c---- Calculate Ce1 coefficient
c
      do i=1,NumberOfElements
c
        tv2=dmax1(dsqrt(dmax1(TurbulentZeta(i),0.)),tiny)
        Ce1Coefficient(i)=ce1*(1.+0.05/tv2)
c
      enddo
c
c---- Calculate time scale
c
      do i=1,NumberOfElements
c
        tke=dmax1(TurbulentKE(i),0.)
        ted=dmax1(TurbulentED(i),tiny)
        tv2=dmax1(TurbulentZeta(i),tiny)
        StRate=dmax1(StrainRate(i),tiny)
c
        term1=tke/ted
        term2=6.*dsqrt(Viscosity(i)/(Density(i)*ted))
        term3=(alpha/(sqrt6*cmu))/(tv2*StRate)
        TScale(i)=dmax1(dmin1(term1,term3),term2)
c        TScale(i)=dmax1(TScale(i),tiny)
        TScale(i)=dmax1(TScale(i),1.d-10)
c
      enddo
c
      do i=1,NumberOfBCSets      
        do j=1,NBFaces(i)     
c      
          tke=dmax1(BTurbulentKE(i,j),0.)
          ted=dmax1(BTurbulentED(i,j),tiny)
          tv2=dmax1(BTurbulentZeta(i,j),tiny)
          StRate=dmax1(BStrainRate(i,j),tiny)
c
          term1=tke/ted
          term2=6.*dsqrt(BViscosity(i,j)/(BDensity(i,j)*ted))
          term3=(alpha/(sqrt6*cmu))/(tv2*StRate)
          BTScale(i,j)=dmax1(dmin1(term1,term3),term2)
c          BTScale(i,j)=dmax1(BTScale(i,j),tiny)
          BTScale(i,j)=dmax1(BTScale(i,j),1.d-10)
c
        enddo
      enddo
c
c---- Calculate length scale
c
      do i=1,NumberOfElements
c
        tke=dmax1(TurbulentKE(i),0.)
        ted=dmax1(TurbulentED(i),tiny)
        tv2=dmax1(TurbulentZeta(i),tiny)
        StRate=dmax1(StrainRate(i),tiny)
c
        term1=tke**1.5/ted
        term2=(1./(sqrt6*cmu))*(tke**0.5)/(tv2*StRate)
        term3=CetaZeta*((Viscosity(i)/Density(i))**3/ted)**0.25
        LScale(i)=CLZeta*dmax1(dmin1(term1,term2),term3)
c        LScale(i)=dmax1(LScale(i),tiny)
        LScale(i)=dmax1(LScale(i),1.d-10)
c
      enddo
c
      do i=1,NumberOfBCSets      
        do j=1,NBFaces(i)     
c
          tke=dmax1(BTurbulentKE(i,j),0.)
          ted=dmax1(BTurbulentED(i,j),tiny)
          tv2=dmax1(BTurbulentZeta(i,j),tiny)
          StRate=dmax1(BStrainRate(i,j),tiny)
c
          term1=tke**1.5/ted
          term2=(1./(sqrt6*cmu))*(tke**0.5)/(tv2*StRate)
          term3=CetaZeta*((BViscosity(i,j)/BDensity(i,j))**3/ted)**0.25
          BLScale(i,j)=CLZeta*dmax1(dmin1(term1,term2),term3)
c          BLScale(i,j)=dmax1(BLScale(i,j),tiny)
          BLScale(i,j)=dmax1(BLScale(i,j),1.d-10)
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateTurbulenceInterpolationFactors
c
C#############################################################################################
c
      use User0, only: TurbulenceModel,RotationCurvatureMethod,
     *                 LKOmegaSSTRotationCurvatureCorrection,LRotation,
     *                 LRough,WallBCModel,GrainSize
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Variables1, only: TurbulentKE,BTurbulentKE,TurbulentOmega,
     *                      BTurbulentOmega,TKEGradx,TKEGrady,TKEGradz,
     *                      TOmegaGradx,TomegaGrady,TomegaGradz,
     *                      BTKEGradx,BTKEGrady,BTKEGradz,
     *                      BTOmegaGradx,BTomegaGrady,BTomegaGradz,
     *                      MaterialDerivative
      use PhysicalProperties1, only: Density,Viscosity,BDensity,
     *                               BViscosity
      use WallDistance1, only: WallDistance,iTau
      use Turbulence1, only: cmu,sigTED2,F1factor,BF1factor,
     *                       F2factor,BF2factor,F3factor,BF3factor,
     *                       F4factor,BF4factor,Crc,StrainRate,
     *                       BStrainRate,Vorticity,BVorticity,
     *                       cr1,cr2,cr3,fr1Coefficient,
     *                       S11,S12,S13,S22,S23,S33,
     *                       BS11,BS12,BS13,BS22,BS23,BS33,
     *                       S11Old,S12Old,S13Old,S22Old,S23Old,S33Old,
     *                       S11OldOld,S12OldOld,S13OldOld,S22OldOld,
     *                       S23OldOld,S33OldOld,DS11Dt,DS12Dt,DS13Dt,
     *                       DS22Dt,DS23Dt,DS33Dt,
     *                       W11,W12,W13,W22,W23,W33
      use Coriolis1, only: AngularVelocityX,AngularVelocityY,
     *                     AngularVelocityZ
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      character*10 Variable
      double precision :: tke1,tomega1,dNorm,gama1,cdKOmega,
     *                    factor1,factor2,factor3,factor4,factor5,
     *                    factor6,gama2,Smagnitude,Wmagnitude,Ratio,
     *                    Richardson,rstar,rtelda,term1,term2,term3,
     *                    D2,D3,a1,a2,a3,a4,a5,a6,a7,a8,a9
c
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE CalculateMaterialDerivative(Variable,FiTemp,
     *                                 BFiTemp,FiTempOld,FiTempOldOld)
c--------------------------------------------------------------------------------
          character*10 Variable
          double precision FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
          double precision, dimension(:) :: FiTemp
          double precision, dimension(:) :: FiTempOld
          double precision, dimension(:) :: FiTempOldOld
          double precision, dimension(:,:) :: BFiTemp
c--------------------------------------------------------------------------------
        end SUBROUTINE CalculateMaterialDerivative
c--------------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------------
c
      do i=1,NumberOfElements
c      
        tke1=dmax1(TurbulentKE(i),0.)
        tomega1=dmax1(TurbulentOmega(i),tiny)
        dNorm=WallDistance(i)
c
        if(LRough) then
          if(TurbulenceModel.eq.'komegasst'.or.
     *                        TurbulenceModel.eq.'komegabsl') then
            if(WallBCModel.eq.'nikuradse'.or.WallBCModel.eq.'colebrook')
     *       dNorm=dNorm+0.03*GrainSize(iTau(i))
          endif
        endif
c
        factor1=dsqrt(tke1)/(cmu*tomega1*dNorm)
        factor2=500.*Viscosity(i)/(Density(i)*(dNorm**2)*tomega1)
        factor4=2.*Density(i)*sigTED2/tomega1
        factor5=TKEGradx(i)*TOmegaGradx(i)+
     *          TKEGrady(i)*TOmegaGrady(i)+TKEGradz(i)*TOmegaGradz(i)
        factor6=factor4*factor5
        cdKOmega=dmax1(factor6,1.d-10)
        factor3=4.*Density(i)*sigTED2*tke1/(cdKOmega*(dNorm**2))
        gama1=dmin1(dmin1(dmax1(factor1,factor2),factor3),10.)
        F1factor(i)=dtanh(gama1**4)
c
      enddo
c
      do i=1,NumberOfBCSets      
c
        do j=1,NBFaces(i)     
c      
          k=NBFaceOwner(i,j)
c      
          tke1=dmax1(BTurbulentKE(i,j),0.)
          tomega1=dmax1(BTurbulentOmega(i,j),tiny)
          dNorm=WallDistance(k)
c
          if(LRough) then
            if(TurbulenceModel.eq.'komegasst'.or.
     *                        TurbulenceModel.eq.'komegabsl') then
              if(WallBCModel.eq.'nikuradse'.or.
     *                             WallBCModel.eq.'colebrook')
     *           dNorm=dNorm+0.03*GrainSize(iTau(i))
            endif
          endif
c
          factor1=dsqrt(tke1)/(cmu*tomega1*dNorm)
          factor2=500.*BViscosity(i,j)/
     *                   (BDensity(i,j)*(dNorm**2)*tomega1)
          factor4=2.*BDensity(i,j)*sigTED2/tomega1
          factor5=BTKEGradx(i,j)*BTOmegaGradx(i,j)+
     *                    BTKEGrady(i,j)*BTOmegaGrady(i,j)+
     *                              BTKEGradz(i,j)*BTOmegaGradz(i,j)
          factor6=factor4*factor5
          cdKOmega=dmax1(factor6,1.d-10)
          factor3=4.*BDensity(i,j)*sigTED2*tke1/(cdKOmega*(dNorm**2))
          gama1=dmin1(dmin1(dmax1(factor1,factor2),factor3),10.)
          BF1factor(i,j)=dtanh(gama1**4)
c            
        enddo
c
      enddo
c
      if(TurbulenceModel.eq.'sstgamaretheta'.or.
     *                         TurbulenceModel.eq.'sstgama') then
c
        do i=1,NumberOfElements
c      
          tke1=dmax1(TurbulentKE(i),0.)
          dNorm=WallDistance(i)
c
          factor1=Density(i)*dNorm*dsqrt(tke1)/Viscosity(i)
          factor2=dexp(-(factor1/120.)**8)
          F1factor(i)=dmax1(F1factor(i),factor2)
c
        enddo
c
        do i=1,NumberOfBCSets      
c
          do j=1,NBFaces(i)     
c      
            k=NBFaceOwner(i,j)
c      
            tke1=dmax1(BTurbulentKE(i,j),0.)
            dNorm=WallDistance(k)
c
            factor1=BDensity(i,j)*dNorm*dsqrt(tke1)/BViscosity(i,j)
            factor2=dexp(-(factor1/120.)**8)
            BF1factor(i,j)=dmax1(BF1factor(i,j),factor2)
c            
          enddo
c
        enddo
c
      endif
c
      if(TurbulenceModel.eq.'komegasst'.or.
     *         TurbulenceModel.eq.'sstgamaretheta'.or.
     *                       TurbulenceModel.eq.'sstgama') then
c
        do i=1,NumberOfElements
c
          tke1=dmax1(TurbulentKE(i),0.)
          tomega1=dmax1(TurbulentOmega(i),tiny)
          dNorm=WallDistance(i)
c
          if(LRough.and.TurbulenceModel.eq.'komegasst') then
            if(WallBCModel.eq.'nikuradse'.or.WallBCModel.eq.'colebrook')
     *         dNorm=dNorm+0.03*GrainSize(iTau(i))
          endif
c
          factor1=2.*dsqrt(tke1)/(cmu*tomega1*dNorm)
          factor2=500.*Viscosity(i)/
     *                 (Density(i)*(dNorm**2)*tomega1)
          gama2=dmin1(dmax1(factor1,factor2),100.)
          F2factor(i)=dtanh(gama2**2)
c
          if(LRough.and.TurbulenceModel.eq.'komegasst') then
c
            if(WallBCModel.eq.'wilcox') then
c
              gama2=150.*Viscosity(i)/(Density(i)*tomega1*(dNorm**2))
              gama2=dmin1(gama2,10.)
              F3factor(i)=1.-dtanh(gama2**4)
c
            else
c
              F3factor(i)=1. 
c
            endif
c
          else
c
            F3factor(i)=1. 
c
          endif
c          
          if(LKOmegaSSTRotationCurvatureCorrection.and.
     *              RotationCurvatureMethod.eq.'hellsten') then
            Smagnitude=StrainRate(i)
            Wmagnitude=Vorticity(i)
            Ratio=Wmagnitude/dmax1(Smagnitude,tiny)
c
            Richardson=Ratio*(Ratio-1.)
            F4factor(i)=1./(1.+Crc*Richardson)
c
          else
c
            F4factor(i)=1.
c
          endif
c          
        enddo
c
        do i=1,NumberOfBCSets      
c
          do j=1,NBFaces(i)     
c      
            k=NBFaceOwner(i,j)
c
            tke1=dmax1(BTurbulentKE(i,j),0.)
            tomega1=dmax1(BTurbulentOmega(i,j),tiny)
            dNorm=WallDistance(k)
c
            if(LRough.and.TurbulenceModel.eq.'komegasst') then
              if(WallBCModel.eq.'nikuradse'.or.
     *                           WallBCModel.eq.'colebrook')
     *         dNorm=dNorm+0.03*GrainSize(iTau(k))
            endif
c
            factor1=2.*dsqrt(tke1)/(cmu*tomega1*dNorm)
            factor2=500.*BViscosity(i,j)/
     *                        (BDensity(i,j)*(dNorm**2)*tomega1)
            gama2=dmin1(dmax1(factor1,factor2),100.)
            BF2factor(i,j)=dtanh(gama2**2)
c
            if(LRough.and.TurbulenceModel.eq.'komegasst') then
c
              if(WallBCModel.eq.'wilcox') then
c
                gama2=150.*BViscosity(i,j)/
     *                      (BDensity(i,j)*tomega1*(dNorm**2))
                gama2=dmin1(gama2,10.)
                BF3factor(i,j)=1.-dtanh(gama2**4)
c
              else
c
                BF3factor(i,j)=1. 
c
              endif
c
            else
c
              BF3factor(i,j)=1. 
c
            endif
c          
            if(LKOmegaSSTRotationCurvatureCorrection.and.
     *              RotationCurvatureMethod.eq.'hellsten') then
              Smagnitude=BStrainRate(i,j)
              Wmagnitude=BVorticity(i,j)
              Ratio=Wmagnitude/dmax1(Smagnitude,tiny)
c
              Richardson=Ratio*(Ratio-1.)
              BF4factor(i,j)=1./(1.+Crc*Richardson)
c
            else
c
              BF4factor(i,j)=1.
c
            endif
          enddo
c
        enddo
c
c--- Calculate the rotation/curvature correction term
c
        if(TurbulenceModel.eq.'komegasst'.and.
     *          LKOmegaSSTRotationCurvatureCorrection) then
c
          if(RotationCurvatureMethod.eq.'spalartshur') then
c
            allocate(MaterialDerivative(NumberOfElements))
            Variable='S11'
            call CalculateMaterialDerivative(Variable,S11,
     *                               BS11,S11Old,S11OldOld)
            DS11Dt=MaterialDerivative
            Variable='S12'
            call CalculateMaterialDerivative(Variable,S12,
     *                               BS12,S12Old,S12OldOld)
            DS12Dt=MaterialDerivative
            Variable='S13'
            call CalculateMaterialDerivative(Variable,S13,
     *                               BS13,S13Old,S13OldOld)
            DS13Dt=MaterialDerivative
            Variable='S22'
            call CalculateMaterialDerivative(Variable,S22,
     *                               BS22,S22Old,S22OldOld)
            DS22Dt=MaterialDerivative
            Variable='S23'
            call CalculateMaterialDerivative(Variable,S23,
     *                               BS23,S23Old,S23OldOld)
            DS23Dt=MaterialDerivative
            Variable='S33'
            call CalculateMaterialDerivative(Variable,S33,
     *                               BS33,S33Old,S33OldOld)
            DS33Dt=MaterialDerivative
            deallocate(MaterialDerivative)
c           
            term1=0.d0
            term2=0.d0
            term3=0.d0
c
            do i=1,NumberOfElements      
c
              rstar=StrainRate(i)/dmax1(Vorticity(i),tiny)         
c
              a1=W12(i)*S12(i)+W13(i)*S13(i)
              a2=W12(i)*S22(i)+W13(i)*S23(i)
              a3=W12(i)*S23(i)+W13(i)*S33(i)
              a4=-W12(i)*S11(i)+W23(i)*S13(i)
              a5=-W12(i)*S12(i)+W23(i)*S23(i)
              a6=-W12(i)*S13(i)+W23(i)*S33(i)
              a7=-W13(i)*S11(i)-W23(i)*S12(i)
              a8=-W13(i)*S12(i)-W23(i)*S22(i)
              a9=-W13(i)*S13(i)-W23(i)*S23(i)
c
              term1=a1*DS11Dt(i)+a2*DS12Dt(i)+a3*DS13Dt(i)+
     *          a4*DS12Dt(i)+a5*DS22Dt(i)+a6*DS23Dt(i)+
     *          a7*DS13Dt(i)+a8*DS23Dt(i)+a9*DS33Dt(i)
c
              if(LRotation) then
c
                term2=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a3*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a4*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a7*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a8*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
                term3=
     *          a1*(S13(i)*AngularVelocityY-S12(i)*AngularVelocityZ)+
     *          a2*(S11(i)*AngularVelocityZ-S13(i)*AngularVelocityX)+
     *          a3*(S12(i)*AngularVelocityX-S11(i)*AngularVelocityY)+
     *          a4*(S23(i)*AngularVelocityY-S22(i)*AngularVelocityZ)+
     *          a5*(S12(i)*AngularVelocityZ-S23(i)*AngularVelocityX)+
     *          a6*(S22(i)*AngularVelocityX-S12(i)*AngularVelocityY)+
     *          a7*(S33(i)*AngularVelocityY-S23(i)*AngularVelocityZ)+
     *          a8*(S13(i)*AngularVelocityZ-S33(i)*AngularVelocityX)+
     *          a9*(S23(i)*AngularVelocityX-S13(i)*AngularVelocityY)
              endif
c
              D2=dmax1(StrainRate(i)**2,0.09*TurbulentOmega(i)**2)
              D3=D2**1.5
c
              rtelda=2.*(term1+term2+term3)/dmax1(Vorticity(i)*D3,tiny)
c
              fr1Coefficient(i)=(1.+cr1)*(2.*rstar/(1+rstar))*
     *                         (1.-cr3*datan(cr2*rtelda))-cr1
c
              fr1Coefficient(i)=dmax1(dmin1(fr1Coefficient(i),1.25),0.)

            enddo
c
          endif
c
        endif
c
      endif
c
      return
      end      