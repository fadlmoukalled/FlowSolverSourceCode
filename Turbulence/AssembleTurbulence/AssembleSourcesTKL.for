c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentKLSources
c
C#############################################################################################
c
      use User0, only: MethodCalcGradientTurbulentKL,
     *                 nIterGradientTurbulentKL,
     *                 LimitGradientTurbulentKL,
     *                 LimitGradientTurbulentKLMethod
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BgDiff,BFaceTx,BFaceTy,BFaceTz
      use WallDistance1, only: WallDistance,BWallDistance
      use Variables1, only: TurbulentKE,TurbulenceProduction,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      TurbulentOmega,BTurbulentKE,
     *                      TKEGradx,TKEGrady,TKEGradz,
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
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BeDiffCoefficient
      use Turbulence1, only: cappa,cmu75,
     *                       c11,c12,xi1,xi2,xi3,cd1,ModelNumber,
     *                       sqrtTurbulentKE,BsqrtTurbulentKE,
     *                       sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz,
     *                       BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz,
     *                       StrainRate,TurbulentViscosityTl,
     *                       ProductionKL,ReT,SourceRbp,SourceRnat
      use Constants1, only: tiny
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j
      character*10 Variable
      double precision :: term1,dNorm,gamf,FluxCflocal,
     *                    FluxFflocal,FluxVflocal,
     *                    dfidxTf,dfidyTf,dfidzTf,
     *                    Uprime1,Uprime2,tke1,Lvk,LvkMin,LvkMax,fp,
     *                    Cphi1,Cphi2,xi,fphi,SourceDL
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
      end interface
c********************************************************************************************
c
      CalculateTurbulentKLSources: select case(ModelNumber)
c
c-----------------------------------------------------------------------------
        case(18) CalculateTurbulentKLSources           !kklmodel
c-----------------------------------------------------------------------------
c
c---- Calculate the second Gradients of velocity components
c
          Variable='velgrad'
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      uVelGradx,uVelGrad2x,uvVelGradxy,uwVelGradxz,BuVelGradx,
     *       BuVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                              LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      uVelGrady,uvVelGradxy,uVelGrad2y,uwVelGradxz,BuVelGrady,
     *       BuvVelGradxy,BuVelGrad2y,BuwVelGradxz,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                              LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      uVelGradz,uvVelGradxy,uwVelGradxz,uVelGrad2z,BuVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BuVelGrad2z,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                              LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      vVelGradx,vVelGrad2x,uvVelGradxy,uwVelGradxz,BvVelGradx,
     *       BvVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                              LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      vVelGrady,uvVelGradxy,vVelGrad2y,uwVelGradxz,BvVelGrady,
     *       BuvVelGradxy,BvVelGrad2y,BuwVelGradxz,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                              LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      vVelGradz,uvVelGradxy,uwVelGradxz,vVelGrad2z,BvVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BvVelGrad2z,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                              LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      wVelGradx,wVelGrad2x,uvVelGradxy,uwVelGradxz,BwVelGradx,
     *       BwVelGrad2x,BuvVelGradxy,BuwVelGradxz,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                             LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      wVelGrady,uvVelGradxy,wVelGrad2y,uwVelGradxz,BwVelGrady,
     *       BuvVelGradxy,BwVelGrad2y,BuwVelGradxz,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                             LimitGradientTurbulentKLMethod)
          call Gradient(Variable,MethodCalcGradientTurbulentKL,
     *      wVelGradz,uvVelGradxy,uwVelGradxz,wVelGrad2z,BwVelGradz,
     *       BuvVelGradxy,BuwVelGradxz,BwVelGrad2z,
     *        nIterGradientTurbulentKL,LimitGradientTurbulentKL,
     *                            LimitGradientTurbulentKLMethod)
c
          do i=1,NumberOfElements    
c
            Uprime1=StrainRate(i)
c
            Uprime2=dmax1(dsqrt(
     *         (uVelGrad2x(i)+uVelGrad2y(i)+uVelGrad2z(i))**2+
     *         (vVelGrad2x(i)+vVelGrad2y(i)+vVelGrad2z(i))**2+
     *         (wVelGrad2x(i)+wVelGrad2y(i)+wVelGrad2z(i))**2),tiny)
c
            Lvk=cappa*Uprime1/Uprime2     
c
            tke1=dmax1(TurbulentKE(i),tiny)       
            LvkMin=dmax1(TurbulentKL(i),0.)/(c11*tke1)
            term1=cmu75*Density(i)*(tke1**2.5)/
     *                              dmax1(TurbulentKL(i),tiny)
            fp=dmin1(dmax1(TurbulenceProduction(i)/term1,0.5),1.0)
            LvkMax=c12*cappa*WallDistance(i)*fp         
c
            if(Lvk.lt.LvkMin) then
              Lvk=LvkMin
            elseif(Lvk.Gt.LvkMax) then
              Lvk=LvkMax
            endif
c
            Cphi1=xi1-xi2*((dmax1(TurbulentKL(i),tiny)/(tke1*LvK))**2)
            Cphi2=xi3
            xi=Density(i)*WallDistance(i)*dsqrt(0.3*tke1)/
     *                                     (20.*Viscosity(i))
            fphi=(1.+cd1*xi)/(1.+xi**4)
c
            term1=Cphi1*dmax1(TurbulentKL(i),0.)*
     *                              TurbulenceProduction(i)/tke1
            FluxTE(i)=FluxTE(i)-term1*Volume(i)
c
            term1=Cphi2*Density(i)*(tke1**1.5)/
     *                    dmax1(TurbulentKL(i),tiny)+
     *                        6*Viscosity(i)*fphi/(WallDistance(i)**2)
            FluxCE(i)=FluxCE(i)+term1*Volume(i)
            FluxTE(i)=FluxTE(i)+term1*Volume(i)*dmax1(TurbulentKL(i),0.)
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
            BTurbulentKL(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTurbulentKLGradx(i2,i3)
            dfidyTf=BTurbulentKLGrady(i2,i3)
            dfidzTf=BTurbulentKLGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKL(i1)+
     *                 FluxFflocal*BTurbulentKL(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------
        case(21) CalculateTurbulentKLSources           !transitional k-kl-w model
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            sqrtTurbulentKE(i)=dsqrt(dmax1(TurbulentKL(i),0.))
c
          enddo          
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BsqrtTurbulentKE(i,j)=dsqrt(dmax1(BTurbulentKL(i,j),0.))
c
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
            ProductionKL(i)=
     *           TurbulentViscosityTl(i)*StrainRate(i)*StrainRate(i)
c
            SourceDL=2.*(Viscosity(i)/Density(i))*
     *        (sqrtTKEGradx(i)**2+sqrtTKEGrady(i)**2+sqrtTKEGradz(i)**2)
c
            FluxCE(i)=FluxCE(i)+Density(i)*Volume(i)*(SourceRbp(i)+
     *                SourceRnat(i)+SourceDL/dmax1(Turbulentkl(i),tiny))
            FluxTE(i)=FluxTE(i)-
     *             Density(i)*Volume(i)*(ProductionKL(i)-
     *                   dmax1(Turbulentkl(i),0.)*(SourceRbp(i)+
     *                                    SourceRnat(i))-SourceDL)
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
            BTurbulentKL(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BTurbulentKLGradx(i2,i3)
            dfidyTf=BTurbulentKLGrady(i2,i3)
            dfidzTf=BTurbulentKLGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *           dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*TurbulentKL(i1)+
     *                FluxFflocal*BTurbulentKL(i2,i3)+FluxVflocal
c
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTurbulentKLSources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end