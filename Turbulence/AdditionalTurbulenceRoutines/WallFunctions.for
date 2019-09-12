c
c#############################################################################################
c
      SUBROUTINE ScalarWallFunctions
c
C#############################################################################################
c
      use User0, only: WallFunctionKind,LRough
c
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      if(LRough) then
c
        call CalculateRoughccScalar
c
        if(WallFunctionKind.eq.'equilibrium') then
c
          call EquilibriumScalarRoughWallFunctions
c
        elseif(WallFunctionKind.eq.'nonequilibrium') then
c
          call NonEquilibriumScalarRoughWallFunctions
c      
        endif
c      
      else
c
        if(WallFunctionKind.eq.'equilibrium') then
c
          call EquilibriumScalarWallFunctions
c
        elseif(WallFunctionKind.eq.'nonequilibrium') then
c
          call NonEquilibriumScalarWallFunctions
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
      SUBROUTINE EquilibriumScalarWallFunctions
c
C#############################################################################################
c
      use User0, only: ScalarWallFunctionType
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use BoundaryFluxes, only: ScalarFlux
      use Scalar1, only: Scalar,BScalar
      use BoundaryConditionsScalar1, only: wallTypeScalar
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,
     *                               BeDiffCoefficient,
     *                               DiffusionCoefficient,
     *                               BDiffusionCoefficient,
     *                               SpecificHeatScalar
      use Turbulence1, only: yplus,WallViscosityS,yplusS,
     *                       cappa,cc,sigScalar,uTau,ctrans
      use Constants1, only: tiny
      use WallDistance1, only: WallDistance
      use Scalar2
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,iScalar
      double precision :: dNorm,Sch,Schp,Schstar,factor,gamma,Sstar
c********************************************************************************************
c
      iScalar=iScalarVariable
c
      if(ScalarWallFunctionType.eq.'standard') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                         DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2
c          
          Schstar=Sch*yplus(i)
          if(Schstar.lt.ctrans) then 
c
            Sstar=Schstar
c
          else
c
            Sstar=2.12*dlog(Schstar)+Schp
c
          endif       
c
          WallViscosityS(i,iScalar)=uTau(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'scalable') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2
c          
          Schstar=dmax1(Sch*yplus(i),ctrans)
c
          Sstar=2.12*dlog(Schstar)+Schp
c
          WallViscosityS(i,iScalar)=uTau(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'automatic') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2
c          
          Schstar=dmax1(Sch*yplus(i),ctrans)
          gamma=0.01*Schstar**4/(1.+5.*Sch*Sch*Schstar)+tiny
c
          Sstar=Schstar*dexp(-gamma)+
     *            (2.12*dlog(Schstar)+Schp)*dexp(-1./gamma)
c
          WallViscosityS(i,iScalar)=uTau(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
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
      SUBROUTINE NonEquilibriumScalarWallFunctions
c
C#############################################################################################
c
      use User0, only: ScalarWallFunctionType
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use BoundaryFluxes, only: ScalarFlux
      use Scalar1, only: Scalar,BScalar
      use BoundaryConditionsScalar1, only: wallTypeScalar
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,
     *                               BeDiffCoefficient,
     *                               DiffusionCoefficient,
     *                               BDiffusionCoefficient,
     *                               SpecificHeatScalar
      use Turbulence1, only: ystar,WallViscosityS,yplusS,
     *                       cappa,cc,sigScalar,ustar,ctrans
      use Constants1, only: tiny
      use WallDistance1, only: WallDistance
      use Scalar2
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,iScalar
      double precision :: dNorm,Sch,Schp,Schstar,factor,gamma,Sstar
c********************************************************************************************
c
      iScalar=iScalarVariable
c
      if(ScalarWallFunctionType.eq.'standard') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2
c          
          Schstar=Sch*ystar(i)
          if(Schstar.lt.ctrans) then 
c
            Sstar=Schstar
c
          else
c
            Sstar=2.12*dlog(Schstar)+Schp
c
          endif       
c
          WallViscosityS(i,iScalar)=ustar(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'scalable') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2
c          
          Schstar=dmax1(Sch*ystar(i),ctrans)
c
          Sstar=2.12*dlog(Schstar)+Schp
c
          WallViscosityS(i,iScalar)=ustar(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'automatic') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2
c          
          Schstar=dmax1(Sch*ystar(i),ctrans)
          gamma=0.01*Schstar**4/(1.+5.*Sch*Sch*Schstar)+tiny
c
          Sstar=Schstar*dexp(-gamma)+
     *            (2.12*dlog(Schstar)+Schp)*dexp(-1./gamma)
c
          WallViscosityS(i,iScalar)=ustar(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
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
      SUBROUTINE EquilibriumScalarRoughWallFunctions
c
C#############################################################################################
c
      use User0, only: ScalarWallFunctionType
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use BoundaryFluxes, only: ScalarFlux
      use Scalar1, only: Scalar,BScalar
      use BoundaryConditionsScalar1, only: wallTypeScalar
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,
     *                               BeDiffCoefficient,
     *                               DiffusionCoefficient,
     *                               BDiffusionCoefficient,
     *                               SpecificHeatScalar
      use Turbulence1, only: yplus,WallViscosityS,yplusS, cappa,cc,
     *                       sigScalar,uTau,ctrans,RoughccS
      use Constants1, only: tiny
      use WallDistance1, only: WallDistance
      use Scalar2
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,iScalar
      double precision :: dNorm,Sch,Schp,Schstar,factor,gamma,Sstar
c********************************************************************************************
c
      iScalar=iScalarVariable
c
      if(ScalarWallFunctionType.eq.'standard') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2-RoughccS(i,iScalar)
c          
          Schstar=Sch*yplus(i)
          if(Schstar.lt.ctrans) then 
c
            Sstar=Schstar
c
          else
c
            Sstar=2.12*dlog(Schstar)+Schp
c
          endif       
c
          WallViscosityS(i,iScalar)=uTau(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'scalable') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2-RoughccS(i,iScalar)
c          
          Schstar=dmax1(Sch*yplus(i),ctrans)
c
          Sstar=2.12*dlog(Schstar)+Schp
c
          WallViscosityS(i,iScalar)=uTau(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'automatic') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2-RoughccS(i,iScalar)
c          
          Schstar=dmax1(Sch*yplus(i),ctrans)
          gamma=0.01*Schstar**4/(1.+5.*Sch*Sch*Schstar)+tiny
c
          Sstar=Schstar*dexp(-gamma)+
     *            (2.12*dlog(Schstar)+Schp)*dexp(-1./gamma)
c
          WallViscosityS(i,iScalar)=uTau(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
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
      SUBROUTINE NonEquilibriumScalarRoughWallFunctions
c
C#############################################################################################
c
      use User0, only: ScalarWallFunctionType
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use BoundaryFluxes, only: ScalarFlux
      use Scalar1, only: Scalar,BScalar
      use BoundaryConditionsScalar1, only: wallTypeScalar
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,
     *                               BeDiffCoefficient,
     *                               DiffusionCoefficient,
     *                               BDiffusionCoefficient,
     *                               SpecificHeatScalar
      use Turbulence1, only: ystar,WallViscosityS,yplusS,cappa,cc,
     *                       sigScalar,ustar,ctrans,RoughccS
      use Constants1, only: tiny
      use WallDistance1, only: WallDistance
      use Scalar2
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,iScalar
      double precision :: dNorm,Sch,Schp,Schstar,factor,gamma,Sstar
c********************************************************************************************
c
      iScalar=iScalarVariable
c
      if(ScalarWallFunctionType.eq.'standard') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2-RoughccS(i,iScalar)
c          
          Schstar=Sch*ystar(i)
          if(Schstar.lt.ctrans) then 
c
            Sstar=Schstar
c
          else
c
            Sstar=2.12*dlog(Schstar)+Schp
c
          endif       
c
          WallViscosityS(i,iScalar)=ustar(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'scalable') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2-RoughccS(i,iScalar)
c          
          Schstar=dmax1(Sch*ystar(i),ctrans)
c
          Sstar=2.12*dlog(Schstar)+Schp
c
          WallViscosityS(i,iScalar)=ustar(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      elseif(ScalarWallFunctionType.eq.'automatic') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
c
          Sch=Viscosity(i1)*SpecificHeatScalar(i1,iScalar)/
     *                      DiffusionCoefficient(i1,iScalar)
          Schp=(3.85*(Sch**(1./3.))-1.3)**2-RoughccS(i,iScalar)
c          
          Schstar=dmax1(Sch*ystar(i),ctrans)
          gamma=0.01*Schstar**4/(1.+5.*Sch*Sch*Schstar)+tiny
c
          Sstar=Schstar*dexp(-gamma)+
     *            (2.12*dlog(Schstar)+Schp)*dexp(-1./gamma)
c
          WallViscosityS(i,iScalar)=ustar(i)*
     *                 SpecificHeatScalar(i1,iScalar)*
     *                       BDensity(i2,i3)*dNorm/dmax1(Sstar,tiny)
          BeDiffCoefficient(i2,i3)=dmax1(WallViscosityS(i,iScalar),
     *                           BDiffusionCoefficient(i2,i3,iScalar))
c
          if(wallTypeScalar(i2,i3,iScalar).eq.'dirichlet') then
c
            ScalarFlux(i2,i3,iScalar)=BeDiffCoefficient(i2,i3)*
     *              (BScalar(i2,i3,iScalar)-Scalar(i1,iScalar))/dNorm
c
          elseif(wallTypeScalar(i2,i3,iScalar).eq.'vonneumann') then
c
            BScalar(i2,i3,iScalar)=Scalar(i1,iScalar)+
     *         ScalarFlux(i2,i3,iScalar)*dNorm/BeDiffCoefficient(i2,i3)
c
          endif
c
        enddo
c
      endif
c
      return
      end
c
c
c#############################################################################################
c
      SUBROUTINE EnergyWallFunctions
c
C#############################################################################################
c
      use User0, only: WallFunctionKind,LRough
c
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      if(LRough) then
c
        call CalculateRoughccEnergy
c
        if(WallFunctionKind.eq.'equilibrium') then
c
          call EquilibriumEnergyRoughWallFunctions
c
        elseif(WallFunctionKind.eq.'nonequilibrium') then
c
          call NonEquilibriumEnergyRoughWallFunctions
c      
        endif
c      
      else
c
        if(WallFunctionKind.eq.'equilibrium') then
c
          call EquilibriumEnergyWallFunctions
c
        elseif(WallFunctionKind.eq.'nonequilibrium') then
c
          call NonEquilibriumEnergyWallFunctions
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
      SUBROUTINE EquilibriumEnergyWallFunctions
c
C#############################################################################################
c
      use User0, only: EnergyWallFunctionType,Lcompressible,
     *                 LSolveTurbulenceDissipationRate
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use BoundaryFluxes, only: HeatFlux
      use Variables1, only: TurbulentKE,BTemperature,Temperature,
     *                      uVelocity,BuVelocity,vVelocity,BvVelocity,
     *                      wVelocity,BwVelocity
      use BoundaryConditions1, only: wallTypeEnergy
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,
     *                               BeDiffCoefficient,SpecificHeat,
     *                               Conductivity,BConductivity,
     *                               BSpecificHeat
      use Turbulence1, only: yplus,WallViscosityT,cmu25,ctrans,cappa,
     *                       cc,yplusT,sigT,uTau,uplus,TauWall,
     *                       WallViscosity
      use Constants1, only: tiny
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,factor,PrYplus,gamma,Tplus,
     *                    Pr,Prp,WallVelocity,
     *                    WallVelocity2,u1plusC,u2plusC,yplus1,
     *                    aCoeff,bCoeff,rCoeff,hCoeff,u1plus,
     *                    u2plus,Tcoefficient,tempdiff,TplusL,
     *                    TplusT,dotproduct
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
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(EnergyWallFunctionType.eq.'standard') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2
c
            PrYplus=Pr*yplus(i)        
c
            if(PrYplus.lt.ctrans) then
c
              Tplus=PrYplus
c
            else
c
              Tplus=2.12*dlog(PrYplus)+Prp
c
            endif
c
            WallViscosityT(i)=utau(i)*BDensity(i2,i3)*
     *                      SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'scalable') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2
c
            PrYplus=dmax1(Pr*yplus(i),ctrans)        
c
            Tplus=2.12*dlog(PrYplus)+Prp
c
            WallViscosityT(i)=uTau(i)*BDensity(i2,i3)*
     *                                 SpecificHeat(i1)*dNorm/Tplus
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2
c
            PrYplus=Pr*yplus(i)        
            gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)+tiny
c
            Tplus=PrYplus*dexp(-gamma)+
     *            (2.12*dlog(PrYplus)+Prp)*dexp(-1./gamma)
c
            WallViscosityT(i)=uTau(i)*BDensity(i2,i3)*
     *                          SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        endif
c
      elseif(Lcompressible) then
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
          Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
          PrYplus=Pr*yplus(i)        
c
          Prp=(3.85*(Pr**(1./3.))-1.3)**2+2.12*dlog(Pr)
          gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)
c           
c--- Calculate the compressible U+ for the log region 
c
          u1plusC=yplus(i)
c
          yplus1=dmax1(yplus(i),dexp(-Prp/2.12))
c
          if(yplus1.gt.0.) then
c            
            u2plusC=(2.12*dlog(yplus1)+Prp)/sigT
c
          else
c
            u2plusC=0.
c
          endif
c
c--- Perform inverse Van Driest velocity transformation
c
          aCoeff=HeatFlux(i2,i3)/dmax1(TauWall(i),tiny)
          bCoeff=2.*BSpecificHeat(i2,i3)*BTemperature(i2,i3)/sigT
c
          rCoeff=uTau(i)/dsqrt(bCoeff)
          hCoeff=aCoeff/dmax1(uTau(i),tiny)
c
          u1plus=(dsin(rCoeff*u1plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u1plusC))         
          u2plus=(dsin(rCoeff*u2plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u2plusC))         
c
          Tplus=Pr*u1plus*dexp(-gamma)+
     *                    sigT*u2plus*dexp(-1./(gamma+tiny))
c
          factor=dmin1(PrYplus/(Tplus+tiny),(Viscosity(i1)+
     *             TurbulentViscosity(i1))/(Viscosity(i1)+tiny)) 
c
          Tcoefficient=uTau(i)*factor*BDensity(i2,i3)*
     *                         SpecificHeat(i1)/(PRYPLUS+tiny)
          tempdiff=BTemperature(i2,i3)-Temperature(i1)
c
          if(tempdiff.eq.0.) tempdiff=tiny
c
          factor=dmax1(0.d0,
     *            1.-sigT*WallVelocity2/(2.*SpecificHeat(i1)*tempdiff))
          WallViscosityT(i)=factor*Tcoefficient*dNorm     !/SpecificHeat(i1)
c
          BeDiffCoefficient(i2,i3)=WallViscosityT(i)
c
          if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
            HeatFlux(i2,i3)=Tcoefficient*
     *              (BTemperature(i2,i3)-Temperature(i1)-
     *                  0.5*sigT*WallVelocity2/SpecificHeat(i1))
c
          elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
            BTemperature(i2,i3)=Temperature(i1)+
     *                       HeatFlux(i2,i3)/Tcoefficient+
     *                          0.5*sigT*WallVelocity2/SpecificHeat(i1)
c
          endif
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
      SUBROUTINE NonEquilibriumEnergyWallFunctions
c
C#############################################################################################
c
      use User0, only: EnergyWallFunctionType,Lcompressible,
     *                 LSolveTurbulenceDissipationRate
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use BoundaryFluxes, only: HeatFlux
      use Variables1, only: TurbulentKE,BTemperature,Temperature,
     *                      uVelocity,BuVelocity,vVelocity,BvVelocity,
     *                      wVelocity,BwVelocity
      use BoundaryConditions1, only: wallTypeEnergy
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,Conductivity,
     *                               BeDiffCoefficient,SpecificHeat,
     *                               BConductivity,BSpecificHeat
      use Turbulence1, only: WallViscosityT,cmu25,ctrans,cappa,
     *                       cc,yplusT,sigT,ustar,TauWall,ystar
      use Constants1, only: tiny
      use MultiGrid2, only: nIter
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,rstar,factor,PrYplus,gamma,Tplus,Tstar,
     *                    Pr,Prp,WallVelocity,
     *                    WallVelocity2,u1plusC,u2plusC,yplus1,
     *                    aCoeff,bCoeff,rCoeff,hCoeff,u1plus,
     *                    u2plus,Tcoefficient,tempdiff,dotproduct
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
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(EnergyWallFunctionType.eq.'standard') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2
c
            PrYplus=Pr*ystar(i)        
c
            if(PrYplus.lt.ctrans) then
c
              Tplus=PrYplus
c
            else
c
              Tplus=2.12*dlog(PrYplus)+Prp
c
            endif
c
            WallViscosityT(i)=ustar(i)*BDensity(i2,i3)*
     *                      SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'scalable') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2
c
            PrYplus=dmax1(Pr*ystar(i),ctrans)        
c
            Tplus=2.12*dlog(PrYplus)+Prp
c
            WallViscosityT(i)=ustar(i)*BDensity(i2,i3)*
     *                                 SpecificHeat(i1)*dNorm/Tplus
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2
c
            PrYplus=Pr*ystar(i)        
            gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)+tiny
c
            Tplus=PrYplus*dexp(-gamma)+
     *            (2.12*dlog(PrYplus)+Prp)*dexp(-1./gamma)
c
            WallViscosityT(i)=ustar(i)*BDensity(i2,i3)*
     *                          SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        endif
c
      elseif(Lcompressible) then
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
          Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
          PrYplus=Pr*ystar(i)        
c
          Prp=(3.85*(Pr**(1./3.))-1.3)**2+2.12*dlog(Pr)
          gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)
c           
c--- Calculate the compressible U+ for the log region 
c
          u1plusC=ystar(i)
c
          yplus1=dmax1(ystar(i),dexp(-Prp/2.12))
c
          if(yplus1.gt.0.) then
c            
            u2plusC=(2.12*dlog(yplus1)+Prp)/sigT
c
          else
c
            u2plusC=0.
c
          endif
c
c--- Perform inverse Van Driest velocity transformation
c
          aCoeff=HeatFlux(i2,i3)/dmax1(TauWall(i),tiny)
          bCoeff=2.*BSpecificHeat(i2,i3)*BTemperature(i2,i3)/sigT
c
          rCoeff=ustar(i)/dsqrt(bCoeff)
          hCoeff=aCoeff/dmax1(ustar(i),tiny)
c
          u1plus=(dsin(rCoeff*u1plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u1plusC))         
          u2plus=(dsin(rCoeff*u2plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u2plusC))         
c
          Tplus=Pr*u1plus*dexp(-gamma)+
     *                    sigT*u2plus*dexp(-1./(gamma+tiny))
c
          factor=dmin1(PrYplus/(Tplus+tiny),(Viscosity(i1)+
     *             TurbulentViscosity(i1))/(Viscosity(i1)+tiny)) 
c
          Tcoefficient=ustar(i)*factor*BDensity(i2,i3)*
     *                         SpecificHeat(i1)/(PRYPLUS+tiny)
          tempdiff=BTemperature(i2,i3)-Temperature(i1)
c
          if(tempdiff.eq.0.) tempdiff=tiny
c
          factor=dmax1(0.d0,
     *            1.-sigT*WallVelocity2/(2.*SpecificHeat(i1)*tempdiff))
          WallViscosityT(i)=factor*Tcoefficient*dNorm     !/SpecificHeat(i1)
c
          BeDiffCoefficient(i2,i3)=WallViscosityT(i)
c
          if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
            HeatFlux(i2,i3)=Tcoefficient*
     *              (BTemperature(i2,i3)-Temperature(i1)-
     *                  0.5*sigT*WallVelocity2/SpecificHeat(i1))
c
          elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
            BTemperature(i2,i3)=Temperature(i1)+
     *                       HeatFlux(i2,i3)/Tcoefficient+
     *                          0.5*sigT*WallVelocity2/SpecificHeat(i1)
c
          endif
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
      SUBROUTINE EquilibriumEnergyRoughWallFunctions
c
C#############################################################################################
c
      use User0, only: EnergyWallFunctionType,Lcompressible,
     *                 LSolveTurbulenceDissipationRate
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use BoundaryFluxes, only: HeatFlux
      use Variables1, only: TurbulentKE,BTemperature,Temperature,
     *                      uVelocity,BuVelocity,vVelocity,BvVelocity,
     *                      wVelocity,BwVelocity
      use BoundaryConditions1, only: wallTypeEnergy
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,
     *                               BeDiffCoefficient,SpecificHeat,
     *                               Conductivity,BConductivity,
     *                               BSpecificHeat
      use Turbulence1, only: yplus,WallViscosityT,cmu25,ctrans,cappa,
     *                       cc,yplusT,sigT,uTau,uplus,TauWall,
     *                       WallViscosity,RoughccE
      use Constants1, only: tiny
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,factor,PrYplus,gamma,Tplus,
     *                    Pr,Prp,WallVelocity,
     *                    WallVelocity2,u1plusC,u2plusC,yplus1,
     *                    aCoeff,bCoeff,rCoeff,hCoeff,u1plus,
     *                    u2plus,Tcoefficient,tempdiff,TplusL,
     *                    TplusT,dotproduct
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
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(EnergyWallFunctionType.eq.'standard') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2-RoughccE(i)
c
            PrYplus=Pr*yplus(i)        
c
            if(PrYplus.lt.ctrans) then
c
              Tplus=PrYplus
c
            else
c
              Tplus=2.12*dlog(PrYplus)+Prp
c
            endif
c
            WallViscosityT(i)=utau(i)*BDensity(i2,i3)*
     *                      SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'scalable') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2-RoughccE(i)
c
            PrYplus=dmax1(Pr*yplus(i),ctrans)        
c
            Tplus=2.12*dlog(PrYplus)+Prp
c
            WallViscosityT(i)=uTau(i)*BDensity(i2,i3)*
     *                                 SpecificHeat(i1)*dNorm/Tplus
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2-RoughccE(i)
c
            PrYplus=Pr*yplus(i)        
            gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)+tiny
c
            Tplus=PrYplus*dexp(-gamma)+
     *            (2.12*dlog(PrYplus)+Prp)*dexp(-1./gamma)
c
            WallViscosityT(i)=uTau(i)*BDensity(i2,i3)*
     *                          SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        endif
c
      elseif(Lcompressible) then
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
          Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
          PrYplus=Pr*yplus(i)        
c
          Prp=(3.85*(Pr**(1./3.))-1.3)**2+2.12*dlog(Pr)-RoughccE(i)
          gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)
c           
c--- Calculate the compressible U+ for the log region 
c
          u1plusC=yplus(i)
c
          yplus1=dmax1(yplus(i),dexp(-Prp/2.12))
c
          if(yplus1.gt.0.) then
c            
            u2plusC=(2.12*dlog(yplus1)+Prp)/sigT
c
          else
c
            u2plusC=0.
c
          endif
c
c--- Perform inverse Van Driest velocity transformation
c
          aCoeff=HeatFlux(i2,i3)/dmax1(TauWall(i),tiny)
          bCoeff=2.*BSpecificHeat(i2,i3)*BTemperature(i2,i3)/sigT
c
          rCoeff=uTau(i)/dsqrt(bCoeff)
          hCoeff=aCoeff/dmax1(uTau(i),tiny)
c
          u1plus=(dsin(rCoeff*u1plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u1plusC))         
          u2plus=(dsin(rCoeff*u2plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u2plusC))         
c
          Tplus=Pr*u1plus*dexp(-gamma)+
     *                    sigT*u2plus*dexp(-1./(gamma+tiny))
c
          factor=dmin1(PrYplus/(Tplus+tiny),(Viscosity(i1)+
     *             TurbulentViscosity(i1))/(Viscosity(i1)+tiny)) 
c
          Tcoefficient=uTau(i)*factor*BDensity(i2,i3)*
     *                         SpecificHeat(i1)/(PRYPLUS+tiny)
          tempdiff=BTemperature(i2,i3)-Temperature(i1)
c
          if(tempdiff.eq.0.) tempdiff=tiny
c
          factor=dmax1(0.d0,
     *            1.-sigT*WallVelocity2/(2.*SpecificHeat(i1)*tempdiff))
          WallViscosityT(i)=factor*Tcoefficient*dNorm     !/SpecificHeat(i1)
c
          BeDiffCoefficient(i2,i3)=WallViscosityT(i)
c
          if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
            HeatFlux(i2,i3)=Tcoefficient*
     *              (BTemperature(i2,i3)-Temperature(i1)-
     *                  0.5*sigT*WallVelocity2/SpecificHeat(i1))
c
          elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
            BTemperature(i2,i3)=Temperature(i1)+
     *                       HeatFlux(i2,i3)/Tcoefficient+
     *                          0.5*sigT*WallVelocity2/SpecificHeat(i1)
c
          endif
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
      SUBROUTINE NonEquilibriumEnergyRoughWallFunctions
c
C#############################################################################################
c
      use User0, only: EnergyWallFunctionType,Lcompressible,
     *                 LSolveTurbulenceDissipationRate
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use BoundaryFluxes, only: HeatFlux
      use Variables1, only: TurbulentKE,BTemperature,Temperature,
     *                      uVelocity,BuVelocity,vVelocity,BvVelocity,
     *                      wVelocity,BwVelocity
      use BoundaryConditions1, only: wallTypeEnergy
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,Conductivity,
     *                               BeDiffCoefficient,SpecificHeat,
     *                               BConductivity,BSpecificHeat
      use Turbulence1, only: WallViscosityT,cmu25,ctrans,cappa,
     *                       cc,yplusT,sigT,ustar,TauWall,ystar,
     *                       RoughccE
      use Constants1, only: tiny
      use MultiGrid2, only: nIter
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,rstar,factor,PrYplus,gamma,Tplus,Tstar,
     *                    Pr,Prp,WallVelocity,
     *                    WallVelocity2,u1plusC,u2plusC,yplus1,
     *                    aCoeff,bCoeff,rCoeff,hCoeff,u1plus,
     *                    u2plus,Tcoefficient,tempdiff,dotproduct
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
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(EnergyWallFunctionType.eq.'standard') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2-RoughccE(i)
c
            PrYplus=Pr*ystar(i)        
c
            if(PrYplus.lt.ctrans) then
c
              Tplus=PrYplus
c
            else
c
              Tplus=2.12*dlog(PrYplus)+Prp
c
            endif
c
            WallViscosityT(i)=ustar(i)*BDensity(i2,i3)*
     *                      SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'scalable') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2-RoughccE(i)
c
            PrYplus=dmax1(Pr*ystar(i),ctrans)        
c
            Tplus=2.12*dlog(PrYplus)+Prp
c
            WallViscosityT(i)=ustar(i)*BDensity(i2,i3)*
     *                                 SpecificHeat(i1)*dNorm/Tplus
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        elseif(EnergyWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
c
            Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
            Prp=(3.85*(Pr**(1./3.))-1.3)**2-RoughccE(i)
c
            PrYplus=Pr*ystar(i)        
            gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)+tiny
c
            Tplus=PrYplus*dexp(-gamma)+
     *            (2.12*dlog(PrYplus)+Prp)*dexp(-1./gamma)
c
            WallViscosityT(i)=ustar(i)*BDensity(i2,i3)*
     *                          SpecificHeat(i1)*dNorm/dmax1(Tplus,tiny)
c
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosityT(i),BConductivity(i2,i3))
c
            if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
              HeatFlux(i2,i3)=BeDiffCoefficient(i2,i3)*
     *              (BTemperature(i2,i3)-Temperature(i1))/dNorm
c
            elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
              BTemperature(i2,i3)=Temperature(i1)+
     *                HeatFlux(i2,i3)*dNorm/BeDiffCoefficient(i2,i3)
c
            endif
c
          enddo
c
        endif
c
      elseif(Lcompressible) then
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
          Pr=Viscosity(i1)*SpecificHeat(i1)/Conductivity(i1)
          PrYplus=Pr*ystar(i)        
c
          Prp=(3.85*(Pr**(1./3.))-1.3)**2+2.12*dlog(Pr)-RoughccE(i)
          gamma=0.01*PrYplus**4/(1.+5.*Pr*Pr*PrYplus)
c           
c--- Calculate the compressible U+ for the log region 
c
          u1plusC=ystar(i)
c
          yplus1=dmax1(ystar(i),dexp(-Prp/2.12))
c
          if(yplus1.gt.0.) then
c            
            u2plusC=(2.12*dlog(yplus1)+Prp)/sigT
c
          else
c
            u2plusC=0.
c
          endif
c
c--- Perform inverse Van Driest velocity transformation
c
          aCoeff=HeatFlux(i2,i3)/dmax1(TauWall(i),tiny)
          bCoeff=2.*BSpecificHeat(i2,i3)*BTemperature(i2,i3)/sigT
c
          rCoeff=ustar(i)/dsqrt(bCoeff)
          hCoeff=aCoeff/dmax1(ustar(i),tiny)
c
          u1plus=(dsin(rCoeff*u1plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u1plusC))         
          u2plus=(dsin(rCoeff*u2plusC))/rCoeff-
     *                 hCoeff*(1.-dcos(rCoeff*u2plusC))         
c
          Tplus=Pr*u1plus*dexp(-gamma)+
     *                    sigT*u2plus*dexp(-1./(gamma+tiny))
c
          factor=dmin1(PrYplus/(Tplus+tiny),(Viscosity(i1)+
     *             TurbulentViscosity(i1))/(Viscosity(i1)+tiny)) 
c
          Tcoefficient=ustar(i)*factor*BDensity(i2,i3)*
     *                         SpecificHeat(i1)/(PRYPLUS+tiny)
          tempdiff=BTemperature(i2,i3)-Temperature(i1)
c
          if(tempdiff.eq.0.) tempdiff=tiny
c
          factor=dmax1(0.d0,
     *            1.-sigT*WallVelocity2/(2.*SpecificHeat(i1)*tempdiff))
          WallViscosityT(i)=factor*Tcoefficient*dNorm     !/SpecificHeat(i1)
c
          BeDiffCoefficient(i2,i3)=WallViscosityT(i)
c
          if(wallTypeEnergy(i2,i3).eq.'dirichlet') then
c
            HeatFlux(i2,i3)=Tcoefficient*
     *              (BTemperature(i2,i3)-Temperature(i1)-
     *                  0.5*sigT*WallVelocity2/SpecificHeat(i1))
c
          elseif(wallTypeEnergy(i2,i3).eq.'vonneumann') then
c
            BTemperature(i2,i3)=Temperature(i1)+
     *                       HeatFlux(i2,i3)/Tcoefficient+
     *                          0.5*sigT*WallVelocity2/SpecificHeat(i1)
c
          endif
c
        enddo
c
      endif
c
      return
      end
c
c
c#############################################################################################
c
      SUBROUTINE MomentumWallFunctions
c
C#############################################################################################
c
      use User0, only: LSolveTurbulenceKineticEnergy,WallFunctionKind,
     *                 LRough,LReadOldSolution
      use MultiGrid2, only: nIter
c
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      if(LRough) then
c
        if(nIter.eq.1.and..not.LReadOldSolution) call CalculateuTau
        call CalculateKsPlus
        call CalculateRoughccMomentum
c
        if(WallFunctionKind.eq.'equilibrium') then
c
          call EquilibriumMomentumRoughWallFunctions
c
        elseif(WallFunctionKind.eq.'nonequilibrium') then
c
          if(.not.LSolveTurbulenceKineticEnergy) call ExtractTurbulentKE
          call NonEquilibriumMomentumRoughWallFunctions
c      
        endif
c      
      else
c
        if(WallFunctionKind.eq.'equilibrium') then
c
          call EquilibriumMomentumWallFunctions
c
        elseif(WallFunctionKind.eq.'nonequilibrium') then
c
          if(.not.LSolveTurbulenceKineticEnergy) call ExtractTurbulentKE
          call NonEquilibriumMomentumWallFunctions
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
      SUBROUTINE EquilibriumMomentumWallFunctions
c
C#############################################################################################
c
      use User0, only: MomentumWallFunctionType,Lcompressible,
     *                 LReadOldSolution,LSolveTurbulenceDissipationRate
      use Multigrid2, only: nIter
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use Variables1, only: TurbulentKE,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      BModifiedED,Temperature
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,SpecificHeat,
     *                               BeDiffCoefficient
      use Constants1, only: tiny
      use BoundaryFluxes, only: HeatFlux
      use Turbulence1, only: yplus,WallViscosity,duplusdyplus,uTau,
     *                       uplus,cmu25,cappa,cc,ctrans,uStar,ystar,
     *                       Bfv1Coefficient,sigT,TauWall,Stelda
      use WallDistance1, only: WallDistance
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,k
      double precision :: nx,ny,nz,dNorm,cmu25k,yplusmin,rstar,ustar1,
     *                    factor,yplus1,uWall1,vWall1,wWall1,dotproduct,
     *                    WallVelocity,uTauViscous,uTauLog,
     *                    WallVelocity2,TauWallStar,Ec,Bstar,TwallOverT,
     *                    Dstar,Ecert,Ucomp
c
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
      if(nIter.eq.1.and..not.LReadOldSolution) then
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
c
          uTau(i)=dsqrt((Viscosity(i1)+TurbulentViscosity(i1))*
     *                              WallVelocity/(dNorm*Density(i1)))
          yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
c
        enddo
c
      endif
c
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(MomentumWallFunctionType.eq.'standard') then
c
          yplusmin=0.2
c        
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            yplus1=dNorm*BDensity(i2,i3)*uTau(i)/BViscosity(i2,i3)
            yplus(i)=dmax1(yplus1,yplusmin)
c
            if(yplus(i).lt.ctrans) then
c
              uplus(i)=yplus(i)
              uTau(i)=WallVelocity/uplus(i)
              duplusdyplus(i)=1.
c          
            else
c          
              uplus(i)=dlog(yplus(i))/cappa+cc
              uTau(i)=WallVelocity/uplus(i)
              duplusdyplus(i)=1./(cappa*yplus(i))
c
            endif
c
            TauWall(i)=BDensity(i2,i3)*uTau(i)*uTau(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *              dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
          ustar=uTau
          ystar=yplus
c
        elseif(MomentumWallFunctionType.eq.'scalable') then
c
          yplusmin=ctrans+tiny
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
c
            yplus1=dNorm*BDensity(i2,i3)*uTau(i)/BViscosity(i2,i3)
            yplus(i)=dmax1(yplus1,yplusmin)
            duplusdyplus(i)=1./(cappa*yplus(i))
            uplus(i)=dlog(yplus(i))/cappa+cc
            uTau(i)=WallVelocity/uplus(i)
c
            TauWall(i)=BDensity(i2,i3)*uTau(i)*uTau(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *             dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
          ustar=uTau
          ystar=yplus
c
        elseif(MomentumWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            uTauViscous=dsqrt(BViscosity(i2,i3)*WallVelocity/
     *                             (BDensity(i2,i3)*dNorm))
            uTauLog=WallVelocity/(dlog(yplus(i))/cappa+cc)
            uTau(i)=(uTauViscous**4+uTauLog**4)**0.25
c
            uplus(i)=WallVelocity/dmax1(uTau(i),tiny)
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
c        
            factor=0.01*(yplus(i)**4)/(1.+5.*yplus(i))+tiny
            duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*yplus(i)))*dexp(-1./factor)
c
            TauWall(i)=BDensity(i2,i3)*uTau(i)*uTau(i)
c
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny) 
c
            BTurbulentViscosity(i2,i3)=
     *                dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
          ustar=uTau
          ystar=yplus
c
        endif
c
      elseif(Lcompressible) then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
          WallVelocity=TangentialVelocity(i1)
c
          uTauViscous=BViscosity(i2,i3)*WallVelocity/
     *                             (BDensity(i2,i3)*dNorm)
          uTauLog=WallVelocity/(dlog(yplus(i))/cappa+cc)
c
          uTau(i)=(uTauViscous**4+uTauLog**4)**0.25
          yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
c
          factor=0.01*(yplus(i)**4)/(1.+5.*yplus(i))+tiny
          duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*yplus(i)))*dexp(-1./factor)
          uplus(i)= yplus(i)*dexp(-factor)+
     *            ((1./cappa)*dlog(yplus(i)+tiny)+cc)*dexp(-1./factor)
          factor=dmin1(yplus(i)/(uplus(i)+tiny),
     *            (Viscosity(i1)+TurbulentViscosity(i1))/Viscosity(i1))
          WallVelocity2=WallVelocity*WallVelocity
          TauWallStar=Density(i1)*uTau(i)*uTau(i)
          Ec=WallVelocity2/(SpecificHeat(i1)*Temperature(i1)+tiny)
          Bstar=HeatFlux(i2,i3)*WallVelocity/
     *             (SpecificHeat(i1)*Temperature(i1)*TauWallStar+tiny)
          TwallOverT=dmax1(1.+sigT*Bstar+sigT*Ec/2.,0.)
          Dstar=dsqrt(Bstar*Bstar+2.*TwallOverT*Ec/sigT)
          Ecert=(dsqrt(2.*TwallOverT/sigT))*
     *                   (dasin((Bstar+Ec)/(Dstar+tiny))-
     *                                      dasin(Bstar/(Dstar+tiny)))
          Ucomp=Ecert*dsqrt(SpecificHeat(i1)*Temperature(i1))
          TauWall(i)=Density(i1)*uTau(i)*Ucomp*factor/(yplus(i)+tiny)
c            
          WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
          BTurbulentViscosity(i2,i3)=
     *             dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
          BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
        enddo
c
        ustar=uTau
        ystar=yplus
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE NonEquilibriumMomentumWallFunctions
c
C#############################################################################################
c
      use User0, only: MomentumWallFunctionType,Lcompressible,
     *                 LSolveTurbulenceDissipationRate
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use Variables1, only: TurbulentKE,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,Temperature
      use BoundaryFluxes, only: HeatFlux
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,SpecificHeat,
     *                               BeDiffCoefficient
      use Constants1, only: tiny

      use Turbulence1, only: yplus,WallViscosity,duplusdyplus,ustar,
     *                       uplus,cmu25,cappa,cc,ctrans,sigT,uTau,
     *                       TauWall,ystar,a1sst
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,yplusmin,rstar,ustar1,factor,yplus1,
     *                    WallVelocity,uTauViscous,uTauLogStar,
     *                    WallVelocity2,TauWallStar,Ec,Bstar,TwallOverT,
     *                    Dstar,Ecert,Ucomp,dotproduct
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
      if(LSolveTurbulenceDissipationRate.and.
     *             MomentumWallFunctionType.eq.'automatic') then
c     
        print*,
     *  'automatic wall functions are not implemented for k-e models'
        pause
        stop
c
      endif
c
      call CalculateuStar
c
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(MomentumWallFunctionType.eq.'standard') then
c
          yplusmin=0.2
c        
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            yplus1=BDensity(i2,i3)*ustar(i)*dNorm/BViscosity(i2,i3)
            ystar(i)=dmax1(yplusmin,yplus1)
c
            if(ystar(i).lt.ctrans) then
c
              uTau(i)=dsqrt(BViscosity(i2,i3)*
     *                  WallVelocity/(BDensity(i2,i3)*dNorm))
              duplusdyplus(i)=1.
              uplus(i)=ystar(i)
c          
            else
c          
              uplus(i)= dlog(ystar(i))/cappa+cc
              uTau(i)=WallVelocity/uplus(i)
              duplusdyplus(i)=1./(cappa*ystar(i))
c
            endif
c
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
            TauWall(i)=BDensity(i2,i3)*uTau(i)*ustar(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
        elseif(MomentumWallFunctionType.eq.'scalable') then
c
c-------------------------------------------------------------------
c     Smooth walls with scalable wall function - Menter wall functions
c-------------------------------------------------------------------
c
          yplusmin=ctrans+tiny
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
            WallVelocity=TangentialVelocity(i1)
c
            yplus1=Density(i1)*ustar(i)*dNorm/(Viscosity(i1)+tiny)
            ystar(i)=dmax1(yplusmin,yplus1)
            duplusdyplus(i)=1./(cappa*ystar(i))
            uplus(i)=dlog(ystar(i))/cappa+cc
            uTau(i)=WallVelocity/uplus(i)
c
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
            TauWall(i)=BDensity(i2,i3)*uTau(i)*ustar(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
        elseif(MomentumWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            uTauViscous=dsqrt(BViscosity(i2,i3)*
     *                  WallVelocity/(BDensity(i2,i3)*dNorm))
            uTauLogStar=ustar(i)
c
            ustar(i)=(uTauViscous**4+uTauLogStar**4)**0.25
            ystar(i)=BDensity(i2,i3)*ustar(i)*dNorm/BViscosity(i2,i3)
            ystar(i)=dmax1(ystar(i),tiny)
            uTauLogStar=WallVelocity/(dlog(ystar(i))/cappa+cc)
            uTau(i)=(uTauViscous**4+uTauLogStar**4)**0.25
c        
            factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
            duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*ystar(i)))*dexp(-1./factor)
            uplus(i)= ystar(i)*dexp(-factor)+
     *               (dlog(ystar(i))/cappa+cc)*dexp(-1./factor)
c
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
            TauWall(i)=BDensity(i2,i3)*uTau(i)*ustar(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *               dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
        endif
c     
      elseif(Lcompressible) then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
          WallVelocity=TangentialVelocity(i1)
c
          uTauViscous=dsqrt(BViscosity(i2,i3)*
     *                  WallVelocity/(BDensity(i2,i3)*dNorm))
          uTauLogStar=ustar(i)
c
          ustar(i)=(uTauViscous**4+uTauLogStar**4)**0.25
          ystar(i)=BDensity(i2,i3)*ustar(i)*dNorm/BViscosity(i2,i3)
          ystar(i)=dmax1(ystar(i),tiny)
c
          factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
          duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*ystar(i)))*dexp(-1./factor)
          uplus(i)= ystar(i)*dexp(-factor)+
     *               (dlog(ystar(i))/cappa+cc)*dexp(-1./factor)
          factor=dmin1(ystar(i)/(uplus(i)+tiny),
     *            (Viscosity(i1)+TurbulentViscosity(i1))/Viscosity(i1))
          WallVelocity2=WallVelocity*WallVelocity
          TauWallStar=Density(i1)*ustar(i)*ustar(i)
          Ec=WallVelocity2/(SpecificHeat(i1)*Temperature(i1)+tiny)
          Bstar=HeatFlux(i2,i3)*WallVelocity/
     *             (SpecificHeat(i1)*Temperature(i1)*TauWallStar+tiny)
          TwallOverT=dmax1(1.+sigT*Bstar+sigT*Ec/2.,0.)
          Dstar=dsqrt(Bstar*Bstar+2.*TwallOverT*Ec/sigT)
          Ecert=(dsqrt(2.*TwallOverT/sigT))*
     *                   (dasin((Bstar+Ec)/(Dstar+tiny))-
     *                                      dasin(Bstar/(Dstar+tiny)))
          Ucomp=Ecert*dsqrt(SpecificHeat(i1)*Temperature(i1))
          TauWall(i)=Density(i1)*ustar(i)*Ucomp*factor/(ystar(i)+tiny)
c            
          WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
          uTau(i)=Ucomp/(dlog(ystar(i)*dNorm*
     *                  BDensity(i2,i3)/BViscosity(i2,i3))/cappa+cc)
          yplus(i)=BDensity(i2,i3)*dNorm*uTau(i)/BViscosity(i2,i3)
c
          BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
          BeDiffCoefficient(i2,i3)=
     *               dmax1(WallViscosity(i),BViscosity(i2,i3))
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
      SUBROUTINE EquilibriumMomentumRoughWallFunctions
c
C#############################################################################################
c
      use User0, only: MomentumWallFunctionType,Lcompressible,
     *                 LReadOldSolution,LSolveTurbulenceDissipationRate
      use Multigrid2, only: nIter
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use Variables1, only: TurbulentKE,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      BModifiedED,Temperature
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,SpecificHeat,
     *                               BeDiffCoefficient
      use Constants1, only: tiny
      use BoundaryFluxes, only: HeatFlux
      use Turbulence1, only: yplus,WallViscosity,duplusdyplus,uTau,
     *                       uplus,cmu25,cappa,cc,ctrans,uStar,ystar,
     *                       Bfv1Coefficient,sigT,TauWall,Stelda,
     *                       RoughccM,KsPlus
      use WallDistance1, only: WallDistance
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,k
      double precision :: nx,ny,nz,dNorm,cmu25k,yplusmin,rstar,ustar1,
     *                    factor,yplus1,uWall1,vWall1,wWall1,dotproduct,
     *                    WallVelocity,uTauViscous,uTauLog,
     *                    WallVelocity2,TauWallStar,Ec,Bstar,TwallOverT,
     *                    Dstar,Ecert,Ucomp
c
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
      if(nIter.eq.1.and..not.LReadOldSolution) then
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
c
          uTau(i)=dsqrt((Viscosity(i1)+TurbulentViscosity(i1))*
     *                              WallVelocity/(dNorm*Density(i1)))
          yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
c
        enddo
c
      endif
c
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(MomentumWallFunctionType.eq.'standard') then
c
          yplusmin=0.2
c        
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            yplus1=dNorm*BDensity(i2,i3)*uTau(i)/BViscosity(i2,i3)
            yplus(i)=dmax1(yplus1,yplusmin)
c
            if(yplus(i).lt.dmin1(ctrans,KsPlus(i)/2.)) then
c
              uplus(i)=yplus(i)
              uTau(i)=WallVelocity/dmax1(uplus(i),tiny)
              duplusdyplus(i)=1.
c          
            else
c          
              uplus(i)=dlog(yplus(i))/cappa+cc-RoughccM(i)
              uTau(i)=WallVelocity/dmax1(uplus(i),tiny)
              duplusdyplus(i)=1./(cappa*yplus(i))
c
            endif
c
            TauWall(i)=BDensity(i2,i3)*uTau(i)*uTau(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *              dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
          ustar=uTau
          ystar=yplus
c
        elseif(MomentumWallFunctionType.eq.'scalable') then
c
          yplusmin=ctrans+tiny
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
c
            yplus1=dNorm*BDensity(i2,i3)*uTau(i)/BViscosity(i2,i3)
            yplus(i)=dmax1(yplus1,yplusmin,Ksplus(i)/2.)
            duplusdyplus(i)=1./(cappa*yplus(i))
            uplus(i)=dlog(yplus(i))/cappa+cc-RoughccM(i)
            uTau(i)=WallVelocity/dmax1(uplus(i),tiny)
c
            TauWall(i)=BDensity(i2,i3)*uTau(i)*uTau(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *             dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
          ustar=uTau
          ystar=yplus
c
        elseif(MomentumWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            uTauViscous=dsqrt(BViscosity(i2,i3)*WallVelocity/
     *                             (BDensity(i2,i3)*dNorm))
            uTauLog=WallVelocity/(dlog(yplus(i))/cappa+cc-RoughccM(i))
            uTau(i)=(uTauViscous**4+uTauLog**4)**0.25
c
            uplus(i)=WallVelocity/dmax1(uTau(i),tiny)
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
c        
            factor=0.01*(yplus(i)**4)/(1.+5.*yplus(i))+tiny
            duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*yplus(i)))*dexp(-1./factor)
c
            TauWall(i)=BDensity(i2,i3)*uTau(i)*uTau(i)
c
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny) 
c
            BTurbulentViscosity(i2,i3)=
     *                dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
          ustar=uTau
          ystar=yplus
c
        endif
c
      elseif(Lcompressible) then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
          WallVelocity=TangentialVelocity(i1)
c
          uTauViscous=BViscosity(i2,i3)*WallVelocity/
     *                             (BDensity(i2,i3)*dNorm)
          uTauLog=WallVelocity/(dlog(yplus(i))/cappa+cc-RoughccM(i))
c
          uTau(i)=(uTauViscous**4+uTauLog**4)**0.25
          yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
c
          factor=0.01*(yplus(i)**4)/(1.+5.*yplus(i))+tiny
          duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*yplus(i)))*dexp(-1./factor)
          uplus(i)= yplus(i)*dexp(-factor)+((1./cappa)*
     *            dlog(yplus(i)+tiny)+cc-RoughccM(i))*dexp(-1./factor)
          factor=dmin1(yplus(i)/(uplus(i)+tiny),
     *            (Viscosity(i1)+TurbulentViscosity(i1))/Viscosity(i1))
          WallVelocity2=WallVelocity*WallVelocity
          TauWallStar=Density(i1)*uTau(i)*uTau(i)
          Ec=WallVelocity2/(SpecificHeat(i1)*Temperature(i1)+tiny)
          Bstar=HeatFlux(i2,i3)*WallVelocity/
     *             (SpecificHeat(i1)*Temperature(i1)*TauWallStar+tiny)
          TwallOverT=dmax1(1.+sigT*Bstar+sigT*Ec/2.,0.)
          Dstar=dsqrt(Bstar*Bstar+2.*TwallOverT*Ec/sigT)
          Ecert=(dsqrt(2.*TwallOverT/sigT))*
     *                   (dasin((Bstar+Ec)/(Dstar+tiny))-
     *                                      dasin(Bstar/(Dstar+tiny)))
          Ucomp=Ecert*dsqrt(SpecificHeat(i1)*Temperature(i1))
          TauWall(i)=Density(i1)*uTau(i)*Ucomp*factor/(yplus(i)+tiny)
c            
          WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
          BTurbulentViscosity(i2,i3)=
     *             dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
          BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
        enddo
c
        ustar=uTau
        ystar=yplus
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE NonEquilibriumMomentumRoughWallFunctions
c
C#############################################################################################
c
      use User0, only: MomentumWallFunctionType,Lcompressible,
     *                 LSolveTurbulenceDissipationRate
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz
      use Variables1, only: TurbulentKE,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,Temperature
      use BoundaryFluxes, only: HeatFlux
      use PhysicalProperties1, only: Density,BDensity,Viscosity,
     *                               BViscosity,BTurbulentViscosity,
     *                               TurbulentViscosity,SpecificHeat,
     *                               BeDiffCoefficient
      use Constants1, only: tiny

      use Turbulence1, only: yplus,WallViscosity,duplusdyplus,ustar,
     *                       uplus,cmu25,cappa,cc,ctrans,sigT,uTau,
     *                       TauWall,ystar,a1sst,RoughccM,KsPlus
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,yplusmin,rstar,ustar1,factor,yplus1,
     *                    WallVelocity,uTauViscous,uTauLogStar,
     *                    WallVelocity2,TauWallStar,Ec,Bstar,TwallOverT,
     *                    Dstar,Ecert,Ucomp,dotproduct
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
      if(LSolveTurbulenceDissipationRate.and.
     *             MomentumWallFunctionType.eq.'automatic') then
c     
        print*,
     *  'automatic wall functions are not implemented for k-e models'
        pause
        stop
c
      endif
c
      call CalculateuStar
c
      if(.not.Lcompressible.or.LSolveTurbulenceDissipationRate) then
c
        if(MomentumWallFunctionType.eq.'standard') then
c
          yplusmin=0.2
c        
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            yplus1=BDensity(i2,i3)*ustar(i)*dNorm/BViscosity(i2,i3)
            ystar(i)=dmax1(yplusmin,yplus1)
c
            if(ystar(i).lt.dmin1(ctrans,KsPlus(i)/2.)) then
c
              uTau(i)=dsqrt(BViscosity(i2,i3)*
     *                  WallVelocity/(BDensity(i2,i3)*dNorm))
              duplusdyplus(i)=1.
              uplus(i)=ystar(i)
c          
            else
c          
              uplus(i)=dlog(ystar(i))/cappa+cc-RoughccM(i)
              uTau(i)=WallVelocity/dmax1(uplus(i),tiny)
              duplusdyplus(i)=1./(cappa*ystar(i))
c
            endif
c
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
            TauWall(i)=BDensity(i2,i3)*uTau(i)*ustar(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
        elseif(MomentumWallFunctionType.eq.'scalable') then
c
c-------------------------------------------------------------------
c     Smooth walls with scalable wall function - Menter wall functions
c-------------------------------------------------------------------
c
          yplusmin=ctrans+tiny
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1) 
            WallVelocity=TangentialVelocity(i1)
c
            yplus1=Density(i1)*ustar(i)*dNorm/(Viscosity(i1)+tiny)
            ystar(i)=dmax1(yplusmin,yplus1,Ksplus(i)/2.)
            duplusdyplus(i)=1./(cappa*ystar(i))
            uplus(i)=dlog(ystar(i))/cappa+cc-RoughccM(i)
            uTau(i)=WallVelocity/dmax1(uplus(i),tiny)
c
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
            TauWall(i)=BDensity(i2,i3)*uTau(i)*ustar(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *             dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
        elseif(MomentumWallFunctionType.eq.'automatic') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
c
            uTauViscous=dsqrt(BViscosity(i2,i3)*
     *                  WallVelocity/(BDensity(i2,i3)*dNorm))
            uTauLogStar=ustar(i)
c
            ustar(i)=(uTauViscous**4+uTauLogStar**4)**0.25
            ystar(i)=BDensity(i2,i3)*ustar(i)*dNorm/BViscosity(i2,i3)
            ystar(i)=dmax1(ystar(i),tiny)
            uTauLogStar=WallVelocity/
     *                       (dlog(ystar(i))/cappa+cc-RoughccM(i))
            uTau(i)=(uTauViscous**4+uTauLogStar**4)**0.25
c        
            factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
            duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*ystar(i)))*dexp(-1./factor)
            uplus(i)= ystar(i)*dexp(-factor)+
     *           (dlog(ystar(i))/cappa+cc-RoughccM(i))*dexp(-1./factor)
c
            yplus(i)=BDensity(i2,i3)*uTau(i)*dNorm/BViscosity(i2,i3)
            TauWall(i)=BDensity(i2,i3)*uTau(i)*ustar(i)
            WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
c
            BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
            BeDiffCoefficient(i2,i3)=
     *               dmax1(WallViscosity(i),BViscosity(i2,i3))
c
          enddo
c
        endif
c     
      elseif(Lcompressible) then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          dNorm=WallDistance(i1)
          WallVelocity=TangentialVelocity(i1)
c
          uTauViscous=dsqrt(BViscosity(i2,i3)*
     *                  WallVelocity/(BDensity(i2,i3)*dNorm))
          uTauLogStar=ustar(i)
c
          ustar(i)=(uTauViscous**4+uTauLogStar**4)**0.25
          ystar(i)=BDensity(i2,i3)*ustar(i)*dNorm/BViscosity(i2,i3)
          ystar(i)=dmax1(ystar(i),tiny)
c
          factor=0.01*(ystar(i)**4)/(1.+5.*ystar(i))+tiny
          duplusdyplus(i)= 1.*dexp(-factor)+
     *                   (1./(cappa*ystar(i)))*dexp(-1./factor)
          uplus(i)= ystar(i)*dexp(-factor)+
     *           (dlog(ystar(i))/cappa+cc-RoughccM(i))*dexp(-1./factor)
          factor=dmin1(ystar(i)/(uplus(i)+tiny),
     *            (Viscosity(i1)+TurbulentViscosity(i1))/Viscosity(i1))
          WallVelocity2=WallVelocity*WallVelocity
          TauWallStar=Density(i1)*ustar(i)*ustar(i)
          Ec=WallVelocity2/(SpecificHeat(i1)*Temperature(i1)+tiny)
          Bstar=HeatFlux(i2,i3)*WallVelocity/
     *             (SpecificHeat(i1)*Temperature(i1)*TauWallStar+tiny)
          TwallOverT=dmax1(1.+sigT*Bstar+sigT*Ec/2.,0.)
          Dstar=dsqrt(Bstar*Bstar+2.*TwallOverT*Ec/sigT)
          Ecert=(dsqrt(2.*TwallOverT/sigT))*
     *                   (dasin((Bstar+Ec)/(Dstar+tiny))-
     *                                      dasin(Bstar/(Dstar+tiny)))
          Ucomp=Ecert*dsqrt(SpecificHeat(i1)*Temperature(i1))
          TauWall(i)=Density(i1)*ustar(i)*Ucomp*factor/(ystar(i)+tiny)
c            
          WallViscosity(i)=TauWall(i)*dNorm/dmax1(WallVelocity,tiny)
          uTau(i)=Ucomp/(dlog(ystar(i)*dNorm*
     *          BDensity(i2,i3)/BViscosity(i2,i3))/cappa+cc-RoughccM(i))
          yplus(i)=BDensity(i2,i3)*dNorm*uTau(i)/BViscosity(i2,i3)
c
          BTurbulentViscosity(i2,i3)=
     *               dmax1(WallViscosity(i)-BViscosity(i2,i3),0.)
          BeDiffCoefficient(i2,i3)=
     *               dmax1(WallViscosity(i),BViscosity(i2,i3))
c
        enddo
c
      endif
c
      return
      end