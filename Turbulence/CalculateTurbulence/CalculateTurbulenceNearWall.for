c
c#############################################################################################
c
      SUBROUTINE CalculateWallTurbulence
c
C#############################################################################################
c
      use User0, only: LSolveTurbulenceKineticEnergy,LRough
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      if(.not.LSolveTurbulenceKineticEnergy) call ExtractTurbulentKE
c
      call CalculateuTau
      call Calculateustar
      call Calculateyplus
      call Calculateystar
      if(LRough) call CalculateKsPlus
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE Calculateystar
c
C#############################################################################################
c
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use Turbulence1, only: cmu25,ustar,ystar
      use Variables1, only: TurbulentKE
      use WallDistance1, only: WallDistance
      use PhysicalProperties1, only: BDensity,BViscosity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
c********************************************************************************************
c
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        ystar(i)=WallDistance(i1)*ustar(i)*
     *               BDensity(i2,i3)/BViscosity(i2,i3)
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE Calculateyplus
c
C#############################################################################################
c
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use Turbulence1, only: uTau,yplus
      use WallDistance1, only: WallDistance
      use PhysicalProperties1, only: BDensity,BViscosity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
c********************************************************************************************
c
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        yplus(i)=WallDistance(i1)*uTau(i)*
     *               BDensity(i2,i3)/BViscosity(i2,i3)
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateuStar
c
C#############################################################################################
c
      use user0, only: LSolveTurbulenceSpecificDissipationRate
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner
      use Turbulence1, only: a1sst,ustar,cmu25
      use Variables1, only: TurbulentKE
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1
c********************************************************************************************
c
      if(LSolveTurbulenceSpecificDissipationRate) then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
c
          ustar(i)=dsqrt(dmax1(a1sst*TurbulentKE(i1),0.))
c
        enddo
c
      else
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
c
          ustar(i)=cmu25*dsqrt(dmax1(TurbulentKE(i1),0.))
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
      SUBROUTINE CalculateuTau
c
C#############################################################################################
c
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use WallDistance1, only: WallDistance
      use PhysicalProperties1, only: BDensity,BViscosity,
     *                               BTurbulentViscosity
      use Turbulence1, only: uTau
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,WallVelocity
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
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        dNorm=WallDistance(i1)
        WallVelocity=TangentialVelocity(i1)
        uTau(i)=dsqrt((BViscosity(i2,i3)+BTurbulentViscosity(i2,i3))*
     *                             WallVelocity/(dNorm*BDensity(i2,i3)))
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateKsPlus
c
C#############################################################################################
c
      use User0, only: GrainSize
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use PhysicalProperties1, only: BDensity,BViscosity
      use Turbulence1, only: uTau,KsPlus
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: dNorm,WallVelocity
c********************************************************************************************
c
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        KsPlus(i)=BDensity(i2,i3)*GrainSize(i2)*
     *                             uTau(i)/BViscosity(i2,i3)
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateRoughccMomentum
c
C#############################################################################################
c
      use BoundaryConditionsTurbulence2, only: IwallTurbulence
      use Turbulence1, only: KsPlus,RoughccM,cappa
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
c********************************************************************************************
c
      do i=1,IwallTurbulence
c
        RoughccM(i)=(dlog(1.+0.3*KsPlus(i)))/cappa
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateRoughccEnergy
c
C#############################################################################################
c
      use User0, only: CRoughEnergy
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use PhysicalProperties1, only: BViscosity,BSpecificHeat,
     *                               BConductivity
      use Turbulence1, only: KsPlus,RoughccE,cappa
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: Pr
c********************************************************************************************
c
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        Pr=BViscosity(i2,i3)*BSpecificHeat(i2,i3)/BConductivity(i2,i3)
        RoughccE(i)=(dlog(1.+CRoughEnergy*0.3*Pr*KsPlus(i)))/cappa
c
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateRoughccScalar
c
C#############################################################################################
c
      use User0, only: CRoughScalar
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use PhysicalProperties1, only: BViscosity,BSpecificHeatScalar,
     *                               BDiffusionCoefficient
      use Turbulence1, only: KsPlus,RoughccS,cappa
      use Scalar2, only: iScalarVariable
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,iScalar
      double precision :: Sch
c********************************************************************************************
c
      iScalar=iScalarVariable
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        Sch=BViscosity(i2,i3)*BSpecificHeatScalar(i2,i3,iScalar)/
     *                             BDiffusionCoefficient(i2,i3,iScalar)
        RoughccS(i,iScalar)=
     *         (dlog(1.+CRoughScalar(iScalar)*0.3*Sch*KsPlus(i)))/cappa
c
      enddo
c
      return
      end