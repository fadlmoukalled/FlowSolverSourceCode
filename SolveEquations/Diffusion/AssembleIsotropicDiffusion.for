c
C#############################################################################################
      SUBROUTINE AssembleDiffusionTerm(Variable,
     *       Gam,BGam,FiT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *                           dfidxfT,dfidyfT,dfidzfT)
C#############################################################################################
      use User0
      use Variables1
      use Variables2
      use Geometry1
      use Geometry3
      use Geometry4
      use PhysicalProperties1
      use BoundaryConditions1
      use BoundaryConditions2
      use BoundaryConditionsScalar1
      use BoundaryConditionsScalar2
      use BoundaryConditionsTurbulence1
      use BoundaryConditionsTurbulence2
      use Scalar1
      use Scalar2
      use BoundaryFluxes
      use Constants1
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k
      integer i1,i2,i3,i4,i5
      character*10 Variable
      double precision Gamf,dfidxTf,dfidyTf,dfidzTf,term1,term2,cpf
      double precision FluxCflocal,FluxFflocal,FluxVflocal,FluxTfLocal
      double precision dNorm,area,nx,ny,nz
      character*16 GradientInterpolationScheme
      integer k1,j1,j2,j3,j4
      double precision xF1,yF1,zF1,distance1,distance2,GFactCF,DistCFx,
     *                 DistCFy,DistCFz,DistCF,DistCFux,DistCFuy,
     *                 DistCFuz,DotProduct,FEx,FEy,FEz,FE,
     *                 fgDiff,FTx,FTy,FTz,FT,dfidxf,dfidyf,dfidzf,
     *                 dfidxf1,dfidyf1,dfidzf1,PhiB,thetaR
c
      double precision, dimension(:) :: FiT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
      double precision, dimension(:) :: dfidxfT
      double precision, dimension(:) :: dfidyfT
      double precision, dimension(:) :: dfidzfT
      double precision, dimension(:) :: Gam
      double precision, dimension(:,:) :: Bgam
c********************************************************************************************
c--- Interfaces
c********************************************************************************************
      interface
c********************************************************************************************
c--------------------------------------------------------------
        SUBROUTINE InterpolateGamaToFace(Variable,Gam,BGam)
c--------------------------------------------------------------
          character*10 Variable
          double precision, dimension(:) :: Gam
          double precision, dimension(:,:) :: Bgam
c--------------------------------------------------------------
        end SUBROUTINE InterpolateGamaToFace
c********************************************************************************************
      end interface
c********************************************************************************************
c--- Interior faces
c----------------------------------------------------------------------
c
      call InterpolateGamaToFace(Variable,Gam,BGam)

      if(Variable.eq.'temp'.and.EnergyEquation.eq.'htotal') then
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gamf=GamaFace(k)
c
          FluxCflocal= gamf*gDiff(k)
          FluxFflocal=-gamf*gDiff(k)
          FluxVflocal=-gamf*(dfidxfT(k)*FaceTx(k)+
     *                   dfidyfT(k)*FaceTy(k)+dfidzfT(k)*FaceTz(k))
c
          FluxVf(k)=FluxVf(k)+FluxVflocal
          FluxTf(k)= FluxTf(k)+
     *             FluxCflocal*FiT(i)+FluxFflocal*FiT(j)+FluxVflocal
c
          cpf=GFactorCF(k)*SpecificHeat(i)+
     *                           (1.-GFactorCF(k))*SpecificHeat(j)
          FluxCf(k)=FluxCf(k)+FluxCflocal/cpf
          FluxFf(k)=FluxFf(k)+FluxFflocal/cpf
c
        enddo
c
      else
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gamf=GamaFace(k)
c
          FluxCflocal= gamf*gDiff(k)
          FluxFflocal=-gamf*gDiff(k)
          FluxVflocal=-gamf*(dfidxfT(k)*FaceTx(k)+
     *                   dfidyfT(k)*FaceTy(k)+dfidzfT(k)*FaceTz(k))
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
          FluxVf(k)=FluxVf(k)+FluxVflocal
          FluxTf(k)= FluxTf(k)+
     *             FluxCflocal*FiT(i)+FluxFflocal*FiT(j)+FluxVflocal
c
        enddo
c
      endif
c----------------------------------------------------------------------
c--- Boundary faces
c----------------------------------------------------------------------
      if(variable.eq.'velx'.or.variable.eq.'vely'.or.
     *                                   variable.eq.'velz') then
c----------------------------------------------------------------------
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dNorm=BDistanceCFx(i2,i3)*BFaceAreanx(i2,i3)+
     *                BDistanceCFy(i2,i3)*BFaceAreany(i2,i3)+
     *                BDistanceCFz(i2,i3)*BFaceAreanz(i2,i3)
          area=BFaceArea(i2,i3)
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
c
          if(Variable.eq.'velx') then
c
            FluxCflocal= gamf*area*(1.-nx**2)/dNorm
            FluxFflocal=0.
            FluxVflocal=-(gamf*area/dNorm)*
     *        (BFiT(i2,i3)*(1.-nx**2)+
     *           (vVelocity(i1)-BvVelocity(i2,i3))*nx*ny+
     *           (wVelocity(i1)-BwVelocity(i2,i3))*nx*nz)
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+
     *                         FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          elseif(Variable.eq.'vely') then
c
            FluxCflocal= gamf*area*(1.-ny**2)/dNorm
            FluxFflocal=0.
            FluxVflocal=-(gamf*area/dNorm)*
     *        ((uVelocity(i1)-BuVelocity(i2,i3))*nx*ny+
     *             BFiT(i2,i3)*(1.-ny**2)+
     *               (wVelocity(i1)-BwVelocity(i2,i3))*nz*ny)
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+
     *                         FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          elseif(Variable.eq.'velz') then
c
            FluxCflocal= gamf*area*(1.-nz**2)/dNorm
            FluxFflocal=0.
            FluxVflocal=-(gamf*area/dNorm)*
     *        ((uVelocity(i1)-BuVelocity(i2,i3))*nx*nz+
     *               (vVelocity(i1)-BvVelocity(i2,i3))*ny*nz+
     *         BFiT(i2,i3)*(1.-nz**2))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+
     *                         FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          endif
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iinletsupersonic
c
          i1=IinletsupersonicOwner(i)
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *             dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i1=IinletSpecifiedVelocityOwner(i)
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i1=IinletSpecifiedStaticPressureOwner(i)
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationPressure
c
          i1=IinletSpecifiedStagnationPressureOwner(i)
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0 
          FluxFf(i4)=FluxFf(i4)+0.D0 
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i1=IoutletspecifiedVelocityOwner(i)
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *            dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedStaticPressure
c
          i1=IoutletSpecifiedStaticPressureOwner(i)
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0 
          FluxFf(i4)=FluxFf(i4)+0.D0 
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i1=IoutletSpecifiedAverageStaticPressureOwner(i)
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0 
          FluxFf(i4)=FluxFf(i4)+0.D0 
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedResistance
c
          i1=IoutletSpecifiedResistanceOwner(i)
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0 
          FluxFf(i4)=FluxFf(i4)+0.D0 
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletFullyDeveloped
c
          i1=IoutletFullyDevelopedOwner(i)
          i2=IoutletFullyDevelopedNumberOfBCSets(i)
          i3=IoutletFullyDevelopedNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
c
          dNorm=BDistanceCFx(i2,i3)*BFaceAreanx(i2,i3)+
     *                BDistanceCFy(i2,i3)*BFaceAreany(i2,i3)+
     *                BDistanceCFz(i2,i3)*BFaceAreanz(i2,i3)
          area=BFaceArea(i2,i3)
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
c
          if(Variable.eq.'velx') then
c
            FluxCflocal= 2.*gamf*area*(nx**2)/dNorm
            FluxFflocal=0.
            FluxVflocal=-(-2.*gamf*area/dNorm)*
     *               (vVelocity(i1)*ny+wVelocity(i1)*nz)*nx
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+
     *                         FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          elseif(Variable.eq.'vely') then
c
            FluxCflocal= 2.*gamf*area*(ny**2)/dNorm
            FluxFflocal=0.
            FluxVflocal=-(-2.*gamf*area/dNorm)*
     *               (uVelocity(i1)*nx+wVelocity(i1)*nz)*ny
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+
     *                         FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          elseif(Variable.eq.'velz') then
c
            FluxCflocal= 2.*gamf*area*(nz**2)/dNorm
            FluxFflocal=0.
            FluxVflocal=-(-2.*gamf*area/dNorm)*
     *               (uVelocity(i1)*nx+vVelocity(i1)*ny)*nz
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+
     *                         FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          endif
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          if(.not.LPeriodicImplicit) then
c
            if(LRotationalPeriodicity) then
            
              if(Variable.eq.'velx') then
c
                PhiB=a1r(j2)*uVelocity(j1)+
     *                  b1r(j2)*vVelocity(j1)+c1r(j2)*wVelocity(j1)
c
              elseif(Variable.eq.'vely') then
c
                PhiB=a2r(j2)*uVelocity(j1)+
     *                  b2r(j2)*vVelocity(j1)+c2r(j2)*wVelocity(j1)
c
              elseif(Variable.eq.'velz') then
c
                PhiB=a3r(j2)*uVelocity(j1)+
     *                  b3r(j2)*vVelocity(j1)+c3r(j2)*wVelocity(j1)
c
              endif
c
            elseif(LTranslationalPeriodicity) then
c 
              PhiB=FiT(j1)
c
            endif
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)= FluxTf(i4)+
     *             FluxCflocal*FiT(i1)+FluxFflocal*PhiB+FluxVflocal
c
          elseif(LPeriodicImplicit) then
c
            if(LRotationalPeriodicity) then
c
              thetaR=theta(j2)*pi/180.
c
              if(Variable.eq.'velx') then
c
                FluxFflocal=FluxFflocal*
     *               (1.+(a2Axis(j2)**2+a3Axis(j2)**2)*
     *                             dmax1(dcos(thetaR),0.))
                FluxVflocal=FluxVflocal-FluxCflocal*
     *           ((a2Axis(j2)**2+a3Axis(j2)**2)*
     *                 (dmax1(-dcos(thetaR),0.)-1.)*uVelocity(j1)+
     *           (a1Axis(j2)*a2Axis(j2)*(1.-dcos(thetaR))-
     *                          a3Axis(j2)*dsin(thetaR))*vVelocity(j1)+
     *           (a1Axis(j2)*a3Axis(j2)*(1.-dcos(thetaR))+
     *                          a2Axis(j2)*dsin(thetaR))*wVelocity(j1))
c              
                FluxCf(i4)=FluxCf(i4)+FluxCflocal
                FluxFf(i4)=FluxFf(i4)+FluxFflocal
                FluxVf(i4)=FluxVf(i4)+FluxVflocal
                FluxTf(i4)= FluxTf(i4)+
     *             FluxCflocal*FiT(i1)+FluxFflocal*FiT(j1)+FluxVflocal
c
              elseif(Variable.eq.'vely') then
c
                FluxFflocal=FluxFflocal*
     *               (1.+(a1Axis(j2)**2+a3Axis(j2)**2)*
     *                             dmax1(dcos(thetaR),0.))
                FluxVflocal=FluxVflocal-FluxCflocal*
     *           ((a1Axis(j2)**2+a3Axis(j2)**2)*
     *                 (dmax1(-dcos(thetaR),0.)-1.)*vVelocity(j1)+
     *           (a1Axis(j2)*a2Axis(j2)*(1.-dcos(thetaR))+
     *                          a3Axis(j2)*dsin(thetaR))*uVelocity(j1)+
     *           (a2Axis(j2)*a3Axis(j2)*(1.-dcos(thetaR))-
     *                          a1Axis(j2)*dsin(thetaR))*wVelocity(j1))
c              
                FluxCf(i4)=FluxCf(i4)+FluxCflocal
                FluxFf(i4)=FluxFf(i4)+FluxFflocal
                FluxVf(i4)=FluxVf(i4)+FluxVflocal
                FluxTf(i4)= FluxTf(i4)+
     *             FluxCflocal*FiT(i1)+FluxFflocal*FiT(j1)+FluxVflocal
c
              elseif(Variable.eq.'velz') then
c
                FluxFflocal=FluxFflocal*
     *               (1.+(a1Axis(j2)**2+a2Axis(j2)**2)*
     *                             dmax1(dcos(thetaR),0.))
                FluxVflocal=FluxVflocal-FluxCflocal*
     *           ((a1Axis(j2)**2+a2Axis(j2)**2)*
     *                 (dmax1(-dcos(thetaR),0.)-1.)*wVelocity(j1)+
     *           (a1Axis(j2)*a3Axis(j2)*(1.-dcos(thetaR))-
     *                          a2Axis(j2)*dsin(thetaR))*uVelocity(j1)+
     *           (a2Axis(j2)*a3Axis(j2)*(1.-dcos(thetaR))+
     *                          a1Axis(j2)*dsin(thetaR))*vVelocity(j1))
c              
                FluxCf(i4)=FluxCf(i4)+FluxCflocal
                FluxFf(i4)=FluxFf(i4)+FluxFflocal
                FluxVf(i4)=FluxVf(i4)+FluxVflocal
                FluxTf(i4)= FluxTf(i4)+
     *             FluxCflocal*FiT(i1)+FluxFflocal*FiT(j1)+FluxVflocal
c
              endif
c
            elseif(LTranslationalPeriodicity) then
c
              FluxCf(i4)=FluxCf(i4)+FluxCflocal
              FluxFf(i4)=FluxFf(i4)+FluxFflocal
              FluxVf(i4)=FluxVf(i4)+FluxVflocal
              FluxTf(i4)= FluxTf(i4)+
     *             FluxCflocal*FiT(i1)+FluxFflocal*FiT(j1)+FluxVflocal
c
            endif
c
          endif
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
      elseif(variable.eq.'frelax') then
c----------------------------------------------------------------------
c
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
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
          BFiT(i2,i3)=0.
c
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
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
          BFiT(i2,i3)=FiT(i1)
c
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
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
          BFiT(i2,i3)=FiT(i1)
c
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c
        do i=1,IinletTurbulence
c
          i1=IinletTurbulenceOwner(i)
          i2=IinletTurbulenceNumberOfBCSets(i)
          i3=IinletTurbulenceNBFaces(i)
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
          BFiT(i2,i3)=FiT(i1)
c
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0 
          FluxFf(i4)=FluxFf(i4)+0.d0 
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
      elseif(variable.eq.'tv2') then
c----------------------------------------------------------------------
c
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
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
          BFiT(i2,i3)=0.
c
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
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
          BFiT(i2,i3)=FiT(i1)
c
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
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
          if(Bmdot(i2,i3).lt.0.) then
c
            gamf=BGamaFace(i2,i3)
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          else  
c
            FluxCf(i4)=FluxCf(i4)+0.d0
            FluxFf(i4)=FluxFf(i4)+0.d0
            FluxVf(i4)=FluxVf(i4)+0.d0
            FluxTf(i4)=FluxTf(i4)+0.d0
c
          endif
c
        enddo
c
        do i=1,IinletTurbulence
c
          i1=IinletTurbulenceOwner(i)
          i2=IinletTurbulenceNumberOfBCSets(i)
          i3=IinletTurbulenceNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0 
          FluxFf(i4)=FluxFf(i4)+0.d0 
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
      elseif(variable.eq.'tzeta') then
c----------------------------------------------------------------------
c
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
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
          BFiT(i2,i3)=0.
c
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
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
          BFiT(i2,i3)=FiT(i1)
c
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
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
          if(Bmdot(i2,i3).lt.0.) then
c
            gamf=BGamaFace(i2,i3)
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          else  
c
            FluxCf(i4)=FluxCf(i4)+0.d0
            FluxFf(i4)=FluxFf(i4)+0.d0
            FluxVf(i4)=FluxVf(i4)+0.d0
            FluxTf(i4)=FluxTf(i4)+0.d0
c
          endif
c
        enddo
c
        do i=1,IinletTurbulence
c
          i1=IinletTurbulenceOwner(i)
          i2=IinletTurbulenceNumberOfBCSets(i)
          i3=IinletTurbulenceNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0 
          FluxFf(i4)=FluxFf(i4)+0.d0 
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
      elseif(variable.eq.'tgamma'.or.variable.eq.'tretheta') then
c----------------------------------------------------------------------
c
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
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
          BFiT(i2,i3)=FiT(i1)
c
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
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
          BFiT(i2,i3)=FiT(i1)
c
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
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
          if(Bmdot(i2,i3).lt.0.) then
c
            gamf=BGamaFace(i2,i3)
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          endif
c
        enddo
c
        do i=1,IinletTurbulence
c
          i1=IinletTurbulenceOwner(i)
          i2=IinletTurbulenceNumberOfBCSets(i)
          i3=IinletTurbulenceNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0 
          FluxFf(i4)=FluxFf(i4)+0.d0 
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
c----------------------------------------------------------------------
      elseif(variable.eq.'tke'.or.variable.eq.'ted'.or.
     *         variable.eq.'tomega'.or.variable.eq.'med'.or.
     *                                variable.eq.'tkl') then
c----------------------------------------------------------------------
c
        if(variable.eq.'med') then
c
          do i=1,IWallSlip
c
            i1=IWallSlipOwner(i)
            i2=IWallSlipNumberOfBCSets(i)
            i3=IWallSlipNBFaces(i)
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
            BFiT(i2,i3)=0.
c
            gamf=BGamaFace(i2,i3)
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
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
            BFiT(i2,i3)=FiT(i1)
c
            FluxCf(i4)=FluxCf(i4)+0.d0
            FluxFf(i4)=FluxFf(i4)+0.d0
            FluxVf(i4)=FluxVf(i4)+0.d0
            FluxTf(i4)=FluxTf(i4)+0.d0
c
          enddo
c
        endif
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
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
          if(Bmdot(i2,i3).lt.0.) then
c
            gamf=BGamaFace(i2,i3)
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+
     *        FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          endif
c
        enddo
c
        do i=1,IinletTurbulence
c
          i1=IinletTurbulenceOwner(i)
          i2=IinletTurbulenceNumberOfBCSets(i)
          i3=IinletTurbulenceNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0 
          FluxFf(i4)=FluxFf(i4)+0.d0 
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
c----------------------------------------------------------------------
      elseif(variable.eq.'temp'.and.
     *                EnergyEquation.eq.'temperature') then
c----------------------------------------------------------------------
        do i=1,IWallDirichlet
c
          i1=IWallDirichletOwner(i)
          i2=IWallDirichletNumberOfBCSets(i)
          i3=IWallDirichletNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallVonNeumann
c
          i1=IWallVonNeumannOwner(i)
          i2=IWallVonNeumannNumberOfBCSets(i)
          i3=IWallVonNeumannNBFaces(i)
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
          FluxVflocal=-HeatFlux(i2,i3)*BFaceArea(i2,i3)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallRobin
c
          i1=IWallRobinOwner(i)
          i2=IWallRobinNumberOfBCSets(i)
          i3=IWallRobinNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          term1=HinfinityRobin(i)*BFaceArea(i2,i3)+gamf*BgDiff(i2,i3)
          FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                                  BgDiff(i2,i3)*gamf/term1
          FluxFflocal=0.
          FluxVflocal=-TinfinityRobin(i)*FluxCflocal-
     *         HinfinityRobin(i)*BFaceArea(i2,i3)*gamf*
     *          (dfidxTf*BFaceTx(i2,i3)+dfidyTf*BFaceTy(i2,i3)
     *                            +dfidzTf*BFaceTz(i2,i3))/term1
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSupersonic
c
          i1=IinletSupersonicOwner(i)
          i2=IinletSupersonicNumberOfBCSets(i)
          i3=IinletSupersonicNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticTemperature
c
          i1=IinletSpecifiedStaticTemperatureOwner(i)
          i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticTemperatureNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationTemperature
c
          i1=IinletSpecifiedStagnationTemperatureOwner(i)
          i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletFullyDevelopedEnergy
c
          i1=IoutletFullyDevelopedEnergyOwner(i)
          i2=IoutletFullyDevelopedEnergyNumberOfBCSets(i)
          i3=IoutletFullyDevelopedEnergyNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSupersonic
c
          i1=IoutletSupersonicOwner(i)
          i2=IoutletSupersonicNumberOfBCSets(i)
          i3=IoutletSupersonicNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
c----------------------------------------------------------------------
      elseif(variable.eq.'temp'.and.EnergyEquation.eq.'htotal') then
c----------------------------------------------------------------------
        do i=1,IWallDirichlet
c
          i1=IWallDirichletOwner(i)
          i2=IWallDirichletNumberOfBCSets(i)
          i3=IWallDirichletNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal/BSpecificHeat(i2,i3)
          FluxFf(i4)=FluxFf(i4)+FluxFflocal/BSpecificHeat(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallVonNeumann
c
          i1=IWallVonNeumannOwner(i)
          i2=IWallVonNeumannNumberOfBCSets(i)
          i3=IWallVonNeumannNBFaces(i)
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
          FluxVflocal=-HeatFlux(i2,i3)*BFaceArea(i2,i3)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallRobin
c
          i1=IWallRobinOwner(i)
          i2=IWallRobinNumberOfBCSets(i)
          i3=IWallRobinNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          term1=HinfinityRobin(i)*BFaceArea(i2,i3)+gamf*BgDiff(i2,i3)
          FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                                  BgDiff(i2,i3)*gamf/term1
          FluxFflocal=0.
          FluxVflocal=-TinfinityRobin(i)*FluxCflocal-
     *         HinfinityRobin(i)*BFaceArea(i2,i3)*gamf*
     *          (dfidxTf*BFaceTx(i2,i3)+dfidyTf*BFaceTy(i2,i3)
     *                            +dfidzTf*BFaceTz(i2,i3))/term1
c
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+FluxVflocal
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal/BSpecificHeat(i2,i3)
          FluxFf(i4)=FluxFf(i4)+FluxFflocal/BSpecificHeat(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSupersonic
c
          i1=IinletSupersonicOwner(i)
          i2=IinletSupersonicNumberOfBCSets(i)
          i3=IinletSupersonicNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal/BSpecificHeat(i2,i3)
          FluxFf(i4)=FluxFf(i4)+FluxFflocal/BSpecificHeat(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticTemperature
c
          i1=IinletSpecifiedStaticTemperatureOwner(i)
          i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticTemperatureNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal/BSpecificHeat(i2,i3)
          FluxFf(i4)=FluxFf(i4)+FluxFflocal/BSpecificHeat(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationTemperature
c
          i1=IinletSpecifiedStagnationTemperatureOwner(i)
          i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                  dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal/BSpecificHeat(i2,i3)
          FluxFf(i4)=FluxFf(i4)+FluxFflocal/BSpecificHeat(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletFullyDevelopedEnergy
c
          i1=IoutletFullyDevelopedEnergyOwner(i)
          i2=IoutletFullyDevelopedEnergyNumberOfBCSets(i)
          i3=IoutletFullyDevelopedEnergyNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSupersonic
c
          i1=IoutletSupersonicOwner(i)
          i2=IoutletSupersonicNumberOfBCSets(i)
          i3=IoutletSupersonicNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
          BFiT(i2,i3)=FiT(i1)
          BHtotal(i2,i3)=Htotal(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
          cpf=GFactCF*SpecificHeat(i1)+(1.-GFactCF)*SpecificHeat(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal/cpf
          FluxFf(i4)=FluxFf(i4)+FluxFflocal/cpf
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
          BFiT(i2,i3)=FiT(i1)
          BHtotal(i2,i3)=Htotal(i1)
c
        enddo
c----------------------------------------------------------------------
      else
c----------------------------------------------------------------------
c
        i5=iScalarVariable
c
        do i=1,IWallDirichletScalar(i5)
c
          i1=IWallDirichletScalarOwner(i,i5)
          i2=IWallDirichletScalarNumberOfBCSets(i,i5)
          i3=IWallDirichletScalarNBFaces(i,i5)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallVonNeumannScalar(i5)
c
          i1=IWallVonNeumannScalarOwner(i,i5)
          i2=IWallVonNeumannScalarNumberOfBCSets(i,i5)
          i3=IWallVonNeumannScalarNBFaces(i,i5)
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
          FluxVflocal=-ScalarFlux(i2,i3,i5)*BFaceArea(i2,i3)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallRobinScalar(i5)
c
          i1=IWallRobinScalarOwner(i,i5)
          i2=IWallRobinScalarNumberOfBCSets(i,i5)
          i3=IWallRobinScalarNBFaces(i,i5)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          term1=ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)+
     *                                              gamf*BgDiff(i2,i3)
          FluxCflocal=ConvectionCoefficientRobin(i,i5)*
     *                      BFaceArea(i2,i3)*BgDiff(i2,i3)*gamf/term1
          FluxFflocal=0.
          FluxVflocal=-PhiinfinityRobin(i,i5)*FluxCflocal-
     *       ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)*gamf*
     *            (dfidxTf*BFaceTx(i2,i3)+dfidyTf*BFaceTy(i2,i3)+
     *                                   dfidzTf*BFaceTz(i2,i3))/term1

c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxCflocal*FiT(i1)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSupersonicScalar(i5)
c
          i1=IinletSupersonicScalarOwner(i,i5)
          i2=IinletSupersonicScalarNumberOfBCSets(i,i5)
          i3=IinletSupersonicScalarNBFaces(i,i5)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedValueScalar(i5)
c
          i1=IinletSpecifiedValueScalarOwner(i,i5)
          i2=IinletSpecifiedValueScalarNumberOfBCSets(i,i5)
          i3=IinletSpecifiedValueScalarNBFaces(i,i5)
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
          gamf=BGamaFace(i2,i3)
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*BgDiff(i2,i3)
          FluxFflocal=-gamf*BgDiff(i2,i3)
          FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                 dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+
     *      FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletFullyDevelopedScalar(i5)
c
          i1=IoutletFullyDevelopedScalarOwner(i,i5)
          i2=IoutletFullyDevelopedScalarNumberOfBCSets(i,i5)
          i3=IoutletFullyDevelopedScalarNBFaces(i,i5)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSupersonicScalar(i5)
c
          i1=IoutletSupersonicScalarOwner(i,i5)
          i2=IoutletSupersonicScalarNumberOfBCSets(i,i5)
          i3=IoutletSupersonicScalarNBFaces(i,i5)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do k1=1,j2-1
c
            j4=j4+NBFaces(k1)
c
          enddo
c
          j4=j4+j3
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
!          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
!     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
!     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
!          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
!     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
!     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
!c
!          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          if(MethodDecomposeS.eq.1) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeS.eq.3) then
c
            DotProduct=DistCFux*BFaceAreax(i2,i3)+
     *                     DistCFuy*BFaceAreay(i2,i3)+
     *                        DistCFuz*BFaceAreaz(i2,i3)
            FEx=(BFaceArea(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceArea(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreax(i2,i3)-FEx
            FTy=BFaceAreay(i2,i3)-FEy
            FTz=BFaceAreaz(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          gamf=BGamaFace(i2,i3) !GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= gamf*fgDiff
          FluxFflocal=-gamf*fgDiff
          FluxVflocal=-gamf*(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.D0
          FluxFf(i4)=FluxFf(i4)+0.D0
          FluxVf(i4)=FluxVf(i4)+0.D0
          FluxTf(i4)=FluxTf(i4)+0.D0
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c----------------------------------------------------------------------
c
      endif
c
      return
      end