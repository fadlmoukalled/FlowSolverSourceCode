c
C#############################################################################################
      SUBROUTINE AssembleAnisotropicDiffusionTerm
     *       (Variable,FiT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *                            dfidxfT,dfidyfT,dfidzfT)
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
      use Scalar1
      use Scalar2
      use BoundaryFluxes
      use Constants1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k
      integer i1,i2,i3,i4,i5
      character*10 Variable
      double precision Gamf,dfidxTf,dfidyTf,dfidzTf,term1,term2,cpf
      double precision FluxCflocal,FluxFflocal,FluxVflocal,FluxTfLocal
      double precision dNorm,area,nx,ny,nz
      character*16 GradientInterpolationScheme
      integer k1,j1,j2,j3,j4
      double precision xF1,yF1,zF1,distance1,distance2,distance3,
     *                 GFactCF,DistCFx,DistCFy,DistCFz,DistCF,
     *                 DistCFux,DistCFuy,DistCFuz,DotProduct,
     *                 FEx,FEy,FEz,FE,fgDiff,FTx,FTy,FTz,FT,
     *                 dfidxf,dfidyf,dfidzf,
     *                 dfidxf1,dfidyf1,dfidzf1,PhiB,thetaR,
     *                 Tfx,Tfy,Tfz,Tf,denom,numer,ratio,
     *                 DotProduct1,DotProduct2
c
      double precision, dimension(:) :: FiT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
      double precision, dimension(:) :: dfidxfT
      double precision, dimension(:) :: dfidyfT
      double precision, dimension(:) :: dfidzfT
c*********************************************************************************************
c--- Interior faces
c----------------------------------------------------------------------
c
      call calculateSprime(Variable)
c
      if(Variable.eq.'temp'.and.EnergyEquation.eq.'htotal') then
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          FluxCflocal= gDiffp(k)
          FluxFflocal=-gDiffp(k)
          FluxVflocal=-(dfidxfT(k)*FaceTxp(k)+
     *                 dfidyfT(k)*FaceTyp(k)+dfidzfT(k)*FaceTzp(k))
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
      elseif(Variable.eq.'lambda') then
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          FluxCflocal= gDiffp(k)
          FluxFflocal=-gDiffp(k)
          FluxVflocal=-(dfidxfT(k)*FaceTxp(k)+
     *                 dfidyfT(k)*FaceTyp(k)+dfidzfT(k)*FaceTzp(k))
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
          FluxVf(k)=FluxVf(k)+FluxVflocal
          FluxTf(k)= FluxTf(k)+
     *             FluxCflocal*FiT(i)+FluxFflocal*FiT(j)+FluxVflocal
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
          FluxCflocal= gDiffp(k)
          FluxFflocal=-gDiffp(k)
          FluxVflocal=-(dfidxfT(k)*FaceTxp(k)+
     *                 dfidyfT(k)*FaceTyp(k)+dfidzfT(k)*FaceTzp(k))
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
c
c----------------------------------------------------------------------
c--- Boundary faces
c----------------------------------------------------------------------
c
      if(variable.eq.'temp'.and.EnergyEquation.eq.'temperature') then
c
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *            dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          term1=HinfinityRobin(i)*BFaceArea(i2,i3)+BgDiffp(i2,i3)
          FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                           BgDiffp(i2,i3)/term1
          FluxFflocal=0.
          FluxVflocal=-TinfinityRobin(i)*FluxCflocal-
     *         HinfinityRobin(i)*BFaceArea(i2,i3)*
     *          (dfidxTf*BFaceTxp(i2,i3)+dfidyTf*BFaceTyp(i2,i3)+
     *                  +dfidzTf*BFaceTzp(i2,i3))/term1

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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *                 dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *                  dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *               dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
c          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
c     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
c     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
c          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
c     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
c     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
c          GFactCF=distance2/(distance1+distance2)
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
          if(MethodDecomposeSprime.eq.1) then
c
            DotProduct=DistCFux*BFaceAreaxp(i2,i3)+
     *                     DistCFuy*BFaceAreayp(i2,i3)+
     *                        DistCFuz*BFaceAreazp(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.2) then
c
            FEx=BFaceAreap(i2,i3)*DistCFux
            FEy=BFaceAreap(i2,i3)*DistCFuy
            FEz=BFaceAreap(i2,i3)*DistCFuz
            FE=BFaceAreap(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.3) then
c
            DotProduct=DistCFux*BFaceAreaxp(i2,i3)+
     *                     DistCFuy*BFaceAreayp(i2,i3)+
     *                        DistCFuz*BFaceAreazp(i2,i3)
            FEx=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.4) then
c
            dotproduct1=BFaceAreaxp(i2,i3)*BFaceAreanx(i2,i3)+
     *                     BFaceAreayp(i2,i3)*BFaceAreany(i2,i3)+
     *                           BFaceAreazp(i2,i3)*BFaceAreanz(i2,i3)
            dotproduct2=DistCFux*BFaceAreanx(i2,i3)+
     *                     DistCFuy*BFaceAreany(i2,i3)+
     *                             DistCFuz*BFaceAreanz(i2,i3)
c
            ratio=dotproduct1/dotproduct2
c
            FEx=ratio*DistCFux
            FEy=ratio*DistCFuy
            FEz=ratio*DistCFuz
            FE=dabs(ratio)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= fgDiff
          FluxFflocal=-fgDiff
          FluxVflocal=-(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
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
      elseif(variable.eq.'temp'.and.EnergyEquation.eq.'htotal') then
c
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *            dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          term1=HinfinityRobin(i)*BFaceArea(i2,i3)+BgDiffp(i2,i3)
          FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                           BgDiffp(i2,i3)/term1
          FluxFflocal=0.
          FluxVflocal=-TinfinityRobin(i)*FluxCflocal-
     *         HinfinityRobin(i)*BFaceArea(i2,i3)*
     *          (dfidxTf*BFaceTxp(i2,i3)+dfidyTf*BFaceTyp(i2,i3)+
     *                  +dfidzTf*BFaceTzp(i2,i3))/term1

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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *                 dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *                  dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *               dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          if(MethodDecomposeSprime.eq.1) then
c
            DotProduct=DistCFux*BFaceAreaxp(i2,i3)+
     *                     DistCFuy*BFaceAreayp(i2,i3)+
     *                        DistCFuz*BFaceAreazp(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.2) then
c
            FEx=BFaceAreap(i2,i3)*DistCFux
            FEy=BFaceAreap(i2,i3)*DistCFuy
            FEz=BFaceAreap(i2,i3)*DistCFuz
            FE=BFaceAreap(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.3) then
c
            DotProduct=DistCFux*BFaceAreaxp(i2,i3)+
     *                     DistCFuy*BFaceAreayp(i2,i3)+
     *                        DistCFuz*BFaceAreazp(i2,i3)
            FEx=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.4) then
c
            dotproduct1=BFaceAreaxp(i2,i3)*BFaceAreanx(i2,i3)+
     *                     BFaceAreayp(i2,i3)*BFaceAreany(i2,i3)+
     *                           BFaceAreazp(i2,i3)*BFaceAreanz(i2,i3)
            dotproduct2=DistCFux*BFaceAreanx(i2,i3)+
     *                     DistCFuy*BFaceAreany(i2,i3)+
     *                             DistCFuz*BFaceAreanz(i2,i3)
c
            ratio=dotproduct1/dotproduct2
c
            FEx=ratio*DistCFux
            FEy=ratio*DistCFuy
            FEz=ratio*DistCFuz
            FE=dabs(ratio)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= fgDiff
          FluxFflocal=-fgDiff
          FluxVflocal=-(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
c
          FluxVf(i4)=FluxVf(i4)+FluxVflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+
     *                FluxFflocal*FiT(j1)+FluxVflocal
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
          cpf=GFactCF*SpecificHeat(i1)+(1.-GFactCF)*SpecificHeat(j1)
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
c
      elseif(variable.eq.'lambda') then
c
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
          BFiT(i2,i3)=0.
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *            dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          BFiT(i2,i3)=FiT(i1)
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

          i4=i4+i3
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *           dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          term1=ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)+
     *                                              BgDiffp(i2,i3)
          FluxCflocal=ConvectionCoefficientRobin(i,i5)*
     *                      BFaceArea(i2,i3)*BgDiffp(i2,i3)/term1
          FluxFflocal=0.
          FluxVflocal=-PhiinfinityRobin(i,i5)*FluxCflocal-
     *       ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)*
     *        (dfidxTf*BFaceTxp(i2,i3)+dfidyTf*BFaceTyp(i2,i3)+
     *                            dfidzTf*BFaceTzp(i2,i3))/term1

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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *                dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= BgDiffp(i2,i3)
          FluxFflocal=-BgDiffp(i2,i3)
          FluxVflocal=-(dfidxTf*BFaceTxp(i2,i3)+
     *               dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
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
c          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
c     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
c     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
c          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
c     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
c     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
c          GFactCF=distance2/(distance1+distance2)
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
          if(MethodDecomposeSprime.eq.1) then
c
            DotProduct=DistCFux*BFaceAreaxp(i2,i3)+
     *                     DistCFuy*BFaceAreayp(i2,i3)+
     *                        DistCFuz*BFaceAreazp(i2,i3)
            FEx=DotProduct*DistCFux
            FEy=DotProduct*DistCFuy
            FEz=DotProduct*DistCFuz
            FE=dabs(DotProduct)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.2) then
c
            FEx=BFaceAreap(i2,i3)*DistCFux
            FEy=BFaceAreap(i2,i3)*DistCFuy
            FEz=BFaceAreap(i2,i3)*DistCFuz
            FE=BFaceAreap(i2,i3)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.3) then
c
            DotProduct=DistCFux*BFaceAreaxp(i2,i3)+
     *                     DistCFuy*BFaceAreayp(i2,i3)+
     *                        DistCFuz*BFaceAreazp(i2,i3)
            FEx=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFux
            FEy=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFuy
            FEz=(BFaceAreap(i2,i3)**2/DotProduct)*DistCFuz
            FE=dsqrt(FEx**2+FEy**2+FEz**2)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          elseif(MethodDecomposeSprime.eq.4) then
c
            dotproduct1=BFaceAreaxp(i2,i3)*BFaceAreanx(i2,i3)+
     *                     BFaceAreayp(i2,i3)*BFaceAreany(i2,i3)+
     *                           BFaceAreazp(i2,i3)*BFaceAreanz(i2,i3)
            dotproduct2=DistCFux*BFaceAreanx(i2,i3)+
     *                     DistCFuy*BFaceAreany(i2,i3)+
     *                             DistCFuz*BFaceAreanz(i2,i3)
c
            ratio=dotproduct1/dotproduct2
c
            FEx=ratio*DistCFux
            FEy=ratio*DistCFuy
            FEz=ratio*DistCFuz
            FE=dabs(ratio)
            fgDiff=FE/DistCF
c
            FTx=BFaceAreaxp(i2,i3)-FEx
            FTy=BFaceAreayp(i2,i3)-FEy
            FTz=BFaceAreazp(i2,i3)-FEz
            FT=dsqrt(FTx**2+FTy**2+FTz**2)
c
          endif
c
          dfidxTf=BdfidxT(i2,i3)
          dfidyTf=BdfidyT(i2,i3)
          dfidzTf=BdfidzT(i2,i3)
c
          FluxCflocal= fgDiff
          FluxFflocal=-fgDiff
          FluxVflocal=-(dfidxTf*FTx+dfidyTf*FTy+dfidzTf*FTz)
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