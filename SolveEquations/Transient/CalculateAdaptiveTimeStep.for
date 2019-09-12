c
C#############################################################################################
c
      SUBROUTINE AdaptiveTimeStep
c
C#############################################################################################
c
      use User0, only: LSolveMomentum,LSolveTurbulenceKineticEnergy,
     *                 LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LSolveTurbulentKL,LSolveModifiedED,LSolveEnergy,
     *                 LsolveScalar,LsolverField,Linviscid,
     *                 EnergyEquation,NumberOfScalarsToSolve,
     *                 LFreeSurfaceFlow,NumberOfrFieldsToSolve,
     *                 minimumdtChangeFactor,maximumdtChangeFactor,
     *                 dt,LanisotropicDiffusion,
     *                 LSolveTurbulenceGammaEquation,
     *                 LSolveTurbulenceReynoldsThetaEquation,
     *                 LSolveTurbulenceV2Equation,
     *                 LSolveTurbulenceZetaEquation
      use Constants1, only: big
      use Scalar2, only: iScalarVariable
      use VolumeOfFluid2, only: irFieldVariable
c--------------------------------------------------------------------------------
      implicit none  
      character*10 :: Variable
      integer :: iScalar,irField
      double precision, save :: dtlocal,dtlocal1
c--------------------------------------------------------------------------------
c
c--- Calculate time step if momentum is solved
c
      dtlocal1=big
      if(LSolveMomentum) then
        Variable='velx'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if turbulence kinetic energy is solved
c
      if(LSolveTurbulenceKineticEnergy) then
        Variable='tke'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if turbulence dissipation rate is solved
c
      if(LSolveTurbulenceDissipationRate) then
        Variable='ted'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if turbulence specific dissipation rate is solved
c
      if(LSolveTurbulenceSpecificDissipationRate) then
        Variable='tomega'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if turbulence Gamma in transition model is solved
c
      if(LSolveTurbulenceGammaEquation) then
        Variable='tgamma'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if turbulence Re Theta in transition model is solved
c
      if(LSolveTurbulenceReynoldsThetaEquation) then
        Variable='retheta'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if turbulent KL equation is solved
c
      if(LSolveTurbulentKL) then
        Variable='tkl'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if the v2 equation is solved
c
      if(LSolveTurbulenceV2Equation) then
        Variable='tv2'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if the v2 equation is solved
c
      if(LSolveTurbulenceZetaEquation) then
        Variable='tzeta'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if the modified eddy diffusivity equation is solved
c
      if(LSolveModifiedED) then
        Variable='med'
        if(.not.Linviscid)
     *       call CalculateEffectiveDiffusionCoefficient(Variable)
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if energy is solved
c
      if(LSolveEnergy.and.EnergyEquation.eq.'temperature') then
        Variable='temp'
        if(.not.Linviscid) then
          if(LanisotropicDiffusion) then
            call calculateSprime(Variable)
          else
            call CalculateEffectiveDiffusionCoefficient(Variable)
          endif
        endif
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
      if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
        if(.not.Linviscid) then
          Variable='temp'
          if(LanisotropicDiffusion) then
            call calculateSprime(Variable)
          else
            call CalculateEffectiveDiffusionCoefficient(Variable)
          endif
        endif
        Variable='htotal'
        call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
        dtlocal1=dmin1(dtlocal1,dtlocal)
      endif
c
c--- Calculate time step if scalars are solved
c
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) then
          iScalarVariable=iScalar
          Variable='scalar'
          if(.not.Linviscid) then
            if(LanisotropicDiffusion) then
              call calculateSprime(Variable)
            else
              call CalculateEffectiveDiffusionCoefficient(Variable)
            endif
          endif
          call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
          dtlocal1=dmin1(dtlocal1,dtlocal)
        endif
      enddo
c
c--- Calculate time step if volume fractions are solved
c
      if(LFreeSurfaceFlow) then
        do irField=1,NumberOfrFieldsToSolve
          if(LsolverField(irField)) then
            irFieldVariable=irField
            Variable='rfield'
            call CalculateAdaptiveTimeStepValue(Variable,dtlocal)
            dtlocal1=dmin1(dtlocal1,dtlocal)
          endif
        enddo
      endif
c
      dtlocal1=dmax1(dtlocal1,minimumdtChangeFactor*dt)
      dt=dmin1(dtlocal1,maximumdtChangeFactor*dt)
c      
      write(*,*) 'calculated time step = ',dt
      write(13,*) 'calculated time step = ',dt
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE CalculateAdaptiveTimeStepValue(Variable,dtlocal)
c
C#############################################################################################
c
      use User0, only: dt,GlobalCourantNumber,
     *                 minimumdtChangeFactor,maximumdtChangeFactor,
     *                 EnergyEquation,MethodDecomposeS,
     *                 LanisotropicDiffusion,MethodDecomposeSprime
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces,
     *                     NIFaceOwner,NIFaceNeighbor,NBFaceOwner
      use Geometry4, only: GFactorCF,gDiff,BgDiff,Volume,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BFaceArea,xc,yc,zc,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz,
     *                     BFaceAreax,BFaceAreay,BFaceAreaz,
     *                     gDiffp,BgDiffp,BFaceAreap,
     *                     BFaceAreaxp,BFaceAreayp,BFaceAreazp
      use Variables1, only: mdot,Bmdot
      use PhysicalProperties1, only: Densityf,Density,BDensity,
     *                               eDiffCoefficient,BeDiffCoefficient,
     *                               SpecificHeat,BSpecificHeat,
     *                               SpecificHeatScalar,
     *                               BSpecificHeatScalar
      use BoundaryConditions1
      use BoundaryConditions2
      use BoundaryConditionsScalar1
      use BoundaryConditionsScalar2
      use BoundaryConditionsTurbulence1
      use BoundaryConditionsTurbulence2
      use scalar2, only: iScalarVariable
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,i1,i2,i3,i5,j,j1,j2,j3,k
      double precision :: dtlocal
      double precision, save :: gamf,cpf,dNorm,FluxCflocal,
     *                          area,nx,ny,nz,xF1,yF1,zF1,term1,
     *                          distance1,distance2,GFactCF,
     *                          DistCFx,DistCFy,DistCFz,
     *                          DistCFux,DistCFuy,DistCFuz,DistCF,
     *                          DotProduct,FEx,FEy,FEz,FE,fgDiff,
     *                          DotProduct1,DotProduct2,ratio
      double precision, save, dimension(:), allocatable :: TimeStep
c********************************************************************************************
c
      allocate(TimeStep(NumberOfElements))
c
      TimeStep=0.
c
c---- Calculate diffusion contribution to time step
c
      if(Variable.eq.'rfield') then
c      
        continue
c
      elseif(Variable.eq.'temp'.and.EnergyEquation.eq.'htotal') then
c
        if(LanisotropicDiffusion) then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            cpf=GFactorCF(k)*SpecificHeat(i)+
     *                 (1.-GFactorCF(k))*SpecificHeat(j)
c
            TimeStep(i)= TimeStep(i)+gDiffp(k)/cpf
            TimeStep(j)= TimeStep(j)+gDiffp(k)/cpf
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
            gamf=GFactorCF(k)*eDiffCoefficient(i)+
     *                 (1.-GFactorCF(k))*eDiffCoefficient(j)
c
            cpf=GFactorCF(k)*SpecificHeat(i)+
     *                 (1.-GFactorCF(k))*SpecificHeat(j)
c
            TimeStep(i)= TimeStep(i)+gamf*gDiff(k)/cpf
            TimeStep(j)= TimeStep(j)+gamf*gDiff(k)/cpf
c
          enddo
c
        endif
c
      elseif(Variable.eq.'temp'.and.EnergyEquation.eq.'temperature')then
c
        if(LanisotropicDiffusion) then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            TimeStep(i)= TimeStep(i)+gDiffp(k)
            TimeStep(j)= TimeStep(j)+gDiffp(k)
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
            gamf=GFactorCF(k)*eDiffCoefficient(i)+
     *                 (1.-GFactorCF(k))*eDiffCoefficient(j)
c
            TimeStep(i)= TimeStep(i)+gamf*gDiff(k)
            TimeStep(j)= TimeStep(j)+gamf*gDiff(k)
c
          enddo
c
        endif
c
      elseif(Variable.eq.'velx'.or.Variable.eq.'tke'.or.
     *        Variable.eq.'ted'.or.Variable.eq.'tomega'.or.
     *         Variable.eq.'med'.or.Variable.eq.'tkl'.or.
     *           Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *              Variable.eq.'tv2'.or.Variable.eq.'tzeta') then
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gamf=GFactorCF(k)*eDiffCoefficient(i)+
     *                 (1.-GFactorCF(k))*eDiffCoefficient(j)
c
          TimeStep(i)= TimeStep(i)+gamf*gDiff(k)
          TimeStep(j)= TimeStep(j)+gamf*gDiff(k)
c
        enddo
c
      elseif(Variable.eq.'scalar') then
c
        if(LanisotropicDiffusion) then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            TimeStep(i)= TimeStep(i)+gDiffp(k)
            TimeStep(j)= TimeStep(j)+gDiffp(k)
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
            gamf=GFactorCF(k)*eDiffCoefficient(i)+
     *                 (1.-GFactorCF(k))*eDiffCoefficient(j)
c
            TimeStep(i)= TimeStep(i)+gamf*gDiff(k)
            TimeStep(j)= TimeStep(j)+gamf*gDiff(k)
c
          enddo
c
        endif
c
      endif
c
c----------------------------------------------------------------------
c--- Boundary faces
c----------------------------------------------------------------------
      if(variable.eq.'velx') then
c----------------------------------------------------------------------
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
c
          gamf=BeDiffCoefficient(i2,i3)
          dNorm=BDistanceCFx(i2,i3)*BFaceAreanx(i2,i3)+
     *                BDistanceCFy(i2,i3)*BFaceAreany(i2,i3)+
     *                BDistanceCFz(i2,i3)*BFaceAreanz(i2,i3)
          area=BFaceArea(i2,i3)
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
c
          FluxCflocal= gamf*area*(1.-nx**2)/dNorm
          FluxCflocal= dmax1(FluxCflocal,gamf*area*(1.-ny**2)/dNorm)
          FluxCflocal= dmax1(FluxCflocal,gamf*area*(1.-nz**2)/dNorm)
c
          TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iinletsupersonic
c
          i1=IinletsupersonicOwner(i)
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
c
          TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i1=IinletSpecifiedVelocityOwner(i)
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
c
          TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i1=IinletSpecifiedStaticPressureOwner(i)
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
c
          TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationPressure
c
          i1=IinletSpecifiedStagnationPressureOwner(i)
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
c
          TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i1=IoutletspecifiedVelocityOwner(i)
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
c
          TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
c
          TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          gamf=BeDiffCoefficient(i2,i3)
c
          dNorm=BDistanceCFx(i2,i3)*BFaceAreanx(i2,i3)+
     *                BDistanceCFy(i2,i3)*BFaceAreany(i2,i3)+
     *                BDistanceCFz(i2,i3)*BFaceAreanz(i2,i3)
          area=BFaceArea(i2,i3)
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
c
          FluxCflocal= 2.*gamf*area*(nx**2)/dNorm
          FluxCflocal= dmax1(FluxCflocal,2.*gamf*area*(ny**2)/dNorm)
          FluxCflocal= dmax1(FluxCflocal,2.*gamf*area*(nz**2)/dNorm)
c
          TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
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
          endif
c
          gamf=GFactCF*eDiffCoefficient(i1)+
     *                   (1.-GFactCF)*eDiffCoefficient(j1)
          TimeStep(i1)=TimeStep(i1)+gamf*fgDiff
c
        enddo
c----------------------------------------------------------------------
      elseif(variable.eq.'tke'.or.variable.eq.'ted'.or.
     *        variable.eq.'tomega'.or.variable.eq.'med'.or.
     *          variable.eq.'tkl'.or.variable.eq.'tgamma'.or.
     *            variable.eq.'tretheta'.or.variable.eq.'tv2'.or.
     *                                    variable.eq.'tzeta') then
c----------------------------------------------------------------------
c
        if(variable.eq.'med') then
c
          do i=1,IWallSlip
c
            i1=IWallSlipOwner(i)
            i2=IWallSlipNumberOfBCSets(i)
            i3=IWallSlipNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *                    BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
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
c
          if(Bmdot(i2,i3).lt.0.) then
c
            TimeStep(i1)=TimeStep(i1)+
     *                    BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
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
c
          TimeStep(i1)=TimeStep(i1)+
     *                    BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          elseif(MethodDecomposeS.eq.2) then
c
            FEx=BFaceArea(i2,i3)*DistCFux
            FEy=BFaceArea(i2,i3)*DistCFuy
            FEz=BFaceArea(i2,i3)*DistCFuz
            FE=BFaceArea(i2,i3)
            fgDiff=FE/DistCF
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
          endif
c
          gamf=GFactCF*eDiffCoefficient(i1)+
     *                    (1.-GFactCF)*eDiffCoefficient(j1)
c
          TimeStep(i1)=TimeStep(i1)+gamf*fgDiff
c
        enddo
c
c----------------------------------------------------------------------
      elseif(variable.eq.'temp'.and.
     *                EnergyEquation.eq.'temperature') then
c----------------------------------------------------------------------
        if(LanisotropicDiffusion) then
c
          do i=1,IWallDirichlet
c
            i1=IWallDirichletOwner(i)
            i2=IWallDirichletNumberOfBCSets(i)
            i3=IWallDirichletNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+BgDiffp(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+BgDiffp(i2,i3)
            FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                                  BgDiffp(i2,i3)/term1
c
            TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSupersonic
c
            i1=IinletSupersonicOwner(i)
            i2=IinletSupersonicNumberOfBCSets(i)
            i3=IinletSupersonicNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+BgDiffp(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStaticTemperature
c
            i1=IinletSpecifiedStaticTemperatureOwner(i)
            i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStaticTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+BgDiffp(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStagnationTemperature
c
            i1=IinletSpecifiedStagnationTemperatureOwner(i)
            i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+BgDiffp(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
c
            j2=PeriodicPair(i2)         
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
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
            elseif(MethodDecomposeSprime.eq.2) then
c
              FEx=BFaceAreap(i2,i3)*DistCFux
              FEy=BFaceAreap(i2,i3)*DistCFuy
              FEz=BFaceAreap(i2,i3)*DistCFuz
              FE=BFaceAreap(i2,i3)
              fgDiff=FE/DistCF
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
            endif
c
            TimeStep(i1)=TimeStep(i1)+fgDiff
c
          enddo
c          
        else       
c
          do i=1,IWallDirichlet
c
            i1=IWallDirichletOwner(i)
            i2=IWallDirichletNumberOfBCSets(i)
            i3=IWallDirichletNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
            gamf=BeDiffCoefficient(i2,i3)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+gamf*BgDiff(i2,i3)
            FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                                  BgDiff(i2,i3)*gamf/term1
c
            TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSupersonic
c
            i1=IinletSupersonicOwner(i)
            i2=IinletSupersonicNumberOfBCSets(i)
            i3=IinletSupersonicNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStaticTemperature
c
            i1=IinletSpecifiedStaticTemperatureOwner(i)
            i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStaticTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStagnationTemperature
c
            i1=IinletSpecifiedStagnationTemperatureOwner(i)
            i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *              BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
c
            j2=PeriodicPair(i2)         
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
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
            elseif(MethodDecomposeS.eq.2) then
c
              FEx=BFaceArea(i2,i3)*DistCFux
              FEy=BFaceArea(i2,i3)*DistCFuy
              FEz=BFaceArea(i2,i3)*DistCFuz
              FE=BFaceArea(i2,i3)
              fgDiff=FE/DistCF
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
            endif
c
            gamf=GFactCF*eDiffCoefficient(i1)+
     *                   (1.-GFactCF)*eDiffCoefficient(j1)
c
            TimeStep(i1)=TimeStep(i1)+gamf*fgDiff
c
          enddo
c          
        endif
c
c----------------------------------------------------------------------
      elseif(variable.eq.'temp'.and.EnergyEquation.eq.'htotal') then
c----------------------------------------------------------------------
        if(LanisotropicDiffusion) then
c
          do i=1,IWallDirichlet
c
            i1=IWallDirichletOwner(i)
            i2=IWallDirichletNumberOfBCSets(i)
            i3=IWallDirichletNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *                        BgDiffp(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+BgDiffp(i2,i3)
            FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                                  BgDiffp(i2,i3)/term1
c
            TimeStep(i1)=TimeStep(i1)+FluxCflocal/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSupersonic
c
            i1=IinletSupersonicOwner(i)
            i2=IinletSupersonicNumberOfBCSets(i)
            i3=IinletSupersonicNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *                      BgDiffp(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStaticTemperature
c
            i1=IinletSpecifiedStaticTemperatureOwner(i)
            i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStaticTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *                       BgDiffp(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStagnationTemperature
c
            i1=IinletSpecifiedStagnationTemperatureOwner(i)
            i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *                       BgDiffp(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
c
            j2=PeriodicPair(i2)         
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
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
            elseif(MethodDecomposeSprime.eq.2) then
c
              FEx=BFaceAreap(i2,i3)*DistCFux
              FEy=BFaceAreap(i2,i3)*DistCFuy
              FEz=BFaceAreap(i2,i3)*DistCFuz
              FE=BFaceAreap(i2,i3)
              fgDiff=FE/DistCF
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
            endif
c
            cpf=GFactCF*SpecificHeat(i1)+(1.-GFactCF)*SpecificHeat(j1)
c
            TimeStep(i1)=TimeStep(i1)+fgDiff/cpf
c
          enddo
c
        else       
c
          do i=1,IWallDirichlet
c
            i1=IWallDirichletOwner(i)
            i2=IWallDirichletNumberOfBCSets(i)
            i3=IWallDirichletNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *      BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
            gamf=BeDiffCoefficient(i2,i3)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+gamf*BgDiff(i2,i3)
            FluxCflocal=HinfinityRobin(i)*BFaceArea(i2,i3)*
     *                                  BgDiff(i2,i3)*gamf/term1
c
            TimeStep(i1)=TimeStep(i1)+FluxCflocal/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSupersonic
c
            i1=IinletSupersonicOwner(i)
            i2=IinletSupersonicNumberOfBCSets(i)
            i3=IinletSupersonicNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *       BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStaticTemperature
c
            i1=IinletSpecifiedStaticTemperatureOwner(i)
            i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStaticTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *       BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedStagnationTemperature
c
            i1=IinletSpecifiedStagnationTemperatureOwner(i)
            i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
            i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
c
            TimeStep(i1)=TimeStep(i1)+
     *      BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)/BSpecificHeat(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
c
            j2=PeriodicPair(i2)         
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
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
            elseif(MethodDecomposeS.eq.2) then
c
              FEx=BFaceArea(i2,i3)*DistCFux
              FEy=BFaceArea(i2,i3)*DistCFuy
              FEz=BFaceArea(i2,i3)*DistCFuz
              FE=BFaceArea(i2,i3)
              fgDiff=FE/DistCF
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
            endif
c
            gamf=GFactCF*eDiffCoefficient(i1)+
     *                     (1.-GFactCF)*eDiffCoefficient(j1)
            cpf=GFactCF*SpecificHeat(i1)+(1.-GFactCF)*SpecificHeat(j1)
c
            TimeStep(i1)=TimeStep(i1)+gamf*fgDiff/cpf
c
          enddo
c
        endif
c----------------------------------------------------------------------
      elseif(Variable.eq.'scalar') then
c----------------------------------------------------------------------
c
        i5=iScalarVariable
c
        if(LanisotropicDiffusion) then
c
          do i=1,IWallDirichletScalar(i5)
c
            i1=IWallDirichletScalarOwner(i,i5)
            i2=IWallDirichletScalarNumberOfBCSets(i,i5)
            i3=IWallDirichletScalarNBFaces(i,i5)
c
            TimeStep(i1)=TimeStep(i1)+BgDiffp(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IWallRobinScalar(i5)
c
            i1=IWallRobinScalarOwner(i,i5)
            i2=IWallRobinScalarNumberOfBCSets(i,i5)
            i3=IWallRobinScalarNBFaces(i,i5)
c
            term1=ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)+
     *                                                 BgDiffp(i2,i3)
            FluxCflocal=ConvectionCoefficientRobin(i,i5)*
     *                      BFaceArea(i2,i3)*BgDiffp(i2,i3)/term1
            TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSupersonicScalar(i5)
c
            i1=IinletSupersonicScalarOwner(i,i5)
            i2=IinletSupersonicScalarNumberOfBCSets(i,i5)
            i3=IinletSupersonicScalarNBFaces(i,i5)
c
            TimeStep(i1)=TimeStep(i1)+BgDiffp(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedValueScalar(i5)
c
            i1=IinletSpecifiedValueScalarOwner(i,i5)
            i2=IinletSpecifiedValueScalarNumberOfBCSets(i,i5)
            i3=IinletSpecifiedValueScalarNBFaces(i,i5)
c
            TimeStep(i1)=TimeStep(i1)+BgDiffp(i2,i3)
c
          enddo
c----------------------------------------------------------------------
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
c
            j2=PeriodicPair(i2)         
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
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
            elseif(MethodDecomposeSprime.eq.2) then
c
              FEx=BFaceAreap(i2,i3)*DistCFux
              FEy=BFaceAreap(i2,i3)*DistCFuy
              FEz=BFaceAreap(i2,i3)*DistCFuz
              FE=BFaceAreap(i2,i3)
              fgDiff=FE/DistCF
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
            endif
c
            TimeStep(i1)=TimeStep(i1)+fgDiff
c
          enddo
c
        else       
c
          do i=1,IWallDirichletScalar(i5)
c
            i1=IWallDirichletScalarOwner(i,i5)
            i2=IWallDirichletScalarNumberOfBCSets(i,i5)
            i3=IWallDirichletScalarNBFaces(i,i5)
c
            FluxCflocal= BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
            TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
          enddo
c----------------------------------------------------------------------
          do i=1,IWallRobinScalar(i5)
c
            i1=IWallRobinScalarOwner(i,i5)
            i2=IWallRobinScalarNumberOfBCSets(i,i5)
            i3=IWallRobinScalarNBFaces(i,i5)
c
            gamf=BeDiffCoefficient(i2,i3)
c
            term1=ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)+
     *                                              gamf*BgDiff(i2,i3)
            FluxCflocal=ConvectionCoefficientRobin(i,i5)*
     *                      BFaceArea(i2,i3)*BgDiff(i2,i3)*gamf/term1
            TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSupersonicScalar(i5)
c
            i1=IinletSupersonicScalarOwner(i,i5)
            i2=IinletSupersonicScalarNumberOfBCSets(i,i5)
            i3=IinletSupersonicScalarNBFaces(i,i5)
c
            FluxCflocal= BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
            TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
          enddo
c----------------------------------------------------------------------
          do i=1,IinletSpecifiedValueScalar(i5)
c
            i1=IinletSpecifiedValueScalarOwner(i,i5)
            i2=IinletSpecifiedValueScalarNumberOfBCSets(i,i5)
            i3=IinletSpecifiedValueScalarNBFaces(i,i5)
c
            FluxCflocal= BeDiffCoefficient(i2,i3)*BgDiff(i2,i3)
c
            TimeStep(i1)=TimeStep(i1)+FluxCflocal
c
          enddo
c----------------------------------------------------------------------
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
c
            j2=PeriodicPair(i2)         
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
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
            elseif(MethodDecomposeS.eq.2) then
c
              FEx=BFaceArea(i2,i3)*DistCFux
              FEy=BFaceArea(i2,i3)*DistCFuy
              FEz=BFaceArea(i2,i3)*DistCFuz
              FE=BFaceArea(i2,i3)
              fgDiff=FE/DistCF
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
            endif
c
            gamf=GFactCF*eDiffCoefficient(i1)+
     *                      (1.-GFactCF)*eDiffCoefficient(j1)
            TimeStep(i1)=TimeStep(i1)+gamf*fgDiff
c
          enddo
c
        endif
c
      endif
c
c---- Calculate convection contribution to time step
c
      if(Variable.eq.'rfield') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          TimeStep(i)=TimeStep(i)+dmax1(mdot(k)/Densityf(k),0.)    
          TimeStep(j)=TimeStep(j)-dmax1(-mdot(k)/Densityf(k),0.)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            TimeStep(k)=TimeStep(k)+dmax1(Bmdot(i,j)/BDensity(i,j),0.)
c
          enddo  
        enddo  
c
      elseif(variable.eq.'temp'.and.
     *                EnergyEquation.eq.'temperature') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          cpf=GFactCF*SpecificHeat(i)+(1.-GFactCF)*SpecificHeat(j)
          TimeStep(i)=TimeStep(i)+dmax1(mdot(k)*cpf,0.)    
          TimeStep(j)=TimeStep(j)-dmax1(-mdot(k)*cpf,0.)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            cpf=BSpecificHeat(i,j)
            TimeStep(k)=TimeStep(k)+dmax1(Bmdot(i,j)*cpf,0.)
c
          enddo  
        enddo  
c
      elseif(variable.eq.'scalar') then
c
        i5=iScalarVariable
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          cpf=GFactCF*SpecificHeatScalar(i,i5)+
     *                 (1.-GFactCF)*SpecificHeatScalar(j,i5)
          TimeStep(i)=TimeStep(i)+dmax1(mdot(k)*cpf,0.)    
          TimeStep(j)=TimeStep(j)-dmax1(-mdot(k)*cpf,0.)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            cpf=BSpecificHeat(i,j)
            TimeStep(k)=TimeStep(k)+dmax1(Bmdot(i,j)*cpf,0.)
c
          enddo  
        enddo  
c
      else
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          TimeStep(i)=TimeStep(i)+dmax1(mdot(k),0.)    
          TimeStep(j)=TimeStep(j)-dmax1(-mdot(k),0.)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            TimeStep(k)=TimeStep(k)+dmax1(Bmdot(i,j),0.)
c
          enddo  
        enddo  
c
      endif
c
c--- Select the minimum time step
c
      dtlocal=-1
c
      if(variable.eq.'temp'.and.
     *                EnergyEquation.eq.'temperature') then
c
        do i=1,NumberOfElements
c
          TimeStep(i)=TimeStep(i)/(Density(i)*SpecificHeat(i)*Volume(i))
          dtlocal=dmax1(dtlocal,TimeStep(i))
c
        enddo
c
      elseif(variable.eq.'scalar') then
c
        i5=iScalarVariable
        do i=1,NumberOfElements
c
          TimeStep(i)=TimeStep(i)/
     *           (Density(i)*SpecificHeatScalar(i,i5)*Volume(i))
          dtlocal=dmax1(dtlocal,TimeStep(i))
c
        enddo
c
      elseif(variable.eq.'rfield') then
c
        do i=1,NumberOfElements
c
          TimeStep(i)=TimeStep(i)/(Volume(i))
          dtlocal=dmax1(dtlocal,TimeStep(i))
c
        enddo

      else
c
        do i=1,NumberOfElements
c
          TimeStep(i)=TimeStep(i)/(Density(i)*Volume(i))
          dtlocal=dmax1(dtlocal,TimeStep(i))
c
        enddo
c
      endif
c
      dtlocal=GlobalCourantNumber/dtlocal
c
      deallocate(TimeStep)
c
      return
      end