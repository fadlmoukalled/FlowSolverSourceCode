c
C#############################################################################################
      SUBROUTINE allocateBCStorage
C#############################################################################################
      use User0
      use Geometry1, only:NumberOfBCSets
      use Geometry3, only:NBFacesMax
      use BoundaryConditions1
      use BoundaryConditions2, only:uVelocityFarField,vVelocityFarField,
     *                            wVelocityFarField,
     *                            SpecificHeatFarField,pressureFarField,
     *                                 TemperatureFarField,MachFarField,
     *                    xFlowDirectionFarField,yFlowDirectionFarField,
     *                    zFlowDirectionFarField,PeriodicPair,theta,
     *                    a1r,b1r,c1r,a2r,b2r,c2r,a3r,b3r,c3r,
     *                    a1Axis,a2Axis,a3Axis,
     *              xTranslation,yTranslation,zTranslation,TKEFarField,
     *              TEDFarField,TOmegaFarField,TurbulentKLFarField,
     *              MEDFarField
      use BoundaryConditionsScalar1
      use BoundaryConditionsrField1
      use BoundaryConditionsTurbulence1, only: inletTypeT,
     *                                    inletTurbulenceIntensity,
     *                                    inletTurbulenceLengthScale,
     *                                    inletTurbulenceViscosityRatio
      use BoundaryFluxes
      use WallStress1
      use Variables1, only: BStagnationPressure,BStagnationTemperature
      use FlowInOut1, only: MassFlowFraction,mdotout,mdot1,
     *                      MassFlowRate,LPrintMassFlowRate
      use ArteryResistance1, only: ArteryResistance,geoDiffSum,
     *                             geoDiffB,PressureCResistance,
     *                             LArteryExplicit,urfPressureResistance
      use AveragePressure1, only: AverageOutletPressure
      use WindKessel1, only: ResistanceToBloodFlow,
     *                       TotalPeripheralResistance,ComplianceC,
     *                       InertiaL,OutletPressure,MassFlowRateOld,
     *                       MassFlowRateOldOld,OutletPressureOld,
     *                       OutletPressureOldOld
      use PhysicalProperties1, only: InterpolationSchemeGamaScalar
c
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i
c*********************************************************************************************
c
      allocate(BCType(NumberOfBCSets,NBFacesMax))
      allocate(BoundaryType(NumberOfBCSets))
      allocate(wallTypeLambda(NumberOfBCSets,NBFacesMax))      
      allocate(wallTypeL(NumberOfBCSets))
c
      allocate(wallTypeMomentum(NumberOfBCSets,NBFacesMax))      
      allocate(wallTypeContinuity(NumberOfBCSets,NBFacesMax))      
      allocate(wallTypeEnergy(NumberOfBCSets,NBFacesMax))      
      allocate(wallTypeM(NumberOfBCSets))      
      allocate(wallTypeC(NumberOfBCSets))      
      allocate(wallTypeE(NumberOfBCSets))      
c
      allocate(inletTypeMomentum(NumberOfBCSets,NBFacesMax))      
      allocate(inletTypeContinuity(NumberOfBCSets,NBFacesMax))      
      allocate(inletTypeEnergy(NumberOfBCSets,NBFacesMax))      
      allocate(inletTypeM(NumberOfBCSets))      
      allocate(inletTypeC(NumberOfBCSets))      
      allocate(inletTypeE(NumberOfBCSets))      
c
      if(LTurbulentFlow) then
c
        allocate(inletTypeT(NumberOfBCSets))
        allocate(inletTurbulenceIntensity(NumberOfBCSets))
        allocate(inletTurbulenceLengthScale(NumberOfBCSets))
        allocate(inletTurbulenceViscosityRatio(NumberOfBCSets))
        allocate(TKEFarField(NumberOfBCSets))
        allocate(TEDFarField(NumberOfBCSets))
        allocate(TOmegaFarField(NumberOfBCSets))
        allocate(TurbulentKLFarField(NumberOfBCSets))
        allocate(MEDFarField(NumberOfBCSets))
c
      endif
c
      allocate(outletTypeMomentum(NumberOfBCSets,NBFacesMax))      
      allocate(outletTypeContinuity(NumberOfBCSets,NBFacesMax))      
      allocate(outletTypeEnergy(NumberOfBCSets,NBFacesMax))      
      allocate(outletTypeM(NumberOfBCSets))
      allocate(outletTypeC(NumberOfBCSets))
      allocate(outletTypeE(NumberOfBCSets))
      allocate(MassFlowFraction(NumberOfBCSets))
      allocate(mdotout(NumberOfBCSets))
      allocate(mdot1(NumberOfBCSets))
      allocate(MassFlowRate(NumberOfBCSets))
      allocate(LPrintMassFlowRate(NumberOfBCSets))
      allocate(ArteryResistance(NumberOfBCSets))
      allocate(geoDiffSum(NumberOfBCSets))
      allocate(geoDiffB(NumberOfBCSets,NBFacesMax))
      allocate(PressureCResistance(NumberOfBCSets))
      allocate(LArteryExplicit(NumberOfBCSets))
      allocate(urfPressureResistance(NumberOfBCSets))
      allocate(AverageOutletPressure(NumberOfBCSets))
      allocate(MaximumShearStress(NumberOfBCSets))
      allocate(xLocationMaxWallShearStress(NumberOfBCSets))
      allocate(yLocationMaxWallShearStress(NumberOfBCSets))
      allocate(zLocationMaxWallShearStress(NumberOfBCSets))
c
      allocate(GrainSize(NumberOfBCSets))
c
      allocate(MassFlowRateOld(NumberOfBCSets))
      allocate(MassFlowRateOldOld(NumberOfBCSets))
      allocate(OutletPressure(NumberOfBCSets))
      allocate(OutletPressureOld(NumberOfBCSets))
      allocate(OutletPressureOldOld(NumberOfBCSets))
      allocate(ResistanceToBloodFlow(NumberOfBCSets))
      allocate(TotalPeripheralResistance(NumberOfBCSets))
      allocate(ComplianceC(NumberOfBCSets))
      allocate(InertiaL(NumberOfBCSets))
      allocate(LWindKessel(NumberOfBCSets))
c
      allocate(HeatFlux(NumberOfBCSets,NBFacesMax))      
      allocate(Tinfinity(NumberOfBCSets,NBFacesMax))      
      allocate(Hinfinity(NumberOfBCSets,NBFacesMax))      
c
      allocate(BStagnationPressure(NumberOfBCSets,NBFacesMax))
      allocate(BStagnationTemperature(NumberOfBCSets,NBFacesMax))
c
      allocate(uVelocityFarField(NumberOfBCSets))
      allocate(vVelocityFarField(NumberOfBCSets))
      allocate(wVelocityFarField(NumberOfBCSets))
      allocate(SpecificHeatFarField(NumberOfBCSets))
      allocate(pressureFarField(NumberOfBCSets))
      allocate(TemperatureFarField(NumberOfBCSets))
      allocate(MachFarField(NumberOfBCSets))
      allocate(xFlowDirectionFarField(NumberOfBCSets))
      allocate(yFlowDirectionFarField(NumberOfBCSets))
      allocate(zFlowDirectionFarField(NumberOfBCSets))
c
      allocate(PeriodicPair(NumberOfBCSets))
      allocate(theta(NumberOfBCSets))
      allocate(a1r(NumberOfBCSets))
      allocate(b1r(NumberOfBCSets))
      allocate(c1r(NumberOfBCSets))
      allocate(a2r(NumberOfBCSets))
      allocate(b2r(NumberOfBCSets))
      allocate(c2r(NumberOfBCSets))
      allocate(a3r(NumberOfBCSets))
      allocate(b3r(NumberOfBCSets))
      allocate(c3r(NumberOfBCSets))
      allocate(a1Axis(NumberOfBCSets))
      allocate(a2Axis(NumberOfBCSets))
      allocate(a3Axis(NumberOfBCSets))
      allocate(xTranslation(NumberOfBCSets))
      allocate(yTranslation(NumberOfBCSets))
      allocate(zTranslation(NumberOfBCSets))
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        i=NumberOfScalarsToSolve
c
        allocate(wallTypeScalar(NumberOfBCSets,NBFacesMax,i))      
        allocate(inletTypeScalar(NumberOfBCSets,NBFacesMax,i))      
        allocate(outletTypeScalar(NumberOfBCSets,NBFacesMax,i))      
        allocate(wallTypeS(NumberOfBCSets,i))      
        allocate(inletTypeS(NumberOfBCSets,i))      
        allocate(outletTypeS(NumberOfBCSets,i))      
c
        allocate(ScalarFlux(NumberOfBCSets,NBFacesMax,i))
        allocate(ScalarConvectionCoefficient
     *                         (NumberOfBCSets,NBFacesMax,i))
        allocate(ScalarPhiInfinity(NumberOfBCSets,NBFacesMax,i))
c
        allocate(ScalarName(i))
        allocate(MethodCalcGradientScalar(i))
        allocate(nIterGradientScalar(i))
        allocate(LimitGradientScalar(i))
        allocate(LimitGradientScalarMethod(i))
        allocate(LSolveScalar(i))
        allocate(LMultigridScalar(i))
        allocate(rrFScalar(i))
        allocate(ASSolverScalar(i))
        allocate(ASIterScalar(i))
        allocate(LRelaxScalar(i))
        allocate(LFalseTransientScalar(i))
        allocate(urfScalar(i))
        allocate(FalseDtScalar(i))
        allocate(LNVFScalar(i))
        allocate(LTVDScalar(i))
        allocate(ConvectionSchemeScalar(i))
        allocate(BleedScalar(i))
        allocate(GradientInterpolationSchemeScalar(i))
        allocate(ConstantDiffusionCoefficientScalar(i))
        allocate(ConstantSpecificHeatScalar(i))
        allocate(InterpolationSchemeGamaScalar(i))
        allocate(CRoughScalar(i))
c
      endif
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        i=NumberOfrFieldsToSolve
c
        allocate(wallTyperField(NumberOfBCSets,NBFacesMax,i))      
        allocate(inletTyperField(NumberOfBCSets,NBFacesMax,i))      
        allocate(outletTyperField(NumberOfBCSets,NBFacesMax,i))      
        allocate(wallTypeR(NumberOfBCSets,i))      
        allocate(inletTypeR(NumberOfBCSets,i))      
        allocate(outletTypeR(NumberOfBCSets,i))      
c
        allocate(rFieldFlux(NumberOfBCSets,NBFacesMax,i))
        allocate(rFieldConvectionCoefficient
     *                         (NumberOfBCSets,NBFacesMax,i))
        allocate(rFieldPhiInfinity(NumberOfBCSets,NBFacesMax,i))
c
        allocate(rFieldName(i))
        allocate(MethodCalcGradientrField(i))
        allocate(nIterGradientrField(i))
        allocate(LimitGradientrField(i))
        allocate(LimitGradientrFieldMethod(i))
        allocate(LSolverField(i))
        allocate(LMultigridrField(i))
        allocate(rrFrField(i))
        allocate(ASSolverrField(i))
        allocate(ASIterrField(i))
        allocate(LRelaxrField(i))
        allocate(LFalseTransientrField(i))
        allocate(urfrField(i))
        allocate(FalseDtrField(i))
        allocate(LNVFrField(i))
        allocate(LTVDrField(i))
        allocate(ConvectionSchemerField(i))
        allocate(BleedrField(i))
        allocate(GradientInterpolationSchemerField(i))
        allocate(ConstantConductivityrField(i))
        allocate(ConstantSpecificHeatrField(i))
        allocate(ConstantDensityrField(i))
        allocate(ConstantViscosityrField(i))
c
      endif
c
c--- Initializze
c
      MassFlowFraction=1.
      HeatFlux=0.     
      Tinfinity=0.   
      Hinfinity=0.
      ScalarFlux=0.
      ScalarConvectionCoefficient=0.
      ScalarPhiInfinity=0.
      BStagnationPressure=constantStagnationPressure
      BStagnationTemperature=constantStagnationTemperature
c
      mdotout=0.
      mdot1=0.
      MassFlowRate=0.
      LPrintMassFlowRate=.false.
      ArteryResistance=0.
      geoDiffSum=0.
      geoDiffB=0.
      PressureCResistance=0.
      LArteryExplicit=.false.
      urfPressureResistance=1.
      AverageOutletPressure=0.
c
      MassFlowRateOld=0.
      MassFlowRateOldOld=0.
      OutletPressure=0.
      OutletPressureOld=0.
      OutletPressureOldOld=0.
      ResistanceToBloodFlow=0.
      TotalPeripheralResistance=0.
      ComplianceC=0.
      InertiaL=0.
      LWindKessel=.false.
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        ScalarName='scalar'
        MethodCalcGradientScalar=2
        nIterGradientScalar=2
        LimitGradientScalar=.false.
        LimitGradientScalarMethod=2
        LSolveScalar=.false.
        LMultigridScalar=.false.
        rrFScalar=0.3
        ASSolverScalar='ilu'
        ASIterScalar=3
        LRelaxScalar=.false.
        LFalseTransientScalar=.false.
        urfScalar=1.
        FalseDtScalar=1.
        LNVFScalar=.false.
        LTVDScalar=.false.
        ConvectionSchemeScalar='upwind'
        BleedScalar=0.
        GradientInterpolationSchemeScalar='average'
        ConstantDiffusionCoefficientScalar=1.
        ConstantSpecificHeatScalar=1.
        InterpolationSchemeGamaScalar='average'
        CRoughScalar=0.2
c
      endif
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        i=NumberOfrFieldsToSolve
c
        rFieldFlux=0.
        rFieldConvectionCoefficient=0.
        rFieldPhiInfinity=0.
        rFieldName='rfield'
        MethodCalcGradientrField=2
        nIterGradientrField=2
        LimitGradientrField=.false.
        LimitGradientrFieldMethod=2
        LSolverField=.false.
        LMultigridrField=.false.
        rrFrField=0.3
        ASSolverrField='ilu'
        ASIterrField=3
        LRelaxrField=.false.
        LFalseTransientrField=.false.
        urfrField=1.
        FalseDtrField=1.
        LNVFrField=.false.
        LTVDrField=.false.
        ConvectionSchemerField='upwind'
        BleedrField=0.
        GradientInterpolationSchemerField='average'
        ConstantConductivityrField=1.
        ConstantSpecificHeatrField=1.
        ConstantDensityrField=1.
        ConstantViscosityrField=1.
c
      endif
c
      return
      end