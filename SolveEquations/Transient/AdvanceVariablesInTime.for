c
C#############################################################################################
c
      SUBROUTINE updateVariablesInTime
c
C#############################################################################################
c
      use User0, only: LsolveMomentum,LsolveContinuity,LsolveEnergy,dt,
     *                 NumberOfScalarsToSolve,LSolveScalar,
     *                 LSolveTurbulenceKineticEnergy,
     *                 LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LSolveTurbulentKL,LSolveModifiedED,
     *                 EnergyEquation,NumberOfrFieldsToSolve,
     *                 LSolverField,LWindKessel,
     *                 LTurbulentFlow,TurbulenceModel,
     *                 LSpalartAllmarasRotationCurvatureCorrection,
     *                 LWrayAgarwalRotationCurvatureCorrection,
     *                 RotationCurvatureMethod,
     *                 LSolveTurbulenceGammaEquation,
     *                 LSolveTurbulenceReynoldsThetaEquation,
     *                 LSolveTurbulenceV2Equation,
     *                 LSolveTurbulenceZetaEquation
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces 
      use Scalar1, only: Scalar,ScalarOld,ScalarOldOld,
     *                   BScalar,BScalarOld,BScalarOldOld
      use VolumeOfFluid1, only: rField,rFieldOld,rFieldOldOld,
     *                          BrField,BrFieldOld,BrFieldOldOld
      use Transient1
      use Variables1, only: uVelocity,uVelocityOld,uVelocityOldOld,
     *                      vVelocity,vVelocityOld,vVelocityOldOld,
     *                      wVelocity,wVelocityOld,wVelocityOldOld,
     *                      Pressure,PressureOld,PressureOldOld,
     *                      Temperature,TemperatureOld,
     *                      TemperatureOldOld,
     *                      BuVelocity,BuVelocityOld,BuVelocityOldOld,
     *                      BvVelocity,BvVelocityOld,BvVelocityOldOld,
     *                      BwVelocity,BwVelocityOld,BwVelocityOldOld,
     *                      BPressure,BPressureOld,BPressureOldOld,
     *                      BTemperature,BTemperatureOld,
     *                      BTemperatureOldOld,
     *                      mdotOldOld,mdotOld,mdot,
     *                      BmdotOldOld,BmdotOld,Bmdot,
     *                      TurbulentKE,TurbulentED,TurbulentOmega,
     *                      BTurbulentKE,BTurbulentED,BTurbulentOmega,
     *                      TurbulentKEOld,TurbulentEDOld,
     *                      TurbulentOmegaOld,BTurbulentKEOld,
     *                      BTurbulentEDOld,BTurbulentOmegaOld,
     *                      TurbulentKEOldOld,TurbulentEDOldOld,
     *                      TurbulentOmegaOldOld,BTurbulentKEOldOld,
     *                      BTurbulentEDOldOld,BTurbulentOmegaOldOld,
     *                      ModifiedED,BModifiedED,ModifiedEDOld,
     *                      BModifiedEDOld,ModifiedEDOldOld,
     *                      BModifiedEDOldOld,
     *                      TurbulentKL,BTurbulentKL,
     *                      TurbulentKLOld,BTurbulentKLOld,
     *                      TurbulentKLOldOld,BTurbulentKLOldOld,
     *                      Htotal,BHtotal,HtotalOld,BHtotalOld,
     *                      HtotalOldOld,BHtotalOldOld,
     *                      TGamma,BTGamma,TReTheta,BTReTheta,
     *                      TGammaOld,BTGammaOld,
     *                      TReThetaOld,BTReThetaOld,
     *                      TGammaOldOld,BTGammaOldOld,
     *                      TReThetaOldOld,BTReThetaOldOld,
     *                      TurbulentV2,BTurbulentV2,
     *                      TurbulentZeta,BTurbulentZeta,
     *                      TurbulentV2Old,BTurbulentV2Old,
     *                      TurbulentZetaOld,BTurbulentZetaOld,
     *                      TurbulentV2OldOld,BTurbulentV2OldOld,
     *                      TurbulentZetaOldOld,BTurbulentZetaOldOld
      use PhysicalProperties1,only: 
     *                Density,DensityOld,DensityOldOld,
     *                BDensity,BDensityOld,BDensityOldOld,
     *                SpecificHeat,SpecificHeatOld,SpecificHeatOldOld,
     *                BSpecificHeat,BSpecificHeatOld,
     *                BSpecificHeatOldOld,SpecificHeatScalar,
     *                SpecificHeatScalarOld,SpecificHeatScalarOldOld,
     *                BSpecificHeatScalar,BSpecificHeatScalarOld,
     *                BSpecificHeatScalarOldOld,SpecificHeatrField,
     *                SpecificHeatrFieldOld,SpecificHeatrFieldOldOld,
     *                BSpecificHeatrField,BSpecificHeatrFieldOld,
     *                BSpecificHeatrFieldOldOld
      use WindKessel1, only: OutletPressure,OutletPressureOld,
     *                       OutletPressureOldOld,MassFlowrateOld,
     *                       MassFlowrateOldOld
      use FlowInOut1, only: MassFlowRate
      use BoundaryConditions1, only: BoundaryType,outletTypeM
      use Turbulence1, only: S11,S12,S13,S22,S23,S33,
     *                       S11Old,S12Old,S13Old,S22Old,S23Old,S33Old,
     *                       S11OldOld,S12OldOld,S13OldOld,S22OldOld,
     *                       S23OldOld,S33OldOld
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,j,iScalar,irField
c
c********************************************************************************************
c
      dtOldOld=dtOld
      dtOld=dt
c
      DensityOldOld=DensityOld
      DensityOld=Density
      BDensityOldOld=BDensityOld
      BDensityOld=BDensity
c
      SpecificHeatOldOld=SpecificHeatOld
      SpecificHeatOld=SpecificHeat
      BSpecificHeatOldOld=BSpecificHeatOld
      BSpecificHeatOld=BSpecificHeat
c
      if(LsolveMomentum) then
c         
        uVelocityOldOld=uVelocityOld
        vVelocityOldOld=vVelocityOld
        wVelocityOldOld=wVelocityOld
        uVelocityOld=uVelocity
        vVelocityOld=vVelocity
        wVelocityOld=wVelocity
c
        BuVelocityOldOld=BuVelocityOld
        BvVelocityOldOld=BvVelocityOld
        BwVelocityOldOld=BwVelocityOld
        BuVelocityOld=BuVelocity
        BvVelocityOld=BvVelocity
        BwVelocityOld=BwVelocity
c          
      endif
c
      if(LsolveContinuity) then
c         
        PressureOldOld=PressureOld
        PressureOld=Pressure
        BPressureOldOld=BPressureOld
        BPressureOld=BPressure
        mdotOldOld=mdotOld
        mdotOld=mdot
        BmdotOldOld=BmdotOld
        BmdotOld=Bmdot
c          
      endif
c
      if(LsolveEnergy) then
c         
        TemperatureOldOld=TemperatureOld
        TemperatureOld=Temperature
        BTemperatureOldOld=BTemperatureOld
        BTemperatureOld=BTemperature
c          
      endif
c
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') then
c         
        HtotalOldOld=HtotalOld
        HtotalOld=Htotal
        BHtotalOldOld=BHtotalOld
        BHtotalOld=BHtotal
c          
      endif
c
      if(LsolveTurbulenceKineticEnergy) then
c         
        TurbulentKEOldOld=TurbulentKEOld
        TurbulentKEOld=TurbulentKE
        BTurbulentKEOldOld=BTurbulentKEOld
        BTurbulentKEOld=BTurbulentKE
c          
      endif
c
      if(LSolveTurbulenceDissipationRate) then
c         
        TurbulentEDOldOld=TurbulentEDOld
        TurbulentEDOld=TurbulentED
        BTurbulentEDOldOld=BTurbulentEDOld
        BTurbulentEDOld=BTurbulentED
c          
      endif
c
      if(LSolveTurbulenceSpecificDissipationRate) then
c         
        TurbulentOmegaOldOld=TurbulentOmegaOld
        TurbulentOmegaOld=TurbulentOmega
        BTurbulentOmegaOldOld=BTurbulentOmegaOld
        BTurbulentOmegaOld=BTurbulentOmega
c          
      endif
c
      if(LSolveTurbulenceGammaEquation) then
c         
        TGammaOldOld=TGammaOld
        TGammaOld=TGamma
        BTGammaOldOld=BTGammaOld
        BTGammaOld=BTGamma
c          
      endif
c
      if(LSolveTurbulenceReynoldsThetaEquation) then
c         
        TReThetaOldOld=TReThetaOld
        TReThetaOld=TReTheta
        BTReThetaOldOld=BTReThetaOld
        BTReThetaOld=BTReTheta
c          
      endif
c
      if(LSolveTurbulentKL) then
c         
        TurbulentKLOldOld=TurbulentKLOld
        TurbulentKLOld=TurbulentKL
        BTurbulentKLOldOld=BTurbulentKLOld
        BTurbulentKLOld=BTurbulentKL
c          
      endif
c
      if(LSolveModifiedED) then
c         
        ModifiedEDOldOld=ModifiedEDOld
        ModifiedEDOld=ModifiedED
        BModifiedEDOldOld=BModifiedEDOld
        BModifiedEDOld=BModifiedED
c          
      endif
c
      if(LSolveTurbulenceV2Equation) then
c         
        TurbulentV2OldOld=TurbulentV2Old
        TurbulentV2Old=TurbulentV2
        BTurbulentV2OldOld=BTurbulentV2Old
        BTurbulentV2Old=BTurbulentV2
c          
      endif
c
      if(LSolveTurbulenceZetaEquation) then
c         
        TurbulentZetaOldOld=TurbulentZetaOld
        TurbulentZetaOld=TurbulentZeta
        BTurbulentZetaOldOld=BTurbulentZetaOld
        BTurbulentZetaOld=BTurbulentZeta
c          
      endif
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.'spalartallmaras'.and.
     *               LSpalartAllmarasRotationCurvatureCorrection.and.
     *                  RotationCurvatureMethod.eq.'spalartshur') then
c
        S11OldOld=S11Old
        S11Old=S11
        S12OldOld=S12Old
        S12Old=S12
        S13OldOld=S13Old
        S13Old=S13
        S22OldOld=S22Old
        S22Old=S22
        S23OldOld=S23Old
        S23Old=S23
        S33OldOld=S33Old
        S33Old=S33
c          
      endif
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.'wrayagarwal'.and.
     *               LWrayAgarwalRotationCurvatureCorrection.and.
     *                  RotationCurvatureMethod.eq.'spalartshur') then
c
        S11OldOld=S11Old
        S11Old=S11
        S12OldOld=S12Old
        S12Old=S12
        S13OldOld=S13Old
        S13Old=S13
        S22OldOld=S22Old
        S22Old=S22
        S23OldOld=S23Old
        S23Old=S23
        S33OldOld=S33Old
        S33Old=S33
c          
      endif
c
      do irField=1,NumberOfrFieldsToSolve
c
        if(LSolverField(irField)) then
c
          do i=1,NumberOfElements
c      
            rFieldOldOld(i,irField)=rFieldOld(i,irField)
            rFieldOld(i,irField)=rField(i,irField)
            SpecificHeatrFieldOldOld(i,irField)=
     *                              SpecificHeatrFieldOld(i,irField)
            SpecificHeatrFieldOld(i,irField)=
     *                                 SpecificHeatrField(i,irField)
c
          enddo      
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BrFieldOldOld(i,j,irField)=BrFieldOld(i,j,irField)
              BrFieldOld(i,j,irField)=BrField(i,j,irField)
              BSpecificHeatrFieldOldOld(i,j,irField)=
     *                              BSpecificHeatrFieldOld(i,j,irField)
              BSpecificHeatrFieldOld(i,j,irField)=
     *                                 BSpecificHeatrField(i,j,irField)
c
            enddo
          enddo      
c      
        endif
c
      enddo
c
      do iScalar=1,NumberOfScalarsToSolve
c
        if(LSolveScalar(iScalar)) then
c
          do i=1,NumberOfElements
c      
            ScalarOldOld(i,iScalar)=ScalarOld(i,iScalar)
            ScalarOld(i,iScalar)=Scalar(i,iScalar)
            SpecificHeatScalarOldOld(i,iScalar)=
     *                              SpecificHeatScalarOld(i,iScalar)
            SpecificHeatScalarOld(i,iScalar)=
     *                                 SpecificHeatScalar(i,iScalar)
c
          enddo      
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BScalarOldOld(i,j,iScalar)=BScalarOld(i,j,iScalar)
              BScalarOld(i,j,iScalar)=BScalar(i,j,iScalar)
              BSpecificHeatScalarOldOld(i,j,iScalar)=
     *                              BSpecificHeatScalarOld(i,j,iScalar)
              BSpecificHeatScalarOld(i,j,iScalar)=
     *                                 BSpecificHeatScalar(i,j,iScalar)
c
            enddo
          enddo      
c      
        endif
c
      enddo
c
      do i=1,NumberOfBCSets
c      
        if(BoundaryType(i).eq.'outlet'.and.LWindKessel(i)) then
c
          if(outletTypeM(i).eq.'specifiedstaticpressure'.or.
     *              outletTypeM(i).eq.'specifiedresistance') then
c          
            MassFlowrateOldOld(i)=MassFlowrateOld(i)
            MassFlowrateOld(i)=MassFlowrate(i)
            OutletPressureOldOld(i)=OutletPressureOld(i)
            OutletPressureOld(i)=OutletPressure(i)
c
          endif      
c
        endif      
c      
      enddo
c
      return
      end