c
C#############################################################################################
c
      SUBROUTINE CrankNicolson
c
C#############################################################################################
c
      use User0, only: LSolveMomentum,LSolveContinuity,LSolveEnergy,
     *                 EnergyEquation,NumberOfScalarsToSolve,
     *                 LSolveScalar,LWindKessel,
     *                 LSolveTurbulenceKineticEnergy,
     *                 LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LSolveTurbulentKL,LSolveModifiedED,
     *                 LTurbulentFlow,TurbulenceModel,
     *                 LSpalartAllmarasRotationCurvatureCorrection,
     *                 LWrayAgarwalRotationCurvatureCorrection,
     *                 RotationCurvatureMethod,
     *                 LSolveTurbulenceGammaEquation,
     *                 LSolveTurbulenceReynoldsThetaEquation
      use Variables1, only: uVelocity,vVelocity,wVelocity,Pressure,
     *                      Temperature,Htotal,BHtotal,
     *                      HtotalOld,BHtotalOld,HtotalOldOld,
     *                      BHtotalOldOld,
     *                      uVelocityOld,vVelocityOld,wVelocityOld,
     *                      PressureOld,TemperatureOld,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      BPressure,BTemperature,BuVelocityOld,
     *                      BvVelocityOld,BwVelocityOld,BPressureOld,
     *                      BTemperatureOld,TurbulentKE,BTurbulentKE,
     *                      TurbulentKEOld,BTurbulentKEOld,
     *                      TurbulentED,BTurbulentED,
     *                      TurbulentEDOld,BTurbulentEDOld,
     *                      TurbulentOmega,BTurbulentOmega,
     *                      TurbulentOmegaOld,BTurbulentOmegaOld,
     *                      TurbulentKL,BTurbulentKL,
     *                      TurbulentKLOld,BTurbulentKLOld,
     *                      ModifiedED,BModifiedED,
     *                      ModifiedEDOld,BModifiedEDOld,
     *                      TGamma,BTGamma,TReTheta,BTReTheta,
     *                      TGammaOld,BTGammaOld,
     *                      TReThetaOld,BTReThetaOld
      use Scalar1, only: Scalar,ScalarOld,BScalar,BScalarOld
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use WindKessel1, only: OutletPressure,OutletPressureOld,
     *                       MassFlowrateOld
      use FlowInOut1, only: MassFlowRate
      use BoundaryConditions1, only: BoundaryType,outletTypeM
      use Turbulence1, only: S11,S12,S13,S22,S23,S33,
     *                       S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,j,iScalar
c********************************************************************************************
c
      if(LSolveMomentum) then
c
        uVelocity=2.*uVelocity-uVelocityOld
        vVelocity=2.*vVelocity-vVelocityOld
        wVelocity=2.*wVelocity-wVelocityOld
        BuVelocity=2.*BuVelocity-BuVelocityOld
        BvVelocity=2.*BvVelocity-BvVelocityOld
        BwVelocity=2.*BwVelocity-BwVelocityOld
c
      endif  
c
      if(LSolveContinuity) then
c
        Pressure=2.*Pressure-PressureOld
        BPressure=2.*BPressure-BPressureOld
c
      endif  
c
      if(LSolveEnergy) then
c
        Temperature=2.*Temperature-TemperatureOld
        BTemperature=2.*BTemperature-BTemperatureOld
c
      endif  
c
      if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
c
        Htotal=2.*Htotal-HtotalOld
        BHtotal=2.*BHtotal-BHtotalOld
c
      endif  
c
      if(LSolveTurbulenceKineticEnergy) then
c
        TurbulentKE=2.*TurbulentKE-TurbulentKEOld
        BTurbulentKE=2.*BTurbulentKE-BTurbulentKEOld
c
      endif  
c
      if(LSolveTurbulenceDissipationRate) then
c
        TurbulentED=2.*TurbulentED-TurbulentEDOld
        BTurbulentED=2.*BTurbulentED-BTurbulentEDOld
c
      endif  
c
      if(LSolveTurbulenceSpecificDissipationRate) then
c
        TurbulentOmega=2.*TurbulentOmega-TurbulentOmegaOld
        BTurbulentOmega=2.*BTurbulentOmega-BTurbulentOmegaOld
c
      endif  
c
      if(LSolveTurbulenceGammaEquation) then
c
        TGamma=2.*TGamma-TGammaOld
        BTGamma=2.*BTGamma-BTGammaOld
c
      endif  
c
      if(LSolveTurbulenceReynoldsThetaEquation) then
c
        TReTheta=2.*TReTheta-TReThetaOld
        BTReTheta=2.*BTReTheta-BTReThetaOld
c
      endif  
c
      if(LSolveTurbulentKL) then
c
        TurbulentKL=2.*TurbulentKL-TurbulentKLOld
        BTurbulentKL=2.*BTurbulentED-BTurbulentKLOld
c
      endif  
c
      if(LSolveModifiedED) then
c
        ModifiedED=2.*ModifiedED-ModifiedEDOld
        BModifiedED=2.*BModifiedED-BModifiedEDOld
c
      endif  
c
      do iScalar=1,NumberOfScalarsToSolve
c
        if(LSolveScalar(iScalar)) then
c        
          do i=1,NumberOfElements
            Scalar(i,iScalar)=2.*Scalar(i,iScalar)-ScalarOld(i,iScalar)
          enddo       
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BScalar(i,j,iScalar)=2.*BScalar(i,j,iScalar)-
     *                                      BScalarOld(i,j,iScalar)
            enddo
          enddo        
c        
        endif
c
      enddo
c
c--- Recalculate the mass flow rate field based on the new fields
c
      if(LSolveContinuity) call AssembleMdot
c
      do i=1,NumberOfBCSets
c      
        if(BoundaryType(i).eq.'outlet'.and.LWindKessel(i)) then
c
          if(outletTypeM(i).eq.'specifiedstaticpressure'.or.
     *              outletTypeM(i).eq.'specifiedresistance') then
c          
            MassFlowrate(i)=2.*MassFlowrate(i)-MassFlowrateOld(i)
            OutletPressure(i)=2.*OutletPressure(i)-OutletPressureOld(i)
c
          endif      
c
        endif      
c      
      enddo
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.'spalartallmaras'
     *          .and.LSpalartAllmarasRotationCurvatureCorrection.and.
     *            RotationCurvatureMethod.eq.'spalartshur') then


        S11=2.*S11-S11Old
        S12=2.*S12-S12Old
        S13=2.*S13-S13Old
        S22=2.*S22-S22Old
        S23=2.*S23-S23Old
        S33=2.*S33-S33Old
c
      endif
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.'wrayagarwal'
     *          .and.LWrayAgarwalRotationCurvatureCorrection.and.
     *            RotationCurvatureMethod.eq.'spalartshur') then


        S11=2.*S11-S11Old
        S12=2.*S12-S12Old
        S13=2.*S13-S13Old
        S22=2.*S22-S22Old
        S23=2.*S23-S23Old
        S33=2.*S33-S33Old
c
      endif
c
      return
      end