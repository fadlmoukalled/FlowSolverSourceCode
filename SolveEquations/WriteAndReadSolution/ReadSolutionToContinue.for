c
C#############################################################################################
      SUBROUTINE readSolution
C#############################################################################################
c
      use Variables1
      use Variables2
      use Variables3
      use Variables4
      use Scalar1
      use TransferVariable1
      use VolumeOfFluid1
      use TransferrField1
      use Transient1
      use PhysicalProperties1
      use BoundaryConditions1
      use BoundaryFluxes
      use MultiGrid1
      use Residuals1
      use Turbulence1
      use FlowInOut1
      use WindKessel1
      use User0
c
c*********************************************************************************************
      implicit none
c*********************************************************************************************
c
      rewind 4
      if(LSolveMomentum) then
c
        read(4) uVelocity,vVelocity,wVelocity
        read(4) uVelGradx,uVelGrady,uVelGradz
        read(4) vVelGradx,vVelGrady,vVelGradz
        read(4) wVelGradx,wVelGrady,wVelGradz
        read(4) BuVelocity,BvVelocity,BwVelocity
        read(4) BuVelGradx,BuVelGrady,BuVelGradz
        read(4) BvVelGradx,BvVelGrady,BvVelGradz
        read(4) BwVelGradx,BwVelGrady,BwVelGradz
        read(4) uVelGradfx,uVelGradfy,uVelGradfz
        read(4) vVelGradfx,vVelGradfy,vVelGradfz
        read(4) wVelGradfx,wVelGradfy,wVelGradfz
        read(4) xVeldirection,yVeldirection,zVeldirection
        read(4) Viscosity
        read(4) BViscosity
        read(4) ScMomentumx,SbMomentumx
        read(4) ScMomentumy,SbMomentumy
        read(4) ScMomentumz,SbMomentumz
        read(4) mdot,Bmdot,effdiv
c
        read(4) Du1Velocity,Du2Velocity,uVelocityStar
        read(4) Dv1Velocity,Dv2Velocity,vVelocityStar
        read(4) Dw1Velocity,Dw2Velocity,wVelocityStar
c          
        read(4) MachNumber,BMachNumber
c
        if(LUnsteady) then
c        
          read(4) uVelocityOld,uVelocityOldOld
          read(4) BuVelocityOld,BuVelocityOldOld
          read(4) vVelocityOld,vVelocityOldOld
          read(4) BvVelocityOld,BvVelocityOldOld
          read(4) wVelocityOld,wVelocityOldOld
          read(4) BwVelocityOld,BwVelocityOldOld
          read(4) mdotOld,mdotOldOld
          read(4) BmdotOld,BmdotOldOld
c        
        endif
c
        if(Lcompressible) then
c
          read(4) drhodp
          read(4) Bdrhodp
c
        endif
c        
        if(LConvectScalar.and..not.LsolveMomentum) then
c        
          read(4) uVelocity,vVelocity,wVelocity
          read(4) BuVelocity,BvVelocity,BwVelocity
          read(4) mdot
          read(4) Bmdot
          read(4) effdiv

c
        endif
c
      endif
c
      if(LSolveContinuity) then
c
        read(4) Pressure
        read(4) BPressure
        read(4) PressGradx,PressGrady,PressGradz
        read(4) BPressGradx,BPressGrady,BPressGradz
        read(4) PressGradfx,PressGradfy,PressGradfz
c
        if(.not.LsolveMomentum) then
c
          read(4) drhodP
          read(4) BdrhodP
c
        endif
c
        if(LUnsteady) then
c
          read(4) PressureOld,PressureOldOld
          read(4) BPressureOld,BPressureOldOld
          read(4) dpdt
c
          if(.not.LsolveMomentum) read(4) mdotOld
          if(.not.LsolveMomentum) read(4) mdotOldOld
c
        endif
c
      endif
c
      if(LSolveEnergy) then
c
        read(4) Temperature,TempGradx,TempGrady,TempGradz
        read(4) BTemperature,BTempGradx,BTempGrady,BTempGradz
        read(4) TempGradfx,TempGradfy,TempGradfz
c
        if(EnergyEquation.eq.'Htotal') then
c
          read(4) Htotal,HtotalGradx,HtotalGrady,HtotalGradz
          read(4) BHtotal,BHtotalGradx,BHtotalGrady,BHtotalGradz
          read(4) HtotalGradfx,HtotalGradfy,HtotalGradfz
c
        endif
c
        read(4) Conductivity,SpecificHeat
        read(4) BConductivity,BSpecificHeat
        read(4) SHeatGradx,SHeatGrady,SHeatGradz
        read(4) BSHeatGradx,BSHeatGrady,BSHeatGradz
        read(4) ScEnergy,SbEnergy
c
        if(LUnsteady) then
c
          read(4) TemperatureOld,TemperatureOldOld
          read(4) BTemperatureOld,BTemperatureOldOld
          read(4) SpecificHeatOld,SpecificHeatOldOld
          read(4) BSpecificHeatOld,BSpecificHeatOldOld
c
          if(EnergyEquation.eq.'Htotal') then
c
            read(4) HtotalOld,HtotalOldOld
            read(4) BHtotalOld,BHtotalOldOld
c
          endif
c
        endif
c
      endif
c
      if(LSolveMomentum.and.LBuoyancy) then
c
        read(4) Buoyancyx,Buoyancyy,Buoyancyz
        read(4) BBuoyancyx,BBuoyancyy,BBuoyancyz
        read(4) Buoyancyfx,Buoyancyfy,Buoyancyfz
c
      endif
c
      if(Lcompressible.and.LsolveMomentum.and..not.LsolveEnergy) then
c      
        read(4) Temperature
        read(4) BTemperature
        read(4) SpecificHeat
        read(4) BSpecificHeat
c      
      endif
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        read(4) LSoutherLandScalar
        read(4) EquationOfStateScalar
        read(4) RGasScalar,GammaGasScalar,PrLaminarScalar
c
        read(4) Scalar,ScalarGradx,ScalarGrady,ScalarGradz
        read(4) BScalar,BScalarGradx,BScalarGrady,BScalarGradz
        read(4) ScalarGradfx,ScalarGradfy,ScalarGradfz
        read(4) ScScalar,SbScalar
        read(4) DiffusionCoefficient,SpecificHeatScalar
        read(4) BDiffusionCoefficient,BSpecificHeatScalar
c
        if(LUnsteady) then
c
          read(4) ScalarOld,ScalarOldOld
          read(4) BScalarOld,BScalarOldOld
          read(4) SpecificHeatScalarOld,SpecificHeatScalarOldOld
          read(4) BSpecificHeatScalarOld,BSpecificHeatScalarOldOld
c
        endif
c
        if(NumberOfPointSources.ne.0) then
c
          if(.not.LSolveEnergy) then
c
            read(4) iElementPointSource,xLocationOfPointSource,
     *                  yLocationOfPointSource,zLocationOfPointSource
c
          endif
c
          read(4) ScPointSourceScalar,SbPointSourceScalar
c
        endif
c
      endif
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        read(4) LSoutherLandrField
        read(4) EquationOfStaterField
        read(4) RGasrField,GammaGasrField,PrLaminarrField
c
        read(4) rField,rFieldGradx,rFieldGrady,rFieldGradz
        read(4) BrField,BrFieldGradx,BrFieldGrady,BrFieldGradz
        read(4) rFieldGradfx,rFieldGradfy,rFieldGradfz
        read(4) ScrField,SbrField
        read(4) cosThetaF,cosThetaE
        read(4) BcosThetaF
        read(4) ConductivityrField,SpecificHeatrField
        read(4) BConductivityrField,BSpecificHeatrField
        read(4) DensityrField,ViscosityrField
        read(4) BDensityrField,BViscosityrField
c
        if(LUnsteady) then
c
          read(4) rFieldOld,rFieldOldOld
          read(4) BrFieldOld,BrFieldOldOld
          read(4) SpecificHeatrFieldOld,SpecificHeatrFieldOldOld
          read(4) BSpecificHeatrFieldOld,BSpecificHeatrFieldOldOld
c
        endif
c
        if(LSurfaceTension) then
c
          read(4) Curvature
          read(4) BCurvature
c
          read(4) delrFieldMagnitude
          read(4) BdelrFieldMagnitude
          read(4) delrFieldMagnitudeGradx
          read(4) BdelrFieldMagnitudeGradx
          read(4) delrFieldMagnitudeGrady
          read(4) BdelrFieldMagnitudeGrady
          read(4) delrFieldMagnitudeGradz
          read(4) BdelrFieldMagnitudeGradz
c
        endif
c
      endif
c
      if(LSolveLambdaELEEquation) then
        read(4) LambdaELE,BLambdaELE
      endif
c
c--- Density is always needed
c 
      read(4) Density,DensGradx,DensGrady,DensGradz    
      read(4) BDensity,BDensGradx,BDensGrady,BDensGradz    
      read(4) Densityf      
c
      if(LUnsteady) then
c
        read(4) DensityOld,DensityOldOld  
        read(4) BDensityOld,BDensityOldOld  
c
      endif
c
      if(LFalseTransientMomentum) then
c
        read(4) DensityStar
        read(4) BDensityStar
c
      endif
c
      if(NumberOfPointSources.ne.0.and.LSolveEnergy) then
c
        read(4) iElementPointSource
        read(4) xLocationOfPointSource,yLocationOfPointSource,
     *           zLocationOfPointSource
        read(4) ScPointSourceEnergy,SbPointSourceEnergy
c
      endif
c
      if(LanisotropicDiffusion) then
c
        read(4) Conductivity11,Conductivity12,Conductivity13
        read(4) Conductivity22,Conductivity23,Conductivity33
        read(4) BConductivity11,BConductivity12,BConductivity13
        read(4) BConductivity22,BConductivity23,BConductivity33
c
        if(NumberOfScalarsToSolve.gt.0) then
c
          read(4) DiffusionCoefficient11,DiffusionCoefficient12,
     *                                       DiffusionCoefficient13
          read(4) DiffusionCoefficient22,DiffusionCoefficient23,
     *                                       DiffusionCoefficient33
          read(4) BDiffusionCoefficient11,BDiffusionCoefficient12,
     *                                       BDiffusionCoefficient13
          read(4) BDiffusionCoefficient22,BDiffusionCoefficient23,
     *                                       BDiffusionCoefficient33
c
        endif
c
      endif
c
      read(4) eDiffCoefficient
      read(4) BeDiffCoefficient
c
      if(LTurbulentFlow.and.LsolveMomentum) then
c
        read(4) rhok,drhokdx,drhokdy,drhokdz
        read(4) Brhok,Bdrhokdx,Bdrhokdy,Bdrhokdz
        read(4) rhoTED,ReT
        read(4) BrhoTED,BReT
c
        read(4) S11,S12,S13,S22,S23,S33,StrainRate
        read(4) BS11,BS12,BS13,BS22,BS23,BS33,BStrainRate
        read(4) W11,W12,W13,W22,W23,W33,Vorticity
        read(4) BW11,BW12,BW13,BW22,BW23,BW33,BVorticity
        read(4) Tau11,Tau12,Tau13,Tau22,Tau23,Tau33
        read(4) TauWall,ustar,ystar,uplus,yplus,uTau
        read(4) duplusdyplus,WallViscosity,KsPlus
c
        if(LCompressible) read(4) XiStarFmt
c
        if(TurbulenceModel.eq.'kepsilonv2f') then
c
          read(4) Ce1Coefficient,LScale,TScale
          read(4) BLScale,BTScale
c
        endif
c
        if(TurbulenceModel.eq.'kepsilonzetaf') then
c
          read(4) Ce1Coefficient,LScale,TScale
          read(4) BLScale,BTScale
c
          read(4) TurbulentZetaGrad2x,TurbulentZetaGrad2y
          read(4) TurbulentZetaGrad2z,TurbulentZetaGradxy
          read(4) TurbulentZetaGradxz
          read(4) BTurbulentZetaGrad2x,BTurbulentZetaGrad2y
          read(4) BTurbulentZetaGrad2z,BTurbulentZetaGradxy
          read(4) BTurbulentZetaGradxz
c
        endif
c
        if(TurbulenceModel.eq.'kepsilonrng') then
c
          read(4) C2eRNG,sigTKERNG,sigTEDRNG,BsigTKERNG,BsigTEDRNG
c
        endif
c
        if(TurbulenceModel.eq.'nut92') then
c
          read(4) TurbulenceProduction,BTurbulenceProduction
          read(4) TurbulentKE,BTurbulentKE
          read(4) ModifiedDensGradx,ModifiedDensGrady,ModifiedDensGradz
          read(4) BModifiedDensGradx,BModifiedDensGrady
          read(4) BModifiedDensGradz
          read(4) ModifiedDensGrad2x,ModifiedDensGrad2y
          read(4) ModifiedDensGrad2z
          read(4) BModifiedDensGrad2x,BModifiedDensGrad2y
          read(4) BModifiedDensGrad2z
c
          read(4) uVelGrad2x,BuVelGrad2x,uVelGrad2y
          read(4) BuVelGrad2y,uVelGrad2z,BuVelGrad2z
          read(4) vVelGrad2x,BvVelGrad2x,vVelGrad2y
          read(4) BvVelGrad2y,vVelGrad2z,BvVelGrad2z
          read(4) wVelGrad2x,BwVelGrad2x,wVelGrad2y
          read(4) BwVelGrad2y,wVelGrad2z,BwVelGrad2z
          read(4) uvVelGradxy,BuvVelGradxy
          read(4) uwVelGradxz,BuwVelGradxz
c
          read(4) ModifiedEDGrad2x,ModifiedEDGrad2y
          read(4) ModifiedEDGrad2z,BModifiedEDGrad2x
          read(4) BModifiedEDGrad2y,BModifiedEDGrad2z
          read(4) ModifiedMut,BModifiedMut
          read(4) ModifiedMutGradx,ModifiedMutGrady
          read(4) ModifiedMutGradz,BModifiedMutGradx
          read(4) BModifiedMutGrady,BModifiedMutGradz
          read(4) G1Nut92,G2Nut92,F1Nut92,F2Nut92
          read(4) N1Nut92,N2Nut92,BN1Nut92,N1Nut92Gradx
          read(4) N1Nut92Grady,N1Nut92Gradz,BN1Nut92Gradx
          read(4) BN1Nut92Grady,BN1Nut92Gradz
c
        endif
c
        if(TurbulenceModel.eq.'sstgama') then
c
          read(4) F1factor,BF1factor,fr1Coefficient,F2factor,BF2factor
          read(4) F3factor,BF3factor,F4factor,BF4factor
c
          read(4) NormalVelocity,BNormalVelocity
          read(4) NormalVelocityGradx,NormalVelocityGrady
          read(4) NormalVelocityGradz,BNormalVelocityGradx
          read(4) BNormalVelocityGrady,BNormalVelocityGradz
c
        endif
c
        if(TurbulenceModel.eq.'kepsilonrt') then
c
          read(4) fr1Coefficient,f2RT,Bf2RT,velRT,BvelRT
          read(4) velRTGradx,BvelRTGradx,velRTGrady,BvelRTGrady
          read(4) velRTGradz,BvelRTGradz
          read(4) fmuKE,BfmuKE,fmuRT,BfmuRT
          read(4) TurbulentViscosityKE,BTurbulentViscosityKE
          read(4) TurbulentViscosityRT,BTurbulentViscosityRT
c
          if(LKEpsilonRtRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
            read(4) DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
            if(LUnsteady) then
              read(4)S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
              read(4)S11OldOld,S12OldOld,S13OldOld,
     *                  S22OldOld,S23OldOld,S33OldOld
            endif            
          endif
        endif
c
        if(TurbulenceModel.eq.'kklomega') then
c
          read(4) TurbulentViscosity
          read(4) BTurbulentViscosity
c
          read(4) TurbulentViscosityTs,TurbulentViscosityTl
          read(4) BTurbulentViscosityTs,BTurbulentViscosityTl
          read(4) AlfaT,BAlfaT
          read(4) LambdaT,BLambdaT
          read(4) LambdaEff,BLambdaEff
          read(4) coefficientFW,BcoefficientFW
          read(4) TurbulentKTs,BTurbulentKTs
          read(4) ProductionKT,ProductionKL
c
          if(LSolveEnergy) then
c
            read(4) AlfaTheta,BAlfaTheta
c
          endif
c
        endif
c
        if(TurbulenceModel.eq.'realizable') then
c
          read(4) cmuR,c1R
          read(4) BcmuR,Bc1R
c
        endif
c
        if(TurbulenceModel.eq.'spalartallmaras'.or.
     *                    TurbulenceModel.eq.'wrayagarwal') then
c
          read(4) TurbulenceProduction
          read(4) BTurbulenceProduction
          read(4) TurbulentKE,BTurbulentKE
c
          if(TurbulenceModel.eq.'spalartallmaras') then
            read(4) TGamma,BTGamma
            read(4) Stelda,fwCoefficient,ft2Coefficient,fr1Coefficient
            read(4) fv1Coefficient
            read(4) Bfv1Coefficient
            if(LNegativeSpalartAllmaras) then
              read(4) fnCoefficient
              read(4) BfnCoefficient
            endif
            if(LSpalartAllmarasRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
              read(4) DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
              if(LUnsteady) then
                read(4)S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
                read(4)S11OldOld,S12OldOld,S13OldOld,
     *                  S22OldOld,S23OldOld,S33OldOld
              endif            
            endif
          elseif(TurbulenceModel.eq.'wrayagarwal') then
            read(4) f1WA,fmuWA,SRateGradx,SRateGrady,SRateGradz
            read(4) Bf1WA,BfmuWA,BSRateGradx,BSRateGrady,BSRateGradz
            read(4) fr1Coefficient
            if(LWrayAgarwalRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
              read(4)DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
              if(LUnsteady) then
                read(4)S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
                read(4)S11OldOld,S12OldOld,S13OldOld,
     *                  S22OldOld,S23OldOld,S33OldOld
              endif            
            endif
          endif
c
        endif
c
        if(TurbulenceModel.eq.'kepsilonchien'.or.
     *       TurbulenceModel.eq.'kepsilonchc'.or.
     *        TurbulenceModel.eq.'kepsilonkasagi'.or.
     *         TurbulenceModel.eq.'kepsilontagawa'.or.
     *          TurbulenceModel.eq.'kepsilonhishida'.or.
     *           TurbulenceModel.eq.'kepsilonsharma'.or.
     *            TurbulenceModel.eq.'kelambremhorst'.or.
     *             TurbulenceModel.eq.'kelambremhorstm') then
c
          read(4) fmuCoefficient,f1Coefficient,f2Coefficient,LTKE,LTED
          read(4) BfmuCoefficient
c
          if(TurbulenceModel.eq.'kepsilonsharma'.or.
     *          TurbulenceModel.eq.'kepsilonhishida') then
c
            read(4) uVelGrad2x,BuVelGrad2x,uVelGrad2y
            read(4) BuVelGrad2y,uVelGrad2z,BuVelGrad2z
            read(4) vVelGrad2x,BvVelGrad2x,vVelGrad2y
            read(4) BvVelGrad2y,vVelGrad2z,BvVelGrad2z
            read(4) wVelGrad2x,BwVelGrad2x,wVelGrad2y
            read(4) BwVelGrad2y,wVelGrad2z,BwVelGrad2z
            read(4) uvVelGradxy,BuvVelGradxy
            read(4) uwVelGradxz,BuwVelGradxz
c
            read(4) sqrtTurbulentKE
            read(4) BsqrtTurbulentKE
            read(4) sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz
            read(4) BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz
            read(4) SRateGradx,SRateGrady,SRateGradz
            read(4) BSRateGradx,BSRateGrady,BSRateGradz
c
          endif
c
        endif
c
        if(TurbulenceModel.eq.'komega2006lrn') then
c
          read(4) alfa,alfaStar,bettaStar
          read(4) Balfa,BalfaStar,BbettaStar
c
        endif
c
        if(TurbulenceModel.eq.'komegabsl'.or.
     *                    TurbulenceModel.eq.'komegasst') then
c
          read(4) F1factor
          read(4) BF1factor
c
          if(TurbulenceModel.eq.'komegasst') then
c
            read(4) F2factor,F3factor,F4factor,fr1Coefficient
            read(4) BF2factor,BF3factor,BF4factor
c
            if(LKOmegaSSTRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
              read(4) DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
              if(LUnsteady) then
                read(4) S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
                read(4) S11OldOld,S12OldOld,S13OldOld,
     *                   S22OldOld,S23OldOld,S33OldOld
              endif            
            endif
c
          endif
c
        endif
c
        if(TurbulenceModel.eq.'sstgamaretheta') then
c
          read(4) TGammaEff,F1factor,fr1Coefficient,F2factor,
     *             F3factor,F4factor
          read(4) BF1factor,BF2factor,BF3factor,BF4factor
c
        endif
c
        read(4) TurbulentViscosity
        read(4) BTurbulentViscosity
c
        if(TurbulenceModel.eq.'komega2006'.or.
     *             TurbulenceModel.eq.'komega2006lrn') then
c
          read(4) TurbulentViscosity1
          read(4) BTurbulentViscosity1
c
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          read(4) TurbulentKE,TKEGradx,TKEGrady,TKEGradz
          read(4) BTurbulentKE,BTKEGradx,BTKEGrady,BTKEGradz
          read(4) TKEGradfx,TKEGradfy,TKEGradfz
          read(4) ScTKE,SbTKE
          read(4) TurbulenceProduction
          read(4) BTurbulenceProduction
c
          if(LBuoyancy) then
            read(4) TurbulenceProductionB
            read(4) BTurbulenceProductionB
          endif
c
          if(LUnsteady) then
            read(4) TurbulentKEOld,TurbulentKEOldOld
            read(4) BTurbulentKEOld,BTurbulentKEOldOld
          endif
c
        endif
c
        if(LSolveTurbulenceDissipationRate) then
c
          read(4) TurbulentED,TEDGradx,TEDGrady,TEDGradz
          read(4) BTurbulentED,BTEDGradx,BTEDGrady,BTEDGradz
          read(4) TEDGradfx,TEDGradfy,TEDGradfz
          read(4) ScTED,SbTED
c
          if(LUnsteady) then
c
            read(4) TurbulentEDOld,TurbulentEDOldOld
            read(4) BTurbulentEDOld,BTurbulentEDOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceSpecificDissipationRate) then
c
          read(4) TurbulentOmega,TOmegaGradx,TOmegaGrady,TOmegaGradZ
          read(4) BTurbulentOmega
          read(4) BTOmegaGradx,BTOmegaGrady,BTOmegaGradZ
          read(4) TOmegaGradfx,TOmegaGradfy,TOmegaGradfy
          read(4) ScTOmega,SbTOmega
c
          if(LUnsteady) then
c
            read(4) TurbulentOmegaOld,TurbulentOmegaOldOld
            read(4) BTurbulentOmegaOld,BTurbulentOmegaOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceGammaEquation) then
c
          read(4) TGamma,TGammaGradx,TGammaGrady,TGammaGradz
          read(4) BTGamma,BTGammaGradx,BTGammaGrady,BTGammaGradz
          read(4) TGammaGradfx,TGammaGradfy,TGammaGradfz
          read(4) ScTGamma,SbTGamma
c
          if(LUnsteady) then
c
            read(4) TGammaOld,TGammaOldOld
            read(4) BTGammaOld,BTGammaOldOld
c
          endif
c
        endif
c
       if(LSolveTurbulenceReynoldsThetaEquation) then
c
          read(4) TReTheta,TReThetaGradx,TReThetaGrady,TReThetaGradz
          read(4) BTReTheta,BTReThetaGradx,
     *                     BTReThetaGrady,BTReThetaGradz
          read(4) TReThetaGradfx,TReThetaGradfy,TReThetaGradfz
          read(4) ScTReTheta,SbTReTheta
c
          if(LUnsteady) then
c
            read(4) TReThetaOld,TReThetaOldOld
            read(4) BTReThetaOld,BTReThetaOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulentKL) then
c
          read(4) TurbulentKL,TurbulentKLGradx,
     *              TurbulentKLGrady,TurbulentKLGradz
          read(4) BTurbulentKL,BTurbulentKLGradx,
     *              BTurbulentKLGrady,BTurbulentKLGradz
          read(4) TurbulentKLGradfx,TurbulentKLGradfy,TurbulentKLGradfz
          read(4) ScTurbulentKL,SbTurbulentKL
          read(4) uVelGrad2x,uVelGrad2y,uVelGrad2z         
          read(4) vVelGrad2x,vVelGrad2y,vVelGrad2z         
          read(4) wVelGrad2x,wVelGrad2y,wVelGrad2z         
          read(4) BuVelGrad2x,BuVelGrad2y,BuVelGrad2z         
          read(4) BvVelGrad2x,BvVelGrad2y,BvVelGrad2z         
          read(4) BwVelGrad2x,BwVelGrad2y,BwVelGrad2z         
          read(4) uvVelGradxy,uwVelGradxz
          read(4) BuvVelGradxy,BuwVelGradxz
c
          if(LUnsteady) then
c
            read(4) TurbulentKLOld,TurbulentKLOldOld
            read(4) BTurbulentKLOld,BTurbulentKLOldOld
c
          endif
c
        endif
c
        if(LSolveModifiedED) then
c
          read(4) ModifiedED,BModifiedED
          read(4) ModifiedEDGradx,BModifiedEDGradx
          read(4) ModifiedEDGrady,BModifiedEDGrady
          read(4) ModifiedEDGradz,BModifiedEDGradz
          read(4) ModifiedEDGradfx
          read(4) ModifiedEDGradfy
          read(4) ModifiedEDGradfz
          read(4) ScModifiedED,SbModifiedED
c
          if(LUnsteady) then
c
            read(4) ModifiedEDOld,BModifiedEDOld
            read(4) ModifiedEDOldOld,BModifiedEDOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceV2Equation) then
c
          read(4) TurbulentV2,BTurbulentV2,TurbulentV2Gradx
          read(4) BTurbulentV2Gradx,TurbulentV2Grady,BTurbulentV2Grady
          read(4) TurbulentV2Gradz,BTurbulentV2Gradz,TurbulentV2Gradfx
          read(4) TurbulentV2Gradfy,TurbulentV2Gradfz
          read(4) ScTurbulentV2,SbTurbulentV2
c
          if(LUnsteady) then
c
            read(4) TurbulentV2Old,BTurbulentV2Old
            read(4) TurbulentV2OldOld,BTurbulentV2OldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceZetaEquation) then
c
          read(4) TurbulentZeta,BTurbulentZeta,TurbulentZetaGradx
          read(4) BTurbulentZetaGradx,TurbulentZetaGrady
          read(4) BTurbulentZetaGrady,TurbulentZetaGradz
          read(4) BTurbulentZetaGradz,TurbulentZetaGradfx
          read(4) TurbulentZetaGradfy,TurbulentZetaGradfz
          read(4) ScTurbulentZeta,SbTurbulentZeta
c
          if(LUnsteady) then
c
            read(4) TurbulentZetaOld,BTurbulentZetaOld
            read(4) TurbulentZetaOldOld,BTurbulentZetaOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulencefRelaxationEquation) then
c
          read(4) TfRelaxation,BTfRelaxation,TfRelaxationGradx
          read(4) BTfRelaxationGradx,TfRelaxationGrady
          read(4) BTfRelaxationGrady,TfRelaxationGradz
          read(4) BTfRelaxationGradz,TfRelaxationGradfx
          read(4) TfRelaxationGradfy,TfRelaxationGradfz
          read(4) ScTfRelaxation,SbTfRelaxation
c
          if(LUnsteady) then
c
            read(4) TfRelaxationOld,BTfRelaxationOld
            read(4) TfRelaxationOldOld,BTfRelaxationOldOld
c
          endif
c
        endif
c
        if(NumberOfScalarsToSolve.ne.0) then
c
          read(4) SigScalar
c
        endif
c
      endif
c
      read(4) MassFlowRate,MassFlowRateOld,MassFlowRateOldOld
      read(4) OutletPressure,OutletPressureOld,OutletPressureOldOld
c      
      if(LUnsteady) then
        read(4) dt,dtOld,dtOldOld
        read(4) time,ndt
      endif
c
      return
      end