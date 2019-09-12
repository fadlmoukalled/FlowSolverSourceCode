c
C#############################################################################################
      SUBROUTINE writeSolution
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
        write(4) uVelocity,vVelocity,wVelocity
        write(4) uVelGradx,uVelGrady,uVelGradz
        write(4) vVelGradx,vVelGrady,vVelGradz
        write(4) wVelGradx,wVelGrady,wVelGradz
        write(4) BuVelocity,BvVelocity,BwVelocity
        write(4) BuVelGradx,BuVelGrady,BuVelGradz
        write(4) BvVelGradx,BvVelGrady,BvVelGradz
        write(4) BwVelGradx,BwVelGrady,BwVelGradz
        write(4) uVelGradfx,uVelGradfy,uVelGradfz
        write(4) vVelGradfx,vVelGradfy,vVelGradfz
        write(4) wVelGradfx,wVelGradfy,wVelGradfz
        write(4) xVeldirection,yVeldirection,zVeldirection
        write(4) Viscosity
        write(4) BViscosity
        write(4) ScMomentumx,SbMomentumx
        write(4) ScMomentumy,SbMomentumy
        write(4) ScMomentumz,SbMomentumz
        write(4) mdot,Bmdot,effdiv
c
        write(4) Du1Velocity,Du2Velocity,uVelocityStar
        write(4) Dv1Velocity,Dv2Velocity,vVelocityStar
        write(4) Dw1Velocity,Dw2Velocity,wVelocityStar
c          
        write(4) MachNumber,BMachNumber
c
        if(LUnsteady) then
c        
          write(4) uVelocityOld,uVelocityOldOld
          write(4) BuVelocityOld,BuVelocityOldOld
          write(4) vVelocityOld,vVelocityOldOld
          write(4) BvVelocityOld,BvVelocityOldOld
          write(4) wVelocityOld,wVelocityOldOld
          write(4) BwVelocityOld,BwVelocityOldOld
          write(4) mdotOld,mdotOldOld
          write(4) BmdotOld,BmdotOldOld
c        
        endif
c
        if(Lcompressible) then
c
          write(4) drhodp
          write(4) Bdrhodp
c
        endif
c        
        if(LConvectScalar.and..not.LsolveMomentum) then
c        
          write(4) uVelocity,vVelocity,wVelocity
          write(4) BuVelocity,BvVelocity,BwVelocity
          write(4) mdot
          write(4) Bmdot
          write(4) effdiv

c
        endif
c
      endif
c
      if(LSolveContinuity) then
c
        write(4) Pressure
        write(4) BPressure
        write(4) PressGradx,PressGrady,PressGradz
        write(4) BPressGradx,BPressGrady,BPressGradz
        write(4) PressGradfx,PressGradfy,PressGradfz
c
        if(.not.LsolveMomentum) then
c
          write(4) drhodP
          write(4) BdrhodP
c
        endif
c
        if(LUnsteady) then
c
          write(4) PressureOld,PressureOldOld
          write(4) BPressureOld,BPressureOldOld
          write(4) dpdt
c
          if(.not.LsolveMomentum) write(4) mdotOld
          if(.not.LsolveMomentum) write(4) mdotOldOld
c
        endif
c
      endif
c
      if(LSolveEnergy) then
c
        write(4) Temperature,TempGradx,TempGrady,TempGradz
        write(4) BTemperature,BTempGradx,BTempGrady,BTempGradz
        write(4) TempGradfx,TempGradfy,TempGradfz
c
        if(EnergyEquation.eq.'Htotal') then
c
          write(4) Htotal,HtotalGradx,HtotalGrady,HtotalGradz
          write(4) BHtotal,BHtotalGradx,BHtotalGrady,BHtotalGradz
          write(4) HtotalGradfx,HtotalGradfy,HtotalGradfz
c
        endif
c
        write(4) Conductivity,SpecificHeat
        write(4) BConductivity,BSpecificHeat
        write(4) SHeatGradx,SHeatGrady,SHeatGradz
        write(4) BSHeatGradx,BSHeatGrady,BSHeatGradz
        write(4) ScEnergy,SbEnergy
c
        if(LUnsteady) then
c
          write(4) TemperatureOld,TemperatureOldOld
          write(4) BTemperatureOld,BTemperatureOldOld
          write(4) SpecificHeatOld,SpecificHeatOldOld
          write(4) BSpecificHeatOld,BSpecificHeatOldOld
c
          if(EnergyEquation.eq.'Htotal') then
c
            write(4) HtotalOld,HtotalOldOld
            write(4) BHtotalOld,BHtotalOldOld
c
          endif
c
        endif
c
      endif
c
      if(LSolveMomentum.and.LBuoyancy) then
c
        write(4) Buoyancyx,Buoyancyy,Buoyancyz
        write(4) BBuoyancyx,BBuoyancyy,BBuoyancyz
        write(4) Buoyancyfx,Buoyancyfy,Buoyancyfz
c
      endif
c
      if(Lcompressible.and.LsolveMomentum.and..not.LsolveEnergy) then
c      
        write(4) Temperature
        write(4) BTemperature
        write(4) SpecificHeat
        write(4) BSpecificHeat
c      
      endif
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        write(4) LSoutherLandScalar
        write(4) EquationOfStateScalar
        write(4) RGasScalar,GammaGasScalar,PrLaminarScalar
c
        write(4) Scalar,ScalarGradx,ScalarGrady,ScalarGradz
        write(4) BScalar,BScalarGradx,BScalarGrady,BScalarGradz
        write(4) ScalarGradfx,ScalarGradfy,ScalarGradfz
        write(4) ScScalar,SbScalar
        write(4) DiffusionCoefficient,SpecificHeatScalar
        write(4) BDiffusionCoefficient,BSpecificHeatScalar
c
        if(LUnsteady) then
c
          write(4) ScalarOld,ScalarOldOld
          write(4) BScalarOld,BScalarOldOld
          write(4) SpecificHeatScalarOld,SpecificHeatScalarOldOld
          write(4) BSpecificHeatScalarOld,BSpecificHeatScalarOldOld
c
        endif
c
        if(NumberOfPointSources.ne.0) then
c
          if(.not.LSolveEnergy) then
c
            write(4) iElementPointSource,xLocationOfPointSource,
     *                  yLocationOfPointSource,zLocationOfPointSource
c
          endif
c
          write(4) ScPointSourceScalar,SbPointSourceScalar
c
        endif
c
      endif
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        write(4) LSoutherLandrField
        write(4) EquationOfStaterField
        write(4) RGasrField,GammaGasrField,PrLaminarrField
c
        write(4) rField,rFieldGradx,rFieldGrady,rFieldGradz
        write(4) BrField,BrFieldGradx,BrFieldGrady,BrFieldGradz
        write(4) rFieldGradfx,rFieldGradfy,rFieldGradfz
        write(4) ScrField,SbrField
        write(4) cosThetaF,cosThetaE
        write(4) BcosThetaF
        write(4) ConductivityrField,SpecificHeatrField
        write(4) BConductivityrField,BSpecificHeatrField
        write(4) DensityrField,ViscosityrField
        write(4) BDensityrField,BViscosityrField
c
        if(LUnsteady) then
c
          write(4) rFieldOld,rFieldOldOld
          write(4) BrFieldOld,BrFieldOldOld
          write(4) SpecificHeatrFieldOld,SpecificHeatrFieldOldOld
          write(4) BSpecificHeatrFieldOld,BSpecificHeatrFieldOldOld
c
        endif
c
        if(LSurfaceTension) then
c
          write(4) Curvature
          write(4) BCurvature
c
          write(4) delrFieldMagnitude
          write(4) BdelrFieldMagnitude
          write(4) delrFieldMagnitudeGradx
          write(4) BdelrFieldMagnitudeGradx
          write(4) delrFieldMagnitudeGrady
          write(4) BdelrFieldMagnitudeGrady
          write(4) delrFieldMagnitudeGradz
          write(4) BdelrFieldMagnitudeGradz
c
        endif
c
      endif
c
      if(LSolveLambdaELEEquation) then
        write(4) LambdaELE,BLambdaELE
      endif
c
c--- Density is always needed
c 
      write(4) Density,DensGradx,DensGrady,DensGradz    
      write(4) BDensity,BDensGradx,BDensGrady,BDensGradz    
      write(4) Densityf      
c
      if(LUnsteady) then
c
        write(4) DensityOld,DensityOldOld  
        write(4) BDensityOld,BDensityOldOld  
c
      endif
c
      if(LFalseTransientMomentum) then
c
        write(4) DensityStar
        write(4) BDensityStar
c
      endif
c
      if(NumberOfPointSources.ne.0.and.LSolveEnergy) then
c
        write(4) iElementPointSource
        write(4) xLocationOfPointSource,yLocationOfPointSource,
     *           zLocationOfPointSource
        write(4) ScPointSourceEnergy,SbPointSourceEnergy
c
      endif
c
      if(LanisotropicDiffusion) then
c
        write(4) Conductivity11,Conductivity12,Conductivity13
        write(4) Conductivity22,Conductivity23,Conductivity33
        write(4) BConductivity11,BConductivity12,BConductivity13
        write(4) BConductivity22,BConductivity23,BConductivity33
c
        if(NumberOfScalarsToSolve.gt.0) then
c
          write(4) DiffusionCoefficient11,DiffusionCoefficient12,
     *                                       DiffusionCoefficient13
          write(4) DiffusionCoefficient22,DiffusionCoefficient23,
     *                                       DiffusionCoefficient33
          write(4) BDiffusionCoefficient11,BDiffusionCoefficient12,
     *                                       BDiffusionCoefficient13
          write(4) BDiffusionCoefficient22,BDiffusionCoefficient23,
     *                                       BDiffusionCoefficient33
c
        endif
c
      endif
c
      write(4) eDiffCoefficient
      write(4) BeDiffCoefficient
c
      if(LTurbulentFlow.and.LsolveMomentum) then
c
        write(4) rhok,drhokdx,drhokdy,drhokdz
        write(4) Brhok,Bdrhokdx,Bdrhokdy,Bdrhokdz
        write(4) rhoTED,ReT
        write(4) BrhoTED,BReT
c
        write(4) S11,S12,S13,S22,S23,S33,StrainRate
        write(4) BS11,BS12,BS13,BS22,BS23,BS33,BStrainRate
        write(4) W11,W12,W13,W22,W23,W33,Vorticity
        write(4) BW11,BW12,BW13,BW22,BW23,BW33,BVorticity
        write(4) Tau11,Tau12,Tau13,Tau22,Tau23,Tau33
        write(4) TauWall,ustar,ystar,uplus,yplus,uTau
        write(4) duplusdyplus,WallViscosity,KsPlus
c
        if(LCompressible) write(4) XiStarFmt
c
        if(TurbulenceModel.eq.'kepsilonv2f') then
c
          write(4) Ce1Coefficient,LScale,TScale
          write(4) BLScale,BTScale
c
        endif
c
        if(TurbulenceModel.eq.'kepsilonzetaf') then
c
          write(4) Ce1Coefficient,LScale,TScale
          write(4) BLScale,BTScale
c
          write(4) TurbulentZetaGrad2x,TurbulentZetaGrad2y
          write(4) TurbulentZetaGrad2z,TurbulentZetaGradxy
          write(4) TurbulentZetaGradxz
          write(4) BTurbulentZetaGrad2x,BTurbulentZetaGrad2y
          write(4) BTurbulentZetaGrad2z,BTurbulentZetaGradxy
          write(4) BTurbulentZetaGradxz
c
        endif
c
        if(TurbulenceModel.eq.'kepsilonrng') then
c
          write(4) C2eRNG,sigTKERNG,sigTEDRNG,BsigTKERNG,BsigTEDRNG
c
        endif
c
        if(TurbulenceModel.eq.'nut92') then
c
          write(4) TurbulenceProduction,BTurbulenceProduction
          write(4) TurbulentKE,BTurbulentKE
          write(4) ModifiedDensGradx,ModifiedDensGrady,ModifiedDensGradz
          write(4) BModifiedDensGradx,BModifiedDensGrady
          write(4) BModifiedDensGradz
          write(4) ModifiedDensGrad2x,ModifiedDensGrad2y
          write(4) ModifiedDensGrad2z
          write(4) BModifiedDensGrad2x,BModifiedDensGrad2y
          write(4) BModifiedDensGrad2z
c
          write(4) uVelGrad2x,BuVelGrad2x,uVelGrad2y
          write(4) BuVelGrad2y,uVelGrad2z,BuVelGrad2z
          write(4) vVelGrad2x,BvVelGrad2x,vVelGrad2y
          write(4) BvVelGrad2y,vVelGrad2z,BvVelGrad2z
          write(4) wVelGrad2x,BwVelGrad2x,wVelGrad2y
          write(4) BwVelGrad2y,wVelGrad2z,BwVelGrad2z
          write(4) uvVelGradxy,BuvVelGradxy
          write(4) uwVelGradxz,BuwVelGradxz
c
          write(4) ModifiedEDGrad2x,ModifiedEDGrad2y
          write(4) ModifiedEDGrad2z,BModifiedEDGrad2x
          write(4) BModifiedEDGrad2y,BModifiedEDGrad2z
          write(4) ModifiedMut,BModifiedMut
          write(4) ModifiedMutGradx,ModifiedMutGrady
          write(4) ModifiedMutGradz,BModifiedMutGradx
          write(4) BModifiedMutGrady,BModifiedMutGradz
          write(4) G1Nut92,G2Nut92,F1Nut92,F2Nut92
          write(4) N1Nut92,N2Nut92,BN1Nut92,N1Nut92Gradx
          write(4) N1Nut92Grady,N1Nut92Gradz,BN1Nut92Gradx
          write(4) BN1Nut92Grady,BN1Nut92Gradz
c
        endif
c
        if(TurbulenceModel.eq.'sstgama') then
c
          write(4) F1factor,BF1factor,fr1Coefficient,F2factor,BF2factor
          write(4) F3factor,BF3factor,F4factor,BF4factor
c
          write(4) NormalVelocity,BNormalVelocity
          write(4) NormalVelocityGradx,NormalVelocityGrady
          write(4) NormalVelocityGradz,BNormalVelocityGradx
          write(4) BNormalVelocityGrady,BNormalVelocityGradz
c
        endif
c
        if(TurbulenceModel.eq.'kepsilonrt') then
c
          write(4) fr1Coefficient,f2RT,Bf2RT,velRT,BvelRT
          write(4) velRTGradx,BvelRTGradx,velRTGrady,BvelRTGrady
          write(4) velRTGradz,BvelRTGradz
          write(4) fmuKE,BfmuKE,fmuRT,BfmuRT
          write(4) TurbulentViscosityKE,BTurbulentViscosityKE
          write(4) TurbulentViscosityRT,BTurbulentViscosityRT
c
          if(LKEpsilonRtRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
            write(4) DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
            if(LUnsteady) then
              write(4)S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
              write(4)S11OldOld,S12OldOld,S13OldOld,
     *                  S22OldOld,S23OldOld,S33OldOld
            endif            
          endif
        endif
c
        if(TurbulenceModel.eq.'kklomega') then
c
          write(4) TurbulentViscosity
          write(4) BTurbulentViscosity
c
          write(4) TurbulentViscosityTs,TurbulentViscosityTl
          write(4) BTurbulentViscosityTs,BTurbulentViscosityTl
          write(4) AlfaT,BAlfaT
          write(4) LambdaT,BLambdaT
          write(4) LambdaEff,BLambdaEff
          write(4) coefficientFW,BcoefficientFW
          write(4) TurbulentKTs,BTurbulentKTs
          write(4) ProductionKT,ProductionKL
c
          if(LSolveEnergy) then
c
            write(4) AlfaTheta,BAlfaTheta
c
          endif
c
        endif
c
        if(TurbulenceModel.eq.'realizable') then
c
          write(4) cmuR,c1R
          write(4) BcmuR,Bc1R
c
        endif
c
        if(TurbulenceModel.eq.'spalartallmaras'.or.
     *                    TurbulenceModel.eq.'wrayagarwal') then
c
          write(4) TurbulenceProduction
          write(4) BTurbulenceProduction
          write(4) TurbulentKE,BTurbulentKE
c
          if(TurbulenceModel.eq.'spalartallmaras') then
            write(4) TGamma,BTGamma
            write(4) Stelda,fwCoefficient,ft2Coefficient,fr1Coefficient
            write(4) fv1Coefficient
            write(4) Bfv1Coefficient
            if(LNegativeSpalartAllmaras) then
              write(4) fnCoefficient
              write(4) BfnCoefficient
            endif
            if(LSpalartAllmarasRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
              write(4) DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
              if(LUnsteady) then
                write(4)S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
                write(4)S11OldOld,S12OldOld,S13OldOld,
     *                  S22OldOld,S23OldOld,S33OldOld
              endif            
            endif
          elseif(TurbulenceModel.eq.'wrayagarwal') then
            write(4) f1WA,fmuWA,SRateGradx,SRateGrady,SRateGradz
            write(4) Bf1WA,BfmuWA,BSRateGradx,BSRateGrady,BSRateGradz
            write(4) fr1Coefficient
            if(LWrayAgarwalRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
              write(4)DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
              if(LUnsteady) then
                write(4)S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
                write(4)S11OldOld,S12OldOld,S13OldOld,
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
          write(4) fmuCoefficient,f1Coefficient,f2Coefficient,LTKE,LTED
          write(4) BfmuCoefficient
c
          if(TurbulenceModel.eq.'kepsilonsharma'.or.
     *          TurbulenceModel.eq.'kepsilonhishida') then
c
            write(4) uVelGrad2x,BuVelGrad2x,uVelGrad2y
            write(4) BuVelGrad2y,uVelGrad2z,BuVelGrad2z
            write(4) vVelGrad2x,BvVelGrad2x,vVelGrad2y
            write(4) BvVelGrad2y,vVelGrad2z,BvVelGrad2z
            write(4) wVelGrad2x,BwVelGrad2x,wVelGrad2y
            write(4) BwVelGrad2y,wVelGrad2z,BwVelGrad2z
            write(4) uvVelGradxy,BuvVelGradxy
            write(4) uwVelGradxz,BuwVelGradxz
c
            write(4) sqrtTurbulentKE
            write(4) BsqrtTurbulentKE
            write(4) sqrtTKEGradx,sqrtTKEGrady,sqrtTKEGradz
            write(4) BsqrtTKEGradx,BsqrtTKEGrady,BsqrtTKEGradz
            write(4) SRateGradx,SRateGrady,SRateGradz
            write(4) BSRateGradx,BSRateGrady,BSRateGradz
c
          endif
c
        endif
c
        if(TurbulenceModel.eq.'komega2006lrn') then
c
          write(4) alfa,alfaStar,bettaStar
          write(4) Balfa,BalfaStar,BbettaStar
c
        endif
c
        if(TurbulenceModel.eq.'komegabsl'.or.
     *                    TurbulenceModel.eq.'komegasst') then
c
          write(4) F1factor
          write(4) BF1factor
c
          if(TurbulenceModel.eq.'komegasst') then
c
            write(4) F2factor,F3factor,F4factor,fr1Coefficient
            write(4) BF2factor,BF3factor,BF4factor
c
            if(LKOmegaSSTRotationCurvatureCorrection.and.
     *                    RotationCurvatureMethod.eq.'spalartshur') then
              write(4) DS11Dt,DS12Dt,DS13Dt,DS22Dt,DS23Dt,DS33Dt
c            
              if(LUnsteady) then
                write(4) S11Old,S12Old,S13Old,S22Old,S23Old,S33Old
                write(4) S11OldOld,S12OldOld,S13OldOld,
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
          write(4) TGammaEff,F1factor,fr1Coefficient,F2factor,
     *             F3factor,F4factor
          write(4) BF1factor,BF2factor,BF3factor,BF4factor
c
        endif
c
        write(4) TurbulentViscosity
        write(4) BTurbulentViscosity
c
        if(TurbulenceModel.eq.'komega2006'.or.
     *             TurbulenceModel.eq.'komega2006lrn') then
c
          write(4) TurbulentViscosity1
          write(4) BTurbulentViscosity1
c
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          write(4) TurbulentKE,TKEGradx,TKEGrady,TKEGradz
          write(4) BTurbulentKE,BTKEGradx,BTKEGrady,BTKEGradz
          write(4) TKEGradfx,TKEGradfy,TKEGradfz
          write(4) ScTKE,SbTKE
          write(4) TurbulenceProduction
          write(4) BTurbulenceProduction
c
          if(LBuoyancy) then
            write(4) TurbulenceProductionB
            write(4) BTurbulenceProductionB
          endif
c
          if(LUnsteady) then
            write(4) TurbulentKEOld,TurbulentKEOldOld
            write(4) BTurbulentKEOld,BTurbulentKEOldOld
          endif
c
        endif
c
        if(LSolveTurbulenceDissipationRate) then
c
          write(4) TurbulentED,TEDGradx,TEDGrady,TEDGradz
          write(4) BTurbulentED,BTEDGradx,BTEDGrady,BTEDGradz
          write(4) TEDGradfx,TEDGradfy,TEDGradfz
          write(4) ScTED,SbTED
c
          if(LUnsteady) then
c
            write(4) TurbulentEDOld,TurbulentEDOldOld
            write(4) BTurbulentEDOld,BTurbulentEDOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceSpecificDissipationRate) then
c
          write(4) TurbulentOmega,TOmegaGradx,TOmegaGrady,TOmegaGradZ
          write(4) BTurbulentOmega
          write(4) BTOmegaGradx,BTOmegaGrady,BTOmegaGradZ
          write(4) TOmegaGradfx,TOmegaGradfy,TOmegaGradfy
          write(4) ScTOmega,SbTOmega
c
          if(LUnsteady) then
c
            write(4) TurbulentOmegaOld,TurbulentOmegaOldOld
            write(4) BTurbulentOmegaOld,BTurbulentOmegaOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceGammaEquation) then
c
          write(4) TGamma,TGammaGradx,TGammaGrady,TGammaGradz
          write(4) BTGamma,BTGammaGradx,BTGammaGrady,BTGammaGradz
          write(4) TGammaGradfx,TGammaGradfy,TGammaGradfz
          write(4) ScTGamma,SbTGamma
c
          if(LUnsteady) then
c
            write(4) TGammaOld,TGammaOldOld
            write(4) BTGammaOld,BTGammaOldOld
c
          endif
c
        endif
c
       if(LSolveTurbulenceReynoldsThetaEquation) then
c
          write(4) TReTheta,TReThetaGradx,TReThetaGrady,TReThetaGradz
          write(4) BTReTheta,BTReThetaGradx,
     *                     BTReThetaGrady,BTReThetaGradz
          write(4) TReThetaGradfx,TReThetaGradfy,TReThetaGradfz
          write(4) ScTReTheta,SbTReTheta
c
          if(LUnsteady) then
c
            write(4) TReThetaOld,TReThetaOldOld
            write(4) BTReThetaOld,BTReThetaOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulentKL) then
c
          write(4) TurbulentKL,TurbulentKLGradx,
     *              TurbulentKLGrady,TurbulentKLGradz
          write(4) BTurbulentKL,BTurbulentKLGradx,
     *              BTurbulentKLGrady,BTurbulentKLGradz
          write(4) TurbulentKLGradfx,TurbulentKLGradfy,TurbulentKLGradfz
          write(4) ScTurbulentKL,SbTurbulentKL
          write(4) uVelGrad2x,uVelGrad2y,uVelGrad2z         
          write(4) vVelGrad2x,vVelGrad2y,vVelGrad2z         
          write(4) wVelGrad2x,wVelGrad2y,wVelGrad2z         
          write(4) BuVelGrad2x,BuVelGrad2y,BuVelGrad2z         
          write(4) BvVelGrad2x,BvVelGrad2y,BvVelGrad2z         
          write(4) BwVelGrad2x,BwVelGrad2y,BwVelGrad2z         
          write(4) uvVelGradxy,uwVelGradxz
          write(4) BuvVelGradxy,BuwVelGradxz
c
          if(LUnsteady) then
c
            write(4) TurbulentKLOld,TurbulentKLOldOld
            write(4) BTurbulentKLOld,BTurbulentKLOldOld
c
          endif
c
        endif
c
        if(LSolveModifiedED) then
c
          write(4) ModifiedED,BModifiedED
          write(4) ModifiedEDGradx,BModifiedEDGradx
          write(4) ModifiedEDGrady,BModifiedEDGrady
          write(4) ModifiedEDGradz,BModifiedEDGradz
          write(4) ModifiedEDGradfx
          write(4) ModifiedEDGradfy
          write(4) ModifiedEDGradfz
          write(4) ScModifiedED,SbModifiedED
c
          if(LUnsteady) then
c
            write(4) ModifiedEDOld,BModifiedEDOld
            write(4) ModifiedEDOldOld,BModifiedEDOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceV2Equation) then
c
          write(4) TurbulentV2,BTurbulentV2,TurbulentV2Gradx
          write(4) BTurbulentV2Gradx,TurbulentV2Grady,BTurbulentV2Grady
          write(4) TurbulentV2Gradz,BTurbulentV2Gradz,TurbulentV2Gradfx
          write(4) TurbulentV2Gradfy,TurbulentV2Gradfz
          write(4) ScTurbulentV2,SbTurbulentV2
c
          if(LUnsteady) then
c
            write(4) TurbulentV2Old,BTurbulentV2Old
            write(4) TurbulentV2OldOld,BTurbulentV2OldOld
c
          endif
c
        endif
c
        if(LSolveTurbulenceZetaEquation) then
c
          write(4) TurbulentZeta,BTurbulentZeta,TurbulentZetaGradx
          write(4) BTurbulentZetaGradx,TurbulentZetaGrady
          write(4) BTurbulentZetaGrady,TurbulentZetaGradz
          write(4) BTurbulentZetaGradz,TurbulentZetaGradfx
          write(4) TurbulentZetaGradfy,TurbulentZetaGradfz
          write(4) ScTurbulentZeta,SbTurbulentZeta
c
          if(LUnsteady) then
c
            write(4) TurbulentZetaOld,BTurbulentZetaOld
            write(4) TurbulentZetaOldOld,BTurbulentZetaOldOld
c
          endif
c
        endif
c
        if(LSolveTurbulencefRelaxationEquation) then
c
          write(4) TfRelaxation,BTfRelaxation,TfRelaxationGradx
          write(4) BTfRelaxationGradx,TfRelaxationGrady
          write(4) BTfRelaxationGrady,TfRelaxationGradz
          write(4) BTfRelaxationGradz,TfRelaxationGradfx
          write(4) TfRelaxationGradfy,TfRelaxationGradfz
          write(4) ScTfRelaxation,SbTfRelaxation
c
          if(LUnsteady) then
c
            write(4) TfRelaxationOld,BTfRelaxationOld
            write(4) TfRelaxationOldOld,BTfRelaxationOldOld
c
          endif
c
        endif
c
        if(NumberOfScalarsToSolve.ne.0) then
c
          write(4) SigScalar
c
        endif
c
      endif
c
      write(4) MassFlowRate,MassFlowRateOld,MassFlowRateOldOld
      write(4) OutletPressure,OutletPressureOld,OutletPressureOldOld
c      
      if(LUnsteady) then
        write(4) dt,dtOld,dtOldOld
        write(4) time,ndt
      endif
c
      return
      end