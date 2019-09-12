c
C#############################################################################################
c
      SUBROUTINE InitializeOldVariables
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
     *                 MethodCalcGradientMomentum,LimitGradientMomentum,
     *                 nIterGradientMomentum,
     *                 LimitGradientMomentumMethod,
     *                 LSpalartAllmarasRotationCurvatureCorrection,
     *                 LTurbulentFlow,TurbulenceModel,
     *                 LWrayAgarwalRotationCurvatureCorrection,
     *                 LKOmegaSSTRotationCurvatureCorrection,
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
     *                      TGamma,TGammaOld,TGammaOldOld,
     *                      BTGamma,BTGammaOld,BTGammaOldOld,
     *                      TReTheta,TReThetaOld,TReThetaOldOld,
     *                      BTReTheta,BTReThetaOld,BTReThetaOldOld,
     *                      ModifiedED,BModifiedED,ModifiedEDOld,
     *                      BModifiedEDOld,ModifiedEDOldOld,
     *                      BModifiedEDOldOld,
     *                      TurbulentKL,BTurbulentKL,
     *                      TurbulentKLOld,BTurbulentKLOld,
     *                      TurbulentKLOldOld,BTurbulentKLOldOld,
     *                      Htotal,BHtotal,HtotalOld,BHtotalOld,
     *                      HtotalOldOld,BHtotalOldOld,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
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
      character*10 Variable
c
c********************************************************************************************
c
c--- Interfaces
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------
          character*10 Variable
          integer :: MethodCalcGradient,nIterGradientPhi
          logical :: LimitGradient
          integer :: LimitGradientMethod
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------
      end interface
C*********************************************************************************************
c
      dtOldOld=dtOld
      dtOld=dt
c
      DensityOld=Density
      DensityOldOld=DensityOld
      BDensityOld=BDensity
      BDensityOldOld=BDensityOld
c
      SpecificHeatOld=SpecificHeat
      SpecificHeatOldOld=SpecificHeatOld
      BSpecificHeatOld=BSpecificHeat
      BSpecificHeatOldOld=BSpecificHeatOld
c
      if(LsolveMomentum) then
c         
        uVelocityOld=uVelocity
        vVelocityOld=vVelocity
        wVelocityOld=wVelocity
        uVelocityOldOld=uVelocityOld
        vVelocityOldOld=vVelocityOld
        wVelocityOldOld=wVelocityOld
c
        BuVelocityOld=BuVelocity
        BvVelocityOld=BvVelocity
        BwVelocityOld=BwVelocity
        BuVelocityOldOld=BuVelocityOld
        BvVelocityOldOld=BvVelocityOld
        BwVelocityOldOld=BwVelocityOld
c          
      endif
c
      if(LsolveContinuity) then
c         
        PressureOld=Pressure
        PressureOldOld=PressureOld
        BPressureOld=BPressure
        BPressureOldOld=BPressureOld
        mdotOld=mdot
        mdotOldOld=mdotOld
        BmdotOld=Bmdot
        BmdotOldOld=BmdotOld
c          
      endif
c
      if(LsolveEnergy) then
c         
        TemperatureOld=Temperature
        TemperatureOldOld=TemperatureOld
        BTemperatureOld=BTemperature
        BTemperatureOldOld=BTemperatureOld
c          
      endif
c
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') then
c         
        HtotalOld=Htotal
        HtotalOldOld=HtotalOld
        BHtotalOld=BHtotal
        BHtotalOldOld=BHtotalOld
c          
      endif
c
      if(LsolveTurbulenceKineticEnergy) then
c         
        TurbulentKEOld=TurbulentKE
        TurbulentKEOldOld=TurbulentKEOld
        BTurbulentKEOld=BTurbulentKE
        BTurbulentKEOldOld=BTurbulentKEOld
c          
      endif
c
      if(LSolveTurbulenceDissipationRate) then
c         
        TurbulentEDOld=TurbulentED
        TurbulentEDOldOld=TurbulentEDOld
        BTurbulentEDOld=BTurbulentED
        BTurbulentEDOldOld=BTurbulentEDOld
c          
      endif
c
      if(LSolveTurbulenceSpecificDissipationRate) then
c         
        TurbulentOmegaOld=TurbulentOmega
        TurbulentOmegaOldOld=TurbulentOmegaOld
        BTurbulentOmegaOld=BTurbulentOmega
        BTurbulentOmegaOldOld=BTurbulentOmegaOld
c          
      endif
c
      if(LSolveTurbulenceGammaEquation) then
c         
        TGammaOld=TGamma
        TGammaOldOld=TGammaOld
        BTGammaOld=BTGamma
        BTGammaOldOld=BTGammaOld
c          
      endif
c
      if(LSolveTurbulenceReynoldsThetaEquation) then
c         
        TReThetaOld=TReTheta
        TReThetaOldOld=TReThetaOld
        BTReThetaOld=BTReTheta
        BTReThetaOldOld=BTReThetaOld
c          
      endif
c
      if(LSolveTurbulentKL) then
c         
        TurbulentKLOld=TurbulentKL
        TurbulentKLOldOld=TurbulentKLOld
        BTurbulentKLOld=BTurbulentKL
        BTurbulentKLOldOld=BTurbulentKLOld
c          
      endif
c
      if(LSolveModifiedED) then
c         
        ModifiedEDOld=ModifiedED
        ModifiedEDOldOld=ModifiedEDOld
        BModifiedEDOld=BModifiedED
        BModifiedEDOldOld=BModifiedEDOld
c          
      endif
c
      if(LSolveTurbulenceV2Equation) then
c         
        TurbulentV2Old=TurbulentV2
        TurbulentV2OldOld=TurbulentV2Old
        BTurbulentV2Old=BTurbulentV2
        BTurbulentV2OldOld=BTurbulentV2Old
c          
      endif
c
      if(LSolveTurbulenceZetaEquation) then
c         
        TurbulentZetaOld=TurbulentZeta
        TurbulentZetaOldOld=TurbulentZetaOld
        BTurbulentZetaOld=BTurbulentZeta
        BTurbulentZetaOldOld=BTurbulentZetaOld
c          
      endif
c
      do irField=1,NumberOfrFieldsToSolve
c
        if(LSolverField(irField)) then
c
          do i=1,NumberOfElements
c      
            rFieldOld(i,irField)=rField(i,irField)
            rFieldOldOld(i,irField)=rFieldOld(i,irField)
            SpecificHeatrFieldOld(i,irField)=
     *                                 SpecificHeatrField(i,irField)
            SpecificHeatrFieldOldOld(i,irField)=
     *                              SpecificHeatrFieldOld(i,irField)
c
          enddo      
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BrFieldOld(i,j,irField)=BrField(i,j,irField)
              BrFieldOldOld(i,j,irField)=BrFieldOld(i,j,irField)
              BSpecificHeatrFieldOld(i,j,irField)=
     *                                 BSpecificHeatrField(i,j,irField)
              BSpecificHeatrFieldOldOld(i,j,irField)=
     *                              BSpecificHeatrFieldOld(i,j,irField)
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
            ScalarOld(i,iScalar)=Scalar(i,iScalar)
            ScalarOldOld(i,iScalar)=ScalarOld(i,iScalar)
            SpecificHeatScalarOld(i,iScalar)=
     *                                 SpecificHeatScalar(i,iScalar)
            SpecificHeatScalarOldOld(i,iScalar)=
     *                              SpecificHeatScalarOld(i,iScalar)
c
          enddo      
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BScalarOld(i,j,iScalar)=BScalar(i,j,iScalar)
              BScalarOldOld(i,j,iScalar)=BScalarOld(i,j,iScalar)
              BSpecificHeatScalarOld(i,j,iScalar)=
     *                                 BSpecificHeatScalar(i,j,iScalar)
              BSpecificHeatScalarOldOld(i,j,iScalar)=
     *                              BSpecificHeatScalarOld(i,j,iScalar)
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
            MassFlowrateOld(i)=MassFlowrate(i)
            MassFlowrateOldOld(i)=MassFlowrateOld(i)
            OutletPressureOld(i)=OutletPressure(i)
            OutletPressureOldOld(i)=OutletPressureOld(i)
c
          endif      
c
        endif      
c      
      enddo
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.'spalartallmaras'.and.
     *               LSpalartAllmarasRotationCurvatureCorrection.and.
     *                  RotationCurvatureMethod.eq.'spalartshur') then
c
c---- Calculate Gradients of the u,v,and p variables
c
        Variable='velx'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelocity,uVelGradx,uVelGrady,uVelGradz,BuVelocity,
     *       BuVelGradx,BuVelGrady,BuVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        Variable='vely'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelocity,vVelGradx,vVelGrady,vVelGradz,BvVelocity,
     *        BvVelGradx,BvVelGrady,BvVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        Variable='velz'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelocity,wVelGradx,wVelGrady,wVelGradz,BwVelocity,
     *        BwVelGradx,BwVelGrady,BwVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        call calculateStrainRatetensor
c
        S11old=S11
        S12old=S12
        S13old=S13
        S22old=S22
        S23old=S23
        S33old=S33
c
        S11oldold=S11old
        S12oldold=S12old
        S13oldold=S13old
        S22oldold=S22old
        S23oldold=S23old
        S33oldold=S33old
c
      endif
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.'wrayagarwal'.and.
     *               LWrayAgarwalRotationCurvatureCorrection.and.
     *                  RotationCurvatureMethod.eq.'spalartshur') then
c
c---- Calculate Gradients of the u,v,and p variable
c
        Variable='velx'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelocity,uVelGradx,uVelGrady,uVelGradz,BuVelocity,
     *       BuVelGradx,BuVelGrady,BuVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        Variable='vely'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelocity,vVelGradx,vVelGrady,vVelGradz,BvVelocity,
     *        BvVelGradx,BvVelGrady,BvVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        Variable='velz'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelocity,wVelGradx,wVelGrady,wVelGradz,BwVelocity,
     *        BwVelGradx,BwVelGrady,BwVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        call calculateStrainRatetensor
c
        S11old=S11
        S12old=S12
        S13old=S13
        S22old=S22
        S23old=S23
        S33old=S33
c
        S11oldold=S11old
        S12oldold=S12old
        S13oldold=S13old
        S22oldold=S22old
        S23oldold=S23old
        S33oldold=S33old
c
      endif
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.'komegasst'.and.
     *               LKOmegaSSTRotationCurvatureCorrection.and.
     *                  RotationCurvatureMethod.eq.'spalartshur') then
c
c---- Calculate Gradients of the u,v,and p variable
c
        Variable='velx'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelocity,uVelGradx,uVelGrady,uVelGradz,BuVelocity,
     *       BuVelGradx,BuVelGrady,BuVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        Variable='vely'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelocity,vVelGradx,vVelGrady,vVelGradz,BvVelocity,
     *        BvVelGradx,BvVelGrady,BvVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        Variable='velz'
c
        call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelocity,wVelGradx,wVelGrady,wVelGradz,BwVelocity,
     *        BwVelGradx,BwVelGrady,BwVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
        call calculateStrainRatetensor
c
        S11old=S11
        S12old=S12
        S13old=S13
        S22old=S22
        S23old=S23
        S33old=S33
c
        S11oldold=S11old
        S12oldold=S12old
        S13oldold=S13old
        S22oldold=S22old
        S23oldold=S23old
        S33oldold=S33old
c
      endif
c
      return
      end        