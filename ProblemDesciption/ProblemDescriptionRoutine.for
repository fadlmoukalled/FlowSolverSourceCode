c
C#############################################################################################
      SUBROUTINE DescribeProblem
C#############################################################################################
c
      use User0
      use Geometry1
      use Geometry2
      use Geometry3
      use Geometry4
      use BoundaryConditions1
      use BoundaryConditions2
      use BoundaryConditionsScalar1
      use BoundaryConditionsScalar2
      use BoundaryConditionsrField1
      use BoundaryConditionsrField2
      use BoundaryConditionsTurbulence1
      use BoundaryConditionsTurbulence2
      use Scalar1
      use VolumeOfFluid1
      use BoundaryFluxes
      use WallStress1
      use PhysicalProperties1
      use MultiGrid2
      use Variables1
      use Variables4
      use Turbulence1
      use constants1
      use Coriolis1
      use ReferenceValues1
      use Transient1, only: ndt
      use FlowInOut1, only: MassFlowFraction,LPrintMassFlowRate,
     *                      MassFlowRate
      use ArteryResistance1, only: ArteryResistance,LArteryExplicit,
     *                             urfPressureResistance
      use AveragePressure1, only: AverageOutletPressure
      use WindKessel1, only: ResistanceToBloodFlow,ComplianceC,InertiaL,
     *                       TotalPeripheralResistance,OutletPressure
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k
      double precision :: Uinlet,Vinlet,Winlet,Pinlet,Tinlet,ActualMdot,
     *                    temp,Reynolds,Minlet,RhoInlet,angle,
     *                    chord,Taw
c
      double precision :: period,localTime,TKEMin,TKEMax
      integer :: timeFloor
c*********************************************************************************************
      interface
c*********************************************************************************************
        SUBROUTINE PhiFieldMinMax(FiT,BFiT,phiMin,phiMax)
c*********************************************************************************************
          double precision :: phiMin,phiMax
          double precision, dimension(:) :: FiT
          double precision, dimension(:,:) :: BFiT
c*********************************************************************************************
        end SUBROUTINE PhiFieldMinMax
c*********************************************************************************************
      end interface
c*********************************************************************************************
c
c--- GradientInterpolationScheme : 'average', 'upwind', 'downwind', 'averagecorrected'
c
C********************************************************************************************
      entry OpenFiles   !start with a unit number 20 and above (below 20 reserved for internal use)
C********************************************************************************************
c       
      PolyMeshDirectory ="C:\Fadl\SourceDecomposed July 20 2019\
     *Console1\Console1\Console1\Testing\polyMesh"
      NeutralMeshDirectory ="C:\Fadl\SourceDecomposed July 20 2019\
     *Console1\Console1\Console1\Testing\NeutralMesh"
      SolutionDirectory="C:\Fadl\SourceDecomposed July 20 2019\
     *Console1\Console1\Console1\Testing\Solution"
c
!      open(unit=20,status='unknown',
!     *   file=trim(SolutionDirectory)//"/MassFlowRate",
!     *              access='sequential',form='formatted')
!      open(unit=21,status='unknown',
!     *   file=trim(SolutionDirectory)//"/Pressure",
!     *              access='sequential',form='formatted')
!      open(unit=22,status='unknown',
!     *   file=trim(SolutionDirectory)//"/Shear",
!     *              access='sequential',form='formatted')
!      open(unit=23,status='unknown',
!     *   file=trim(SolutionDirectory)//"/TKE",
!     *              access='sequential',form='formatted')
c
      return      
c
C********************************************************************************************
      entry SetVariablesToSolve
C********************************************************************************************
      LReadOldSolution=.false.
      LStructured=.false.
      LReadSavedGrid=.false.
      LReadSavedWallDistance=.false.
      LimitTemperature=.false.
      LTurbulentFlow=.false.
      TurbulenceModel='spalartallmaras' 
      WallTreatment='lowreynoldsnumber' 
c
      MeshType='polymesh'
      LprintParaviewNodeBased=.true.
      LprintParaviewCellBased=.true.
      LprintParaviewFile=.true.
c
c--- Set the type of flow to solve
c
      Lcompressible=.false.
      Linviscid=.false.
      LUnsteady=.false.
      LCoriolis=.false.
      LFreeSurfaceFlow=.false.
      LanisotropicDiffusion=.true.  !should always be true for LambdaELE
      LBuoyancy=.false.
      LSurfaceTension=.false.
c
c---- Set variables to solve
c
      LSolveMomentum=.false.
      LSolveContinuity=.false.
      LSolveEnergy=.false.
      LSolveTurbulenceKineticEnergy=.false.
      LSolveTurbulenceDissipationRate=.false.
      LSolveTurbulenceSpecificDissipationRate=.false.
      LSolveTurbulentKL=.false.
      LSolveModifiedED=.false.
      LSolveTurbulenceGammaEquation=.false.
      LSolveTurbulenceReynoldsThetaEquation=.false.
      LSolveTurbulencefRelaxationEquation=.false.
      LSolveTurbulenceV2Equation=.false.
      LSolveTurbulenceZetaEquation=.false.
      LSolveLambdaELEEquation=.true.
c
      EnergyEquation='temperature'
c
c--- In case momentum is not solve and only the convection of a scalar is sought,
c    set the below parameter to true
c
      LConvectScalar=.false.
c
c--- Set the number of rFields to solve for free surface flows
c
      NumberOfrFieldsToSolve=0
c
c--- Set the number of additional scalar variables to solve
c
      NumberOfScalarsToSolve=0
c
c--- Set the number of point sources in the domain
c
      NumberofPointSources=0
c
c--- Declare whether false transient underrelaxation will be used
c
      LFalseTransientMomentum=.false.
      LFalseTransientEnergy=.false.
      LFalseTransientScalar=.false.
      LFalseTransientTKE=.false.
      LFalseTransientTED=.false.
      LFalseTransientTOmega=.false.
      LFalseTransientTurbulentKL=.false.
      LFalseTransientModifiedED=.false.
      LFalseTransientTGamma=.false.
      LFalseTransientTReTheta=.false.
      LFalseTransientTfRelaxation=.false.
      LFalseTransientTurbulentV2=.false.
      LFalseTransientTurbulentZeta=.false.
      LFalseTransientLambdaELE=.false.
      LFalseTransientrField=.false.
c
c--- Assign the false time step values
c
      FalseDtMomentum=1.e-2
      FalseDtEnergy=1.e-3
      FalseDtTKE=1.e-5
      FalseDtTOmega=1.e-5
      FalseDtModifiedED=1.e-3
      FalseDtTED=1.e-5
      FalseDtTGamma=1.e-5
      FalseDtTReTheta=1.e-5
      FalseDtTfRelaxation=1.e-5
      FalseDtTurbulentV2=1.e-5
      FalseDtTurbulentZeta=1.e-5
      FalseDtLambdaELE=1.e-5
c      FalseDtScalar=2.e-4
      FalseDtrField=1.e-2
c
      return
c
C********************************************************************************************
      entry SetGlobalVariables
C********************************************************************************************
c
      if(LsolveMomentum) then
c
c--- Set the algorithm type
c
        Lsimple=.true.
        Lsimplec=.false.
c
        AngularVelocityX=0.
        AngularVelocityY=0.
        AngularVelocityZ=0.
        AxisOfRotationOriginX=0.
        AxisOfRotationOriginY=0.
        AxisOfRotationOriginZ=0.
c
c--- Maximum number of global iterations
c
        IterMax=200
c
c--- Assign the monitoring location

        xMonitor=3.
        yMonitor=0.5
        zMonitor=0.0
c
c--- Reference pressure location for incompressible flows
c
        xRefPressure=3.
        yRefPressure=0.5
        zRefPressure=0.
c
c--- Stop criteria  1: sum of absolute residuals; 2: maximum absolute residual;
c                   3: rms residual             ; 4: normalized residual;
c
      NstopType=2
      maximumResidual=5.e-6
c
c--- Set the whether the pressure will be fixed at a point
c
      LfixPressure=.false.
c
c--- Set the location where the pressure will be fixed and its value
c
      FixedPressureValue=0.
      xFixPressure=0.
      yFixPressure=0.0499
      zFixPressure=0.
c
c--- Algebraic solvers of basic variables (pbcg, ilu, sor, direct)
c
      ASSolverMomentum='ilu'
      ASSolverContinuity='ilu'
      ASSolverTKE='ilu'
      ASSolverTED='ilu'
      ASSolverTOmega='ilu'
      ASSolverTurbulentKL='ilu'
      ASSolverEnergy='ilu'
      ASSolverModifiedED='ilu'
      ASSolverTGamma='ilu'
      ASSolverTReTheta='ilu'
      ASSolverTfRelaxation='ilu'
      ASSolverTurbulentV2='ilu'
      ASSolverTurbulentZeta='ilu'
c
c--- Residual reduction factor of variables while solving their algebraic equations
c
      rrfMomentum=0.1
      rrfContinuity=0.001
      rrfTKE=0.1
      rrfTED=0.1
      rrfTOmega=0.1
      rrfTurbulentKL=0.1
      rrfEnergy=0.01
      rrfModifiedED=0.1
      rrfTGamma=0.1
      rrfTReTheta=0.1
      rrfTfRelaxation=0.1
      rrfTurbulentV2=0.1
      rrfTurbulentZeta=0.1
c
c--- maximum number of algebraic solver iterations
c
      ASIterMomentum=10
      ASIterContinuity=30
      ASIterTKE=3
      ASIterTED=3
      ASIterTOmega=3
      ASIterTurbulentKL=3
      ASIterEnergy=3
      ASIterModifiedED=3
      ASIterTGamma=3
      ASIterTReTheta=3
      ASIterTfRelaxation=3
      ASIterTurbulentV2=3
      ASIterTurbulentZeta=3
c
c--- Set the relaxation parameters
c
      LRelaxMomentum=.true.
      LRelaxPressure=.true.
      LRelaxEnergy=.true.
      LRelaxTKE=.true.
      LRelaxTED=.true.
      LRelaxTOmega=.true.
      LRelaxTurbulentKL=.true.
      LRelaxModifiedED=.true.
      LRelaxTGamma=.true.
      LRelaxTReTheta=.true.
      LRelaxTfRelaxation=.true.
      LRelaxTurbulentV2=.true.
      LRelaxTurbulentZeta=.true.
c
c--- Assign the underrelaxation factors values
c
      urfMomentum=0.7 !0.4 !0.3
      urfPressure=0.7 !0.6  !0.7
      urfTKE=0.5 !0.3
      urfTED=0.5
      urfTOmega=0.5 !0.3
      urfTurbulentKL=0.8 
      urfEnergy=0.3
      urfTViscosity=0.5
      urfModifiedED=0.8
      urfTGamma=0.8
      urfTReTheta=0.8
      urfTfRelaxation=0.1
      urfTurbulentV2=0.3
      urfTurbulentZeta=0.1
c
c--- Set whether to use upwind, NVF or TVD convection schemes
c
      LNVFMomentum=.false.
      LNVFTKE=.false.
      LNVFTED=.false.
      LNVFTOmega=.false.
      LNVFTurbulentKL=.false.
      LNVFEnergy=.false.
      LNVFModifiedED=.false.
      LNVFDensity=.false.
      LNVFTGamma=.false.
      LNVFTReTheta=.false.
      LNVFTfRelaxation=.false.
      LNVFTurbulentV2=.false.
      LNVFTurbulentZeta=.false.
c
      LTVDMomentum=.true.
      LTVDTKE=.true.
      LTVDTED=.true.
      LTVDTOmega=.true.
      LTVDTurbulentKL=.true.
      LTVDEnergy=.true.
      LTVDModifiedED=.true.
      LTVDDensity=.true.
      LTVDTGamma=.true.
      LTVDTReTheta=.true.
      LTVDTfRelaxation=.true.
      LTVDTurbulentV2=.true.
      LTVDTurbulentZeta=.true.
c
c--- Set when to start applying HR schemes
c
      nIterStartApplyingHR=1
c
c--- Set the name of scheme to use for variables
c
      ConvectionSchemeMomentum='minmod'
      ConvectionSchemeTKE='minmod'
      ConvectionSchemeTOmega='minmod'
      ConvectionSchemeTED='minmod'
      ConvectionSchemeTurbulentKL='minmod'
      ConvectionSchemeModifiedED='minmod'
      ConvectionSchemeEnergy='smart'
      ConvectionSchemeDensity='minmod'
      ConvectionSchemeTGamma='minmod'
      ConvectionSchemeTReTheta='minmod'
      ConvectionSchemeTfRelaxation='minmod'
      ConvectionSchemeTurbulentV2='minmod'
      ConvectionSchemeTurbulentZeta='minmod'
c
c--- Set the value of the coefficient by which to bleed the HR scheme with upwind scheme 
c    (0= no bleeding, 1= upwind) 
c
      BleedMomentum=0.
      BleedTKE=0.0
      BleedTED=0.
      BleedTOmega=0.0
      BleedTurbulentKL=0.
      BleedModifiedED=0.
      BleedEnergy=0.
      BleedDensity=0.
      BleedTGamma=0.0
      BleedTReTheta=0.0
      BleedTfRelaxation=0.0
      BleedTurbulentV2=0.0
      BleedTurbulentZeta=0.0
c----------------------------------------------------------------------------------
c--- Set the multigrid variables
c----------------------------------------------------------------------------------
      MGType='algebraic' !'geometricelement'
c
      LMultigridMomentum=.false.
      LMultigridTKE=.false.
      LMultigridTED=.false.
      LMultigridTOMega=.false.
      LMultigridTurbulentKL=.false.
      LMultigridModifiedED=.false.
      LMultigridContinuity=.false.
      LMultigridEnergy=.false.
      LMultigridTGamma=.false.
      LMultigridTReTheta=.false.
      LMultigridTfRelaxation=.false.
      LMultigridTurbulentV2=.false.
      LMultigridTurbulentZeta=.false.
c
c--- Set whether to print data during multigrid iterations (for testing)
c
      LPrintMultiGridResiduals=.false.
      LTestMultiGrid=.true.
      nprintMG=2000
c
c--- Set the variable on which to base the grid agglomoration (could be a scalar)
c
      MGVariable='pressurecorrection'
c
c--- Set the multigrid cycle to use (V, F, or W)
c
      MultiGridCycleType = 'wcycle'
c
c--- Set the number of pre and post sweep number of iterations and the number of MG cycles
c
      MultiGridpreSweep=1
      MultiGridpostSweep=3
      MaxMultiGridCycles=30
c
c--- Set when to re-agglomorate fine grids to create coarse grid levels
c
      reAgglomorate=50000000
c
c--- Set the minimum number of elements of a coarse level and the maximum number of coarse levels
c
      minNumberofParents=10
      maxNumberofCoarseLevels=10
c
c--- Set when to start applying multigrid (a number >=1)
c
      nIterStartApplyingMG=1
c
c--- Set the multigrid Residual reduction factor
c
      MultiGridrrf=0.1
c----------------------------------------------------------------------------------
c--- Set the gradient variables
c----------------------------------------------------------------------------------
c
c--- Method to calculate the gradient (1 to 6)
c
      MethodCalcGradientMomentum =2
      MethodCalcGradientTKE =2
      MethodCalcGradientTED =2
      MethodCalcGradientTOmega =2
      MethodCalcGradientTurbulentKL =2
      MethodCalcGradientModifiedED =2
      MethodCalcGradientContinuity=2
      MethodCalcGradientEnergy=2
      MethodCalcGradientDensity=2
      MethodCalcGradientTGamma=2
      MethodCalcGradientTReTheta=2
      MethodCalcGradientTfRelaxation=2
      MethodCalcGradientTurbulentV2=2
      MethodCalcGradientTurbulentZeta=2
c
c--- number of iterations to perform with isterative methods (2 to 4)
c
      nIterGradientMomentum=2
      nIterGradientTKE=2
      nIterGradientTED=2
      nIterGradientTOmega=2
      nIterGradientTurbulentKL=2
      nIterGradientModifiedED=2
      nIterGradientContinuity=2 
      nIterGradientEnergy=2
      nIterGradientDensity=2
      nIterGradientTGamma=2
      nIterGradientTReTheta=2
      nIterGradientTfRelaxation=2
      nIterGradientTurbulentV2=2
      nIterGradientTurbulentZeta=2
c
c--- The power to use for the inverse distance with the Least Square Gradient (method 6) 
c
      InvDistancePower=1.
c
c--- Set whether or not to limit the gradient 
c
      LimitGradientMomentum=.false.
      LimitGradientTKE=.false.
      LimitGradientTED=.false.
      LimitGradientTOmega=.false.
      LimitGradientTurbulentKL=.false.
      LimitGradientModifiedED=.false.
      LimitGradientContinuity=.false.
      LimitGradientEnergy=.false.
      LimitGradientDensity=.false.
      LimitGradientTGamma=.false.
      LimitGradientTReTheta=.false.
      LimitGradientTfRelaxation=.false.
      LimitGradientTurbulentV2=.false.
      LimitGradientTurbulentZeta=.false.
c
c--- Set the method to use to limit the gradient 
c
      LimitGradientMomentumMethod=2
      LimitGradientTKEMethod=2
      LimitGradientTEDMethod=2
      LimitGradientTOmegaMethod=2
      LimitGradientTurbulentKLMethod=2
      LimitGradientModifiedEDMethod=2
      LimitGradientContinuityMethod=2
      LimitGradientEnergyMethod=2
      LimitGradientDensityMethod=2
      LimitGradientTGammaMethod=2
      LimitGradientTReThetaMethod=2
      LimitGradientTfRelaxationMethod=2
      LimitGradientTurbulentV2Method=2
      LimitGradientTurbulentZetaMethod=2
c
c--- Set the type of gradient interpolation to face
c    available types: upwind, downwind, average, and averagecorrected 
c
      GradientInterpolationSchemeMomentum='average'
      GradientInterpolationSchemeTKE='average'
      GradientInterpolationSchemeTED='average'
      GradientInterpolationSchemeTOmega='average'
      GradientInterpolationSchemeTurbulentKL='average'
      GradientInterpolationSchemeModifiedED='average'
      GradientInterpolationSchemeContinuity='averagecorrected'
      GradientInterpolationSchemeEnergy ='average'
      GradientInterpolationSchemeTGamma ='average'
      GradientInterpolationSchemeTReTheta ='average'
      GradientInterpolationSchemeTfRelaxation ='average'
      GradientInterpolationSchemeTurbulentV2 ='average'
      GradientInterpolationSchemeTurbulentZeta ='average'
c
c----------------------------------------------------------------------------------
c--- Set the transient variables
c----------------------------------------------------------------------------------
c
c--- Set the scheme to use (euler, cranknicolson1 (2 steps), 
c                                  cranknicolson2, and adamsmoulton)
c
      TransientScheme='cranknicolson1'
c
c--- Set the step size
c
      LadaptiveTimeStep=.false.
      GlobalCourantNumber=2.
      minimumdtChangeFactor=0.5
      maximumdtChangeFactor=5
      dt=1.d-2
c
c--- Set the maximum time for which the solution to be computed
c
      timemax=30.
c
c--- Set the time to start computations
c
      time=0.
c
c--- Set the number of time steps tp print results
c
      nTimePrint=10
      nTimeSave=10
c
c----------------------------------------------------------------------------------
c--- Set buoyancy variables
c----------------------------------------------------------------------------------
c
c--- Two models are implemented: boussinesq and rhog
c
      LBuoyancy=.true.
      BuoyancyModel='rhog' !'boussinesq'
c
      GravityX=0.
      GravityY=-9.81
      GravityZ=0.
c----------------------------------------------------------------------------------
c--- Start describing the problem
c----------------------------------------------------------------------------------
c
c--- Assign the types of boundary conditions
c
c      MassFlowFraction(3)=1.
c      MassFlowFraction(3)=0.2
c      MassFlowFraction(4)=0.2
c      MassFlowFraction(5)=0.1
c      MassFlowFraction(6)=0.1
c
c       LPrintMassFlowRate(1)=.true.
c       LPrintMassFlowRate(2)=.true.
c       LPrintMassFlowRate(3)=.true.
c       LPrintMassFlowRate(4)=.true.
c       LPrintMassFlowRate(5)=.true.
c       LPrintMassFlowRate(6)=.true.
c       LPrintMassFlowRate(7)=.true.
c       LPrintMassFlowRate(8)=.true.
c       LPrintMassFlowRate(9)=.true.
c       LPrintMassFlowRate(10)=.true.
c       LPrintMassFlowRate(11)=.true.
c       LPrintMassFlowRate(12)=.true.
c       LPrintMassFlowRate(13)=.true.
c       LPrintMassFlowRate(14)=.true.
c
!       LArteryExplicit(1)=.false.
!       LArteryExplicit(2)=.false.
!       LArteryExplicit(3)=.false.
!       LArteryExplicit(4)=.false.
!       LArteryExplicit(5)=.false.
!       LArteryExplicit(6)=.false.
!       LArteryExplicit(7)=.false.
!c
!       urfPressureResistance(1)=1. !0.1
!       urfPressureResistance(2)=1. !0.1
!       urfPressureResistance(3)=1. !0.1
!       urfPressureResistance(4)=1. !0.1
!       urfPressureResistance(5)=1. !0.1
!       urfPressureResistance(6)=1. !0.1
!       urfPressureResistance(7)=1. !0.1
!c
!      ArteryResistance(1)=6452.8e5
!      ArteryResistance(2)=6452.8e5
!      ArteryResistance(3)=10056.3e5
!      ArteryResistance(4)=10615.0e5
!      ArteryResistance(5)=9553.5e5
!      ArteryResistance(6)=9553.5e5
!      ArteryResistance(7)=9553.5e5
c      
c---- WindKessel Parameters
c      
!      LWindKessel=.true.
!      WindKesselType=3
!      urfOutletPressure=0.1
!c
!      ResistanceToBloodFlow(2)=5.084e7
!      TotalPeripheralResistance(2)=8.57e8
!      ComplianceC(2)=1.193e-9
!c      
!      ResistanceToBloodFlow(3)=8.343e7
!      TotalPeripheralResistance(3)=1.406e9
!      ComplianceC(3)=8.95e-10
!c      
!      ResistanceToBloodFlow(4)=7.821e7
!      TotalPeripheralResistance(4)=1.318e9
!      ComplianceC(4)=5.967e-10
!c      
!      ResistanceToBloodFlow(5)=2.933e8
!      TotalPeripheralResistance(5)=7.542e8
!      ComplianceC(5)=5.967e-10
!c      
!      ResistanceToBloodFlow(6)=2.933e8
!      TotalPeripheralResistance(6)=7.542e8
!      ComplianceC(6)=5.967e-10
      do i=1,NumberOfBCSets
c
        BoundaryType(i)='wall'
        wallTypeM(i)='slip'
        wallTypeC(i)='slip'
        if(i.eq.1) then
           BoundaryType(i)='inlet'
           inletTypeM(i)='specifiedvelocity'
           inletTypeC(i)='specifiedvelocity'
        elseif(i.eq.2) then
           BoundaryType(i)='outlet'
           outletTypeM(i)='specifiedstaticpressure'
           outletTypeC(i)='specifiedstaticpressure'
        elseif(i.eq.3) then
           BoundaryType(i)='symmetry'
        elseif(i.eq.4) then
           BoundaryType(i)='wall'
           wallTypeM(i)='slip'
           wallTypeC(i)='slip'
        endif
!
c        wallTypeM(i)='noslip'
c        wallTypeC(i)='noslip'
c         if(BoundaryName(i).eq.'leftWall') BoundaryType(i)='wall'
c         if(BoundaryName(i).eq.'rightWall') BoundaryType(i)='wall'
c         if(BoundaryName(i).eq.'lowerWall') BoundaryType(i)='wall'
c         if(BoundaryName(i).eq.'atmosphere') BoundaryType(i)='outlet'
c         if(BoundaryName(i).eq.'defaultFaces')BoundaryType(i)='symmetry'
c         if(i.eq.4) BoundaryType(i)='wall'
        !if(i.eq.1) BoundaryType(i)='wall'
        !if(i.eq.2) BoundaryType(i)='wall'
        !if(i.eq.3) BoundaryType(i)='wall'
        !if(i.eq.4) BoundaryType(i)='wall'
        !if(i.eq.5) BoundaryType(i)='pressurefarfield'
        !if(i.eq.6) BoundaryType(i)='pressurefarfield'
        !if(i.eq.7) BoundaryType(i)='pressurefarfield'
        !if(i.eq.8) BoundaryType(i)='pressurefarfield'
        !if(i.eq.9) BoundaryType(i)='pressurefarfield'
        !if(i.eq.10) BoundaryType(i)='pressurefarfield'
c        if(i.eq.11) BoundaryType(i)='pressurefarfield'
c        if(i.eq.12) BoundaryType(i)='wall'
c        if(i.eq.13) BoundaryType(i)='wall'
c        if(i.eq.14) BoundaryType(i)='symmetry'
c        if(i.eq.4) BoundaryType(i)='wall'
c        if(i.eq.5) BoundaryType(i)='symmetry'
c        if(i.eq.5) BoundaryType(i)='pressurefarfield'
c        if(i.eq.6) BoundaryType(i)='outlet'
c        if(i.eq.7) BoundaryType(i)='wall'
c        if(i.eq.4) BoundaryType(i)='wall'
c        if(i.eq.5) BoundaryType(i)='symmetry'
c        if(i.eq.6) BoundaryType(i)='symmetry'
c         if(BoundaryName(i).eq.'leftWall') then
c           wallTypeM(i)='noslip'
c           wallTypeC(i)='noslip'
c         elseif(BoundaryName(i).eq.'rightWall') then
c           wallTypeM(i)='noslip'
c           wallTypeC(i)='noslip'
c         elseif(BoundaryName(i).eq.'lowerWall') then
c           wallTypeM(i)='noslip'
c           wallTypeC(i)='noslip'
c         elseif(BoundaryName(i).eq.'atmosphere') then
c           outletTypeM(i)='specifiedstaticpressure'
c           outletTypeC(i)='specifiedstaticpressure'
c
         !elseif(i.eq.9.or.i.eq.10) then
         !  inletTypeM(i)='supersonic'
         !  inletTypeC(i)='supersonic'
         !  inletTypeE(i)='supersonic'
         !elseif(i.eq.7.or.i.eq.8) then
         !  outletTypeM(i)='supersonic'
         !  outletTypeC(i)='supersonic'
         !  outletTypeE(i)='supersonic'
c
c        elseif(i.eq.3.or.i.eq.8.or.i.eq.9.or.i.eq.
c     *                      10.or.i.eq.12.or.i.eq.13) then
c
c          wallTypeM(i)='slip'
c          wallTypeC(i)='slip'
c          wallTypeE(i)='vonneumann'
          !inletTypeM(i)='supersonic'
          !inletTypeC(i)='supersonic'
          !inletTypeE(i)='supersonic'
c          inletTypeM(i)='specifiedvelocity'
c          inletTypeC(i)='specifiedvelocity'
c          inletTypeE(i)='specifiedstatictemperature'
c          inletTypeT(i)='specifiedil'
c          inletTurbulenceIntensity(i)=0.05
c          inletTurbulenceLengthScale(i)=2.3 !0.12
c
c        elseif(i.eq.6.or.i.eq.7) then
c
c          inletTypeM(i)='specifiedvelocity'
c          inletTypeC(i)='specifiedvelocity'
c          inletTypeE(i)='specifiedstatictemperature'
!          outletTypeM(i)='supersonic'
!          outletTypeC(i)='supersonic'
!          outletTypeE(i)='supersonic'
c          AverageOutletPressure(i)=0.
c          outletTypeM(i)='fullydeveloped'
c          outletTypeC(i)='fullydeveloped'
c          outletTypeM(i)='specifiedresistance' 
c
c        elseif(i.eq.3) then
c
c          wallTypeM(i)='noslip'
c          wallTypeC(i)='noslip'
c          wallTypeE(i)='dirichlet'
c
!        elseif(i.eq.4) then
!c
!          wallTypeM(i)='noslip'
!          wallTypeC(i)='noslip'
!          wallTypeE(i)='dirichlet'
c
c        elseif(i.eq.4) then
c
c          wallTypeM(i)='noslip'
c          wallTypeC(i)='noslip'
c          wallTypeE(i)='vonneumann'
c
!        elseif(i.eq.5) then
!c
!          wallTypeM(i)='slip'
!          wallTypeC(i)='noslip'
!c
!        elseif(i.eq.6) then
!c
!          wallTypeM(i)='noslip'
!          wallTypeC(i)='noslip'
c
c        endif
c
c----------------------------------------------------------------------------------
c---  Available types for the momentum, continuity, and energy equations
c----------------------------------------------------------------------------------
c
c        BoundaryType(i)='wall'
c        BoundaryType(i)='inlet'
c        BoundaryType(i)='outlet'
c        BoundaryType(i)='symmetry'
c        BoundaryType(i)='periodic'
c        BoundaryType(i)='axis'
c        BoundaryType(i)='pressurefarfield'

c        wallTypeM(i)='slip'
c        wallTypeM(i)='noslip'
c
c        wallTypeC(i)='slip'
c        wallTypeC(i)='noslip'
c
c        wallTypeE(i)='dirichlet'
c        wallTypeE(i)='vonneumann'
c        wallTypeE(i)='robin'
c
c        inletTypeM(i)='specifiedvelocity'
c        inletTypeM(i)='specifiedstaticpressure'      (need to specify velocity direction)(xVeldirection,yVeldirection,zVeldirection)
c        inletTypeM(i)='specifiedstagnationpressure'  (need to specify velocity direction)(xVeldirection,yVeldirection,zVeldirection)
c        inletTypeM(i)='specifiedmassflowrate'        (need to specify velocity direction)(xVeldirection,yVeldirection,zVeldirection)
c        inletTypeM(i)='supersonic'

c        inletTypeC(i)='specifiedvelocity'
c        inletTypeC(i)='specifiedstaticpressure'
c        inletTypeC(i)='specifiedstagnationpressure'
c        inletTypeC(i)='specifiedmassflowrate'
c        inletTypeC(i)='supersonic'
c
c        inletTypeE(i)='specifiedstatictemperature'
c        inletTypeE(i)='specifiedstagnationtemperature'
c        inletTypeE(i)='supersonic'
c
c        outletTypeM(i)='supersonic'
c        outletTypeM(i)='specifiedvelocity'
c        outletTypeM(i)='specifiedmassflowrate'
c        outletTypeM(i)='specifiedstaticpressure'
c        outletTypeM(i)='fullydeveloped'
c        outletTypeM(i)='specifiedaveragestaticpressure'
c        outletTypeM(i)='specifiedresistance'

c        outletTypeC(i)='supersonic'
c        outletTypeC(i)='specifiedvelocity'
c        outletTypeC(i)='specifiedmassflowRate'
c        outletTypeC(i)='specifiedstaticpressure'
c        outletTypeC(i)='fullydeveloped'
c
c        outletTypeE(i)='fullydevelopedenergy'
c        outletTypeE(i)='supersonic'
c
      enddo
c
c--- Set the boundary heat fluxes, Tinfinity, and Hinfinity
c    for use with von Neumann and Robin conditions
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)      
c
          HeatFlux(i,j)=0.
          Hinfinity(i,j)=0.
          Tinfinity(i,j)=0.
c
        enddo
      enddo
c
c--- Set the value of the gas constant for compressible flows
c
      EquationOfState='constant'  !'constant'
      BulkModulus=2.2d9
      RGas=287.05
      GammaGas=1.4
      PrLaminar=0.72
c      Reynolds=5.e6
c      Reynolds=15.e6
      Tinlet=300. !275.3872
      Minlet=0.9
c      Taw=Tinlet*(1.+0.178*(Minlet**2))
      angle=0.*pi/180.
      Uinfinity=1.e-5 !1.05 !Minlet*dsqrt(GammaGas*RGas*Tinlet) !5.e-3 !0.174814656 !0.1123641468 !0.104998   !1.e-4 !*  !0.14142588
      Uinlet=Uinfinity !*dcos(angle)
      Vinlet=-1. !Uinfinity*dsin(angle)
c      Minlet=Uinlet/dsqrt(GammaGas*RGas*Tinlet)
c      Uinfinity=dsqrt(Uinlet**2+Vinlet**2)
c      Minlet=Uinfinity/dsqrt(GammaGas*RGas*Tinlet)
      SoutherlandLambda=1.512041288e-6 !for air 1.512041288e-6
      SoutherlandC=120.  !for air 120.
c
c--- Set the constant density,specific heat, viscosity, and conductivity values
c
      LSoutherLand=.false.
      ConstantSpecificHeat=GammaGas*RGas/(GammaGas-1.)
      ConstantViscosity=SoutherlandLambda*(Tinlet)**1.5/
     *                                (Tinlet+SoutherlandC) !0.0035 !1.7894e-5 !0.001003 !1.7894e-5 !1.8e-5
      ConstantConductivity=ConstantViscosity*
     *                 ConstantSpecificHeat/PrLaminar   !0.0242
      Pinlet=0. !29765.
      RhoInlet=1000. !Pinlet/(RGas*Tinlet) !Reynolds*ConstantViscosity/(1.*Uinlet)
      Rhoinfinity=RhoInlet
c      Pinlet=46040 !RhoInlet*Rgas*Tinlet                !1000. !0.*110.*133.3224
      ConstantDensity=RhoInlet  !for compressible flows it should be set to zero
c      print*,RhoInlet,Pinlet
c      pause
c
c--- Set values for anisotropic diffusion
c
      LanisotropicDiffusion=.false.
      MethodDecomposeSprime=4
c
c--- Set the constant stagnation pressure and temperature (if needed) 
c
      SpecificHeat=ConstantSpecificHeat
      BSpecificHeat=ConstantSpecificHeat
      Conductivity=ConstantConductivity
      BConductivity=ConstantConductivity
c
c--- Set the variable values at the boundaries
c
!      Vinlet=Uinlet*dsin(angle)
!      angle=dasin(Vinlet/Uinfinity)
      constantStagnationPressure=Pinlet/0.59126007
      constantStagnationTemperature=Tinlet/0.86058519
c
c--- Set the far field conditions for pressure far field boundary condition
c
!      call CalculateLengthScale
      do i=1,NumberOfBCSets
          MachFarField(i)=Minlet
          xFlowDirectionFarField(i)=dcos(angle)
          yFlowDirectionFarField(i)=dsin(angle)
          zFlowDirectionFarField(i)=0.
          PressureFarField(i)=Pinlet
          TemperatureFarField(i)=Tinlet
        !if(i.eq.1) then
        !  MachFarField(i)=58.25/dsqrt(GammaGas*RGas*250.19)
        !  xFlowDirectionFarField(i)=1.
        !  yFlowDirectionFarField(i)=0. 
        !  zFlowDirectionFarField(i)=0.
        !  PressureFarField(i)=47485.
        !  TemperatureFarField(i)=250.19
        !elseif(i.eq.2) then
        !  MachFarField(i)=307.13/dsqrt(GammaGas*RGas*240.5)
        !  xFlowDirectionFarField(i)=1.
        !  yFlowDirectionFarField(i)=0.
        !  zFlowDirectionFarField(i)=0.
        !  PressureFarField(i)=32429.
        !  TemperatureFarField(i)=240.5
        !else
        !  MachFarField(i)=Minlet
        !  xFlowDirectionFarField(i)=dcos(angle)
        !  yFlowDirectionFarField(i)=0. !dsin(angle)
        !  zFlowDirectionFarField(i)=dsin(angle)
        !  PressureFarField(i)=Pinlet
        !  TemperatureFarField(i)=Tinlet
        !endif
        SpecificHeatFarField(i)=GammaGas*RGas/(GammaGas-1.)
c
        if(LTurbulentFlow) then
c
          if(LSolveTurbulenceDissipationRate) then 
c
            TKEFarField(i)=1.5*(1.*Uinlet/100.)**2
            TEDFarField(i)=0.09*RhoInlet*TKEFarField(i)**2/(20.*
     *                                  ConstantViscosity)
c
          endif
c
          if(LSolveTurbulenceSpecificDissipationRate) then 
c
            TKEFarField(i)=1.5*(0.9*Uinlet/100.)**2
            TOmegaFarField(i)=RhoInlet*TKEFarField(i)/(8.*
     *                                  ConstantViscosity)
c            TKEFarField(i)=9.e-9*(1.4*RGas*Tinlet)
c            TOmegaFarField(i)=1.e-6*RhoInlet*
c     *               (1.4*RGas*Tinlet)/constantViscosity      !125.*Uinlet/2. 
c
          endif
c
          if(LSolveTurbulentKL) then 
c
            TKEFarField(i)=9.e-9*(1.4*RGas*Tinlet)
            TurbulentKLFarField(i)=0.*1.5589e-6*ConstantViscosity*
     *                             dsqrt(1.4*RGas*Tinlet)/RhoInlet

          endif
c
          if(LSolveModifiedED) then 
c
            MEDFarField(i)=3.*ConstantViscosity/RhoInlet
c
          endif
c
        endif
c
      enddo
c
      GrainSize=5.e-4
c
      uVelocity=Uinlet
      BuVelocity=0.
      vVelocity=0.
      BvVelocity=0.
      wVelocity=0.
      BwVelocity=0.
c      
c      do i=1,NumberOfElements

c        if(zc(i).gt.0.0093) wVelocity(i)=Uinlet
c        if(zc(i).lt.-0.0093) wVelocity(i)=-Uinlet
c        if(xc(i).lt.-0.0093) uvelocity(i)=-Uinlet
c
c      enddo
c
      if(LSolveTurbulenceDissipationRate) then 
c
        TurbulentKE=1.5*(1.*Uinlet/100.)**2
        BTurbulentKE=1.5*(1.*Uinlet/100.)**2
        TurbulentED=0.09*RhoInlet*TurbulentKE**2/(10.*
     *                                  ConstantViscosity)
        BTurbulentED=0.09*RhoInlet*BTurbulentKE**2/(10.*
     *                                  ConstantViscosity)
c        TurbulentED=((0.09)**0.75)*TurbulentKE**1.5/(0.07*2.)
c        BTurbulentED=((0.09)**0.75)*BTurbulentKE**1.5/(0.07*2.)
        TurbulentViscosity=0.09*RhoInlet*TurbulentKE**2/TurbulentED
        BTurbulentViscosity=0.09*RhoInlet*BTurbulentKE**2/BTurbulentED
c
      endif
c
      if(LSolveTurbulenceSpecificDissipationRate) then 
c
        TurbulentKE=1.5*(1.*Uinlet/100.)**2
        TurbulentOmega=RhoInlet*TurbulentKE/(9.*
     *                                  ConstantViscosity)
        BTurbulentKE=1.5*(1.*Uinlet/100.)**2
        BTurbulentOmega=RhoInlet*BTurbulentKE/(9.*
     *                                  ConstantViscosity)
c        TurbulentKE=9.e-9*(1.4*RGas*Tinlet)
c        BTurbulentKE=9.e-9*(1.4*RGas*Tinlet)
c        TurbulentOmega=1.e-6*RhoInlet*
c     *               (1.4*RGas*Tinlet)/constantViscosity      !125.*Uinlet/2. 
c        BTurbulentOmega=1.e-6*RhoInlet*
c     *               (1.4*RGas*Tinlet)/constantViscosity      !125.*Uinlet/2. 
        TurbulentViscosity=RhoInlet*TurbulentKE/
     *                                dmax1(TurbulentOmega,tiny)
        BTurbulentViscosity=RhoInlet*BTurbulentKE/
     *                                dmax1(BTurbulentOmega,tiny)
c
      endif
c
      if(LSolveTurbulentKL) then 
c
        TurbulentKE=9.e-9*(1.4*RGas*Tinlet)
        BTurbulentKE=9.e-9*(1.4*RGas*Tinlet)
        TurbulentKL=1.5589e-6*ConstantViscosity*
     *                             dsqrt(1.4*RGas*Tinlet)/RhoInlet
        BTurbulentKL=1.5589e-6*ConstantViscosity*
     *                             dsqrt(1.4*RGas*Tinlet)/RhoInlet
        TurbulentViscosity=RhoInlet*((0.09)**0.25)*
     *                           TurbulentKL/dsqrt(TurbulentKE)
        BTurbulentViscosity=RhoInlet*((0.09)**0.25)*
     *                    BTurbulentKL/dsqrt(BTurbulentKE)
c
      endif
c
      if(LSolveModifiedED) then 
c
        ModifiedED=3.*ConstantViscosity/RhoInlet
        BModifiedED=3.*ConstantViscosity/RhoInlet
        TurbulentViscosity=ModifiedED*RhoInlet
        BTurbulentViscosity=BModifiedED*RhoInlet
c
      endif
c
c--- Set the temperature
c
      Temperature=Tinlet
      BTemperature=Tinlet
      Pressure=Pinlet  !-0.5*1050.*Uinlet*Uinlet
      BPressure=Pinlet
      Density=Pinlet/(RGas*Tinlet)
      BDensity=Pinlet/(RGas*Tinlet)
      TGamma=0.5
      TReTheta=0.01
      BTGamma=1.
      BTReTheta=0.01
      TfRelaxation=1.
      BTfRelaxation=0.5
      TurbulentV2=twothird*TurbulentKE
      BTurbulentV2=twothird*BTurbulentKE
      TurbulentZeta=twothird
      BTurbulentZeta=twothird
!c
!c--- Set the pressure 
!c
!      do i=1,NumberOfElements
!c
!        Pressure=
!c
!      enddo
!c
!      do i=1,NumberOfBCSets
!        do j=1,NBFaces(i)
!          BPressure=
!        enddo
!      enddo
c
c--- Set the u and v velocity components
c
      AxisOfRotationOriginX=0.
      AxisOfRotationOriginY=0.
      AxisOfRotationOriginZ=0.
c
!      i=2    
!      do j=1,NBFaces(i)      
!         BuVelocity(i,j)=0.
!         BvVelocity(i,j)=0.
!         BwVelocity(i,j)=0.
!c
!      enddo
!      i=2     
!      do j=1,NBFaces(i)      
!         BuVelocity(i,j)=0.
!         BvVelocity(i,j)=0.
!         BwVelocity(i,j)=0.
!c
!      enddo
!      i=3     
!      do j=1,NBFaces(i)      
!         BuVelocity(i,j)=0.
!         BvVelocity(i,j)=0.
!         BwVelocity(i,j)=0.
!c
!      enddo
!      i=4     
!      do j=1,NBFaces(i)      
!         BuVelocity(i,j)=0.
!         BvVelocity(i,j)=0.
!         BwVelocity(i,j)=0.
!c
!      enddo
      
      
      
      
      
      
      
      
      
      
      
      
      
!      i=9
!      do j=1,NBFaces(i)      
!
!        BstagnationPressure(i,j)=constantStagnationPressure
!        xVeldirection(i,j)=dcos(angle)
!        yVeldirection(i,j)=dsin(angle)
!        zVeldirection(i,j)=0.
!        BstagnationTemperature(i,j)=constantStagnationTemperature
!c
!      enddo
!c
!      i=10
!      do j=1,NBFaces(i)      
!
!        BstagnationPressure(i,j)=constantStagnationPressure
!        xVeldirection(i,j)=dcos(angle)
!        yVeldirection(i,j)=dsin(angle)
!        zVeldirection(i,j)=0.
!        BstagnationTemperature(i,j)=constantStagnationTemperature
!c
!      enddo
!c
!      i=1
!      do j=1,NBFaces(i)      
!c
!        BPressure(i,j)=47485.
!        BTemperature(i,j)=250.19
!        BuVelocity(i,j)=58.25
!        BvVelocity(i,j)=0.
!        BwVelocity(i,j)=0.
!        BDensity(i,j)=BPressure(i,j)/(RGas*BTemperature(i,j))
!c
!      enddo
!      i=2
!      do j=1,NBFaces(i)      
!c
!        BPressure(i,j)=32429.
!        BTemperature(i,j)=240.5
!        BuVelocity(i,j)=307.13 
!        BvVelocity(i,j)=0.
!        BwVelocity(i,j)=0.
!        BDensity(i,j)=BPressure(i,j)/(RGas*BTemperature(i,j))
!c
!c
!      enddo
!      i=3
!      do j=1,NBFaces(i)      
c
!        BPressure(i,j)=Pinlet
!        BTemperature(i,j)=600.
!        BuVelocity(i,j)=0.
!        BvVelocity(i,j)=0.
!        BwVelocity(i,j)=0.
c
!      enddo
c
c      do j=1,NBFaces(2)      
c
c        BuVelocity(2,j)=-(BFaceCentroidy(2,j)-AxisOfRotationOriginY)*
c     *             AngularVelocityZ
c        BvVelocity(2,j)=(BFaceCentroidx(2,j)-AxisOfRotationOriginx)*
c     *             AngularVelocityZ
c        BwVelocity(2,j)=0.
c        print*,BFaceAreanx(2,j),BFaceAreany(2,j),BFaceAreanz(2,j)
c        pause
c
c      enddo
c
c      do j=1,NBFaces(3)      
c
c        BuVelocity(3,j)=-(BFaceCentroidy(3,j)-AxisOfRotationOriginY)*
c     *             AngularVelocityZ
c        BvVelocity(3,j)=(BFaceCentroidx(3,j)-AxisOfRotationOriginx)*
c     *             AngularVelocityZ
c        BwVelocity(3,j)=0.
c
c        BPressure(3,j)=Pinlet
c        BuVelocity(3,j)=Uinlet
c        BvVelocity(3,j)=0.
c        BwVelocity(3,j)=0.
c        print*,BFaceAreanx(3,j),BFaceAreany(3,j),BFaceAreanz(3,j)
c        pause
c        BPressure(3,j)=Pinlet
c        BStagnationPressure(1,j)=constantStagnationPressure
c        xVeldirection(1,j)=1.
c        yVeldirection(1,j)=0.
c        zVeldirection(1,j)=0.
c        BTurbulentViscosity(1,j)=0.
c        BModifiedED(1,j)=0.
c        BTurbulentKE(1,j)=0.
c        BTurbulentED(1,j)=0.
c        BTurbulentKL(1,j)=0.
c        BTurbulentOmega(1,j)=5.*Uinlet/0.4322 
c        BTemperature(3,j)=Tinlet
c
c      enddo
c
c      do j=1,NBFaces(3)      
c
c        BuVelocity(3,j)=Uinlet
c        BvVelocity(3,j)=0.
c        BwVelocity(3,j)=0.
c        BTurbulentViscosity(3,j)=0.
c        BModifiedED(3,j)=0.
c        BTurbulentKE(3,j)=0.
c        BTurbulentED(3,j)=0.
c        BTurbulentKL(3,j)=0.
c        BTurbulentOmega(3,j)=5.*Uinlet/0.4322 
c        BTemperature(3,j)=Tinlet
c        BPressure(3,j)=Pinlet
c
c      enddo
c
!      do j=1,NBFaces(4)      
!c
!        xVeldirection(4,j)=1.
!        yVeldirection(4,j)=0.
!        zVeldirection(4,j)=0.
!        BStagnationPressure(4,j)=constantStagnationPressure
!        BStagnationTemperature(4,j)=constantStagnationTemperature
!        BuVelocity(4,j)=Uinlet
!        BvVelocity(4,j)=0.
!        BwVelocity(4,j)=0.
!c
!      enddo
c
c
c--- Set the viscosity, conductivity, density, and specific heat
c
      Viscosity=ConstantViscosity
      BViscosity=ConstantViscosity
c      TurbulentViscosity=ConstantViscosity
c      BTurbulentViscosity=ConstantViscosity
      Conductivity=ConstantConductivity
      BConductivity=ConstantConductivity
      SpecificHeat=ConstantSpecificHeat
      BSpecificHeat=ConstantSpecificHeat
c
c--- Overwrite/set values at the boundaries
c
c      do i=1,NumberOfBCSets
c        do j=1,NBFaces(i)
c          if(i.eq.1)Bpressure(i,j)=0.99303138*constantStagnationPressure
c          if(i.eq.3)Bpressure(i,j)=0.98562651*constantStagnationPressure
c           Bdensity=Bpressure(i,j)/(RGas*Tinlet)
c          if(i.eq.2) BvVelocity(i,j)=0.
c          if(i.eq.1) BStagnationPressure(i,j)=
c     *                     constantStagnationPressure/ 0.93946969
c          if(i.eq.1) xVeldirection(i,j)=-BFaceAreanx(i,j)
c          if(i.eq.1) yVeldirection(i,j)=-BFaceAreany(i,j)
c        enddo
c      enddo
c--- In two dimensional situations a1 and a2 are zero ==> A3=B3=C1=C2=0 C3=1 
c    for every periodic face the user has to specify the rotation vector
c    and rotation angle (positive in the clockwise direction)
c    for every face the user has to specify the corresponding transformed face     
c
c
      LRotationalPeriodicity=.false.
      LTranslationalPeriodicity=.false.
      LcalculateBeta=.false.
      relaxBeta=1.
      BetaIterations=1
      nIterCorrectBeta=30
      LPeriodicImplicit=.true.
c
      if(LRotationalPeriodicity) then
c
        a1Axis(1)=0
        a2Axis(1)=0
        a3Axis(1)=1.
        a1Axis(2)=0
        a2Axis(2)=0
        a3Axis(2)=1.
        PeriodicPair(2)=1
        PeriodicPair(1)=2
        theta(1)=90.
        theta(2)=-90.
c
      elseif(LTranslationalPeriodicity) then
c
        xTranslation(2)=-2.
        yTranslation(2)=0.
        zTranslation(2)=0.
        xTranslation(4)=2.
        yTranslation(4)=0.
        zTranslation(4)=0.
        PeriodicPair(2)=4
        PeriodicPair(4)=2
        periodicBeta=0. !27.015 !0.001028303 ! !1.98679623  ! gives periodicMdot=0.05
        periodicMdot=150.
c
      endif
c
      endif
c----------------------------------------------------------------------------------------
c---  Set parameters for Lambda Euler Lagrange Equatioin
c----------------------------------------------------------------------------------------
c
      if(LSolveLambdaELEEquation) then
c
c--- Maximum number of global iterations
c
        IterMax=5000
c
c--- Assign the monitoring location
c
        xMonitor=2500.
        yMonitor=2500.
        zMonitor=1000.
c
c--- Stop criteria  1: sum of absolute residuals; 2: maximum absolute residual;
c                   3: rms residual             ; 4: normalized residual;
c
        NstopType=2
        maximumResidual=1.e-9
c
c--- Algebraic solvers of basic variables (pbcg, ilu, sor, direct)
c
        ASSolverLambdaELE='ilu'
c
c--- Residual reduction factor of variables while solving their algebraic equations
c
        rrfLambdaELE=0.1
c
c--- maximum number of algebraic solver iterations
c
        ASIterLambdaELE=3
c
c--- Set the relaxation parameters
c
        LRelaxLambdaELE=.true.
c
c--- Assign the underrelaxation factors values
c
        urfLambdaELE=1.
c
c----------------------------------------------------------------------------------
c--- Set the multigrid variables
c----------------------------------------------------------------------------------
        MGType='geometricelement' !'geometricelement' 'algebraic'
c
        LMultigridLambdaELE=.true.
        LPrintMultiGridResiduals=.false.
c
c--- Set whether to print data during multigrid iterations (for testing)
c
        LTestMultiGrid=.true.
        nprintMG=2000
c
c--- Set the variable on which to base the grid agglomoration (could be a scalar)
c
        MGVariable='lambdaele'
c
c--- Set the multigrid cycle to use (v, f, or w)
c
        MultiGridCycleType = 'wcycle'
c
c--- Set the number of pre and post sweep number of iterations and the number of MG cycles
c
        MultiGridpreSweep=1
        MultiGridpostSweep=3
        MaxMultiGridCycles=30
c
c--- Set when to re-agglomorate fine grids to create coarse grid levels
c
        reAgglomorate=50000000
c
c--- Set the minimum number of elements of a coarse level and the maximum number of coarse levels
c
        minNumberofParents=10
        maxNumberofCoarseLevels=10
c
c--- Set when to start applying multigrid (a number >=1)
c
        nIterStartApplyingMG=1
c
c--- Set the multigrid Residual reduction factor
c
        MultiGridrrf=0.1
c----------------------------------------------------------------------------------
c--- Set the gradient variables
c----------------------------------------------------------------------------------
c
c--- Method to calculate the gradient (1 to 6)
c
        MethodCalcGradientLambdaELE=2
        MethodCalcGradientMomentum=2
c
c--- number of iterations to perform with isterative methods (2 to 4)
c
        nIterGradientLambdaELE=2
        nIterGradientMomentum=2
c
c--- The power to use for the inverse distance with the Least Square Gradient (method 6) 
c
        InvDistancePower=1.
c
c--- Set whether or not to limit the gradient 
c
        LimitGradientLambdaELE=.false.
        LimitGradientMomentum=.false.
c
c--- Set the method to use to limit the gradient 
c
        LimitGradientLambdaELEMethod=2
        LimitGradientMomentumMethod=2
c
c--- Set the type of gradient interpolation to face
c    available types: upwind, downwind, average, and averagecorrected 
c
        GradientInterpolationSchemeLambdaELE ='average'
        alpha1Lambda=0.7
        alpha2Lambda=10.
c
c        wallTypeL(i)='dirichlet'     'vonneumann'
c
        do i=1,NumberOfBCSets
c
          if(i.eq.5.or.i.eq.7) then
            BoundaryType(i)='closed'  !ground
          else
            wallTypeL(i)='permeable'  !air around
          endif
c
        enddo
c
c--- Set values for anisotropic diffusion
c
        MethodDecomposeSprime=4
c
c--- Set the constant stagnation pressure and temperature (if needed) 
c
        LambdaELE=0.5
        BLambdaELE=0.
c
      endif      
c
c----------------------------------------------------------------------------------------
c---  Set parameters for the rFields for free surface flows
c----------------------------------------------------------------------------------------
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
c--- Set physical properties related to the fluid used
c
      ConstantDensityrField(1)=998.2
      ConstantViscosityrField(1)=0.001003
      ConstantSpecificHeatrField(1)=4182.
      ConstantConductivityrField(1)=0.6
c
c      ConstantDensityrField(2)=500.
c      ConstantViscosityrField(2)=0.0007003
c      ConstantSpecificHeatrField(2)=1006.43
c      ConstantConductivityrField(2)=0.0242
c
      ConstantDensityrField(2)=1.225
      ConstantViscosityrField(2)=1.7894E-05
      ConstantSpecificHeatrField(2)=1006.43
      ConstantConductivityrField(2)=0.0242
c
c--- Assign the flow field
c
      uVelocity=0.
      vVelocity=0.
      wVelocity=0.
      BuVelocity=0.
      BvVelocity=0.
      BwVelocity=0.
c
c--- Assign name to rFields (fluids)
c
      rFieldName(1)='water'
c      rFieldName(2)='oil'
      rFieldName(2)='air'
c
c--- Declare rFields to solve out of those created
c
      LsolverField(1)=.true.
      LsolverField(2)=.true.
c      LsolverField(3)=.true.
c
c--- Declare which rFields to be solved using the algebraic multigrid solver
c
      LMultigridrField(1)=.false.
      LMultigridrField(2)=.false.
c      LMultigridrField(3)=.false.
c----------------------------------------------------------------------------------
c--- Set the multigrid variables
c----------------------------------------------------------------------------------
c      MGType='geometricelement'
c
c--- Set whether to print data during multigrid iterations (for testing)
c
      LTestMultiGrid=.true.
      LPrintMultiGridResiduals=.false.
c      nprintMG=20000
c
c--- Set the variable on which to base the grid agglomoration (could be a scalar)
c
c      MGVariable=rFieldName(1)
c
c--- Set the multigrid cycle to use (V, F, or W)
c
c      MultiGridCycleType = 'vcycle'
c
c--- Set the number of pre and post sweep number of iterations and the number of MG cycles
c
c      MultiGridpreSweep=1
c      MultiGridpostSweep=3
c      MaxMultiGridCycles=30
c
c--- Set when to re-agglomorate fine grids to create coarse grid levels
c
c      reAgglomorate=50000000
c
c--- Set the minimum number of elements of a coarse level and the maximum number of coarse levels
c
c      minNumberofParents=10
c      maxNumberofCoarseLevels=10
c
c--- Set when to start applying multigrid (a number >=1)
c
c      nIterStartApplyingMG=1
c
c--- Set the multigrid Residual reduction factor
c
c      MultiGridrrf=0.1
c
c--- Algebraic solvers of rFields variables (pbcg, ilu, sor, direct)
c
        ASSolverrField(1)='ilu'
        ASSolverrField(2)='ilu'
c        ASSolverrField(3)='ilu'
c
c--- Residual reduction factor of rFields while solving their algebraic equations
c
        rrfrField(1)=0.001
        rrfrField(2)=0.001
c        rrfrField(3)=0.001
c
c--- maximum number of algebraic solver iterations for rFields
c
        ASIterrField(1)=3
        ASIterrField(2)=3
c        ASIterrField(3)=3
c
c--- Set the relaxation parameters for rFields
c
        LFalseTransientrField(1)=.true.
        LFalseTransientrField(2)=.true.
c        LFalseTransientrField(3)=.true.
        LRelaxrField(1)=.true.
        LRelaxrField(2)=.true.
c        LRelaxrField(3)=.true.
        FalseDtrField(1)=1.e-3
        FalseDtrField(2)=1.e-3
c        FalseDtrField(3)=1.e-2
c
c--- Assign the underrelaxation factor values for rFields
c
        urfrField(1)=0.8
        urfrField(2)=0.8
c        urfrField(3)=0.8
c
c--- Set whether to use upwind, NVF or TVD convection schemes
c
        LNVFrField(1)=.true.
        LNVFrField(2)=.true.
c        LNVFrField(3)=.true.
        LTVDrField(1)=.false.
        LTVDrField(2)=.false.
c        LTVDrField(3)=.false.
c
c--- Set the name of scheme to use for variables
c
        ConvectionSchemerField(1)='stacs'   !sicsam  !hric !stacs
        ConvectionSchemerField(2)='stacs'
c        ConvectionSchemerField(3)='stacs'
c
c      LTVDDensity=.true.
c
c--- Set when to start applying HR schemes
c
c      nIterStartApplyingHR=1
c
c--- Set the name of scheme to use for variables
c
c      ConvectionSchemeDensity='minmod'
c
c--- Set the value of the coefficient by which to bleed the HR scheme with upwind scheme 
c    (0= no bleeding, 1= upwind) 
c
        BleedrField(1)=0.
        BleedrField(2)=0.
c        BleedrField(3)=0.
c
c----------------------------------------------------------------------------------
c--- Set the transient variables
c----------------------------------------------------------------------------------
c
c--- Set the scheme to use (euler, cranknicolson1 (2 steps), 
c                cranknicolson2,  adamsmoulton,tics1.75,tics2.5)
c
      TransientSchemerField='cranknicolson1'
c
c--- Set the step size
c
      dt=1.d-2
c
c--- Set the maximum time for which the solution to be computed
c
      timemax=30.
c
c--- Set the time to start computations
c
      time=0.
c
c--- Set the number of time steps tp print results
c
c      nTimePrint=10
c
c----------------------------------------------------------------------------------
c--- Set the gradient variables for the rFields
c----------------------------------------------------------------------------------
c
c--- Method to calculate the gradient (1 to 6)
c
        MethodCalcGradientrField(1)=2
        MethodCalcGradientrField(2)=2
c        MethodCalcGradientrField(3)=2
        InvDistancePower=1.
c
c--- number of iterations to perform with isterative methods (2 to 4)
c
        nIterGradientrField(1)=2
        nIterGradientrField(2)=2
c        nIterGradientrField(3)=2
c
c--- Set whether or not to limit the gradient 
c
        LimitGradientrField(1)=.false.
        LimitGradientrField(2)=.false.
c        LimitGradientrField(3)=.false.
c
c--- Set the method to use to limit the gradient 
c
        LimitGradientrFieldMethod(1)=2
        LimitGradientrFieldMethod(2)=2
c        LimitGradientrFieldMethod(3)=2
c
c--- Set the type of gradient interpolation to face
c    should be averagecorrected 
c
        GradientInterpolationSchemerField(1)='average'
        GradientInterpolationSchemerField(2)='average'
c        GradientInterpolationSchemerField(3)='average'
c
c--- Assign the types of boundary conditions
c
c
         do i=1,NumberOfBCSets
c
           if(i.eq.1) then
c              
             inletTypeR(i,1)='specifiedvaluerfield' !'dirichletrfield'
             inletTypeR(i,2)='specifiedvaluerfield' !'dirichletrfield'
c            wallTypeR(i,3)='vonneumannrfield' !'dirichletrfield'
           elseif(i.eq.2) then
            outletTypeR(i,1)='fullydevelopedrfield' !'dirichletrfield'
            outletTypeR(i,2)='fullydevelopedrfield' !'dirichletrfield'
c            outletTypeR(i,3)='fullydevelopedrfield' !'dirichletrfield'
           elseif(i.eq.4) then
             wallTypeR(i,1)='vonneumannrfield' !'dirichletrfield'
             wallTypeR(i,2)='vonneumannrfield' !'dirichletrfield'
c            wallTypeR(i,3)='vonneumannrfield' !'dirichletrfield'
          endif
          
c          wallTypeR(i,3)='vonneumannrfield' !'dirichletrfield'
c         if(i.eq.1) inletTypeR(i,1)='specifiedvaluerfield'
c         if(i.eq.6) inletTypeR(i,1)='specifiedvaluerfield'
c         if(i.eq.2) outletTypeR(i,1)='fullydevelopedrfield'
c         if(i.eq.3) outletTypeR(i,1)='fullydevelopedrfield'
c         if(i.eq.4) wallTypeR(i,1)='vonneumannrfield'
c         if(i.eq.5) wallTypeR(i,1)='vonneumannrfield'
c
c           if(BoundaryName(i).eq.'leftWall')
c     *                           wallTypeR(i,1)='vonneumannrfield'
c           if(BoundaryName(i).eq.'rightWall')
c     *                           wallTypeR(i,1)='vonneumannrfield'
c           if(BoundaryName(i).eq.'lowerWall')
c     *                           wallTypeR(i,1)='vonneumannrfield'
c           if(BoundaryName(i).eq.'atmosphere') 
c     *                           outletTypeR(i,1)='fullydevelopedrfield'
c
         enddo
c
c--- Set the boundary heat fluxes, Tinfinity, and Hinfinity
c    for use with von Neumann and Robin conditions for rFields
c
c
c--- Set the boundary and initial values of rFields
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)      
c
          if(i.eq.1) then
            BrField(i,j,1)=1.
            BrField(i,j,2)=0.
          else
            BrField(i,j,1)=0.
            BrField(i,j,2)=1.
          endif
c            print*,BFaceCentroidx(i,j),BFaceCentroidy(i,j),
c     *                             BFaceCentroidz(i,j)
!          if(BFaceCentroidz(i,j).lt.-onethird/20.) then
!c          if(BFaceCentroidz(i,j).lt.0.) then
!            BrField(i,j,1)=0.
!            BrField(i,j,2)=1.
!            BrField(i,j,3)=0.
!          elseif(BFaceCentroidz(i,j).ge.-onethird/20..and.
!     *                    BFaceCentroidz(i,j).le.onethird/20.) then
!            BrField(i,j,1)=0.
!            BrField(i,j,2)=0.
!            BrField(i,j,3)=1.
!          else
!            BrField(i,j,1)=1.
!            BrField(i,j,2)=0.
!            BrField(i,j,3)=0.
!          endif
          
c          if(BFaceCentroidx(i,j).le.(0.1461).and.
c     *             BFaceCentroidy(i,j).le.(0.292)) then
c            BrField(i,j,1)=1.
c          endif
c
        enddo
      enddo
        do i=1,NumberOfElements
c
            rField(i,1)=0.
            rField(i,2)=1.
c
!c          if(zc(i).lt.0.) then
!          if(zc(i).lt.-onethird/20.) then
!            rField(i,1)=0.
!            rField(i,2)=1.
!            rField(i,3)=0.
!          elseif(zc(i).ge.-onethird/20..and.
!     *                          zc(i).le.onethird/20.) then
!            rField(i,1)=0.
!            rField(i,2)=0.
!            rField(i,3)=1.
!          else
!            rField(i,1)=1.
!            rField(i,2)=0.
!            rField(i,3)=0.
!          endif
!c          rField(i,1)=0.
!c          if(xc(i).le.(0.1461).and.yc(i).le.(0.292)) then
!c            rField(i,1)=1.
!c          endif
c
        enddo
!c
!        do i=1,NumberOfBCSets
!          do j=1,NBFaces(i)
!c
!            BrField(i,j,2)=1.-BrField(i,j,1)
!c
!          enddo
!        enddo
!c
!        do i=1,NumberOfElements
!c
!            rField(i,2)=1.-rField(i,1)
!c
!        enddo
!c
        uVelocity=0.
        vVelocity=0.
        wVelocity=0.
        BuVelocity=0.
        BvVelocity=0.
        BwVelocity=0.
c
c--- Set the pressure 
c
      do i=1,NumberOfElements
c
          Pressure(i)=-(rField(i,1)*ConstantDensityrField(1)+
     *                   rField(i,2)*ConstantDensityrField(2))*
     *                   Gravityy*(2.-yc(i))
c          Pressure(i)=-ConstantDensityrField(2)*GravityY*(0.584-yc(i))
c          if(yc(i).lt.0.292.and.xc(i).le.0.1461)
c     *      Pressure(i)=-0.292*ConstantDensityrField(2)*GravityY-  
c     *                   ConstantDensityrField(1)*GravityY*(0.292-yc(i))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BPressure(i,j)=-(BrField(i,j,1)*ConstantDensityrField(1)+
     *                   BrField(i,j,2)*ConstantDensityrField(2))*
     *                   Gravityy*(2.-BFaceCentroidy(i,j))
            if(i.eq.2) BPressure(i,j)=0.
            if(i.eq.1) BvVelocity(i,j)=-1.
!          BPressure(i,j)=-gravityY*
!     *          ConstantDensityrField(2)*(0.584-BFaceCentroidy(i,j))
!          if(BFaceCentroidx(i,j).le.(0.1461).and.
!     *             BFaceCentroidy(i,j).le.(0.292))  
!     *    BPressure(i,j)=-0.292*ConstantDensityrField(2)*GravityY-
!     *     gravityY*ConstantDensityrField(1)*(0.292-BFaceCentroidy(i,j))
!c
        enddo
      enddo
c


!        do i=1,NumberOfBCSets
!          do j=1,NBFaces(i)      
!c
!            rFieldFlux(i,j,1)=0.
!            rFieldFlux(i,j,2)=0.
!            rFieldConvectionCoefficient(i,j,1)=0.
!            rFieldConvectionCoefficient(i,j,2)=0.
!            rFieldPhiInfinity(i,j,1)=0.
!            rFieldPhiInfinity(i,j,2)=0.
!c
!          enddo
!        enddo
!c
!c--- Set the boundary and initial values of rFields
!c
!      do i=1,NumberOfElements
!          if(yc(i).gt.0.5) then
!            rField(i,1)=0.
!            rField(i,2)=1.
!          else
!            rField(i,1)=1.
!            rField(i,2)=0.
!          endif
!      enddo
!c      
!      do i=1,NumberOfBCSets
!        do j=1,NBFaces(i)      
!c
!          k=NBFaceOwner(i,j)
!c
!          if(i.eq.1) then
!            BrField(i,j,1)=1.
!            BrField(i,j,2)=0.
!          elseif(i.eq.6) then
!            BrField(i,j,1)=0.
!            BrField(i,j,2)=1.
!          elseif(i.eq.2) then
!            if(BFaceCentroidy(i,j).gt.0.5) then
!              BrField(i,j,1)=0.
!              BrField(i,j,2)=1.
!            else
!              BrField(i,j,1)=1.
!              BrField(i,j,2)=0.
!            endif
!          elseif(i.eq.3) then
!            BrField(i,j,1)=0.
!            BrField(i,j,2)=1.
!          elseif(i.eq.4.or.i.eq.5) then
!            BrField(i,j,1)=1.
!            BrField(i,j,2)=0.
!          elseif(i.eq.7) then
!            BrField(i,j,1)=rField(k,1)
!            BrField(i,j,2)=rField(k,2)
!          endif
!c
!        enddo
!      enddo
!c
!        do i=1,NumberOfElements
!c
!          if(yc(i).gt.0.5) then
!            Pressure(i)=-ConstantDensityrField(2)*GravityY*(1.-yc(i))
!          else
!            Pressure(i)=-0.5*ConstantDensityrField(2)*GravityY-
!     *                     ConstantDensityrField(1)*GravityY*(0.5-yc(i))
!          endif
!c
!        enddo
!c
!        uVelocity=1.05
!        vVelocity=1.e-5
!        wVelocity=1.e-5
!c        
!        do i=1,NumberOfBCSets
!          do j=1,NBFaces(i)
!c
!            k=NBFaceOwner(i,j)
!c
!            if(i.eq.1) then
!              BPressure(i,j)=Pressure(k)
!              BuVelocity(i,j)=1.05
!              BvVelocity(i,j)=0.
!              BwVelocity(i,j)=0.
!              xVeldirection(i,j)=1.
!              yVeldirection(i,j)=0.
!              zVeldirection(i,j)=0.
!              BstagnationPressure(i,j)=BPressure(i,j)+
!     *               .5*ConstantDensityrField(1)*(BuVelocity(i,j)**2)
!            elseif(i.eq.6) then
!              BPressure(i,j)=Pressure(k)
!              BuVelocity(i,j)=1.05
!              BvVelocity(i,j)=0.
!              BwVelocity(i,j)=0.
!c
!            elseif(i.eq.2) then
!              BPressure(i,j)=Pressure(k)
!              BuVelocity(i,j)=1.05
!              BvVelocity(i,j)=0.
!              BwVelocity(i,j)=0.
!c
!            elseif(i.eq.3) then
!              BPressure(i,j)=0.
!              BuVelocity(i,j)=1.05
!              BvVelocity(i,j)=0.
!              BwVelocity(i,j)=0.
!c
!            elseif(i.eq.4.or.i.eq.5) then
!              BPressure(i,j)=Pressure(k)
!              BuVelocity(i,j)=uVelocity(k)
!              BvVelocity(i,j)=0.
!              BwVelocity(i,j)=0.
!c
!            elseif(i.eq.7) then
!              BPressure(i,j)=Pressure(k)
!              BuVelocity(i,j)=uVelocity(k)
!              BvVelocity(i,j)=vVelocity(k)
!              BwVelocity(i,j)=wVelocity(k)
!c
!            endif
!c
!          enddo
!        enddo
!c
      endif
c
c----------------------------------------------------------------------------------------
c---  Set parameters for the additional scalars to be solved
c----------------------------------------------------------------------------------------
c
      if(NumberOfScalarsToSolve.gt.0) then

      IterMax=200

      IMonitor=NumberOfElements/2
      NstopType=2
      maximumResidual=1.e-10
c
      xMonitor=0.
      yMonitor=0.
      zMonitor=0.
c
      Density=7000. !Pinlet/(RGas*Tinlet)
      Bdensity=7000. !Pinlet/(RGas*Tinlet)
c
c--- Assign name to scalars
c
      ScalarName(1)='scalar1'
c      ScalarName(2)='scalar2'
c
c--- Declare scalars to solve out of those created
c
      LsolveScalar(1)=.true.
c      LsolveScalar(2)=.false.
c
c--- Declare which scalars to be solved using the algebraic multigrid solver
c
      LMultigridScalar(1)=.true.
c      LMultigridScalar(2)=.true.
c
c--- Set the variable on which to base the grid agglomoration (could be a scalar)
c
      MGVariable='scalar1'

      LTestMultiGrid=.true.
      LPrintMultiGridResiduals=.false.
      nprintMG=2000
c
c--- Set the multigrid cycle to use (V, F, or W)
c
      MultiGridCycleType = 'wcycle'
c
c--- Set the number of pre and post sweep number of iterations and the number of MG cycles
c
      MultiGridpreSweep=1
      MultiGridpostSweep=3
      MaxMultiGridCycles=10
c
c--- Set when to re-agglomorate fine grids to create coarse grid levels
c
      reAgglomorate=50000
c
c--- Set the minimum number of elements of a coarse level and the maximum number of coarse levels
c
      minNumberofParents=10
      maxNumberofCoarseLevels=20
c
c--- Set when to start applying multigrid (a number >=1)
c
      nIterStartApplyingMG=1
c
c--- Set the multigrid Residual reduction factor
c
      MultiGridrrf=0.1
      MGType='algebraic'
c
c--- Algebraic solvers of scalar variables (pbcg, ilu, sor, direct)
c
        ASSolverScalar(1)='ilu'
c        ASSolverScalar(2)='ilu'
c
c--- Residual reduction factor of scalar while solving their algebraic equations
c
        rrfScalar(1)=0.3
c        rrfScalar(2)=0.3
c
c--- maximum number of algebraic solver iterations for additional scalars
c
        ASIterScalar(1)=10
c        ASIterScalar(2)=10
c
c--- Set the relaxation parameters for the additional scalars
c
        LRelaxScalar(1)=.true.
c        LRelaxScalar(2)=.true.
c
c--- Assign the underrelaxation factor values for the additional scalars
c
        urfScalar(1)=1.
c        urfScalar(2)=1.
c
c--- Set whether to use upwind, NVF or TVD convection schemes
c
        LNVFScalar(1)=.false.
c        LNVFScalar(2)=.false.
        LTVDScalar(1)=.false.
c        LTVDScalar(2)=.false.
c
c--- Set the name of scheme to use for variables
c
        ConvectionSchemeScalar(1)='minmod'
c        ConvectionSchemeScalar(2)='minmod'
c
c--- Set the value of the coefficient by which to bleed the HR scheme with upwind scheme 
c    (0= no bleeding, 1= upwind) 
c
        BleedScalar(1)=0.
c        BleedScalar(2)=0.
c
c----------------------------------------------------------------------------------
c--- Set the gradient variables for the additional scalars
c----------------------------------------------------------------------------------
c
c--- Method to calculate the gradient (1 to 6)
c
        MethodCalcGradientScalar(1)=2
c        MethodCalcGradientScalar(2)=2
        InvDistancePower=1.
c
c--- number of iterations to perform with isterative methods (2 to 4)
c
        nIterGradientScalar(1)=2
c        nIterGradientScalar(2)=2
c
c--- Set whether or not to limit the gradient 
c
        LimitGradientScalar(1)=.false.
c        LimitGradientScalar(2)=.false.
c
c--- Set the method to use to limit the gradient 
c
        LimitGradientScalarMethod(1)=1
c        LimitGradientScalarMethod(2)=2
c
c--- Set the type of gradient interpolation to face
c    available types: upwind, downwind, average, and averagecorrected 
c
        GradientInterpolationSchemeScalar(1)='average'
c        GradientInterpolationSchemeScalar(2)='average'
c
c--- Assign the types of boundary conditions
c
c
         do i=1,NumberOfBCSets
c
c          if(i.eq.1) then
            BoundaryType(i)='wall'
            wallTypeS(i,1)='dirichletscalar'
c          else
c            BoundaryType(i)='wall'
c            wallTypeS(i,1)='dirichletscalar'
c          endif
c          if(i.eq.2.or.i.eq.4) BoundaryType(i)='wall'
c            
c           BoundaryType(i)='wall'
c          if(i.eq.1) then
c            wallTypeS(i,1)='dirichletscalar'
c            wallTypeS(i,2)='dirichletscalar'
c            wallTypeScalar(i,j,2)='vonneumannscalar'
c          elseif(i.eq.2) then
c            wallTypeS(i,1)='dirichletscalar'
c            wallTypeScalar(i,j,2)='dirichletscalar'
c          elseif(i.eq.3) then
c            wallTypeS(i,1)='dirichletscalar'
c            wallTypeScalar(i,j,2)='vonneumannscalar'
c          elseif(i.eq.4) then
c            wallTypeS(i,1)='dirichletscalar'
c            wallTypeScalar(i,j,2)='dirichletscalar'
c          endif
c
        enddo
c
c---- Location of point sources and source values
c
!      xLocationOfPointSource(1)=2
!      yLocationOfPointSource(1)=2
!      zLocationOfPointSource(1)=2
!c
!      xLocationOfPointSource(2)=3
!      yLocationOfPointSource(2)=3
!      zLocationOfPointSource(2)=2
!c
!      xLocationOfPointSource(3)=2
!      yLocationOfPointSource(3)=2
!      zLocationOfPointSource(3)=3
!c
!      xLocationOfPointSource(4)=3
!      yLocationOfPointSource(4)=3
!      zLocationOfPointSource(4)=3
!c
!      xLocationOfPointSource(5)=.5
!      yLocationOfPointSource(5)=.25
!      zLocationOfPointSource(5)=.75
!c
!      xLocationOfPointSource(6)=.25
!      yLocationOfPointSource(6)=.5
!      zLocationOfPointSource(6)=.75
!c
!      xLocationOfPointSource(7)=.75
!      yLocationOfPointSource(7)=.75
!      zLocationOfPointSource(7)=.75
!c
!      xLocationOfPointSource(8)=.75
!      yLocationOfPointSource(8)=.75
!      zLocationOfPointSource(8)=.25
!c
      !SbPointSourceScalar(1,1)=10000.
      !SbPointSourceScalar(2,1)=10000.
      !SbPointSourceScalar(3,1)=10000.
      !SbPointSourceScalar(4,1)=10000.
!      SbPointSourceScalar(5,1)=1.
!      SbPointSourceScalar(6,1)=1.
!      SbPointSourceScalar(7,1)=10.
!      SbPointSourceScalar(8,1)=10.
!      SbPointSourceScalar(1,2)=1000.
!      SbPointSourceScalar(2,2)=1000.
!      SbPointSourceScalar(3,2)=1000.
!      SbPointSourceScalar(4,2)=1000.
c
c--- Set values for anisotropic diffusion
c
      LanisotropicDiffusion=.false.
      MethodDecomposeSprime=4
c
!      do i=1,NumberOfElements
!c
!        DiffusionCoefficient11(i,1)=yc(i)**2+zc(i)**2+1.
!        DiffusionCoefficient12(i,1)=-xc(i)*yc(i)
!        DiffusionCoefficient13(i,1)=-xc(i)*zc(i)
!        DiffusionCoefficient22(i,1)=xc(i)**2+zc(i)**2+1.
!        DiffusionCoefficient23(i,1)=-yc(i)*zc(i)
!        DiffusionCoefficient33(i,1)=xc(i)**2+yc(i)**2+1.
!c
!      enddo
c
!      do i=1,NumberOfBCSets
!        do j=1,NBFaces(i)      
!
!          BDiffusionCoefficient11(i,j,1)=
!     *        BFacecentroidy(i,j)**2+BFacecentroidz(i,j)**2+1.
!          BDiffusionCoefficient12(i,j,1)=
!     *       -BFacecentroidx(i,j)*BFacecentroidy(i,j)
!          BDiffusionCoefficient13(i,j,1)=
!     *       -BFacecentroidx(i,j)*BFacecentroidz(i,j)
!          BDiffusionCoefficient22(i,j,1)=
!     *        BFacecentroidx(i,j)**2+BFacecentroidz(i,j)**2+1.
!          BDiffusionCoefficient23(i,j,1)=
!     *       -BFacecentroidy(i,j)*BFacecentroidz(i,j)
!          BDiffusionCoefficient33(i,j,1)=
!     *        BFacecentroidx(i,j)**2+BFacecentroidy(i,j)**2+1.
!c
!        enddo
!      enddo
c
      Scalar=300.
      BScalar=300.
c
c--- Set boundary values for scalars
c            if(i.le.4.or.i.eq.6) then
      do i=1,NumberOfBCSets
        if(BoundaryType(i).eq.'wall') then
          do j=1,NBFaces(i)      
c
            if(i.eq.1) then
              BScalar(i,j,1)=1600.
            elseif(i.eq.2) then
              BScalar(i,j,1)=800.
            elseif(i.eq.3) then
              BScalar(i,j,1)=400.
            elseif(i.eq.4) then
              BScalar(i,j,1)=3200.
            elseif(i.eq.5) then
              BScalar(i,j,1)=5000.
            elseif(i.eq.6) then
              BScalar(i,j,1)=6000.
            elseif(i.eq.7) then
              BScalar(i,j,1)=7000.
            endif
c
          enddo
        endif
      enddo
c
c
c--- In two dimensional situations a1 and a2 are zero ==> A3=B3=C1=C2=0 C3=1 
c    for every periodic face the user has to specify the rotation vector
c    and rotation angle (positive in the clockwise direction)
c    for every face the user has to specify the corresponding transformed face     
c
      LRotationalPeriodicity=.false.
      LTranslationalPeriodicity=.false.
      LcalculateBeta=.false.
      relaxBeta=1.
c
      if(LRotationalPeriodicity) then
c
        a1Axis(2)=0
        a2Axis(2)=0
        a3Axis(2)=1.
        a1Axis(4)=0
        a2Axis(4)=0
        a3Axis(4)=1.
        PeriodicPair(2)=4
        PeriodicPair(4)=2
        theta(2)=90.
        theta(4)=-90.
c
      elseif(LTranslationalPeriodicity) then
c
        xTranslation(2)=-10.
        yTranslation(2)=0.
        zTranslation(2)=0.
        xTranslation(4)=10.
        yTranslation(4)=0.
        zTranslation(4)=0.
        PeriodicPair(2)=4
        PeriodicPair(4)=2
        periodicBeta=12.28 
        periodicMdot=0.
c
      endif
c
c--- Set the constant diffusion coefficient values for the additional scalars
c
        ConstantDiffusionCoefficientScalar(1)=400.
c        ConstantDiffusionCoefficientScalar(2)=1
c
c--- Set the constant specific-heat like coefficient for the additional scalars
c
        ConstantSpecificHeatScalar(1)=4000.
c        ConstantSpecificHeatScalar(2)=1.
c
c--- Set the boundary heat fluxes, Tinfinity, and Hinfinity
c    for use with von Neumann and Robin conditions for scalar
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)      
c
            ScalarFlux(i,j,1)=0.
c            ScalarFlux(i,j,2)=0.
c            ScalarConvectionCoefficient(i,j,1)=0.
c            ScalarConvectionCoefficient(i,j,2)=0.
c            ScalarPhiInfinity(i,j,1)=0.
c            ScalarPhiInfinity(i,j,2)=0.
c
          enddo
        enddo
c
c--- Set the diffusion and specific heat like parameters for scalars
c
        do i=1,NumberOfElements
c
          DiffusionCoefficient(i,1)=
     *                  ConstantDiffusionCoefficientScalar(1)
c          DiffusionCoefficient(i,2)=
c     *                  ConstantDiffusionCoefficientScalar(2)
c
          SpecificHeatScalar(i,1)=ConstantSpecificHeatScalar(1)
c          SpecificHeatScalar(i,2)=ConstantSpecificHeatScalar(2)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BDiffusionCoefficient(i,j,1)=
     *                 ConstantDiffusionCoefficientScalar(1)
c            BDiffusionCoefficient(i,j,2)=
c     *                 ConstantDiffusionCoefficientScalar(2)
c
            BSpecificHeatScalar(i,j,1)=
     *                 ConstantSpecificHeatScalar(1)
c            BSpecificHeatScalar(i,j,2)=
c     *                 ConstantSpecificHeatScalar(2)
c
          enddo
        enddo
c
c--- Set the boundary and initial values of scalars
c
!        do i=1,NumberOfElements
!c
!          scalar(i,1)=1.+yc(i)+zc(i)
!c
!        enddo
!        write(99,*) scalar
!        do i=1,NumberOfElements
!c
!          scalar(i,1)=1.+xc(i)+zc(i)
!c
!        enddo
!        write(99,*) scalar
!        do i=1,NumberOfElements
!c
!          scalar(i,1)=1.+xc(i)+yc(i)
!c
!        enddo
!        write(99,*) scalar
!        do i=1,NumberOfElements
!c
!          Scalar(i,1)=xc(i)+yc(i)+zc(i)+
!     *               xc(i)*yc(i)+xc(i)*zc(i)+yc(i)*zc(i)
!c          Scalar(i,2)=1.
!c
!        enddo
!        write(99,*) scalar     
!        scalar=0.5
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            if(i.eq.1) then 
!              BScalar(i,j,1)=BFaceCentroidx(i,j)+
!     *        BFaceCentroidy(i,j)+
!     *        BFaceCentroidz(i,j)+
!     *               BFaceCentroidx(i,j)*BFaceCentroidy(i,j)+
!     *               BFaceCentroidx(i,j)*BFaceCentroidz(i,j)+
!     *               BFaceCentroidy(i,j)*BFaceCentroidz(i,j)
c              BScalar(i,j,2)=3000.
c            elseif(i.eq.2) then 
c              BScalar(i,j,1)=273.
c              BScalar(i,j,2)=3000.
            elseif(i.eq.3) then 
c              BScalar(i,j,1)=1273.
c              BScalar(i,j,2)=3000.
c            elseif(i.eq.4) then 
c              BScalar(i,j,1)=273.
c              BScalar(i,j,2)=4000.
            endif
c
          enddo
        enddo
c
      endif
c
      return
c
c***************************************************************************************
      entry MakeChangesDuringComputations
c***************************************************************************************
c
      
      
      
      
      if(nIter.eq.5000000000) then
!
        LYapCorrection=.true.
!!
!          allocate(fmuCoefficient(NumberOfElements))
!          allocate(BfmuCoefficient(NumberOfBCSets,NBFacesMax))
!          allocate(f1Coefficient(NumberOfElements))
!          allocate(f2Coefficient(NumberOfElements))
!          allocate(LTKE(NumberOfElements))
!          allocate(LTED(NumberOfElements))
! 
!          fmuCoefficient=1.
!          BfmuCoefficient=1.
!          f1Coefficient=1.
!          f2Coefficient=1.
!          LTKE=0.
!          LTED=0.
!
        allocate(fmuCoefficient(NumberOfElements))
        allocate(BfmuCoefficient(NumberOfBCSets,NBFacesMax))
        allocate(f1Coefficient(NumberOfElements))
        allocate(f2Coefficient(NumberOfElements))
        allocate(LTKE(NumberOfElements))
        allocate(LTED(NumberOfElements))
        f1Coefficient=1.
c
        allocate(sqrtTurbulentKE(NumberOfElements))
        allocate(BsqrtTurbulentKE(NumberOfBCSets,NBFacesMax))
        allocate(sqrtTKEGradx(NumberOfElements))
        allocate(sqrtTKEGrady(NumberOfElements))
        allocate(sqrtTKEGradz(NumberOfElements))
        allocate(BsqrtTKEGradx(NumberOfBCSets,NBFacesMax))
        allocate(BsqrtTKEGrady(NumberOfBCSets,NBFacesMax))
        allocate(BsqrtTKEGradz(NumberOfBCSets,NBFacesMax))
        allocate(SRateGradx(NumberOfElements))
        allocate(SRateGrady(NumberOfElements))
        allocate(SRateGradz(NumberOfElements))
        allocate(BSRateGradx(NumberOfBCSets,NBFacesMax))
        allocate(BSRateGrady(NumberOfBCSets,NBFacesMax))
        allocate(BSRateGradz(NumberOfBCSets,NBFacesMax))
c
        allocate(uVelGrad2x(NumberOfElements))         
        allocate(BuVelGrad2x(NumberOfBCSets,NBFacesMax))
        allocate(uVelGrad2y(NumberOfElements))         
        allocate(BuVelGrad2y(NumberOfBCSets,NBFacesMax))
        allocate(uVelGrad2z(NumberOfElements))         
        allocate(BuVelGrad2z(NumberOfBCSets,NBFacesMax))
        allocate(vVelGrad2x(NumberOfElements))         
        allocate(BvVelGrad2x(NumberOfBCSets,NBFacesMax))
        allocate(vVelGrad2y(NumberOfElements))         
        allocate(BvVelGrad2y(NumberOfBCSets,NBFacesMax))
        allocate(vVelGrad2z(NumberOfElements))         
        allocate(BvVelGrad2z(NumberOfBCSets,NBFacesMax))
        allocate(wVelGrad2x(NumberOfElements))         
        allocate(BwVelGrad2x(NumberOfBCSets,NBFacesMax))
        allocate(wVelGrad2y(NumberOfElements))         
        allocate(BwVelGrad2y(NumberOfBCSets,NBFacesMax))
        allocate(wVelGrad2z(NumberOfElements))         
        allocate(BwVelGrad2z(NumberOfBCSets,NBFacesMax))
        allocate(uvVelGradxy(NumberOfElements))         
        allocate(BuvVelGradxy(NumberOfBCSets,NBFacesMax))
        allocate(uwVelGradxz(NumberOfElements))         
        allocate(BuwVelGradxz(NumberOfBCSets,NBFacesMax))
c
!        allocate(TurbulentED(NumberOfElements))
!        allocate(BTurbulentED(NumberOfBCSets,NBFacesMax))
!        allocate(TEDGradx(NumberOfElements))
!        allocate(BTEDGradx(NumberOfBCSets,NBFacesMax))
!        allocate(TEDGrady(NumberOfElements))
!        allocate(BTEDGrady(NumberOfBCSets,NBFacesMax))
!        allocate(TEDGradz(NumberOfElements))
!        allocate(BTEDGradz(NumberOfBCSets,NBFacesMax))
!        allocate(TEDGradfx(NIFaces))
!        allocate(TEDGradfy(NIFaces))
!        allocate(TEDGradfz(NIFaces))
!        allocate(ScTED(NumberOfElements))
!        allocate(SbTED(NumberOfElements))
c
!        if(LUnsteady) then
!c
!          allocate(TurbulentEDOld(NumberOfElements))
!          allocate(BTurbulentEDOld(NumberOfBCSets,NBFacesMax))
!          allocate(TurbulentEDOldOld(NumberOfElements))
!          allocate(BTurbulentEDOldOld(NumberOfBCSets,NBFacesMax))
!c
!        endif
!c
        TurbulenceModel='kepsilonhishida'
!        LSolveTurbulenceDissipationRate=.true.
!        LSolveTurbulenceSpecificDissipationRate=.false.
!c
        call SetTurbulenceModelCoefficients
        WallTreatment='lowreynoldsnumber'
!        do i=1,NumberOfElements
!c
!          TurbulentED(i)=0.09*TurbulentKE(i)*TurbulentOmega(i)
!c
!        enddo
!c
!        do i=1,NumberOfBCSets
!          do j=1,NBFaces(i)
!c
!            BTurbulentED(i,j)=
!     *           0.09*BTurbulentKE(i,j)*BTurbulentOmega(i,j)
!c
!          enddo
!        enddo
!c
      endif
c
!        if(nIter.eq.1) then
!        TurbulenceModel='kepsilon'
!        WallTreatment='wallfunctions'
!        endif
!        if(nIter.eq.1000) then
!        TurbulenceModel='kepsilonsharma'
!        WallTreatment='lowreynoldsnumber'
!        endif
c       if(nIter.eq.100) WallTreatment='lowreynoldsnumber'
c       if(nIter.gt.3000) FalseDtMomentum=5.e-4
c       if(nIter.gt.4000) FalseDtMomentum=1.e-3
c       if(nIter.gt.4500) FalseDtMomentum=5.e-3
c       if(nIter.gt.5000) FalseDtMomentum=1.e-2
c       if(nIter.gt.6000) FalseDtMomentum=5.e-2
c       if(nIter.gt.7000) FalseDtMomentum=0.1

c      if(time.gt.1d-3) then
c        LTVDMomentum=.true.
c        LTVDDensity=.true.
c        TransientScheme='cranknicolson1'
c        LTVDrField(1)=.true.
c        TransientSchemerField='cranknicolson1'
c      endif


      return
c
c***************************************************************************************
      entry SetTimeStepValue
c***************************************************************************************
c
      dt=1.d-3
c
      return
c
c***************************************************************************************
      entry PrintVariablesInTime
c***************************************************************************************
c
c      write(20,*) ndt,time,(MassFlowRate(i)*57.1428571,i=1,6)
c      write(21,*) ndt,time,(OutletPressure(i)*0.007500617,i=1,6)
c
      return
c
c***************************************************************************************
      entry UpdateUnsteadyBoundaryConditions
c***************************************************************************************
c
!      period=0.75
!      timeFloor=int(time/period)
!      localTime=time-float(timeFloor)*period
!c
!      i=1
!      if(localTime.le.0.28) then
!c
!        do j=1,NBFaces(i)
!c
!          BuVelocity(i,j)=0.
!          BvVelocity(i,j)=0.717*dsin(50.*pi*localTime/14.)+0.1
!          BwVelocity(i,j)=0.
!c
!        enddo
!c
!      else
!c
!        do j=1,NBFaces(i)
!c
!          BuVelocity(i,j)=0.
!          BvVelocity(i,j)=0.1
!          BwVelocity(i,j)=0.
!c
!        enddo
!c
!      endif
!c
      return
c
c***************************************************************************************
      entry updateBoundaryConditions
c***************************************************************************************
c
       return
        i=2
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
          BPressure(i,j)=Pressure(k)
c     *          -BDensity(i,j)*GravityY*(1.-BFaceCentroidy(i,j))
c      
        enddo     
c
      return
c***************************************************************************************
      entry updateDensity
c***************************************************************************************
c 



c
      return
c***************************************************************************************
      entry updateViscosity
c***************************************************************************************
c






c
      return
c***************************************************************************************
      entry updateConductivity
c***************************************************************************************
c





c
      return
c***************************************************************************************
      entry updateSpecificHeat
c***************************************************************************************
c

c
      return
c***************************************************************************************
      entry updateSources
c***************************************************************************************
c

c
      return
c***************************************************************************************
      entry updateScalars
c***************************************************************************************
c

c
      return
      end
c
c#############################################################################################