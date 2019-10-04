      MODULE User0
 
      Implicit none
!
!--- Type of mesh used
!
      character*20, save :: MeshType='polymesh' !'neutral' 'polymesh'
      logical, save :: LprintFoamCase=.true.
      logical, save :: LprintParaviewNodeBased=.false.
      logical, save :: LprintParaviewCellBased=.false.
      logical, save :: LprintParaviewFile=.false.
      logical, save :: LprintResultsFile=.false.
      integer, save :: nInterPoints=4     !needed for Lambda Euler Lagrange
      character*500, save :: directory,NeutralMeshdirectory,PolyMeshdirectory,SolutionDirectory,FoamCasedirectory
      double precision, save :: GridScalex=1.d0
      double precision, save :: GridScaley=1.d0
      double precision, save :: GridScalez=1.d0
!
!--- Read old solution to continue
!
      logical, save :: LReadOldSolution=.false.
      integer, save :: nTimeSave=1  ! save solution every 'nTimeSave' step   
!
!--- Flow type
!
      logical, save :: Lcompressible=.false.  !.true. if flow is compressible
      logical, save :: Linviscid=.false.      !.true. if viscosity is zero
      logical, save :: LCoriolis=.false.      !coriolis acceleration
!
!--- Solver
!
      character*10, save :: Algorithm='simple' !'simplec'
!
!--- Monitoring location
!
      integer, save :: IMonitor=1  !Monitoring element determined by the code from the set location
      double precision, save :: xMonitor=0. !x location of monitoring element
      double precision, save :: yMonitor=0. !y location of monitoring element
      double precision, save :: zMonitor=0. !z location of monitoring element
!
!--- Fixing pressure at a location (for incompressible flow when no pressure is specified)
!
      logical, save :: LfixPressure=.false.
      integer, save :: FixAtLocation=1   !Element number where pressure is fixed. Determined by the code from the set location 
      double precision, save :: xFixPressure=0. !x location of the element where pressure is fixed
      double precision, save :: yFixPressure=0. !y location of the element where pressure is fixed
      double precision, save :: zFixPressure=0. !z location of the element where pressure is fixed
      double precision, save :: FixedPressureValue !Fixed value of pressure
!
!--- Reference pressure location
!
      integer, save :: ReferencePressureLocation=1
      double precision, save :: xRefPressure=0.
      double precision, save :: yRefPressure=0.
      double precision, save :: zRefPressure=0.
!
!--- Grid variables
!
      character*80, save :: name
      character*80, save :: FilNeu,FilGrd,Filin,Filout,Filres
      character*80, save :: Filprint,FilMG,Filwdist,Filold,Filgnuplot
      character*100, save :: Filcf
      logical, save :: LStructured=.false.
      logical, save :: LReadSavedGrid=.false.
      logical, save :: LReadSavedWallDistance=.false.
!
!---  Diffusion variables
!
      integer, save :: MethodDecomposeS=3
      logical, save :: LanisotropicDiffusion=.false.
      integer, save :: MethodDecomposeSprime=4
!
!---  Lambda Euler Lagrange variables
!
      double precision, save :: alpha1Lambda,alpha2Lambda  !input
      double precision, save :: IAbsDivergence,IMaxDivergence,IrmsDivergence    !output
      double precision, save :: FAbsDivergence,FMaxDivergence,FrmsDivergence    !output
!
!     Interpolation Scheme: Diffusion Coefficient ('average', 'harmonic', 'upwind', 'downwind')
!
      character*16, save :: InterpolationSchemeGamaMomentum='average'  
      character*16, save :: InterpolationSchemeGamaEnergy='harmonic'
      character*16, save :: InterpolationSchemeGamaTKE='average'
      character*16, save :: InterpolationSchemeGamaTED='average'
      character*16, save :: InterpolationSchemeGamaTOmega='average'
      character*16, save :: InterpolationSchemeGamaTKL='average'
      character*16, save :: InterpolationSchemeGamaMED='average'
      character*16, save :: InterpolationSchemeGamaTGamma='average'
      character*16, save :: InterpolationSchemeGamaTReTheta='average'
      character*16, save :: InterpolationSchemeGamaTfRelaxation='average'
      character*16, save :: InterpolationSchemeGamaTurbulentV2='average'
      character*16, save :: InterpolationSchemeGamaTurbulentZeta='average'
      character*16, save :: InterpolationSchemeGamaLambdaELE='average'
!
!--- Gradient calculation
!
      integer, save :: MethodCalcGradientMomentum =2
      integer, save :: MethodCalcGradientContinuity=2
      integer, save :: MethodCalcGradientEnergy=2
      integer, save :: MethodCalcGradientDensity=2
      integer, save :: MethodCalcGradientTKE=2
      integer, save :: MethodCalcGradientTED=2
      integer, save :: MethodCalcGradientTOmega=2
      integer, save :: MethodCalcGradientModifiedED=2
      integer, save :: MethodCalcGradientTurbulentKL=2
      integer, save :: MethodCalcGradientTGamma=2
      integer, save :: MethodCalcGradientTReTheta=2
      integer, save :: MethodCalcGradientTfRelaxation=2
      integer, save :: MethodCalcGradientTurbulentV2=2
      integer, save :: MethodCalcGradientTurbulentZeta=2
      integer, save :: MethodCalcGradientLambdaELE=2
!
      integer, save :: nIterGradientMomentum=2
      integer, save :: nIterGradientContinuity=2 
      integer, save :: nIterGradientEnergy=2
      integer, save :: nIterGradientDensity=2
      integer, save :: nIterGradientTKE=2
      integer, save :: nIterGradientTED=2
      integer, save :: nIterGradientTOmega=2
      integer, save :: nIterGradientModifiedED=2
      integer, save :: nIterGradientTurbulentKL=2
      integer, save :: nIterGradientTGamma=2
      integer, save :: nIterGradientTReTheta=2
      integer, save :: nIterGradientTfRelaxation=2
      integer, save :: nIterGradientTurbulentV2=2
      integer, save :: nIterGradientTurbulentZeta=2
      integer, save :: nIterGradientLambdaELE=2
!
      double precision, save :: InvDistancePower=1.
!
      logical, save :: LimitGradientMomentum=.false.
      logical, save :: LimitGradientContinuity=.false.
      logical, save :: LimitGradientEnergy=.false.
      logical, save :: LimitGradientDensity=.false.
      logical, save :: LimitGradientTKE=.false.
      logical, save :: LimitGradientTED=.false.
      logical, save :: LimitGradientTOmega=.false.
      logical, save :: LimitGradientModifiedED=.false.
      logical, save :: LimitGradientTurbulentKL=.false.
      logical, save :: LimitGradientTGamma=.false.
      logical, save :: LimitGradientTReTheta=.false.
      logical, save :: LimitGradientTfRelaxation=.false.
      logical, save :: LimitGradientTurbulentV2=.false.
      logical, save :: LimitGradientTurbulentZeta=.false.
      logical, save :: LimitGradientLambdaELE=.false.
!
      integer, save :: LimitGradientMomentumMethod=2
      integer, save :: LimitGradientContinuityMethod=2
      integer, save :: LimitGradientEnergyMethod=2
      integer, save :: LimitGradientTKEMethod=2
      integer, save :: LimitGradientTEDMethod=2
      integer, save :: LimitGradientTOmegaMethod=2
      integer, save :: LimitGradientTurbulentKLMethod=2
      integer, save :: LimitGradientModifiedEDMethod=2
      integer, save :: LimitGradientDensityMethod=2
      integer, save :: LimitGradientTGammaMethod=2
      integer, save :: LimitGradientTfRelaxationMethod=2
      integer, save :: LimitGradientTReThetaMethod=2
      integer, save :: LimitGradientTurbulentV2Method=2
      integer, save :: LimitGradientTurbulentZetaMethod=2
      integer, save :: LimitGradientLambdaELEMethod=2
!
      logical, save :: LrelaxGradientMomentum=.true.
      logical, save :: LrelaxGradientContinuity=.true.
      logical, save :: LrelaxGradientEnergy=.true.
      logical, save :: LrelaxGradientTKE=.true.
      logical, save :: LrelaxGradientTED=.true.
      logical, save :: LrelaxGradientTOmega=.true.
      logical, save :: LrelaxGradientMED=.true.
      logical, save :: LrelaxGradientTKL=.true.
      logical, save :: LrelaxGradientTGamma=.true.
      logical, save :: LrelaxGradientTReTheta=.true.
      logical, save :: LrelaxGradientTv2=.true.
      logical, save :: LrelaxGradientTZeta=.true.
      logical, save :: LrelaxGradientTfRelaxation=.true.
      logical, save :: LrelaxGradientLambdaELE=.true.
      logical, save :: LrelaxGradientDensity=.true.
      logical, save :: LrelaxGradientOthers=.true.
!      
      double precision, save :: urfGradientMomentum=0.75
      double precision, save :: urfGradientContinuity=0.75
      double precision, save :: urfGradientEnergy=0.75
      double precision, save :: urfGradientTKE=0.75
      double precision, save :: urfGradientTED=0.75
      double precision, save :: urfGradientTOmega=0.75
      double precision, save :: urfGradientMED=0.75
      double precision, save :: urfGradientTKL=0.75
      double precision, save :: urfGradientTGamma=0.75
      double precision, save :: urfGradientTReTheta=0.75
      double precision, save :: urfGradientTv2=0.75
      double precision, save :: urfGradientTZeta=0.75
      double precision, save :: urfGradientTfRelaxation=0.75
      double precision, save :: urfGradientLambdaELE=0.75
      double precision, save :: urfGradientDensity=0.75
      double precision, save :: urfGradientOthers=0.75
!
!--- Limit temperature
!
     logical, save :: LimitTemperature=.false.
     double precision, save :: tempmin=100.
     double precision, save :: tempmax=1.e6
!
!--- Transient variables
!
      logical, save :: LUnsteady=.false.
      character*14, save :: TransientScheme='cranknicolson2'  !euler !adamsmoulton !cranknicolson2
      double precision, save :: dt=10.
      double precision, save :: timemax=10000.
      double precision, save :: time=0.
      integer, save :: nTimePrint=1000
      logical, save :: LadaptiveTimeStep=.false.
      double precision, save :: GlobalCourantNumber=2.
      double precision, save :: minimumdtChangeFactor=0.5
      double precision, save :: maximumdtChangeFactor=5.
!
!--- Solving variables
!
      logical, save :: LSolveMomentum=.false.
      logical, save :: LSolveContinuity=.false.
      logical, save :: LSolveEnergy=.false.
      logical, save :: LSolveTurbulenceKineticEnergy=.false.
      logical, save :: LSolveTurbulenceDissipationRate=.false.
      logical, save :: LSolveTurbulenceSpecificDissipationRate=.false.
      logical, save :: LSolveTurbulentKL=.false.
      logical, save :: LSolveModifiedED=.false.
      logical, save :: LSolveTurbulenceGammaEquation=.false.
      logical, save :: LSolveTurbulenceReynoldsThetaEquation=.false.
      logical, save :: LSolveTurbulencefRelaxationEquation=.false.
      logical, save :: LSolveTurbulenceV2Equation=.false.
      logical, save :: LSolveTurbulenceZetaEquation=.false.
      logical, save :: LSolveLambdaELEEquation=.false.
!      
!--- Energy equation to be solve      
!   
      character*13, save :: EnergyEquation='temperature' !'temperature' !'htotal'    
!
!--- Multigrid variables
!
      logical, save :: LTestMultiGrid=.true.
      logical, save :: LPrintMultiGridResiduals=.true.
      logical, save :: LMultigridMomentum=.true.
      logical, save :: LMultigridContinuity=.true.
      logical, save :: LMultigridEnergy=.true.
      logical, save :: LMultigridTKE=.false.
      logical, save :: LMultigridTED=.false.
      logical, save :: LMultigridTOmega=.false.
      logical, save :: LMultigridTurbulentKL=.false.
      logical, save :: LMultigridModifiedED=.false.
      logical, save :: LMultigridTGamma=.false.
      logical, save :: LMultigridTReTheta=.false.
      logical, save :: LMultigridTfRelaxation=.false.
      logical, save :: LMultigridTurbulentV2=.false.
      logical, save :: LMultigridTurbulentZeta=.false.
      logical, save :: LMultigridLambdaELE=.false.
      character*16, save :: MGType='geometricnode' !'algebraic'/'geometricelement'/'geometricnode'
      character*35, save :: MGVariable='scalar1'
      character*6, save :: MultiGridCycleType = 'fcycle'  !'fcycle'  'wcycle' 'vcycle'
      integer, save :: nprintMG=1000
      integer, save :: MultiGridpreSweep=3
      integer, save :: MultiGridpostSweep=5
      integer, save :: MaxMultiGridCycles=30
      integer, save :: reAgglomorate=300
      integer, save :: minNumberofParents=10
      integer, save :: maxNumberofCoarseLevels=10
      integer, save :: nIterStartApplyingMG=5
      double precision, save :: MultiGridrrf=0.1
!
!--- Rate of reduction for algebraic solver 
!
      double precision, save :: rrFMomentum=0.3
      double precision, save :: rrFContinuity=0.1
      double precision, save :: rrFEnergy=0.3
      double precision, save :: rrFTKE=0.3
      double precision, save :: rrFTED=0.3
      double precision, save :: rrFTOmega=0.3
      double precision, save :: rrFTurbulentKL=0.3
      double precision, save :: rrFModifiedED=0.3
      double precision, save :: rrFTGamma=0.3
      double precision, save :: rrFTReTheta=0.3
      double precision, save :: rrFTfRelaxation=0.3
      double precision, save :: rrFTurbulentV2=0.3
      double precision, save :: rrFTurbulentZeta=0.3
      double precision, save :: rrFLambdaELE=0.3
!
!--- Algebraic solver type for different variables
!
      character*6, save :: ASSolverMomentum='ilu'
      character*6, save :: ASSolverContinuity='ilu'
      character*6, save :: ASSolverEnergy='ilu'
      character*6, save :: ASSolverTKE='ilu'
      character*6, save :: ASSolverTED='ilu'
      character*6, save :: ASSolverTOmega='ilu'
      character*6, save :: ASSolverTurbulentKL='ilu'
      character*6, save :: ASSolverModifiedED='ilu'
      character*6, save :: ASSolverTGamma='ilu'
      character*6, save :: ASSolverTReTheta='ilu'
      character*6, save :: ASSolverTfRelaxation='ilu'
      character*6, save :: ASSolverTurbulentV2='ilu'
      character*6, save :: ASSolverTurbulentZeta='ilu'
      character*6, save :: ASSolverLambdaELE='ilu'
!
!--- maximum iterations for the algebraic solver for difeerent variables
!
      integer, save :: ASIterMomentum=5
      integer, save :: ASIterContinuity=20
      integer, save :: ASIterEnergy=5
      integer, save :: ASIterTKE=5
      integer, save :: ASIterTED=5
      integer, save :: ASIterTOmega=5
      integer, save :: ASIterTurbulentKL=5
      integer, save :: ASIterModifiedED=5
      integer, save :: ASIterTGamma=5
      integer, save :: ASIterTReTheta=5
      integer, save :: ASIterTfRelaxation=5
      integer, save :: ASIterTurbulentV2=5
      integer, save :: ASIterTurbulentZeta=5
      integer, save :: ASIterLambdaELE=5
!
!--- Global variables (maximum outer number of iterations, 
!                      stoping criteria type and value)
!
      integer, save :: IterMax=2000
      integer, save :: NstopType=2
      double precision, save :: maximumResidual=1.e-7
!
!--- underrelaxation
!
      logical, save :: LRelaxMomentum=.false.
      logical, save :: LRelaxPressure=.false.
      logical, save :: LRelaxEnergy=.false.
      logical, save :: LRelaxTKE=.false.
      logical, save :: LRelaxTED=.false.
      logical, save :: LRelaxTOmega=.false.
      logical, save :: LRelaxTurbulentKL=.false.
      logical, save :: LRelaxModifiedED=.false.
      logical, save :: LRelaxTGamma=.false.
      logical, save :: LRelaxTReTheta=.false.
      logical, save :: LRelaxTfRelaxation=.false.
      logical, save :: LRelaxTurbulentV2=.false.
      logical, save :: LRelaxTurbulentZeta=.false.
      logical, save :: LRelaxLambdaELE=.false.
      
      logical, save :: LFalseTransientMomentum=.false.
      logical, save :: LFalseTransientEnergy=.false.
      logical, save :: LFalseTransientTKE=.false.
      logical, save :: LFalseTransientTED=.false.
      logical, save :: LFalseTransientTOmega=.false.
      logical, save :: LFalseTransientTurbulentKL=.false.
      logical, save :: LFalseTransientModifiedED=.false.
      logical, save :: LFalseTransientTGamma=.false.
      logical, save :: LFalseTransientTReTheta=.false.
      logical, save :: LFalseTransientTfRelaxation=.false.
      logical, save :: LFalseTransientTurbulentV2=.false.
      logical, save :: LFalseTransientTurbulentZeta=.false.
      logical, save :: LFalseTransientLambdaELE=.false.
      
      double precision, save :: urfMomentum,urfEnergy,urfPressure,urfTGamma,urfTReTheta
      double precision, save :: urfTKE,urfTED,urfTOmega,urfTurbulentKL,urfModifiedED,urfTViscosity
      double precision, save :: urfTfRelaxation,urfTurbulentV2,urfTurbulentZeta,urfLambdaELE
      double precision, save :: FalseDtTKE,FalseDtTED,FalseDtTOmega,FalseDtTurbulentKL,FalseDtModifiedED
      double precision, save :: FalseDtMomentum,FalseDtEnergy,FalseDtTGamma,FalseDtTReTheta
      double precision, save :: FalseDtTfRelaxation,FalseDtTurbulentV2,FalseDtTurbulentZeta,FalseDtLambdaELE
!
!--- High resolution convection schemes
!
      character*4, save :: HRFrameworkMomentum='tvd'        !'nvf' 'none'
      character*4, save :: HRFrameworkEnergy='tvd'          !'nvf' 'none'
      character*4, save :: HRFrameworkDensity='tvd'         !'nvf' 'none'
      character*4, save :: HRFrameworkTKE='tvd'             !'nvf' 'none'
      character*4, save :: HRFrameworkTED='tvd'             !'nvf' 'none'
      character*4, save :: HRFrameworkTOmega='tvd'          !'nvf' 'none'
      character*4, save :: HRFrameworkTurbulentKL='tvd'     !'nvf' 'none'
      character*4, save :: HRFrameworkModifiedED='tvd'      !'nvf' 'none'
      character*4, save :: HRFrameworkTGamma='tvd'          !'nvf' 'none'
      character*4, save :: HRFrameworkTReTheta='tvd'        !'nvf' 'none'
      character*4, save :: HRFrameworkTfRelaxation='tvd'    !'nvf' 'none'
      character*4, save :: HRFrameworkTurbulentV2='tvd'     !'nvf' 'none'
      character*4, save :: HRFrameworkTurbulentZeta='tvd'   !'nvf' 'none'
!      
      character*20, save :: ConvectionSchemeMomentum='minmod'
      character*20, save :: ConvectionSchemeEnergy='minmod'
      character*20, save :: ConvectionSchemeDensity='minmod'
      character*20, save :: ConvectionSchemeTKE='minmod'
      character*20, save :: ConvectionSchemeTED='minmod'
      character*20, save :: ConvectionSchemeTOmega='minmod'
      character*20, save :: ConvectionSchemeTurbulentKL='minmod'
      character*20, save :: ConvectionSchemeModifiedED='minmod'
      character*20, save :: ConvectionSchemeTGamma='minmod'
      character*20, save :: ConvectionSchemeTReTheta='minmod'
      character*20, save :: ConvectionSchemeTfRelaxation='minmod'
      character*20, save :: ConvectionSchemeTurbulentV2='minmod'
      character*20, save :: ConvectionSchemeTurbulentZeta='minmod'

      double precision, save :: BleedMomentum=0.
      double precision, save :: BleedEnergy=0.
      double precision, save :: BleedDensity=0.
      double precision, save :: BleedTKE=0.
      double precision, save :: BleedTED=0.
      double precision, save :: BleedTOmega=0.
      double precision, save :: BleedTurbulentKL=0.
      double precision, save :: BleedModifiedED=0.
      double precision, save :: BleedTGamma=0.
      double precision, save :: BleedTReTheta=0.
      double precision, save :: BleedTfRelaxation=0.
      double precision, save :: BleedTurbulentV2=0.
      double precision, save :: BleedTurbulentZeta=0.
      double precision, save :: BettaConvection=0.1          !      1/10 <= bettam <= 1/2
!
      integer, save :: nIterStartApplyingHR=0
!
!--- Indicator to solve convection without solving the flow problem
!
      logical, save :: LConvectScalar=.false.
!
!--- Type of gradient interpolation to face
!
      logical, save :: LPrintGradients=.false.
      character*16, save :: GradientInterpolationSchemeMomentum='average'
      character*16, save :: GradientInterpolationSchemeContinuity='averagecorrected'
      character*16, save :: GradientInterpolationSchemeEnergy ='average'
      character*16, save :: GradientInterpolationSchemeTKE='average'
      character*16, save :: GradientInterpolationSchemeTED='average'
      character*16, save :: GradientInterpolationSchemeTOmega='average'
      character*16, save :: GradientInterpolationSchemeTurbulentKL='average'
      character*16, save :: GradientInterpolationSchemeModifiedED='average'
      character*16, save :: GradientInterpolationSchemeTGamma='average'
      character*16, save :: GradientInterpolationSchemeTReTheta='average'
      character*16, save :: GradientInterpolationSchemeTfRelaxation='average'
      character*16, save :: GradientInterpolationSchemeTurbulentV2='average'
      character*16, save :: GradientInterpolationSchemeTurbulentZeta='average'
      character*16, save :: GradientInterpolationSchemeLambdaELE='average'
!
!--- Constant values
!
      double precision, save :: ConstantDensity=1.
      double precision, save :: ConstantSpecificHeat=1.
      double precision, save :: ConstantViscosity=0.
      double precision, save :: ConstantConductivity=0.
      double precision, save :: constantStagnationTemperature=1.
      double precision, save :: constantStagnationPressure=1.
!
!--- Bulks modulus for calculating speed of sound in liquids
!
      double precision, save :: BulkModulus=144648.
!
!           Bulk Modulus Values of some fluids
!
!       Acetone	                0.92d9          (PaN/m2)
!       Benzene	                1.05d9          (PaN/m2)
!       Carbon Tetrachloride	1.32d9          (PaN/m2)
!       Ethyl Alcohol	        1.06d9          (PaN/m2)
!       Gasoline	            1.3d9           (PaN/m2)
!       Glycerin	            4.35d9          (PaN/m2)
!       ISO 32 mineral oil	    1.8d9           (PaN/m2)
!       Kerosene	            1.3d9           (PaN/m2)
!       Mercury	                28.5d9          (PaN/m2)
!       Paraffin Oil	        1.66d9          (PaN/m2)
!       Petrol	                1.07d9-1.49d9   (PaN/m2)
!       Phosphate ester	        3d9             (PaN/m2)
!       SAE 30 Oil	            1.5d9           (PaN/m2)
!       Seawater	            2.34d9          (PaN/m2)
!       Sulfuric Acid	        3.0d9           (PaN/m2)
!       Water	                2.15d9          (PaN/m2)
!       Water - glycol	        3.4d9           (PaN/m2)
!       Water in oil emulsion	2.3d9           (PaN/m2)
! 
!   Equivalent for air (for density=1.2, to get a=330 ---> BulkModulus=130680)
!
!--- Periodic boundary condition parameters (should be deleted from here they are in BC2
!
!      logical, save :: TranslationalPeriodicity=.false.
!      logical, save :: RotationalPeriodicity=.false.
!      double precision, save :: xdirTranslate,ydirTranslate,zdirTranslate
!      double precision, save :: AngleOfRotation
!
!--- Buoyancy
!
     character*10, save :: BuoyancyModel='rhog'    !rhog  !boussinesq
     logical, save :: LBuoyancy=.false.
!
!---- Data reporting
!
      double precision, save :: urfViscosity=1.
      double precision, save :: urfConductivity=1.
!
!---  Turbulence Variables
!
     logical, save :: LTurbulentFlow=.false.
     logical, save :: LimitTurbulenceProduction=.true.
     logical, save :: LstagnationCorrectionKE=.false.
     logical, save :: LtestTurbulenceModel=.true.
     logical, save :: LcalculateReTheta=.true.   !for flat plate only
     logical, save :: LRotation=.false.
!                                             ! smooth walls: 'wilcox' 'menter'
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    Smooth and rough wall functions
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Parameters related to rough and smooth walls
!
     logical, save :: LRough=.false.     !make sure LFullyRough is set to false
     double precision, save :: CRoughEnergy=0.2
     double precision, save, dimension(:), allocatable :: CRoughScalar
     double precision, save, dimension(:), allocatable :: GrainSize
     character*20, save :: WallBCModel='menter' ! rough walls: 'wilcox' 'knopp' 'nikuradse'  'colebrook'
!                                               ! smooth walls: 'wilcox' 'menter'
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    End of parameters related to smooth and rough wall functions
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    Fully rough wall functions for atmospheric boundary layers
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
     logical, save :: LFullyRough=.false.   !make sure LRough is set to false
     double precision, save, dimension(:), allocatable :: RoughnessHeight
     double precision, save, dimension(:), allocatable :: DisplacementHeight
!    The roughness height y0 (RoughnessHeight) is related to the size of the roughness elements 
!    on the surface, and is typically between 1/10 and 1/30 of the average height of the roughness  
!    elements.Some commonly-used values of y0 are listed below:
!
!           Surface type	                                Roughness height y0 (m)
!           Calm open water -------------------------------------     0.0002
!           Rough open sea	-------------------------------------     0.001
!           Open flat terrain, grass, few isolated obstacles ----	  0.03
!           Low crops, occasional large obstacles ---------------	  0.10
!           High crops, scattered obstacles ---------------------	  0.25
!           Parkland, bushes, numerous obstacles ----------------	  0.50
!           Suburb, forest, regular large obstacle coverage -----	  0.50 to 1.0
!
!           y0 Values greater than 1m are rare and indicate excessively rough terrain.
!
!    The zero-plane displacement d (DisplacementHeight) is the height above the ground at which 
!    zero wind speed is achieved as a result of flow obstacles such as trees or buildings, 
!    i.e. U=0.0 at y=y0+d. Often the displacement height is zero, but for flow over an array of 
!    densely packed objects (e.g. a forest, cropland or buildings), an offset in height is introduced 
!    into the log law to allow for the upward displacement of the flow by the surface objects. This 
!    displacement height d is usually estimated as 2/3 of the average height of the obstacles.
!    The displacement height d is a positive quantity, although one can set d=-y0 to use wall 
!    functions that are consistent with the inlet wind profiles proposed by Richards & Hoxey (1993), 
!    which are often used in wind-engineering simulations.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    end of Fully rough wall functions for atmospheric boundary layers
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Parameters related to turbulent buoyancy correction
     logical, save :: LBuoyancyCorrection=.false.
! End of parameters related to turbulent buoyancy correction
!
! Parameters related to Komega-SST model
     logical, save :: LKOmegaSSTRotationCurvatureCorrection=.false.
!     character*10, save :: RotationCurvatureMethod='hellsten'   ! 'spalartshur'  'hellsten'
! Parameters related to k-spsilon-Rt model
     logical, save :: LKEpsilonRtRotationCurvatureCorrection=.false.
! End of parameters related to k-spsilon-Rt model
!
! Compressibility correction methods (dilatational dissipation correction)
     logical, save :: LCompressibilityCorrection=.false.
     character*20, save :: CompressibilityCorrectionMethod='wilcox' !'sarkar' 'zeman'  'wilcox' 'nicoetala' 'nicoetalb' 'abdolhamid'
     double precision, save :: Zeman1=0.25   ! 0.25 b.l. flows  0.1 free shear flows
     double precision, save :: Zeman2=0.66   ! 0.66 b.l. flows  0.6 free shear flows
! End of compressibility correction methods (dilatational dissipation correction)
!
! Compressibility correction methods (pressure dilatational correction)
     ! not implemented yet (check Wilcox book)
! End of compressibility correction methods for k-omega based models


! Additional compressibility correction to turbulent viscosity by Temperature correction method
!  applicable to standard k-epsilon model, Wray-Agarwal, and omega-based models
     logical, save :: LCompressibilityTemperatureCorrection=.false.
! end of additional compressibility correction to turbulent viscosity by Temperature correction method
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    Parameter related to the YAP correction with low Reynolds number k-epsilon models 
!             (it can be added to high Reynolds number k-epsilon models)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!   The Yap correction is active in nonequilibrium flows and tends to reduce the departure
!   of the turbulence length scale from its local equilibrium level. It is an ad-hoc fix 
!   which seldom causes any problems and often improves the predictions.
!   Yap showed strongly improved results with the k-epsilon model in separated flows when 
!   using this extra source term. The Yap correction has also been shown to improve results
!   in a stagnation region. Launder [Launder, B. E. (1993)] recommends that the Yap correction
!   should always be used with the epsilon equation.
!
     logical, save :: LYapCorrection=.false.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    end of parameter related to the YAP correction 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Parameters related to Spalart-Allmaras model
     logical, save :: noft2=.true.
     logical, save :: LTransitionalSA=.false.
     logical, save :: LNegativeSpalartAllmaras=.true.
     logical, save :: LSpalartAllmarasRotationCorrection=.false.
     logical, save :: LSpalartAllmarasRotationCurvatureCorrection=.false.
     logical, save :: LSpalartCompressibilityCorrection=.false.
! End of parameters related to Spalart-Allmaras model
!
! Parameters related to Wray-Agarwal model
     logical, save :: LWrayAgarwalRotationCurvatureCorrection=.false.
     logical, save :: LWallDistanceFreeWA=.true.
! End of parameters related to Wray-Agarwal model
!
! Parameters related to Spalart-Allmaras, Wray-Agarwal, and K-Epsilon-Rt models
     character*20, save :: RotationCurvatureMethod='zhangyang'   ! 'spalartshur'  'zhangyang'
! End of parameters related to Spalart-Allmaras, Wray-Agarwal, and K-Epsilon-Rt models
!
     character*20, save :: TurbulenceModel='komega' !'kepsilonchc'  'kepsilon' !'komega' 'komegaepsilon' 'komegasst'
     character*20, save :: WallTreatment='wallfunctions' !'wallfunctions'   !'lowreynoldsnumber'
     character*20, save :: WallFunctionKind='nonequilibrium' !'equilibrium' 'nonequilibrium'
     character*20, save :: MomentumWallFunctionType='scalable' !'scalable' !'standard'  !'automatic'
     character*20, save :: EnergyWallFunctionType='scalable'  !'scalable' !'standard'  !'automatic'
     character*20, save :: ScalarWallFunctionType='scalable' !'scalable' !'standard'  !'automatic'
     double precision, save :: urff2Coefficient=1.
     double precision, save :: urffmuCoefficient=1.  
     double precision, save :: urfWallEpsilon=1.  
     double precision, save :: MaximumViscosityRatio=1.e6
!
!---- Two-equation turbulence models
! 
!                ---------k-epsilon based---------------
!    kepsilon           the standard HR k-e model
!    kepsilonsharma     the k-e low reynolds number launder-sharma model
!    kepsilonchc        the k-e low reynolds number chang-hsieh-chen model
!    kepsilonchien      the k-e low reynolds number chien model
!    kepsilonkasagi     the k-e low reynolds number Myong-Kasagi model
!    kepsilontagawa     the k-e low reynolds number Nagano-Tagawa model
!    kepsilonhishida    the k-e low reynolds number Nagano-Hishida model
!    kelambremhorst     the k-e low reynolds number model of Lam and Bremhorst
!    kelambremhorstm    Lam and Bremhorst with the modifications of Lars Davidson 
!    realizable         Realizable k-e model of shih (High Reynolds number only)
!    kepsilonrng        the renormalization group HR k-e model
!    kepsilonv2f        k-epsilon V2-f model
!    kepsilonzetaf      k-epsilon zeta-f model
!
!                ---------k-omega based---------------
!
!    komega             the k-w model of wilcox (low and high reynolds number)
!    komegaepsilon      the k-w model derived from the k-e model (low and high reynolds number)
!    komega2006         the k-w 2006 model of wilcox (low and high reynolds number)
!    komega2006lrn      the k-w 2006 model of wilcox (low reynolds number)
!    komegabsl          the k-w-BSL model of menter (low and high reynolds number)
!    komegasst          the k-w-SST model of menter (low and high reynolds number)
!
!                ---------the k-kl model---------------
!
!    kklmodel           the k-kl model (low Reynolds number only)
!
!---- One-equation turbulence models
! 
!    spalartallmaras    the one-equation model of spalart-allmaras
!    wrayagarwal        the one-equation model of Wray and Agarwal
!    nut92              the one-equation model of Kovasznay 1967. Based on Shur, M., Streets, M.,.....
!     
!---- Three-equation turbulence model
! 
!    kklomega            K-KL-W three-equation Transitional model (low Reynolds number)
!    kepsilonrt          K-e-Rt Three equation low Reynolds number model
!    sstgama             Menter 3-equation transitional turbulence model
!-------------------------------------------------------------------------------------------
!
!---- Four-equation turbulence model
! 
!    sstgamaretheta      Langtry-Menter 4-equation Transitional SST model (low Reynolds number)
!-------------------------------------------------------------------------------------------
!--- Scalars to solve
!-------------------------------------------------------------------------------------------
!
     integer, save :: NumberOfScalarsToSolve=0
!
     character*10, save, dimension(:), allocatable :: ScalarName
!
!--- Gradient calculation
!
      Integer, save, dimension(:), allocatable :: MethodCalcGradientScalar
      Integer, save, dimension(:), allocatable :: nIterGradientScalar
      logical, save, dimension(:), allocatable :: LimitGradientScalar
      Integer, save, dimension(:), allocatable :: LimitGradientScalarMethod
      logical, save, dimension(:), allocatable :: LrelaxGradientScalar
      double precision, save, dimension(:), allocatable :: urfGradientScalar
!
!--- Solving variables
!
      logical, save, dimension(:), allocatable :: LSolveScalar
!
!--- Multigrid variables
!
      logical, save, dimension(:), allocatable :: LMultigridScalar
!
!--- Rate of reduction for algebraic solver 
!
      double precision, save, dimension(:), allocatable :: rrFScalar
!
!--- Algebraic solver type for different variables
!
      character*6, save, dimension(:), allocatable :: ASSolverScalar
!
!--- maximum iterations for the algebraic solver for difeerent variables
!
      integer, save, dimension(:), allocatable :: ASIterScalar
!
!--- underrelaxation
!
      logical, save, dimension(:), allocatable :: LRelaxScalar
      logical, save, dimension(:), allocatable :: LFalseTransientScalar
      double precision, save, dimension(:), allocatable :: urfScalar
      double precision, save, dimension(:), allocatable :: FalseDtScalar
!
!--- High resolution convection schemes
!
      character*4, save, dimension(:), allocatable :: HRFrameworkScalar   !'tvd' 'nvf' 'none'
      character*20, save, dimension(:), allocatable :: ConvectionSchemeScalar
      double precision, save, dimension(:), allocatable :: BleedScalar
!
!--- Type of gradient interpolation to face
!
      character*16, save, dimension(:), allocatable :: GradientInterpolationSchemeScalar
!
!--- Constant values
!
      double precision, save, dimension(:), allocatable :: ConstantDiffusionCoefficientScalar
      double precision, save, dimension(:), allocatable :: ConstantSpecificHeatScalar
!
!-------------------------------------------------------------------------------------------
!--- rFields to solve
!-------------------------------------------------------------------------------------------
!
     integer, save :: NumberOfrFieldsToSolve=0
     logical, save :: LFreeSurfaceFlow=.false.
     logical, save :: LSurfaceTension=.false.
     logical, save :: LHarmonic=.true.
     double precision, save :: SurfaceTension=0.072
!
     character*10, save, dimension(:), allocatable :: rFieldName
!
!--- Gradient calculation
!
      Integer, save, dimension(:), allocatable :: MethodCalcGradientrField
      Integer, save, dimension(:), allocatable :: nIterGradientrField
      logical, save, dimension(:), allocatable :: LimitGradientrField
      Integer, save, dimension(:), allocatable :: LimitGradientrFieldMethod
      logical, save, dimension(:), allocatable :: LrelaxGradientrField
      double precision, save, dimension(:), allocatable :: urfGradientrField
!
!--- Solving variables
!
      logical, save, dimension(:), allocatable :: LSolverField
!
!--- Multigrid variables
!
      logical, save, dimension(:), allocatable :: LMultigridrField
!
!--- Rate of reduction for algebraic solver 
!
      double precision, save, dimension(:), allocatable :: rrFrField
!
!--- Algebraic solver type for different variables
!
      character*6, save, dimension(:), allocatable :: ASSolverrField
!
!--- maximum iterations for the algebraic solver for difeerent variables
!
      integer, save, dimension(:), allocatable :: ASIterrField
!
      character*14, save :: TransientSchemerField
           !'cranknicolson1'  !euler !adamsmoulton !cranknicolson2 !tics1.75 !tics2.5
!
!--- underrelaxation
!
      logical, save, dimension(:), allocatable :: LRelaxrField
      logical, save, dimension(:), allocatable :: LFalseTransientrField
      double precision, save, dimension(:), allocatable :: urfrField
      double precision, save, dimension(:), allocatable :: FalseDtrField
!
!--- High resolution convection schemes
!
      character*4, save, dimension(:), allocatable :: HRFrameworkrField    !'tvd'  'nvf' 'none'
      character*20, save, dimension(:), allocatable :: ConvectionSchemerField   !stacs  !sicsam !hric
      double precision, save, dimension(:), allocatable :: BleedrField
!
!--- Type of gradient interpolation to face
!
      character*16, save, dimension(:), allocatable :: GradientInterpolationSchemerField
!
!--- Constant values
!
      double precision, save, dimension(:), allocatable :: ConstantConductivityrField
      double precision, save, dimension(:), allocatable :: ConstantSpecificHeatrField
      double precision, save, dimension(:), allocatable :: ConstantDensityrField
      double precision, save, dimension(:), allocatable :: ConstantViscosityrField
!
!--- WindKessel Model input parameter ----------------------------------------------------
!
      logical, save, dimension(:), allocatable :: LWindKessel
      integer, save :: WindKesselType=3 !(2, 3, or 4 elements)
      double precision, save :: urfOutletPressure=1.
!




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                     Fluid Element Stress Field
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                          Notation & Sign Convention
!    1. Txy – the first subscript denotes the face the stress is acting
!       on and the 2nd subscript the direction.
!      Example:
!       Tzy is a shear stress acting on the z face in the y direction
!    2. Planes – sign convention: if the outwardly pointing normal
!       is in the positive direction, the plane is positive.
!    3. Stress Sign Convention – a stress component is positive if
!       both the stress component and plane on which it acts are positive.
!       A component of a stress tensor (e.g., Txy, etc.) is positive if the force 
!       vector component and the area normal are either both positive or both negative.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! End of fluid Element Stress Field sign convention
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$












      end MODULE User0