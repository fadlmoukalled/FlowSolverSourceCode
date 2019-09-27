c
C********************************************************************************************
      SUBROUTINE WriteInputDataToFile
C********************************************************************************************
c
      use User0
      use Geometry1, only: NumberOfBCSets
      use Variables1
      use BoundaryConditions1
      use BoundaryConditions2
      use BoundaryConditionsScalar1
      use BoundaryConditionsScalar2
      use BoundaryConditionsrField1
      use BoundaryConditionsrField2
      use Scalar1
      use PhysicalProperties1
      use FlowInOut1, only: MassFlowFraction,LPrintMassFlowRate,
     *                      MassFlowRate
      use ArteryResistance1, only: ArteryResistance,LArteryExplicit,
     *                             urfPressureResistance
      use AveragePressure1, only: AverageOutletPressure
      use WindKessel1, only: ResistanceToBloodFlow,ComplianceC,InertiaL,
     *                       TotalPeripheralResistance,OutletPressure
      use BoundaryFluxes
      use VolumeOfFluid1
      use Coriolis1
      use BoundaryConditionsTurbulence1
      use ReferenceValues1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j
c*********************************************************************************************
c      
      write(10,*) 'name = ',name
      write(10,*) 'NeutralMeshdirectory=',NeutralMeshdirectory
      write(10,*) 'PolyMeshdirectory=',PolyMeshdirectory
      write(10,*) 'SolutionDirectory=',SolutionDirectory
      write(10,*) 'MeshType =',MeshType
      write(10,*) 'LprintParaviewFile=',LprintParaviewFile
      write(10,*) 'LprintParaviewNodeBased=',LprintParaviewNodeBased
      write(10,*) 'LprintParaviewCellBased=',LprintParaviewCellBased
      write(10,*) 'GridScalex=',GridScalex
      write(10,*) 'GridScaley=',GridScaley
      write(10,*) 'GridScalez=',GridScalez
      write(10,*) 'LReadOldSolution=',LReadOldSolution
      write(10,*) 'nTimeSave=',nTimeSave
      write(10,*) 'Lcompressible=',Lcompressible
      write(10,*) 'Linviscid=',Linviscid
      write(10,*) 'LUnsteady=',LUnsteady
      write(10,*) 'LCoriolis=',LCoriolis
      write(10,*) 'LFreeSurfaceFlow=',LFreeSurfaceFlow
      write(10,*) 'LanisotropicDiffusion=',LanisotropicDiffusion
      write(10,*) 'LSolveMomentum=',LSolveMomentum
      write(10,*) 'LSolveContinuity=',LSolveContinuity
      write(10,*) 'LSolveEnergy=',LSolveEnergy
      write(10,*) 'LSolveTurbulenceKineticEnergy=',
     *                        LSolveTurbulenceKineticEnergy
      write(10,*) 'LSolveTurbulenceDissipationRate=',
     *                        LSolveTurbulenceDissipationRate 
      write(10,*) 'LSolveTurbulenceSpecificDissipationRate=',
     *                        LSolveTurbulenceSpecificDissipationRate
      write(10,*) 'LSolveTurbulentKL=',LSolveTurbulentKL
      write(10,*) 'LSolveModifiedED=',LSolveModifiedED
      write(10,*) 'LSolveTurbulenceGammaEquation=',
     *                        LSolveTurbulenceGammaEquation
      write(10,*) 'LSolveTurbulenceReynoldsThetaEquation=',
     *                        LSolveTurbulenceReynoldsThetaEquation
      write(10,*) 'LSolveTurbulencefRelaxationEquation=',
     *                        LSolveTurbulencefRelaxationEquation
      write(10,*) 'LSolveTurbulenceV2Equation=',
     *                        LSolveTurbulenceV2Equation
      write(10,*) 'LSolveTurbulenceZetaEquation=',
     *                        LSolveTurbulenceZetaEquation
      write(10,*) 'LSolveLambdaELEEquation=',LSolveLambdaELEEquation
      write(10,*) 'EnergyEquation=',EnergyEquation
      write(10,*) 'LConvectScalar=',LConvectScalar
      write(10,*) 'NumberOfrFieldsToSolve=',NumberOfrFieldsToSolve
      write(10,*) 'NumberOfScalarsToSolve=',NumberOfScalarsToSolve
      write(10,*) 'NumberofPointSources=',NumberofPointSources
      write(10,*) 'LFalseTransientMomentum=',LFalseTransientMomentum
      write(10,*) 'LFalseTransientEnergy=',LFalseTransientEnergy
      write(10,*) 'LFalseTransientTKE=',LFalseTransientTKE
      write(10,*) 'LFalseTransientTED=',LFalseTransientTED
      write(10,*) 'LFalseTransientTOmega=',LFalseTransientTOmega
      write(10,*) 'LFalseTransientTurbulentKL=',
     *                                 LFalseTransientTurbulentKL
      write(10,*) 'LFalseTransientModifiedED=',LFalseTransientModifiedED
      write(10,*) 'LFalseTransientTGamma=',LFalseTransientTGamma
      write(10,*) 'LFalseTransientTReTheta=',LFalseTransientTReTheta
      write(10,*) 'LFalseTransientTfRelaxation=',
     *                                 LFalseTransientTfRelaxation
      write(10,*) 'LFalseTransientTurbulentV2=',
     *                                 LFalseTransientTurbulentV2
      write(10,*) 'LFalseTransientTurbulentZeta=',
     *                                 LFalseTransientTurbulentZeta
      write(10,*) 'LFalseTransientLambdaELE=',LFalseTransientLambdaELE
      do i=1,NumberOfScalarsToSolve
        write(10,*) 'LFalseTransientScalar(',i,')=',
     *                                 LFalseTransientScalar
      enddo      
      do i=1,NumberOfrFieldsToSolve
        write(10,*) 'LFalseTransientrField(',i,')=',
     *                                 LFalseTransientrField
      enddo      
      write(10,*) 'FalseDtMomentum=',FalseDtMomentum
      write(10,*) 'FalseDtEnergy=',FalseDtEnergy
      write(10,*) 'FalseDtTKE=',FalseDtTKE
      write(10,*) 'FalseDtTOmega=',FalseDtTOmega
      write(10,*) 'FalseDtModifiedED=',FalseDtModifiedED
      write(10,*) 'FalseDtTED=',FalseDtTED
      write(10,*) 'FalseDtTGamma=',FalseDtTGamma
      write(10,*) 'FalseDtTReTheta=',FalseDtTReTheta
      write(10,*) 'FalseDtTfRelaxation=',FalseDtTfRelaxation
      write(10,*) 'FalseDtTurbulentV2=',FalseDtTurbulentV2
      write(10,*) 'FalseDtTurbulentZeta=',FalseDtTurbulentZeta
      write(10,*) 'FalseDtLambdaELE=',FalseDtLambdaELE
      do i=1,NumberOfScalarsToSolve
        write(10,*) 'FalseDtScalar(',i,')=',FalseDtScalar
      enddo      
      do i=1,NumberOfrFieldsToSolve
        write(10,*) 'FalseDtrField(',i,')=',FalseDtrField
      enddo      
      write(10,*) 'Algorithm=',Algorithm
c
      write(10,*) 'AngularVelocityX=',AngularVelocityX
      write(10,*) 'AngularVelocityY=',AngularVelocityY
      write(10,*) 'AngularVelocityZ=',AngularVelocityZ
      write(10,*) 'AxisOfRotationOriginX=',AxisOfRotationOriginX
      write(10,*) 'AxisOfRotationOriginY=',AxisOfRotationOriginY
      write(10,*) 'AxisOfRotationOriginZ=',AxisOfRotationOriginZ
      write(10,*) 'IterMax=',IterMax
      write(10,*) 'xMonitor=',xMonitor
      write(10,*) 'yMonitor=',yMonitor
      write(10,*) 'zMonitor=',zMonitor
      write(10,*) 'xRefPressure=',xRefPressure
      write(10,*) 'yRefPressure=',yRefPressure
      write(10,*) 'zRefPressure=',zRefPressure
      write(10,*) 'LfixPressure=',LfixPressure
      write(10,*) 'FixedPressureValue=',FixedPressureValue
      write(10,*) 'xFixPressure=',xFixPressure
      write(10,*) 'yFixPressure=',yFixPressure
      write(10,*) 'zFixPressure=',zFixPressure
      write(10,*) 'NstopType=',NstopType
      write(10,*) 'maximumResidual=',maximumResidual
      write(10,*) 'ASSolverMomentum=',ASSolverMomentum
      write(10,*) 'ASSolverContinuity=',ASSolverContinuity
      write(10,*) 'ASSolverTKE=',ASSolverTKE
      write(10,*) 'ASSolverTED=',ASSolverTED
      write(10,*) 'ASSolverTOmega=',ASSolverTOmega
      write(10,*) 'ASSolverTurbulentKL=',ASSolverTurbulentKL
      write(10,*) 'ASSolverEnergy=',ASSolverEnergy
      write(10,*) 'ASSolverModifiedED=',ASSolverModifiedED
      write(10,*) 'ASSolverTGamma=',ASSolverTGamma
      write(10,*) 'ASSolverTReTheta=',ASSolverTReTheta
      write(10,*) 'ASSolverTfRelaxation=',ASSolverTfRelaxation
      write(10,*) 'ASSolverTurbulentV2=',ASSolverTurbulentV2
      write(10,*) 'ASSolverTurbulentZeta=',ASSolverTurbulentZeta  
      write(10,*) 'rrfMomentum=',rrfMomentum
      write(10,*) 'rrfContinuity=',rrfContinuity
      write(10,*) 'rrfTKE=',rrfTKE
      write(10,*) 'rrfTED=',rrfTED
      write(10,*) 'rrfTOmega=',rrfTOmega
      write(10,*) 'rrfTurbulentKL=',rrfTurbulentKL
      write(10,*) 'rrfEnergy=',rrfEnergy
      write(10,*) 'rrfModifiedED=',rrfModifiedED
      write(10,*) 'rrfTGamma=',rrfTGamma
      write(10,*) 'rrfTReTheta=',rrfTReTheta
      write(10,*) 'rrfTfRelaxation=',rrfTfRelaxation
      write(10,*) 'rrfTurbulentV2=',rrfTurbulentV2
      write(10,*) 'rrfTurbulentZeta=',rrfTurbulentZeta
c
c--- maximum number of algebraic solver iterations
c
      write(10,*) 'ASIterMomentum=',ASIterMomentum
      write(10,*) 'ASIterContinuity=',ASIterContinuity
      write(10,*) 'ASIterTKE=',ASIterTKE
      write(10,*) 'ASIterTED=',ASIterTED
      write(10,*) 'ASIterTOmega=',ASIterTOmega
      write(10,*) 'ASIterTurbulentKL=',ASIterTurbulentKL
      write(10,*) 'ASIterEnergy=',ASIterEnergy
      write(10,*) 'ASIterModifiedED=',ASIterModifiedED
      write(10,*) 'ASIterTGamma=',ASIterTGamma
      write(10,*) 'ASIterTReTheta=',ASIterTReTheta
      write(10,*) 'ASIterTfRelaxation=',ASIterTfRelaxation
      write(10,*) 'ASIterTurbulentV2=',ASIterTurbulentV2
      write(10,*) 'ASIterTurbulentZeta=',ASIterTurbulentZeta
c
c--- Set the relaxation parameters
c
      write(10,*) 'LRelaxMomentum=',LRelaxMomentum
      write(10,*) 'LRelaxPressure=',LRelaxPressure
      write(10,*) 'LRelaxEnergy=',LRelaxEnergy
      write(10,*) 'LRelaxTKE=',LRelaxTKE
      write(10,*) 'LRelaxTED=',LRelaxTED
      write(10,*) 'LRelaxTOmega=',LRelaxTOmega
      write(10,*) 'LRelaxTurbulentKL=',LRelaxTurbulentKL
      write(10,*) 'LRelaxModifiedED=',LRelaxModifiedED
      write(10,*) 'LRelaxTGamma=',LRelaxTGamma
      write(10,*) 'LRelaxTReTheta=',LRelaxTReTheta
      write(10,*) 'LRelaxTfRelaxation=',LRelaxTfRelaxation
      write(10,*) 'LRelaxTurbulentV2=',LRelaxTurbulentV2
      write(10,*) 'LRelaxTurbulentZeta=',LRelaxTurbulentZeta
c
c--- Assign the underrelaxation factors values
c
      write(10,*) 'urfMomentum=',urfMomentum
      write(10,*) 'urfPressure=',urfPressure
      write(10,*) 'urfTKE=',urfTKE
      write(10,*) 'urfTED=',urfTED
      write(10,*) 'urfTOmega=',urfTOmega
      write(10,*) 'urfTurbulentKL=',urfTurbulentKL
      write(10,*) 'urfEnergy=',urfEnergy
      write(10,*) 'urfTViscosity=',urfTViscosity
      write(10,*) 'urfModifiedED=',urfModifiedED
      write(10,*) 'urfTGamma=',urfTGamma
      write(10,*) 'urfTReTheta=',urfTReTheta
      write(10,*) 'urfTfRelaxation=',urfTfRelaxation
      write(10,*) 'urfTurbulentV2=',urfTurbulentV2
      write(10,*) 'urfTurbulentZeta=',urfTurbulentZeta
c
c--- Set whether to use upwind, NVF or TVD convection schemes
c
      write(10,*) 'HRFrameworkMomentum=',HRFrameworkMomentum
      write(10,*) 'HRFrameworkTKE=',HRFrameworkTKE
      write(10,*) 'HRFrameworkTED=',HRFrameworkTED
      write(10,*) 'HRFrameworkTOmega=',HRFrameworkTOmega
      write(10,*) 'HRFrameworkTurbulentKL=',HRFrameworkTurbulentKL
      write(10,*) 'HRFrameworkEnergy=',HRFrameworkEnergy
      write(10,*) 'HRFrameworkModifiedED=',HRFrameworkModifiedED
      write(10,*) 'HRFrameworkDensity=',HRFrameworkDensity
      write(10,*) 'HRFrameworkTGamma=',HRFrameworkTGamma
      write(10,*) 'HRFrameworkTReTheta=',HRFrameworkTReTheta
      write(10,*) 'HRFrameworkTfRelaxation=',HRFrameworkTfRelaxation
      write(10,*) 'HRFrameworkTurbulentV2=',HRFrameworkTurbulentV2
      write(10,*) 'HRFrameworkTurbulentZeta=',HRFrameworkTurbulentZeta
c
c--- Set when to start applying HR schemes
c
      write(10,*) 'nIterStartApplyingHR=',nIterStartApplyingHR
c
c--- Set the name of scheme to use for variables
c
      write(10,*) 'ConvectionSchemeMomentum=',ConvectionSchemeMomentum
      write(10,*) 'ConvectionSchemeTKE=',ConvectionSchemeTKE
      write(10,*) 'ConvectionSchemeTOmega=',ConvectionSchemeTOmega
      write(10,*) 'ConvectionSchemeTED=',ConvectionSchemeTED
      write(10,*) 'ConvectionSchemeTurbulentKL=',
     *                                    ConvectionSchemeTurbulentKL
      write(10,*) 'ConvectionSchemeModifiedED=',
     *                                    ConvectionSchemeModifiedED
      write(10,*) 'ConvectionSchemeEnergy=',ConvectionSchemeEnergy
      write(10,*) 'ConvectionSchemeDensity=',ConvectionSchemeDensity
      write(10,*) 'ConvectionSchemeTGamma=',ConvectionSchemeTGamma
      write(10,*) 'ConvectionSchemeTReTheta=',ConvectionSchemeTReTheta
      write(10,*) 'ConvectionSchemeTfRelaxation=',
     *                                    ConvectionSchemeTfRelaxation
      write(10,*) 'ConvectionSchemeTurbulentV2=',
     *                                    ConvectionSchemeTurbulentV2
      write(10,*) 'ConvectionSchemeTurbulentZeta=',
     *                                    ConvectionSchemeTurbulentZeta
c
c--- Set the value of the coefficient by which to bleed the HR scheme with upwind scheme 
c    (0= no bleeding, 1= upwind) 
c
      write(10,*) 'BleedMomentum=',BleedMomentum
      write(10,*) 'BleedTKE=',BleedTKE
      write(10,*) 'BleedTED=',BleedTED
      write(10,*) 'BleedTOmega=',BleedTOmega
      write(10,*) 'BleedTurbulentKL=',BleedTurbulentKL
      write(10,*) 'BleedModifiedED=',BleedModifiedED
      write(10,*) 'BleedEnergy=',BleedEnergy
      write(10,*) 'BleedDensity=',BleedDensity
      write(10,*) 'BleedTGamma=',BleedTGamma
      write(10,*) 'BleedTReTheta=',BleedTReTheta
      write(10,*) 'BleedTfRelaxation=',BleedTfRelaxation
      write(10,*) 'BleedTurbulentV2=',BleedTurbulentV2
      write(10,*) 'BleedTurbulentZeta=',BleedTurbulentZeta
      
      
c----------------------------------------------------------------------------------
c--- Set the multigrid variables
c----------------------------------------------------------------------------------
      write(10,*) 'MGType=',MGType
c
      write(10,*) 'LMultigridMomentum=',LMultigridMomentum
      write(10,*) 'LMultigridTKE=',LMultigridTKE
      write(10,*) 'LMultigridTED=',LMultigridTED
      write(10,*) 'LMultigridTOMega=',LMultigridTOMega
      write(10,*) 'LMultigridTurbulentKL=',LMultigridTurbulentKL
      write(10,*) 'LMultigridModifiedED=',LMultigridModifiedED
      write(10,*) 'LMultigridContinuity=',LMultigridContinuity
      write(10,*) 'LMultigridEnergy=',LMultigridEnergy
      write(10,*) 'LMultigridTGamma=',LMultigridTGamma
      write(10,*) 'LMultigridTReTheta=',LMultigridTReTheta
      write(10,*) 'LMultigridTfRelaxation=',LMultigridTfRelaxation
      write(10,*) 'LMultigridTurbulentV2=',LMultigridTurbulentV2
      write(10,*) 'LMultigridTurbulentZeta=',LMultigridTurbulentZeta
c
c--- Set whether to print data during multigrid iterations (for testing)
c
      write(10,*) 'LTestMultiGrid=',LTestMultiGrid
      write(10,*) 'nprintMG=',nprintMG
c
c--- Set the variable on which to base the grid agglomoration (could be a scalar)
c
      write(10,*) 'MGVariable=',MGVariable
c
c--- Set the multigrid cycle to use (V, F, or W)
c
      write(10,*) 'MultiGridCycleType=',MultiGridCycleType
c
c--- Set the number of pre and post sweep number of iterations and the number of MG cycles
c
      write(10,*) 'MultiGridpreSweep=',MultiGridpreSweep
      write(10,*) 'MultiGridpostSweep=',MultiGridpostSweep
      write(10,*) 'MaxMultiGridCycles=',MaxMultiGridCycles
c
c--- Set when to re-agglomorate fine grids to create coarse grid levels
c
      write(10,*) 'reAgglomorate=',reAgglomorate
c
c--- Set the minimum number of elements of a coarse level and the maximum number of coarse levels
c
      write(10,*) 'minNumberofParents=',minNumberofParents
      write(10,*) 'maxNumberofCoarseLevels=',maxNumberofCoarseLevels
c
c--- Set when to start applying multigrid (a number >=1)
c
      write(10,*) 'nIterStartApplyingMG=',nIterStartApplyingMG
c
c--- Set the multigrid Residual reduction factor
c
      write(10,*) 'MultiGridrrf=',MultiGridrrf
c----------------------------------------------------------------------------------
c--- Set the gradient variables
c----------------------------------------------------------------------------------
c
c--- Method to calculate the gradient (1 to 6)
c
      write(10,*) 'MethodCalcGradientMomentum=',
     *                                MethodCalcGradientMomentum
      write(10,*) 'MethodCalcGradientTKE=',MethodCalcGradientTKE
      write(10,*) 'MethodCalcGradientTED=',MethodCalcGradientTED
      write(10,*) 'MethodCalcGradientTOmega=',MethodCalcGradientTOmega
      write(10,*) 'MethodCalcGradientTurbulentKL=',
     *                                MethodCalcGradientTurbulentKL
      write(10,*) 'MethodCalcGradientModifiedED=',
     *                                MethodCalcGradientModifiedED
      write(10,*) 'MethodCalcGradientContinuity=',
     *                                MethodCalcGradientContinuity
      write(10,*) 'MethodCalcGradientEnergy=',MethodCalcGradientEnergy
      write(10,*) 'MethodCalcGradientDensity=',MethodCalcGradientDensity
      write(10,*) 'MethodCalcGradientTGamma=',MethodCalcGradientTGamma
      write(10,*) 'MethodCalcGradientTReTheta=',
     *                                MethodCalcGradientTReTheta
      write(10,*) 'MethodCalcGradientTfRelaxation=',
     *                                MethodCalcGradientTfRelaxation
      write(10,*) 'MethodCalcGradientTurbulentV2=',
     *                                MethodCalcGradientTurbulentV2
      write(10,*) 'MethodCalcGradientTurbulentZeta=',
     *                                MethodCalcGradientTurbulentZeta
c
c--- number of iterations to perform with isterative methods (2 to 4)
c
      write(10,*) 'nIterGradientMomentum=',nIterGradientMomentum
      write(10,*) 'nIterGradientTKE=',nIterGradientTKE
      write(10,*) 'nIterGradientTED=',nIterGradientTED
      write(10,*) 'nIterGradientTOmega=',nIterGradientTOmega
      write(10,*) 'nIterGradientTurbulentKL=',nIterGradientTurbulentKL
      write(10,*) 'nIterGradientModifiedED=',nIterGradientModifiedED
      write(10,*) 'nIterGradientContinuity=',nIterGradientContinuity
      write(10,*) 'nIterGradientEnergy=',nIterGradientEnergy
      write(10,*) 'nIterGradientDensity=',nIterGradientDensity
      write(10,*) 'nIterGradientTGamma=',nIterGradientTGamma
      write(10,*) 'nIterGradientTReTheta=',nIterGradientTReTheta
      write(10,*) 'nIterGradientTfRelaxation=',nIterGradientTfRelaxation
      write(10,*) 'nIterGradientTurbulentV2=',nIterGradientTurbulentV2
      write(10,*) 'nIterGradientTurbulentZeta=',
     *                                nIterGradientTurbulentZeta
c
c--- The power to use for the inverse distance with the Least Square Gradient (method 6) 
c
      write(10,*) 'InvDistancePower=',InvDistancePower
c
c--- Set whether or not to limit the gradient 
c
      write(10,*) 'LimitGradientMomentum=',LimitGradientMomentum
      write(10,*) 'LimitGradientTKE=',LimitGradientTKE
      write(10,*) 'LimitGradientTED=',LimitGradientTED
      write(10,*) 'LimitGradientTOmega=',LimitGradientTOmega
      write(10,*) 'LimitGradientTurbulentKL=',LimitGradientTurbulentKL
      write(10,*) 'LimitGradientModifiedED=',LimitGradientModifiedED
      write(10,*) 'LimitGradientContinuity=',LimitGradientContinuity
      write(10,*) 'LimitGradientEnergy=',LimitGradientEnergy
      write(10,*) 'LimitGradientDensity=',LimitGradientDensity
      write(10,*) 'LimitGradientTGamma=',LimitGradientTGamma
      write(10,*) 'LimitGradientTReTheta=',LimitGradientTReTheta
      write(10,*) 'LimitGradientTfRelaxation=',LimitGradientTfRelaxation
      write(10,*) 'LimitGradientTurbulentV2=',LimitGradientTurbulentV2
      write(10,*) 'LimitGradientTurbulentZeta=',
     *                                LimitGradientTurbulentZeta
c
c--- Set the method to use to limit the gradient 
c
      write(10,*) 'LimitGradientMomentumMethod=',
     *                                LimitGradientMomentumMethod
      write(10,*) 'LimitGradientTKEMethod=',LimitGradientTKEMethod
      write(10,*) 'LimitGradientTEDMethod=',LimitGradientTEDMethod
      write(10,*) 'LimitGradientTOmegaMethod=',LimitGradientTOmegaMethod
      write(10,*) 'LimitGradientTurbulentKLMethod=',
     *                                LimitGradientTurbulentKLMethod
      write(10,*) 'LimitGradientModifiedEDMethod=',
     *                                LimitGradientModifiedEDMethod
      write(10,*) 'LimitGradientContinuityMethod=',
     *                                LimitGradientContinuityMethod
      write(10,*) 'LimitGradientEnergyMethod=',LimitGradientEnergyMethod
      write(10,*) 'LimitGradientDensityMethod=',
     *                                LimitGradientDensityMethod
      write(10,*) 'LimitGradientTGammaMethod=',LimitGradientTGammaMethod
      write(10,*) 'LimitGradientTReThetaMethod=',
     *                                LimitGradientTReThetaMethod
      write(10,*) 'LimitGradientTfRelaxationMethod=',
     *                                LimitGradientTfRelaxationMethod
      write(10,*) 'LimitGradientTurbulentV2Method=',
     *                                LimitGradientTurbulentV2Method
      write(10,*) 'LimitGradientTurbulentZetaMethod=',
     *                                LimitGradientTurbulentZetaMethod
c
c--- Set the type of gradient interpolation to face
c    available types: upwind, downwind, average, and averagecorrected 
c
      write(10,*) 'GradientInterpolationSchemeMomentum=',
     *                        GradientInterpolationSchemeMomentum
      write(10,*) 'GradientInterpolationSchemeTKE=',
     *                        GradientInterpolationSchemeTKE
      write(10,*) 'GradientInterpolationSchemeTED=',
     *                        GradientInterpolationSchemeTED
      write(10,*) 'GradientInterpolationSchemeTOmega=',
     *                        GradientInterpolationSchemeTOmega
      write(10,*) 'GradientInterpolationSchemeTurbulentKL=',
     *                        GradientInterpolationSchemeTurbulentKL
      write(10,*) 'GradientInterpolationSchemeModifiedED=',
     *                        GradientInterpolationSchemeModifiedED
      write(10,*) 'GradientInterpolationSchemeContinuity=',
     *                        GradientInterpolationSchemeContinuity
      write(10,*) 'GradientInterpolationSchemeEnergy=',
     *                        GradientInterpolationSchemeEnergy
      write(10,*) 'GradientInterpolationSchemeTGamma=',
     *                        GradientInterpolationSchemeTGamma
      write(10,*) 'GradientInterpolationSchemeTReTheta=',
     *                        GradientInterpolationSchemeTReTheta
      write(10,*) 'GradientInterpolationSchemeTfRelaxation=',
     *                        GradientInterpolationSchemeTfRelaxation
      write(10,*) 'GradientInterpolationSchemeTurbulentV2=',
     *                        GradientInterpolationSchemeTurbulentV2
      write(10,*) 'GradientInterpolationSchemeTurbulentZeta=',
     *                        GradientInterpolationSchemeTurbulentZeta
c
c----------------------------------------------------------------------------------
c--- Set the transient variables
c----------------------------------------------------------------------------------
c
c--- Set the scheme to use (Euler, CrankNicolson1 (2 steps), 
c                                  CrankNicolson2, and adamsmoulton)
c
      write(10,*) 'TransientScheme=',TransientScheme
c
c--- Set the step size
c
      write(10,*) 'LadaptiveTimeStep=',LadaptiveTimeStep
      write(10,*) 'GlobalCourantNumber=',GlobalCourantNumber
      write(10,*) 'minimumdtChangeFactor=',minimumdtChangeFactor
      write(10,*) 'maximumdtChangeFactor=',maximumdtChangeFactor
      write(10,*) 'dt=',dt
c
c--- Set the maximum time for which the solution to be computed
c
      write(10,*) 'timemax=',timemax
c
c--- Set the time to start computations
c
      write(10,*) 'time=',time
c
c--- Set the number of time steps tp print results
c
      write(10,*) 'nTimePrint=',nTimePrint
c
c----------------------------------------------------------------------------------
c--- Set buoyancy variables
c----------------------------------------------------------------------------------
c
c--- Two models are implemented: boussinesq and rhog
c
      write(10,*) 'LBuoyancy=',LBuoyancy
      write(10,*) 'BuoyancyModel=',BuoyancyModel
c
      write(10,*) 'GravityX=',GravityX
      write(10,*) 'GravityY=',GravityY
      write(10,*) 'GravityZ=',GravityZ
c----------------------------------------------------------------------------------
c--- Start describing the problem
c----------------------------------------------------------------------------------
      do i=1,NumberOfBCSets
c
c--- Assign the types of boundary conditions
c
        write(10,*) 'MassFlowFraction(',i,')=',MassFlowFraction(i)
        write(10,*) 'LPrintMassFlowRate(',i,')=',LPrintMassFlowRate(i)
        write(10,*) 'LArteryExplicit(',i,')=',LArteryExplicit(i)
        write(10,*) 'urfPressureResistance(',i,')=',
     *                                        urfPressureResistance(i)
        write(10,*) 'ArteryResistance(',i,')=',ArteryResistance(i)
c
      enddo
c      
c---- WindKessel Parameters
c      
      write(10,*) 'LWindKessel=',LWindKessel
      write(10,*) 'WindKesselType=',WindKesselType
      write(10,*) 'urfOutletPressure=',urfOutletPressure
c
      do i=1,NumberOfBCSets
c
        write(10,*) 'ResistanceToBloodFlow(',i,')=',
     *                                 ResistanceToBloodFlow(i)
        write(10,*) 'TotalPeripheralResistance(',i,')=',
     *                              TotalPeripheralResistance(i)
        write(10,*) 'ComplianceC(',i,')=',ComplianceC(i)
c      
      enddo      
c      
      do i=1,NumberOfBCSets
c
        write(10,*) 'BoundaryType(',i,')=',BoundaryType(i)
        if(BoundaryType(i).eq.'wall') then
          write(10,*) 'wallTypeM(',i,')=',wallTypeM(i)
          write(10,*) 'wallTypeC(',i,')=',wallTypeC(i)
          write(10,*) 'wallTypeE(',i,')=',wallTypeE(i)
        elseif(BoundaryType(i).eq.'inlet') then
          write(10,*) 'inletTypeM(',i,')=',inletTypeM(i)
          if(inletTypeM(i).eq.'specifiedstaticpressure'.or.
     *         inletTypeM(i).eq.'specifiedstagnationpressure') then
            write(10,*) 'xVeldirection(',i,':)='
            write(10,*)  xVeldirection(i,:)
            write(10,*) 'yVeldirection(',i,':)='
            write(10,*)  yVeldirection(i,:)
            write(10,*) 'zVeldirection(',i,':)='
            write(10,*)  zVeldirection(i,:)
          endif
          write(10,*) 'inletTypeC(',i,')=',inletTypeC(i)
          write(10,*) 'inletTypeE(',i,')=',inletTypeE(i)
        elseif(BoundaryType(i).eq.'outlet') then    
          write(10,*) 'outletTypeM(',i,')=',outletTypeM(i)
          write(10,*) 'outletTypeC(',i,')=',outletTypeC(i)
          if(outletTypeC(i).eq.'specifiedaveragestaticpressure') then
            write(10,*) 'AverageOutletPressure',i,')=',
     *                                 AverageOutletPressure(i)
          endif
          write(10,*) 'outletTypeE(',i,')=',outletTypeE(i)
        elseif(BoundaryType(i).eq.'pressurefarfield') then    
          write(10,*) 'MachFarField(',i,')=',MachFarField(i)
          write(10,*) 'xFlowDirectionFarField(',i,')=',
     *                                xFlowDirectionFarField(i)
          write(10,*) 'yFlowDirectionFarField(',i,')=',
     *                                yFlowDirectionFarField(i)
          write(10,*) 'zFlowDirectionFarField(',i,')=',
     *                                zFlowDirectionFarField(i)
          write(10,*) 'PressureFarField(',i,')=',PressureFarField(i)
          write(10,*) 'TemperatureFarField(',i,')=',
     *                                TemperatureFarField(i)
          write(10,*) 'SpecificHeatFarField(',i,')=',
     *                                SpecificHeatFarField(i)
          if(LTurbulentFlow) then
            if(LSolveTurbulenceDissipationRate) then 
              write(10,*) 'TKEFarField(',i,')=',TKEFarField(i)
              write(10,*) 'TEDFarField(',i,')=',TEDFarField(i)
            endif
            if(LSolveTurbulenceSpecificDissipationRate) then 
              write(10,*) 'TKEFarField(',i,')=',TKEFarField(i)
              write(10,*) 'TOmegaFarField(',i,')=',TOmegaFarField(i)
            endif
            if(LSolveTurbulentKL) then 
              write(10,*) 'TKEFarField(',i,')=',TKEFarField(i)
              write(10,*) 'TurbulentKLFarField(',i,')=',
     *                                  TurbulentKLFarField(i)
            endif
            if(LSolveModifiedED) then 
              write(10,*) 'MEDFarField(',i,')=',MEDFarField(i)
            endif
          endif
        endif
      enddo      
c
c--- Set the boundary heat fluxes, Tinfinity, and Hinfinity
c    for use with von Neumann and Robin conditions
c
      if(LSolveEnergy) then
        do i=1,NumberOfBCSets
          write(10,*) 'HeatFlux(',i,':)='
          write(10,*) HeatFlux(i,:)
          write(10,*) 'Hinfinity(',i,':)='
          write(10,*) Hinfinity(i,:)
          write(10,*) 'Tinfinity(',i,':)='
          write(10,*) Tinfinity(i,:)
        enddo
      endif
c
c--- Set the value of the gas constant for compressible flows
c
      write(10,*) 'EquationOfState=',EquationOfState
      write(10,*) 'RGas=',RGas
      write(10,*) 'GammaGas=',GammaGas
      write(10,*) 'PrLaminar=',PrLaminar
c
c--- Set the constant density,specific heat, viscosity, and conductivity values
c
      write(10,*) 'Uinfinity=',Uinfinity
      write(10,*) 'Rhoinfinity=',Rhoinfinity
      write(10,*) 'LSoutherLand=',LSoutherLand
      write(10,*) 'SoutherlandLambda=',SoutherlandLambda
      write(10,*) 'SoutherlandC=',SoutherlandC
      write(10,*) 'ConstantSpecificHeat=',ConstantSpecificHeat
      write(10,*) 'ConstantViscosity=',ConstantViscosity
      write(10,*) 'ConstantConductivity=',ConstantConductivity
      write(10,*) 'ConstantDensity=',ConstantDensity
c
c--- Set values for anisotropic diffusion
c
      write(10,*) 'LanisotropicDiffusion=',LanisotropicDiffusion
      write(10,*) 'MethodDecomposeSprime=',MethodDecomposeSprime
      write(10,*) 'constantStagnationPressure=',
     *                                constantStagnationPressure
      write(10,*) 'constantStagnationTemperature=',
     *                                constantStagnationTemperature
c
      if(LSolveMomentum) then 
        do i=1,NumberOfBCSets
          write(10,*) 'BuVelocity(',i,':)=' 
          write(10,*) BuVelocity(i,:)
          write(10,'(/)')
        enddo
        do i=1,NumberOfBCSets
          write(10,*) 'BvVelocity(',i,':)=' 
          write(10,*) BvVelocity(i,:)
          write(10,'(/)')
        enddo
        do i=1,NumberOfBCSets
          write(10,*) 'BwVelocity(',i,':)=' 
          write(10,*) BwVelocity(i,:)
          write(10,'(/)')
        enddo
      endif
      if(LSolveContinuity) then 
        do i=1,NumberOfBCSets
          write(10,*) 'BPressure(',i,':)=' 
          write(10,*) BPressure(i,:)
          write(10,'(/)')
        enddo
      endif
      if(LSolveEnergy) then 
        do i=1,NumberOfBCSets
          write(10,*) 'BTemperature(',i,':)=' 
          write(10,*) BTemperature(i,:)
          write(10,'(/)')
        enddo
      endif
      if(LTurbulentFlow) then
        do i=1,NumberOfBCSets
          write(10,*) 'inletTypeT(',i,':)=',inletTypeT(i)        
          if(inletTypeT(i).eq.'specifiedil') then
            write(10,*) 'inletTurbulenceIntensity(',i,')=',
     *                             inletTurbulenceIntensity(i)
            write(10,*) 'inletTurbulenceLengthScale(',i,')=',
     *                             inletTurbulenceLengthScale(i)
          endif          
        enddo
        do i=1,NumberOfBCSets
          if(LSolveTurbulenceKineticEnergy) then 
            write(10,*) 'BTurbulentKE(',i,':)=' 
            write(10,*) BTurbulentKE(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulenceDissipationRate) then 
            write(10,*) 'BTurbulentED(',i,':)=' 
            write(10,*) BTurbulentED(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulenceSpecificDissipationRate) then 
            write(10,*) 'BTurbulentOmega(',i,':)=' 
            write(10,*) BTurbulentOmega(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulentKL) then 
            write(10,*) 'BTurbulentKL(',i,':)=' 
            write(10,*) BTurbulentKL(i,:)
            write(10,'(/)')
          endif
          if(LSolveModifiedED) then 
            write(10,*) 'BModifiedED(',i,':)=' 
            write(10,*) BModifiedED(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulenceGammaEquation) then 
            write(10,*) 'BTGamma(',i,':)=' 
            write(10,*) BTGamma(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulenceReynoldsThetaEquation) then 
            write(10,*) 'BTReTheta(',i,':)=' 
            write(10,*) BTReTheta(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulencefRelaxationEquation) then 
            write(10,*) 'BTfRelaxation(',i,':)=' 
            write(10,*) BTfRelaxation(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulenceV2Equation) then 
            write(10,*) 'BTurbulentV2(',i,':)=' 
            write(10,*) BTurbulentV2(i,:)
            write(10,'(/)')
          endif
          if(LSolveTurbulenceZetaEquation) then 
            write(10,*) 'BTurbulentZeta(',i,':)=' 
            write(10,*) BTurbulentZeta(i,:)
            write(10,'(/)')
          endif
          write(10,*) 'BTurbulentViscosity(',i,':)=' 
          write(10,*) BTurbulentViscosity(i,:)
          write(10,'(/)')
        enddo
      endif
      if(LSolveLambdaELEEquation) then 
        do i=1,NumberOfBCSets
          write(10,*) 'BLambdaELE(',i,':)=' 
          write(10,*) BLambdaELE(i,:)
          write(10,'(/)')
        enddo
      endif
c 
c--- Periodic boundary    
c      
      write(10,*) 'LRotationalPeriodicity=',LRotationalPeriodicity
      write(10,*) 'LTranslationalPeriodicity=',LTranslationalPeriodicity
      write(10,*) 'LcalculateBeta=',LcalculateBeta
      write(10,*) 'relaxBeta=',relaxBeta
      write(10,*) 'BetaIterations=',BetaIterations
      write(10,*) 'nIterCorrectBeta=',nIterCorrectBeta
      write(10,*) 'LPeriodicImplicit=',LPeriodicImplicit
      if(LRotationalPeriodicity) then
        do i=1,NumberOfBCSets
          write(10,*) 'a1Axis(',i,')=',a1Axis(i)
          write(10,*) 'a2Axis(',i,')=',a2Axis(i)
          write(10,*) 'a3Axis(',i,')=',a3Axis(i)
          write(10,*) 'PeriodicPair(',i,')=',PeriodicPair(i)
          write(10,*) 'theta(',i,')=',theta(1)
        enddo
      elseif(LTranslationalPeriodicity) then
          write(10,*) 'periodicBeta=',periodicBeta
          write(10,*) 'periodicMdot=',periodicMdot
        do i=1,NumberOfBCSets
          write(10,*) 'xTranslation(',i,')=',xTranslation(i)
          write(10,*) 'yTranslation(',i,')=',yTranslation(i)
          write(10,*) 'zTranslation(',i,')=',zTranslation(i)
          write(10,*) 'PeriodicPair(',i,')=',PeriodicPair(i)
        enddo
      endif
c
c----------------------------------------------------------------------------------------
c---  Set parameters for Lambda Euler Lagrange equation
c----------------------------------------------------------------------------------------
c
      if(LSolveLambdaELEEquation) then
        write(10,*) 'ASSolverLambdaELE=',ASSolverLambdaELE
        write(10,*) 'rrfLambdaELE=',rrfLambdaELE
        write(10,*) 'ASIterLambdaELE=',ASIterLambdaELE
        write(10,*) 'LRelaxLambdaELE=',LRelaxLambdaELE
        write(10,*) 'urfLambdaELE=',urfLambdaELE
        write(10,*) 'LMultigridLambdaELE=',LMultigridLambdaELE
        write(10,*) 'MethodCalcGradientLambdaELE=',
     *                            MethodCalcGradientLambdaELE
        write(10,*) 'MethodCalcGradientMomentum=',
     *                            MethodCalcGradientMomentum
        write(10,*) 'nIterGradientLambdaELE=',nIterGradientLambdaELE
        write(10,*) 'nIterGradientMomentum=',nIterGradientMomentum
        write(10,*) 'LimitGradientLambdaELE=',LimitGradientLambdaELE
        write(10,*) 'LimitGradientMomentum=',LimitGradientMomentum
        write(10,*) 'LimitGradientLambdaELEMethod=',
     *                            LimitGradientLambdaELEMethod
        write(10,*) 'LimitGradientMomentumMethod=',
     *                            LimitGradientMomentumMethod
        write(10,*) 'GradientInterpolationSchemeLambdaELE=',
     *                            GradientInterpolationSchemeLambdaELE
        write(10,*) 'alpha1Lambda=',alpha1Lambda
        write(10,*) 'alpha2Lambda=',alpha2Lambda
        do i=1,NumberOfBCSets
          write(10,*) 'BoundaryType(',i,')=',BoundaryType(i)
          write(10,*) 'wallTypeL(',i,')=',wallTypeL(i)
        enddo
        do i=1,NumberOfBCSets
          write(10,*) 'BLambdaELE(',i,':)=' 
          write(10,*) BLambdaELE(i,:)
          write(10,'(/)')
        enddo
      endif            
c
c----------------------------------------------------------------------------------------
c---  Set parameters for the rFields for free surface flows
c----------------------------------------------------------------------------------------
c
      if(NumberOfrFieldsToSolve.gt.0) then
        do i=1,NumberOfrFieldsToSolve
c
c--- Set physical properties related to the fluid used
c
          write(10,*) 'ConstantDensityrField(',i,')=',
     *                            ConstantDensityrField(i)
          write(10,*) 'ConstantSpecificHeatrField(',i,')=',
     *                            ConstantSpecificHeatrField(i)
          write(10,*) 'ConstantConductivityrField(',i,')=',
     *                            ConstantConductivityrField(i)
          write(10,*) 'rFieldName(',i,')=',rFieldName(i)
          write(10,*) 'LsolverField(',i,')=',LsolverField(i)
          write(10,*) 'LMultigridrField(',i,')=',LMultigridrField(i)
          write(10,*) 'ASSolverrField(',i,')=',ASSolverrField(i)
          write(10,*) 'rrfrField(',i,')=',rrfrField(i)
          write(10,*) 'ASIterrField(',i,')=',ASIterrField(i)
          write(10,*) 'LFalseTransientrField(',i,')=',
     *                            LFalseTransientrField(i)
          write(10,*) 'LRelaxrField(',i,')=',LRelaxrField(i)
          write(10,*) 'FalseDtrField(',i,')=',FalseDtrField(i)
          write(10,*) 'urfrField(',i,')=',urfrField(i)
          write(10,*) 'HRFrameworkrField(',i,')=',HRFrameworkrField(i)
          write(10,*) 'ConvectionSchemerField(',i,')=',
     *                            ConvectionSchemerField(i)
          write(10,*) 'BleedrField(',i,')=',BleedrField(i)
          write(10,*) 'MethodCalcGradientrField(',i,')=',
     *                            MethodCalcGradientrField(i)
          write(10,*) 'nIterGradientrField(',i,')=',
     *                            nIterGradientrField(i)
          write(10,*) 'LimitGradientrField(',i,')=',
     *                            LimitGradientrField(i)
          write(10,*) 'LimitGradientrFieldMethod(',i,')=',
     *                            LimitGradientrFieldMethod(i)
          write(10,*) 'GradientInterpolationSchemerField(',i,')=',
     *                            GradientInterpolationSchemerField(i)
          do j=1,NumberOfBCSets
            write(10,*) 'BoundaryType(',j,')=',BoundaryType(j)
            if(BoundaryType(j).eq.'wall') then
              write(10,*) 'wallTypeR(',j,',',i,')=',wallTypeR(j,i)
            endif
            if(BoundaryType(j).eq.'inlet') then
              write(10,*) 'BrField(',j,',:,',i,')=',BrField(j,:,i)
            endif
          enddo
c
c--- Set the boundary heat fluxes, Tinfinity, and Hinfinity
c    for use with von Neumann and Robin conditions for rFields
c
          do j=1,NumberOfBCSets
            write(10,*) 'rFieldFlux(',j,',:,',i,')='
            write(10,*) rFieldFlux(j,:,i)
          enddo
c
c--- Set the boundary and initial values of rFields
c
          do j=1,NumberOfBCSets
            write(10,*) 'BrField(',j,',:,',i,')='
            write(10,*) BrField(j,:,i)
          enddo
        enddo
      endif
c
c----------------------------------------------------------------------------------------
c---  Set parameters for the additional scalars to be solved
c----------------------------------------------------------------------------------------
c
      if(NumberOfScalarsToSolve.gt.0) then
        do i=1,NumberOfScalarsToSolve
c
c--- Set physical properties related to the fluid used
c
          write(10,*) 'ConstantSpecificHeatScalar(',i,')=',
     *                            ConstantSpecificHeatScalar(i)
          write(10,*) 'ConstantDiffusionCoefficientScalar(',i,')=',
     *                            ConstantDiffusionCoefficientScalar(i)
          write(10,*) 'ScalarName(',i,')=',ScalarName(i)
          write(10,*) 'LsolveScalar(',i,')=',LsolveScalar(i)
          write(10,*) 'LMultigridScalar(',i,')=',LMultigridScalar(i)
          write(10,*) 'ASSolverScalar(',i,')=',ASSolverScalar(i)
          write(10,*) 'rrfScalar(',i,')=',rrfScalar(i)
          write(10,*) 'ASIterScalar(',i,')=',ASIterScalar(i)
          write(10,*) 'LFalseTransientScalar(',i,')=',
     *                            LFalseTransientScalar(i)
          write(10,*) 'LRelaxScalar(',i,')=',LRelaxScalar(i)
          write(10,*) 'FalseDtScalar(',i,')=',FalseDtScalar(i)
          write(10,*) 'urfScalar(',i,')=',urfScalar(i)
          write(10,*) 'HRFrameworkScalar(',i,')=',HRFrameworkScalar(i)
          write(10,*) 'ConvectionSchemeScalar(',i,')=',
     *                            ConvectionSchemeScalar(i)
          write(10,*) 'BleedScalar(',i,')=',BleedScalar(i)
          write(10,*) 'MethodCalcGradientScalar(',i,')=',
     *                            MethodCalcGradientScalar(i)
          write(10,*) 'nIterGradientScalar(',i,')=',
     *                            nIterGradientScalar(i)
          write(10,*) 'LimitGradientScalar(',i,')=',
     *                            LimitGradientScalar(i)
          write(10,*) 'LimitGradientScalarMethod(',i,')=',
     *                            LimitGradientScalarMethod(i)
          write(10,*) 'GradientInterpolationSchemeScalar(',i,')=',
     *                            GradientInterpolationSchemeScalar(i)
          do j=1,NumberOfBCSets
            write(10,*) 'BoundaryType(',j,')=',BoundaryType(j)
            if(BoundaryType(j).eq.'wall') then
              write(10,*) 'wallTypeS(',j,',',i,')=',wallTypeS(j,i)
            endif
            if(BoundaryType(j).eq.'inlet') then
              write(10,*) 'BScalar(',j,',:,',i,')=',BScalar(j,:,i)
            endif
          enddo
          do j=1,NumberofPointSources          
            write(10,*) 'xLocationOfPointSource(',j,')=',
     *                                xLocationOfPointSource(j)
            write(10,*) 'yLocationOfPointSource(',j,')=',
     *                                yLocationOfPointSource(j)
            write(10,*) 'zLocationOfPointSource(',j,')=',
     *                                zLocationOfPointSource(j)
            write(10,*) 'SbPointSourceScalar(',j,')=',
     *                                SbPointSourceScalar(j,i)
            write(10,*) 'ScPointSourceScalar(',j,')=',
     *                                ScPointSourceScalar(j,i)
          enddo
c
c--- Set the boundary heat fluxes, Tinfinity, and Hinfinity
c    for use with von Neumann and Robin conditions for rFields
c
          do j=1,NumberOfBCSets
            write(10,*) 'ScalarFlux(',j,',:,',i,')='
            write(10,*) ScalarFlux(j,:,i)
          enddo
c
c--- Set the boundary and initial values of rFields
c
          do j=1,NumberOfBCSets
            write(10,*) 'BScalar(',j,',:,',i,')='
            write(10,*) BScalar(j,:,i)
          enddo
        enddo
      endif
c          
      close(10)
c
      return
      end