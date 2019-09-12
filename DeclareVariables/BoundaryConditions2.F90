MODULE BoundaryConditions2

implicit none

integer, save  :: Iwall,IWallSlip,IWallnoSlip,IWallDirichlet,IWallVonNeumann,IWallRobin

integer, save  :: Iinlet,Iinletsupersonic,IinletSpecifiedVelocity,IinletSpecifiedMassFlowRate
integer, save  :: IinletSpecifiedStaticPressure,IinletSpecifiedStagnationPressure
integer, save  :: IinletSpecifiedStaticTemperature,IinletSpecifiedStagnationTemperature

integer, save  :: Ioutlet,Ioutletsupersonic,IoutletspecifiedVelocity
integer, save  :: IoutletSpecifiedStaticPressure,IoutletSpecifiedMassFlowRate
integer, save  :: IoutletSpecifiedAverageStaticPressure,IoutletSpecifiedResistance
integer, save  :: IoutletFullyDeveloped,IoutletFullyDevelopedEnergy

integer, save  :: IpressureFarField
integer, save  :: IoutletTransmissive

integer, save  :: Isymmetry,Iperiodic,Iaxis


integer, save, dimension(:), allocatable :: IWallSlipOwner
integer, save, dimension(:), allocatable :: IWallSlipNumberOfBCSets
integer, save, dimension(:), allocatable :: IWallSlipNBFaces

integer, save, dimension(:), allocatable :: IWallnoSlipOwner
integer, save, dimension(:), allocatable :: IWallnoSlipNumberOfBCSets
integer, save, dimension(:), allocatable :: IWallnoSlipNBFaces

integer, save, dimension(:), allocatable :: IWallDirichletOwner
integer, save, dimension(:), allocatable :: IWallDirichletNumberOfBCSets
integer, save, dimension(:), allocatable :: IWallDirichletNBFaces

integer, save, dimension(:), allocatable :: IWallVonNeumannOwner
integer, save, dimension(:), allocatable :: IWallVonNeumannNumberOfBCSets
integer, save, dimension(:), allocatable :: IWallVonNeumannNBFaces

integer, save, dimension(:), allocatable :: IWallRobinOwner
integer, save, dimension(:), allocatable :: IWallRobinNumberOfBCSets
integer, save, dimension(:), allocatable :: IWallRobinNBFaces
double precision, save, dimension(:), allocatable :: TinfinityRobin
double precision, save, dimension(:), allocatable :: HinfinityRobin

integer, save, dimension(:), allocatable :: IinletsupersonicOwner
integer, save, dimension(:), allocatable :: IinletsupersonicNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletsupersonicNBFaces

integer, save, dimension(:), allocatable :: IinletSpecifiedVelocityOwner
integer, save, dimension(:), allocatable :: IinletSpecifiedVelocityNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletSpecifiedVelocityNBFaces

integer, save, dimension(:), allocatable :: IinletSpecifiedMassFlowRateOwner
integer, save, dimension(:), allocatable :: IinletSpecifiedMassFlowRateNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletSpecifiedMassFlowRateNBFaces

integer, save, dimension(:), allocatable :: IinletSpecifiedStaticPressureOwner
integer, save, dimension(:), allocatable :: IinletSpecifiedStaticPressureNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletSpecifiedStaticPressureNBFaces

integer, save, dimension(:), allocatable :: IinletSpecifiedStagnationPressureOwner
integer, save, dimension(:), allocatable :: IinletSpecifiedStagnationPressureNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletSpecifiedStagnationPressureNBFaces

integer, save, dimension(:), allocatable :: IinletSpecifiedStaticTemperatureOwner
integer, save, dimension(:), allocatable :: IinletSpecifiedStaticTemperatureNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletSpecifiedStaticTemperatureNBFaces

integer, save, dimension(:), allocatable :: IinletSpecifiedStagnationTemperatureOwner
integer, save, dimension(:), allocatable :: IinletSpecifiedStagnationTemperatureNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletSpecifiedStagnationTemperatureNBFaces

integer, save, dimension(:), allocatable :: IoutletsupersonicOwner
integer, save, dimension(:), allocatable :: IoutletsupersonicNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletsupersonicNBFaces

integer, save, dimension(:), allocatable :: IoutletspecifiedVelocityOwner
integer, save, dimension(:), allocatable :: IoutletspecifiedVelocityNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletspecifiedVelocityNBFaces

integer, save, dimension(:), allocatable :: IoutletSpecifiedStaticPressureOwner
integer, save, dimension(:), allocatable :: IoutletSpecifiedStaticPressureNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletSpecifiedStaticPressureNBFaces

integer, save, dimension(:), allocatable :: IoutletSpecifiedAverageStaticPressureOwner
integer, save, dimension(:), allocatable :: IoutletSpecifiedAverageStaticPressureNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletSpecifiedAverageStaticPressureNBFaces

integer, save, dimension(:), allocatable :: IoutletSpecifiedResistanceOwner
integer, save, dimension(:), allocatable :: IoutletSpecifiedResistanceNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletSpecifiedResistanceNBFaces

integer, save, dimension(:), allocatable :: IoutletSpecifiedMassFlowRateOwner
integer, save, dimension(:), allocatable :: IoutletSpecifiedMassFlowRateNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletSpecifiedMassFlowRateNBFaces

integer, save, dimension(:), allocatable :: IoutletFullyDevelopedOwner
integer, save, dimension(:), allocatable :: IoutletFullyDevelopedNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletFullyDevelopedNBFaces

integer, save, dimension(:), allocatable :: IoutletFullyDevelopedEnergyOwner
integer, save, dimension(:), allocatable :: IoutletFullyDevelopedEnergyNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletFullyDevelopedEnergyNBFaces

integer, save, dimension(:), allocatable :: IoutletTransmissiveOwner
integer, save, dimension(:), allocatable :: IoutletTransmissiveNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletTransmissiveNBFaces

integer, save, dimension(:), allocatable :: IpressurefarfieldOwner
integer, save, dimension(:), allocatable :: IpressurefarfieldNumberOfBCSets
integer, save, dimension(:), allocatable :: IpressurefarfieldNBFaces
double precision, save, dimension(:), allocatable :: uVelocityFarField
double precision, save, dimension(:), allocatable :: vVelocityFarField
double precision, save, dimension(:), allocatable :: wVelocityFarField
double precision, save, dimension(:), allocatable :: SpecificHeatFarField
double precision, save, dimension(:), allocatable :: pressureFarField
double precision, save, dimension(:), allocatable :: TemperatureFarField
double precision, save, dimension(:), allocatable :: MachFarField
double precision, save, dimension(:), allocatable :: xFlowDirectionFarField
double precision, save, dimension(:), allocatable :: yFlowDirectionFarField
double precision, save, dimension(:), allocatable :: zFlowDirectionFarField
double precision, save, dimension(:), allocatable :: TKEFarField
double precision, save, dimension(:), allocatable :: TEDFarField
double precision, save, dimension(:), allocatable :: TOmegaFarField
double precision, save, dimension(:), allocatable :: TurbulentKLFarField
double precision, save, dimension(:), allocatable :: MEDFarField

integer, save, dimension(:), allocatable :: IsymmetryOwner
integer, save, dimension(:), allocatable :: IsymmetryNumberOfBCSets
integer, save, dimension(:), allocatable :: IsymmetryNBFaces


integer, save, dimension(:), allocatable :: IperiodicOwner
integer, save, dimension(:), allocatable :: IperiodicNumberOfBCSets
integer, save, dimension(:), allocatable :: IperiodicNBFaces
integer, save, dimension(:), allocatable :: ElementPeriodicFace
              
logical, save :: LRotationalPeriodicity=.false.
logical, save :: LTranslationalPeriodicity=.false.
logical, save :: LcalculateBeta=.false.
logical, save :: LPeriodicImplicit=.false.
double precision, save :: periodicBeta
integer, save :: BetaIterations
integer, save :: nIterCorrectBeta
double precision, save :: relaxBeta
double precision, save :: periodicMdot
double precision, save, dimension(:), allocatable :: bcSave

integer, save, dimension(:), allocatable :: PeriodicPair
integer, save, dimension(:,:), allocatable :: Icorrespondingface
double precision, save, dimension(:), allocatable :: theta
double precision, save, dimension(:), allocatable :: a1r
double precision, save, dimension(:), allocatable :: b1r
double precision, save, dimension(:), allocatable :: c1r
double precision, save, dimension(:), allocatable :: a2r
double precision, save, dimension(:), allocatable :: b2r
double precision, save, dimension(:), allocatable :: c2r
double precision, save, dimension(:), allocatable :: a3r
double precision, save, dimension(:), allocatable :: b3r
double precision, save, dimension(:), allocatable :: c3r
double precision, save, dimension(:), allocatable :: a1Axis
double precision, save, dimension(:), allocatable :: a2Axis
double precision, save, dimension(:), allocatable :: a3Axis
double precision, save, dimension(:), allocatable :: xTranslation
double precision, save, dimension(:), allocatable :: yTranslation
double precision, save, dimension(:), allocatable :: zTranslation
double precision, save ::  xTranslationUse
double precision, save ::  yTranslationUse
double precision, save ::  zTranslationUse


integer, save, dimension(:), allocatable :: IaxisOwner
integer, save, dimension(:), allocatable :: IaxisNumberOfBCSets
integer, save, dimension(:), allocatable :: IaxisNBFaces




end MODULE BoundaryConditions2
