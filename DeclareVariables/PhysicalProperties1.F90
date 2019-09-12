MODULE PhysicalProperties1
implicit none

character*20, save :: EquationOfState  !'idealgas'  'tait'  'constant'
double precision, save :: ThetaEOS=7.15  ! for water=7.15
double precision, save :: RGas,GammaGas,PrLaminar
double precision, save :: CoefficientOfThermalExpansion=0.003333333333333
double precision, save :: ReferenceDensity=1000.
double precision, save :: ReferenceTemperature=300.
double precision, save :: ReferencePressure=100000.
double precision, save :: GravityX,GravityY,GravityZ

logical, save :: LSoutherLand
double precision, save :: SoutherlandLambda,SoutherlandC

double precision, save, dimension(:), allocatable :: Density
double precision, save, dimension(:), allocatable :: SpecificHeat
double precision, save, dimension(:), allocatable :: Conductivity
double precision, save, dimension(:), allocatable :: Viscosity
double precision, save, dimension(:), allocatable :: Densityf
double precision, save, dimension(:), allocatable :: drhodP
double precision, save, dimension(:,:), allocatable :: BdrhodP

double precision, save, dimension(:), allocatable :: Conductivity11
double precision, save, dimension(:), allocatable :: Conductivity12
double precision, save, dimension(:), allocatable :: Conductivity13
double precision, save, dimension(:), allocatable :: Conductivity22
double precision, save, dimension(:), allocatable :: Conductivity23
double precision, save, dimension(:), allocatable :: Conductivity33
double precision, save, dimension(:,:), allocatable :: BConductivity11
double precision, save, dimension(:,:), allocatable :: BConductivity12
double precision, save, dimension(:,:), allocatable :: BConductivity13
double precision, save, dimension(:,:), allocatable :: BConductivity22
double precision, save, dimension(:,:), allocatable :: BConductivity23
double precision, save, dimension(:,:), allocatable :: BConductivity33

double precision, save, dimension(:,:), allocatable :: BDensity
double precision, save, dimension(:,:), allocatable :: BSpecificHeat
double precision, save, dimension(:,:), allocatable :: BConductivity
double precision, save, dimension(:,:), allocatable :: BViscosity

double precision, save, dimension(:), allocatable :: DensityStar
double precision, save, dimension(:,:), allocatable :: BDensityStar

double precision, save, dimension(:), allocatable :: DensityOld
double precision, save, dimension(:,:), allocatable :: BDensityOld

double precision, save, dimension(:), allocatable :: DensityOldOld
double precision, save, dimension(:,:), allocatable :: BDensityOldOld

double precision, save, dimension(:), allocatable :: SpecificHeatOld
double precision, save, dimension(:,:), allocatable :: BSpecificHeatOld

double precision, save, dimension(:), allocatable :: SpecificHeatOldOld
double precision, save, dimension(:,:), allocatable :: BSpecificHeatOldOld

double precision, save, dimension(:), allocatable :: DensGradx
double precision, save, dimension(:), allocatable :: DensGrady
double precision, save, dimension(:), allocatable :: DensGradz
double precision, save, dimension(:), allocatable :: SHeatGradx
double precision, save, dimension(:), allocatable :: SHeatGrady
double precision, save, dimension(:), allocatable :: SHeatGradz
double precision, save, dimension(:), allocatable :: cpterm

double precision, save, dimension(:,:), allocatable :: BDensGradx
double precision, save, dimension(:,:), allocatable :: BDensGrady
double precision, save, dimension(:,:), allocatable :: BDensGradz
double precision, save, dimension(:,:), allocatable :: BSHeatGradx
double precision, save, dimension(:,:), allocatable :: BSHeatGrady
double precision, save, dimension(:,:), allocatable :: BSHeatGradz


character*16, save, dimension(:), allocatable :: InterpolationSchemeGamaScalar
logical, save, dimension(:), allocatable :: LSoutherLandScalar
character*20, save, dimension(:), allocatable :: EquationOfStateScalar
double precision, save, dimension(:), allocatable :: RGasScalar
double precision, save, dimension(:), allocatable :: GammaGasScalar
double precision, save, dimension(:), allocatable :: PrLaminarScalar


double precision, save, dimension(:,:), allocatable :: DiffusionCoefficient
double precision, save, dimension(:,:,:), allocatable :: BDiffusionCoefficient

double precision, save, dimension(:,:), allocatable :: DiffusionCoefficient11
double precision, save, dimension(:,:), allocatable :: DiffusionCoefficient12
double precision, save, dimension(:,:), allocatable :: DiffusionCoefficient13
double precision, save, dimension(:,:), allocatable :: DiffusionCoefficient22
double precision, save, dimension(:,:), allocatable :: DiffusionCoefficient23
double precision, save, dimension(:,:), allocatable :: DiffusionCoefficient33
double precision, save, dimension(:,:,:), allocatable :: BDiffusionCoefficient11
double precision, save, dimension(:,:,:), allocatable :: BDiffusionCoefficient12
double precision, save, dimension(:,:,:), allocatable :: BDiffusionCoefficient13
double precision, save, dimension(:,:,:), allocatable :: BDiffusionCoefficient22
double precision, save, dimension(:,:,:), allocatable :: BDiffusionCoefficient23
double precision, save, dimension(:,:,:), allocatable :: BDiffusionCoefficient33


double precision, save, dimension(:,:), allocatable :: SpecificHeatScalar
double precision, save, dimension(:,:,:), allocatable :: BSpecificHeatScalar
double precision, save, dimension(:,:), allocatable :: SpecificHeatScalarOld
double precision, save, dimension(:,:,:), allocatable :: BSpecificHeatScalarOld
double precision, save, dimension(:,:), allocatable :: SpecificHeatScalarOldOld
double precision, save, dimension(:,:,:), allocatable :: BSpecificHeatScalarOldOld

double precision, save, dimension(:), allocatable :: eDiffCoefficient
double precision, save, dimension(:,:), allocatable :: BeDiffCoefficient
double precision, save, dimension(:), allocatable :: TurbulentViscosity
double precision, save, dimension(:,:), allocatable :: BTurbulentViscosity
double precision, save, dimension(:), allocatable :: TurbulentViscosity1
double precision, save, dimension(:,:), allocatable :: BTurbulentViscosity1
!
double precision, save, dimension(:), allocatable :: GamaFace
double precision, save, dimension(:,:), allocatable :: BGamaFace
!
!Volume of fluid variables
!
double precision, save, dimension(:,:), allocatable :: DensityrField
double precision, save, dimension(:,:), allocatable :: SpecificHeatrField
double precision, save, dimension(:,:), allocatable :: ConductivityrField
double precision, save, dimension(:,:), allocatable :: ViscosityrField

double precision, save, dimension(:,:,:), allocatable :: BDensityrField
double precision, save, dimension(:,:,:), allocatable :: BSpecificHeatrField
double precision, save, dimension(:,:,:), allocatable :: BConductivityrField
double precision, save, dimension(:,:,:), allocatable :: BViscosityrField

double precision, save, dimension(:,:), allocatable :: SpecificHeatrFieldOld
double precision, save, dimension(:,:,:), allocatable :: BSpecificHeatrFieldOld
double precision, save, dimension(:,:), allocatable :: SpecificHeatrFieldOldOld
double precision, save, dimension(:,:,:), allocatable :: BSpecificHeatrFieldOldOld

logical, save, dimension(:), allocatable :: LSoutherLandrField
character*20, save, dimension(:), allocatable :: EquationOfStaterField
double precision, save, dimension(:), allocatable :: RGasrField
double precision, save, dimension(:), allocatable :: GammaGasrField
double precision, save, dimension(:), allocatable :: PrLaminarrField
double precision, save, dimension(:), allocatable :: SoutherlandLambdarField
double precision, save, dimension(:), allocatable :: SoutherlandCrField

end MODULE PhysicalProperties1