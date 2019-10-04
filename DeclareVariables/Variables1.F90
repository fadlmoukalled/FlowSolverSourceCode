MODULE Variables1
implicit none

double precision, save, dimension(:), allocatable :: uVelocity
double precision, save, dimension(:), allocatable :: vVelocity
double precision, save, dimension(:), allocatable :: wVelocity
double precision, save, dimension(:), allocatable :: PressureC
double precision, save, dimension(:), allocatable :: Pressure
double precision, save, dimension(:), allocatable :: Temperature
double precision, save, dimension(:), allocatable :: Htotal
double precision, save, dimension(:), allocatable :: mdot
double precision, save, dimension(:), allocatable :: mdotOld
double precision, save, dimension(:), allocatable :: mdotOldOld
double precision, save, dimension(:), allocatable :: effdiv

double precision, save, dimension(:), allocatable :: uVelocityOld
double precision, save, dimension(:), allocatable :: vVelocityOld
double precision, save, dimension(:), allocatable :: wVelocityOld
double precision, save, dimension(:), allocatable :: PressureOld
double precision, save, dimension(:), allocatable :: TemperatureOld
double precision, save, dimension(:), allocatable :: HtotalOld

double precision, save, dimension(:), allocatable :: uVelocityOldOld
double precision, save, dimension(:), allocatable :: vVelocityOldOld
double precision, save, dimension(:), allocatable :: wVelocityOldOld
double precision, save, dimension(:), allocatable :: PressureOldOld
double precision, save, dimension(:), allocatable :: TemperatureOldOld
double precision, save, dimension(:), allocatable :: HtotalOldOld

double precision, save, dimension(:,:), allocatable :: BuVelocity
double precision, save, dimension(:,:), allocatable :: BvVelocity
double precision, save, dimension(:,:), allocatable :: BwVelocity
double precision, save, dimension(:,:), allocatable :: BPressureC
double precision, save, dimension(:,:), allocatable :: BPressure
double precision, save, dimension(:,:), allocatable :: BTemperature
double precision, save, dimension(:,:), allocatable :: BHtotal
double precision, save, dimension(:,:), allocatable :: BStagnationPressure
double precision, save, dimension(:,:), allocatable :: BStagnationTemperature
double precision, save, dimension(:,:), allocatable :: Bmdot
double precision, save, dimension(:,:), allocatable :: BmdotOld
double precision, save, dimension(:,:), allocatable :: BmdotOldOld

double precision, save, dimension(:,:), allocatable :: BuVelocityOld
double precision, save, dimension(:,:), allocatable :: BvVelocityOld
double precision, save, dimension(:,:), allocatable :: BwVelocityOld
double precision, save, dimension(:,:), allocatable :: BPressureOld
double precision, save, dimension(:,:), allocatable :: BTemperatureOld
double precision, save, dimension(:,:), allocatable :: BHtotalOld

double precision, save, dimension(:,:), allocatable :: BuVelocityOldOld
double precision, save, dimension(:,:), allocatable :: BvVelocityOldOld
double precision, save, dimension(:,:), allocatable :: BwVelocityOldOld
double precision, save, dimension(:,:), allocatable :: BPressureOldOld
double precision, save, dimension(:,:), allocatable :: BTemperatureOldOld
double precision, save, dimension(:,:), allocatable :: BHtotalOldOld

double precision, save, dimension(:), allocatable :: uVelGradx
double precision, save, dimension(:), allocatable :: uVelGrady
double precision, save, dimension(:), allocatable :: uVelGradz
double precision, save, dimension(:), allocatable :: vVelGradx
double precision, save, dimension(:), allocatable :: vVelGrady
double precision, save, dimension(:), allocatable :: vVelGradz
double precision, save, dimension(:), allocatable :: wVelGradx
double precision, save, dimension(:), allocatable :: wVelGrady
double precision, save, dimension(:), allocatable :: wVelGradz
double precision, save, dimension(:), allocatable :: PressGradx
double precision, save, dimension(:), allocatable :: PressGrady
double precision, save, dimension(:), allocatable :: PressGradz
double precision, save, dimension(:), allocatable :: PressCGradx
double precision, save, dimension(:), allocatable :: PressCGrady
double precision, save, dimension(:), allocatable :: PressCGradz
double precision, save, dimension(:), allocatable :: TempGradx
double precision, save, dimension(:), allocatable :: TempGrady
double precision, save, dimension(:), allocatable :: TempGradz
double precision, save, dimension(:), allocatable :: HtotalGradx
double precision, save, dimension(:), allocatable :: HtotalGrady
double precision, save, dimension(:), allocatable :: HtotalGradz

double precision, save, dimension(:,:), allocatable :: BuVelGradx
double precision, save, dimension(:,:), allocatable :: BuVelGrady
double precision, save, dimension(:,:), allocatable :: BuVelGradz
double precision, save, dimension(:,:), allocatable :: BvVelGradx
double precision, save, dimension(:,:), allocatable :: BvVelGrady
double precision, save, dimension(:,:), allocatable :: BvVelGradz
double precision, save, dimension(:,:), allocatable :: BwVelGradx
double precision, save, dimension(:,:), allocatable :: BwVelGrady
double precision, save, dimension(:,:), allocatable :: BwVelGradz
double precision, save, dimension(:,:), allocatable :: BPressGradx
double precision, save, dimension(:,:), allocatable :: BPressGrady
double precision, save, dimension(:,:), allocatable :: BPressGradz
double precision, save, dimension(:,:), allocatable :: BPressCGradx
double precision, save, dimension(:,:), allocatable :: BPressCGrady
double precision, save, dimension(:,:), allocatable :: BPressCGradz
double precision, save, dimension(:,:), allocatable :: BTempGradx
double precision, save, dimension(:,:), allocatable :: BTempGrady
double precision, save, dimension(:,:), allocatable :: BTempGradz
double precision, save, dimension(:,:), allocatable :: BHtotalGradx
double precision, save, dimension(:,:), allocatable :: BHtotalGrady
double precision, save, dimension(:,:), allocatable :: BHtotalGradz

double precision, save, dimension(:), allocatable :: uVelGradfx
double precision, save, dimension(:), allocatable :: uVelGradfy
double precision, save, dimension(:), allocatable :: uVelGradfz
double precision, save, dimension(:), allocatable :: vVelGradfx
double precision, save, dimension(:), allocatable :: vVelGradfy
double precision, save, dimension(:), allocatable :: vVelGradfz
double precision, save, dimension(:), allocatable :: wVelGradfx
double precision, save, dimension(:), allocatable :: wVelGradfy
double precision, save, dimension(:), allocatable :: wVelGradfz
double precision, save, dimension(:), allocatable :: PressGradfx
double precision, save, dimension(:), allocatable :: PressGradfy
double precision, save, dimension(:), allocatable :: PressGradfz
double precision, save, dimension(:), allocatable :: PressCGradfx
double precision, save, dimension(:), allocatable :: PressCGradfy
double precision, save, dimension(:), allocatable :: PressCGradfz
double precision, save, dimension(:), allocatable :: TempGradfx
double precision, save, dimension(:), allocatable :: TempGradfy
double precision, save, dimension(:), allocatable :: TempGradfz
double precision, save, dimension(:), allocatable :: HtotalGradfx
double precision, save, dimension(:), allocatable :: HtotalGradfy
double precision, save, dimension(:), allocatable :: HtotalGradfz
double precision, save, dimension(:), allocatable :: PressureWork

double precision, save, dimension(:), allocatable :: Du1Velocity
double precision, save, dimension(:), allocatable :: Du2Velocity
double precision, save, dimension(:), allocatable :: uVelocityStar
double precision, save, dimension(:), allocatable :: Dv1Velocity
double precision, save, dimension(:), allocatable :: Dv2Velocity
double precision, save, dimension(:), allocatable :: vVelocityStar
double precision, save, dimension(:), allocatable :: Dw1Velocity
double precision, save, dimension(:), allocatable :: Dw2Velocity
double precision, save, dimension(:), allocatable :: wVelocityStar

double precision, save, dimension(:,:), allocatable :: xVeldirection
double precision, save, dimension(:,:), allocatable :: yVeldirection
double precision, save, dimension(:,:), allocatable :: zVeldirection

double precision, save, dimension(:), allocatable :: MachNumber
double precision, save, dimension(:,:), allocatable :: BMachNumber

double precision, save, dimension(:), allocatable :: Buoyancyx
double precision, save, dimension(:), allocatable :: Buoyancyy
double precision, save, dimension(:), allocatable :: Buoyancyz
double precision, save, dimension(:), allocatable :: Buoyancyfx
double precision, save, dimension(:), allocatable :: Buoyancyfy
double precision, save, dimension(:), allocatable :: Buoyancyfz
double precision, save, dimension(:,:), allocatable :: BBuoyancyx
double precision, save, dimension(:,:), allocatable :: BBuoyancyy
double precision, save, dimension(:,:), allocatable :: BBuoyancyz
!
!--- Lambda Euler Lagrange equation variables
!
double precision, save, dimension(:), allocatable :: LambdaELE
double precision, save, dimension(:,:), allocatable :: BLambdaELE
double precision, save, dimension(:), allocatable :: LambdaELEGradx
double precision, save, dimension(:), allocatable :: LambdaELEGrady
double precision, save, dimension(:), allocatable :: LambdaELEGradz
double precision, save, dimension(:,:), allocatable :: BLambdaELEGradx
double precision, save, dimension(:,:), allocatable :: BLambdaELEGrady
double precision, save, dimension(:,:), allocatable :: BLambdaELEGradz
double precision, save, dimension(:), allocatable :: LambdaELEGradfx
double precision, save, dimension(:), allocatable :: LambdaELEGradfy
double precision, save, dimension(:), allocatable :: LambdaELEGradfz
double precision, save, dimension(:), allocatable :: InitialVelDivergence
double precision, save, dimension(:), allocatable :: FinalVelDivergence
double precision, save, dimension(:,:), allocatable :: BInitialVelDivergence
double precision, save, dimension(:,:), allocatable :: BFinalVelDivergence

!
!--- Turbulence variables
!
double precision, save, dimension(:), allocatable :: TurbulentKE
double precision, save, dimension(:), allocatable :: TurbulentED
double precision, save, dimension(:), allocatable :: TurbulentOmega
double precision, save, dimension(:,:), allocatable :: BTurbulentKE
double precision, save, dimension(:,:), allocatable :: BTurbulentED
double precision, save, dimension(:,:), allocatable :: BTurbulentOmega

double precision, save, dimension(:), allocatable :: TurbulentKEOld
double precision, save, dimension(:), allocatable :: TurbulentEDOld
double precision, save, dimension(:), allocatable :: TurbulentOmegaOld
double precision, save, dimension(:,:), allocatable :: BTurbulentKEOld
double precision, save, dimension(:,:), allocatable :: BTurbulentEDOld
double precision, save, dimension(:,:), allocatable :: BTurbulentOmegaOld

double precision, save, dimension(:), allocatable :: TurbulentKEOldOld
double precision, save, dimension(:), allocatable :: TurbulentEDOldOld
double precision, save, dimension(:), allocatable :: TurbulentOmegaOldOld
double precision, save, dimension(:,:), allocatable :: BTurbulentKEOldOld
double precision, save, dimension(:,:), allocatable :: BTurbulentEDOldOld
double precision, save, dimension(:,:), allocatable :: BTurbulentOmegaOldOld

double precision, save, dimension(:), allocatable :: TKEGradx
double precision, save, dimension(:), allocatable :: TKEGrady
double precision, save, dimension(:), allocatable :: TKEGradz
double precision, save, dimension(:), allocatable :: TEDGradx
double precision, save, dimension(:), allocatable :: TEDGrady
double precision, save, dimension(:), allocatable :: TEDGradz
double precision, save, dimension(:), allocatable :: TOmegaGradx
double precision, save, dimension(:), allocatable :: TOmegaGrady
double precision, save, dimension(:), allocatable :: TOmegaGradz
double precision, save, dimension(:,:), allocatable :: BTKEGradx
double precision, save, dimension(:,:), allocatable :: BTKEGrady
double precision, save, dimension(:,:), allocatable :: BTKEGradz
double precision, save, dimension(:,:), allocatable :: BTEDGradx
double precision, save, dimension(:,:), allocatable :: BTEDGrady
double precision, save, dimension(:,:), allocatable :: BTEDGradz
double precision, save, dimension(:,:), allocatable :: BTOmegaGradx
double precision, save, dimension(:,:), allocatable :: BTOmegaGrady
double precision, save, dimension(:,:), allocatable :: BTOmegaGradz

double precision, save, dimension(:), allocatable :: TKEGradfx
double precision, save, dimension(:), allocatable :: TKEGradfy
double precision, save, dimension(:), allocatable :: TKEGradfz
double precision, save, dimension(:), allocatable :: TEDGradfx
double precision, save, dimension(:), allocatable :: TEDGradfy
double precision, save, dimension(:), allocatable :: TEDGradfz
double precision, save, dimension(:), allocatable :: TOmegaGradfx
double precision, save, dimension(:), allocatable :: TOmegaGradfy
double precision, save, dimension(:), allocatable :: TOmegaGradfz

double precision, save, dimension(:), allocatable :: TurbulenceProduction
double precision, save, dimension(:,:), allocatable :: BTurbulenceProduction
double precision, save, dimension(:), allocatable :: TurbulenceProductionB
double precision, save, dimension(:,:), allocatable :: BTurbulenceProductionB

double precision, save, dimension(:), allocatable :: TGammaEff
double precision, save, dimension(:,:), allocatable :: BTGammaEff
double precision, save, dimension(:), allocatable :: TGamma
double precision, save, dimension(:,:), allocatable :: BTGamma
double precision, save, dimension(:), allocatable :: TGammaOld
double precision, save, dimension(:,:), allocatable :: BTGammaOld
double precision, save, dimension(:), allocatable :: TGammaOldOld
double precision, save, dimension(:,:), allocatable :: BTGammaOldOld
double precision, save, dimension(:), allocatable :: TGammaGradx
double precision, save, dimension(:,:), allocatable :: BTGammaGradx
double precision, save, dimension(:), allocatable :: TGammaGrady
double precision, save, dimension(:,:), allocatable :: BTGammaGrady
double precision, save, dimension(:), allocatable :: TGammaGradz
double precision, save, dimension(:,:), allocatable :: BTGammaGradz
double precision, save, dimension(:), allocatable :: TGammaGradfx
double precision, save, dimension(:), allocatable :: TGammaGradfy
double precision, save, dimension(:), allocatable :: TGammaGradfz

double precision, save, dimension(:), allocatable :: TReTheta
double precision, save, dimension(:,:), allocatable :: BTReTheta
double precision, save, dimension(:), allocatable :: TReThetaOld
double precision, save, dimension(:,:), allocatable :: BTReThetaOld
double precision, save, dimension(:), allocatable :: TReThetaOldOld
double precision, save, dimension(:,:), allocatable :: BTReThetaOldOld
double precision, save, dimension(:), allocatable :: TReThetaGradx
double precision, save, dimension(:,:), allocatable :: BTReThetaGradx
double precision, save, dimension(:), allocatable :: TReThetaGrady
double precision, save, dimension(:,:), allocatable :: BTReThetaGrady
double precision, save, dimension(:), allocatable :: TReThetaGradz
double precision, save, dimension(:,:), allocatable :: BTReThetaGradz
double precision, save, dimension(:), allocatable :: TReThetaGradfx
double precision, save, dimension(:), allocatable :: TReThetaGradfy
double precision, save, dimension(:), allocatable :: TReThetaGradfz

double precision, save, dimension(:), allocatable :: ModifiedED
double precision, save, dimension(:,:), allocatable :: BModifiedED
double precision, save, dimension(:), allocatable :: ModifiedEDOld
double precision, save, dimension(:,:), allocatable :: BModifiedEDOld
double precision, save, dimension(:), allocatable :: ModifiedEDOldOld
double precision, save, dimension(:,:), allocatable :: BModifiedEDOldOld
double precision, save, dimension(:), allocatable :: ModifiedEDGradx
double precision, save, dimension(:,:), allocatable :: BModifiedEDGradx
double precision, save, dimension(:), allocatable :: ModifiedEDGrady
double precision, save, dimension(:,:), allocatable :: BModifiedEDGrady
double precision, save, dimension(:), allocatable :: ModifiedEDGradz
double precision, save, dimension(:,:), allocatable :: BModifiedEDGradz
double precision, save, dimension(:), allocatable :: ModifiedEDGradfx
double precision, save, dimension(:), allocatable :: ModifiedEDGradfy
double precision, save, dimension(:), allocatable :: ModifiedEDGradfz

double precision, save, dimension(:), allocatable :: TurbulentKL
double precision, save, dimension(:,:), allocatable :: BTurbulentKL
double precision, save, dimension(:), allocatable :: TurbulentKLOld
double precision, save, dimension(:,:), allocatable :: BTurbulentKLOld
double precision, save, dimension(:), allocatable :: TurbulentKLOldOld
double precision, save, dimension(:,:), allocatable :: BTurbulentKLOldOld
double precision, save, dimension(:), allocatable :: TurbulentKLGradx
double precision, save, dimension(:,:), allocatable :: BTurbulentKLGradx
double precision, save, dimension(:), allocatable :: TurbulentKLGrady
double precision, save, dimension(:,:), allocatable :: BTurbulentKLGrady
double precision, save, dimension(:), allocatable :: TurbulentKLGradz
double precision, save, dimension(:,:), allocatable :: BTurbulentKLGradz
double precision, save, dimension(:), allocatable :: TurbulentKLGradfx
double precision, save, dimension(:), allocatable :: TurbulentKLGradfy
double precision, save, dimension(:), allocatable :: TurbulentKLGradfz

double precision, save, dimension(:), allocatable :: TurbulentV2
double precision, save, dimension(:,:), allocatable :: BTurbulentV2
double precision, save, dimension(:), allocatable :: TurbulentV2Old
double precision, save, dimension(:,:), allocatable :: BTurbulentV2Old
double precision, save, dimension(:), allocatable :: TurbulentV2OldOld
double precision, save, dimension(:,:), allocatable :: BTurbulentV2OldOld
double precision, save, dimension(:), allocatable :: TurbulentV2Gradx
double precision, save, dimension(:,:), allocatable :: BTurbulentV2Gradx
double precision, save, dimension(:), allocatable :: TurbulentV2Grady
double precision, save, dimension(:,:), allocatable :: BTurbulentV2Grady
double precision, save, dimension(:), allocatable :: TurbulentV2Gradz
double precision, save, dimension(:,:), allocatable :: BTurbulentV2Gradz
double precision, save, dimension(:), allocatable :: TurbulentV2Gradfx
double precision, save, dimension(:), allocatable :: TurbulentV2Gradfy
double precision, save, dimension(:), allocatable :: TurbulentV2Gradfz

double precision, save, dimension(:), allocatable :: TurbulentZeta
double precision, save, dimension(:,:), allocatable :: BTurbulentZeta
double precision, save, dimension(:), allocatable :: TurbulentZetaOld
double precision, save, dimension(:,:), allocatable :: BTurbulentZetaOld
double precision, save, dimension(:), allocatable :: TurbulentZetaOldOld
double precision, save, dimension(:,:), allocatable :: BTurbulentZetaOldOld
double precision, save, dimension(:), allocatable :: TurbulentZetaGradx
double precision, save, dimension(:,:), allocatable :: BTurbulentZetaGradx
double precision, save, dimension(:), allocatable :: TurbulentZetaGrady
double precision, save, dimension(:,:), allocatable :: BTurbulentZetaGrady
double precision, save, dimension(:), allocatable :: TurbulentZetaGradz
double precision, save, dimension(:,:), allocatable :: BTurbulentZetaGradz
double precision, save, dimension(:), allocatable :: TurbulentZetaGradfx
double precision, save, dimension(:), allocatable :: TurbulentZetaGradfy
double precision, save, dimension(:), allocatable :: TurbulentZetaGradfz

double precision, save, dimension(:), allocatable :: TfRelaxation
double precision, save, dimension(:,:), allocatable :: BTfRelaxation
double precision, save, dimension(:), allocatable :: TfRelaxationOld
double precision, save, dimension(:,:), allocatable :: BTfRelaxationOld
double precision, save, dimension(:), allocatable :: TfRelaxationOldOld
double precision, save, dimension(:,:), allocatable :: BTfRelaxationOldOld
double precision, save, dimension(:), allocatable :: TfRelaxationGradx
double precision, save, dimension(:,:), allocatable :: BTfRelaxationGradx
double precision, save, dimension(:), allocatable :: TfRelaxationGrady
double precision, save, dimension(:,:), allocatable :: BTfRelaxationGrady
double precision, save, dimension(:), allocatable :: TfRelaxationGradz
double precision, save, dimension(:,:), allocatable :: BTfRelaxationGradz
double precision, save, dimension(:), allocatable :: TfRelaxationGradfx
double precision, save, dimension(:), allocatable :: TfRelaxationGradfy
double precision, save, dimension(:), allocatable :: TfRelaxationGradfz

double precision, save, dimension(:), allocatable :: uVelGrad2x
double precision, save, dimension(:,:), allocatable :: BuVelGrad2x
double precision, save, dimension(:), allocatable :: uVelGrad2y
double precision, save, dimension(:,:), allocatable :: BuVelGrad2y
double precision, save, dimension(:), allocatable :: uVelGrad2z
double precision, save, dimension(:,:), allocatable :: BuVelGrad2z
double precision, save, dimension(:), allocatable :: vVelGrad2x
double precision, save, dimension(:,:), allocatable :: BvVelGrad2x
double precision, save, dimension(:), allocatable :: vVelGrad2y
double precision, save, dimension(:,:), allocatable :: BvVelGrad2y
double precision, save, dimension(:), allocatable :: vVelGrad2z
double precision, save, dimension(:,:), allocatable :: BvVelGrad2z
double precision, save, dimension(:), allocatable :: wVelGrad2x
double precision, save, dimension(:,:), allocatable :: BwVelGrad2x
double precision, save, dimension(:), allocatable :: wVelGrad2y
double precision, save, dimension(:,:), allocatable :: BwVelGrad2y
double precision, save, dimension(:), allocatable :: wVelGrad2z
double precision, save, dimension(:,:), allocatable :: BwVelGrad2z


double precision, save, dimension(:), allocatable :: uvVelGradxy
double precision, save, dimension(:,:), allocatable :: BuvVelGradxy
double precision, save, dimension(:), allocatable :: uwVelGradxz
double precision, save, dimension(:,:), allocatable :: BuwVelGradxz

double precision, save, dimension(:), allocatable :: MaterialDerivative

double precision, save, dimension(:), allocatable :: NormalVelocity
double precision, save, dimension(:,:), allocatable :: BNormalVelocity
double precision, save, dimension(:), allocatable :: NormalVelocityGradx
double precision, save, dimension(:), allocatable :: NormalVelocityGrady
double precision, save, dimension(:), allocatable :: NormalVelocityGradz
double precision, save, dimension(:,:), allocatable :: BNormalVelocityGradx
double precision, save, dimension(:,:), allocatable :: BNormalVelocityGrady
double precision, save, dimension(:,:), allocatable :: BNormalVelocityGradz

double precision, save, dimension(:), allocatable :: dfidxTstar
double precision, save, dimension(:), allocatable :: dfidyTstar
double precision, save, dimension(:), allocatable :: dfidzTstar


end MODULE Variables1