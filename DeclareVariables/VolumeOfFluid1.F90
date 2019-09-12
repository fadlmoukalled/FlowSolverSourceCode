MODULE VolumeOfFluid1


  double precision, save, dimension(:,:), allocatable :: rField
  double precision, save, dimension(:,:), allocatable :: rFieldOld
  double precision, save, dimension(:,:), allocatable :: rFieldOldOld

  double precision, save, dimension(:,:,:), allocatable :: BrField
  double precision, save, dimension(:,:,:), allocatable :: BrFieldOld
  double precision, save, dimension(:,:,:), allocatable :: BrFieldOldOld

  double precision, save, dimension(:,:), allocatable :: rFieldGradx
  double precision, save, dimension(:,:), allocatable :: rFieldGrady
  double precision, save, dimension(:,:), allocatable :: rFieldGradz
  double precision, save, dimension(:,:,:), allocatable :: BrFieldGradx
  double precision, save, dimension(:,:,:), allocatable :: BrFieldGrady
  double precision, save, dimension(:,:,:), allocatable :: BrFieldGradz

  double precision, save, dimension(:,:), allocatable :: rFieldGradfx
  double precision, save, dimension(:,:), allocatable :: rFieldGradfy
  double precision, save, dimension(:,:), allocatable :: rFieldGradfz

  double precision, save, dimension(:,:), allocatable :: ScrField
  double precision, save, dimension(:,:), allocatable :: SbrField

  double precision, save, dimension(:), allocatable :: cosThetaF
  double precision, save, dimension(:), allocatable :: cosThetaE
  double precision, save, dimension(:,:), allocatable :: BcosThetaF
!
  double precision, save, dimension(:), allocatable :: Curvature
  double precision, save, dimension(:,:), allocatable :: BCurvature
!
  double precision, save, dimension(:), allocatable :: delrFieldMagnitude    
  double precision, save, dimension(:,:), allocatable :: BdelrFieldMagnitude    
  double precision, save, dimension(:), allocatable :: delrFieldMagnitudeGradx    
  double precision, save, dimension(:,:), allocatable :: BdelrFieldMagnitudeGradx
  double precision, save, dimension(:), allocatable :: delrFieldMagnitudeGrady    
  double precision, save, dimension(:,:), allocatable :: BdelrFieldMagnitudeGrady
  double precision, save, dimension(:), allocatable :: delrFieldMagnitudeGradz    
  double precision, save, dimension(:,:), allocatable :: BdelrFieldMagnitudeGradz
!
  double precision, save, dimension(:), allocatable :: rFieldGradxxT
  double precision, save, dimension(:,:), allocatable :: BrFieldGradxxT
  double precision, save, dimension(:), allocatable :: rFieldGradxyT
  double precision, save, dimension(:,:), allocatable :: BrFieldGradxyT
  double precision, save, dimension(:), allocatable :: rFieldGradxzT
  double precision, save, dimension(:,:), allocatable :: BrFieldGradxzT
  double precision, save, dimension(:), allocatable :: rFieldGradyyT
  double precision, save, dimension(:,:), allocatable :: BrFieldGradyyT
  double precision, save, dimension(:), allocatable :: rFieldGradyzT
  double precision, save, dimension(:,:), allocatable :: BrFieldGradyzT
  double precision, save, dimension(:), allocatable :: rFieldGradzzT
  double precision, save, dimension(:,:), allocatable ::BrFieldGradzzT
  
end MODULE VolumeOfFluid1