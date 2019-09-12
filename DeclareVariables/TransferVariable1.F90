MODULE TransferVariable1

implicit none

double precision, save, dimension(:), allocatable :: ScalarT
double precision, save, dimension(:), allocatable :: ScalarOldT
double precision, save, dimension(:), allocatable :: ScalarOldOldT

double precision, save, dimension(:,:), allocatable :: BScalarT
double precision, save, dimension(:,:), allocatable :: BScalarOldT
double precision, save, dimension(:,:), allocatable :: BScalarOldOldT

double precision, save, dimension(:), allocatable :: ScalarGradxT
double precision, save, dimension(:), allocatable :: ScalarGradyT
double precision, save, dimension(:), allocatable :: ScalarGradzT
double precision, save, dimension(:,:), allocatable :: BScalarGradxT
double precision, save, dimension(:,:), allocatable :: BScalarGradyT
double precision, save, dimension(:,:), allocatable :: BScalarGradzT

double precision, save, dimension(:), allocatable :: ScalarGradfxT
double precision, save, dimension(:), allocatable :: ScalarGradfyT
double precision, save, dimension(:), allocatable :: ScalarGradfzT

double precision, save, dimension(:), allocatable :: ScScalarT
double precision, save, dimension(:), allocatable :: SbScalarT

double precision, save, dimension(:), allocatable :: DiffusionCoefficientT
double precision, save, dimension(:,:), allocatable :: BDiffusionCoefficientT


double precision, save, dimension(:), allocatable :: ScPointSourceScalarT
double precision, save, dimension (:), allocatable :: SbPointSourceScalarT

end MODULE TransferVariable1


