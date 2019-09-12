MODULE Scalar1
implicit none

double precision, save, dimension(:,:), allocatable :: Scalar
double precision, save, dimension(:,:), allocatable :: ScalarOld
double precision, save, dimension(:,:), allocatable :: ScalarOldOld

double precision, save, dimension(:,:,:), allocatable :: BScalar
double precision, save, dimension(:,:,:), allocatable :: BScalarOld
double precision, save, dimension(:,:,:), allocatable :: BScalarOldOld

double precision, save, dimension(:,:), allocatable :: ScalarGradx
double precision, save, dimension(:,:), allocatable :: ScalarGrady
double precision, save, dimension(:,:), allocatable :: ScalarGradz
double precision, save, dimension(:,:,:), allocatable :: BScalarGradx
double precision, save, dimension(:,:,:), allocatable :: BScalarGrady
double precision, save, dimension(:,:,:), allocatable :: BScalarGradz

double precision, save, dimension(:,:), allocatable :: ScalarGradfx
double precision, save, dimension(:,:), allocatable :: ScalarGradfy
double precision, save, dimension(:,:), allocatable :: ScalarGradfz

double precision, save, dimension(:,:), allocatable :: ScScalar
double precision, save, dimension(:,:), allocatable :: SbScalar


double precision, save, dimension(:,:), allocatable :: ScPointSourceScalar
double precision, save, dimension (:,:), allocatable :: SbPointSourceScalar


end Module Scalar1