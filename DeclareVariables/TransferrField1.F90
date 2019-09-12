MODULE TransferrField1

implicit none

double precision, save, dimension(:), allocatable :: rFieldT
double precision, save, dimension(:), allocatable :: rFieldOldT
double precision, save, dimension(:), allocatable :: rFieldOldOldT

double precision, save, dimension(:,:), allocatable :: BrFieldT
double precision, save, dimension(:,:), allocatable :: BrFieldOldT
double precision, save, dimension(:,:), allocatable :: BrFieldOldOldT

double precision, save, dimension(:), allocatable :: rFieldGradxT
double precision, save, dimension(:), allocatable :: rFieldGradyT
double precision, save, dimension(:), allocatable :: rFieldGradzT
double precision, save, dimension(:,:), allocatable :: BrFieldGradxT
double precision, save, dimension(:,:), allocatable :: BrFieldGradyT
double precision, save, dimension(:,:), allocatable :: BrFieldGradzT

double precision, save, dimension(:), allocatable :: rFieldGradfxT
double precision, save, dimension(:), allocatable :: rFieldGradfyT
double precision, save, dimension(:), allocatable :: rFieldGradfzT

end MODULE TransferrField1


