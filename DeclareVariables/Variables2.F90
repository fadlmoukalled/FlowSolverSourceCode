MODULE Variables2
implicit none

double precision, save, dimension(:), allocatable :: FluxCf
double precision, save, dimension(:), allocatable :: FluxFf
double precision, save, dimension(:), allocatable :: FluxVf
double precision, save, dimension(:), allocatable :: FluxTf


double precision, save, dimension(:), allocatable :: ac
double precision, save, dimension(:), allocatable :: acStar
double precision, save, dimension(:), allocatable :: acold
double precision, save, dimension(:), allocatable :: acoldold
double precision, save, dimension(:,:), allocatable :: anb
double precision, save, dimension(:), allocatable :: bc
double precision, save, dimension(:), allocatable :: bcOriginal
double precision, save, dimension(:), allocatable :: dphi

double precision, save, dimension(:), allocatable :: dc
double precision, save, dimension(:), allocatable :: rc

      
end MODULE Variables2