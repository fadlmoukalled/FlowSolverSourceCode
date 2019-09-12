MODULE BoundaryConditionsTurbulence1
!use User0
implicit none

character*30, save, dimension(:,:), allocatable :: inletTypeTurbulence
character*30, save, dimension(:), allocatable :: inletTypeT
double precision, save, dimension(:), allocatable :: inletTurbulenceIntensity
double precision, save, dimension(:), allocatable :: inletTurbulenceLengthScale
double precision, save, dimension(:), allocatable :: inletTurbulenceViscosityRatio
end MODULE BoundaryConditionsTurbulence1
