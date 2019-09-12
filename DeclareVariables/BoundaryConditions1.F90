MODULE BoundaryConditions1

implicit none

character*20, save, dimension(:,:), allocatable :: BCType
character*20, save, dimension(:), allocatable :: BoundaryType

character*6, save, dimension(:,:), allocatable :: wallTypeMomentum
character*6, save, dimension(:,:), allocatable :: wallTypeContinuity
character*10, save, dimension(:,:), allocatable :: wallTypeEnergy
character*6, save, dimension(:), allocatable :: wallTypeM
character*6, save, dimension(:), allocatable :: wallTypeC
character*10, save, dimension(:), allocatable :: wallTypeE
character*10, save, dimension(:,:), allocatable :: wallTypeLambda
character*10, save, dimension(:), allocatable :: wallTypeL

character*30, save, dimension(:,:), allocatable :: inletTypeMomentum
character*30, save, dimension(:,:), allocatable :: inletTypeContinuity
character*30, save, dimension(:,:), allocatable :: inletTypeEnergy
character*30, save, dimension(:), allocatable :: inletTypeM
character*30, save, dimension(:), allocatable :: inletTypeC
character*30, save, dimension(:), allocatable :: inletTypeE

character*30, save, dimension(:,:), allocatable :: outletTypeMomentum
character*30, save, dimension(:,:), allocatable :: outletTypeContinuity
character*30, save, dimension(:,:), allocatable :: outletTypeEnergy
character*30, save, dimension(:), allocatable :: outletTypeM
character*30, save, dimension(:), allocatable :: outletTypeC
character*30, save, dimension(:), allocatable :: outletTypeE

integer, save :: NumberOfPointSources=0
integer, save, dimension (:), allocatable :: iElementPointSource
double precision, save, dimension(:), allocatable :: xLocationOfPointSource
double precision, save, dimension (:), allocatable :: yLocationOfPointSource
double precision, save, dimension (:), allocatable :: zLocationOfPointSource
double precision, save, dimension(:), allocatable :: ScPointSourceEnergy
double precision, save, dimension (:), allocatable :: SbPointSourceEnergy


end MODULE BoundaryConditions1
