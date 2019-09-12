MODULE BoundaryConditionsScalar2
!use User0
implicit none

integer, save, dimension(:), allocatable :: IWallDirichletScalar
integer, save, dimension(:), allocatable :: IWallVonNeumannScalar
integer, save, dimension(:), allocatable :: IWallRobinScalar
integer, save, dimension(:), allocatable :: IinletSpecifiedValueScalar
integer, save, dimension(:), allocatable :: IinletSupersonicScalar
integer, save, dimension(:), allocatable :: IoutletFullyDevelopedScalar
integer, save, dimension(:), allocatable :: IoutletSupersonicScalar

integer, save, dimension(:,:), allocatable :: IWallDirichletScalarOwner
integer, save, dimension(:,:), allocatable :: IWallDirichletScalarNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IWallDirichletScalarNBFaces

integer, save, dimension(:,:), allocatable :: IWallVonNeumannScalarOwner
integer, save, dimension(:,:), allocatable :: IWallVonNeumannScalarNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IWallVonNeumannScalarNBFaces

integer, save, dimension(:,:), allocatable :: IWallRobinScalarOwner
integer, save, dimension(:,:), allocatable :: IWallRobinScalarNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IWallRobinScalarNBFaces
double precision, save, dimension(:,:), allocatable :: PhiinfinityRobin
double precision, save, dimension(:,:), allocatable :: ConvectionCoefficientRobin

integer, save, dimension(:,:), allocatable :: IinletsupersonicScalarOwner
integer, save, dimension(:,:), allocatable :: IinletsupersonicScalarNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IinletsupersonicScalarNBFaces

integer, save, dimension(:,:), allocatable :: IinletSpecifiedValueScalarOwner
integer, save, dimension(:,:), allocatable :: IinletSpecifiedValueScalarNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IinletSpecifiedValueScalarNBFaces

integer, save, dimension(:,:), allocatable :: IoutletsupersonicScalarOwner
integer, save, dimension(:,:), allocatable :: IoutletsupersonicScalarNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IoutletsupersonicScalarNBFaces

integer, save, dimension(:,:), allocatable :: IoutletFullyDevelopedScalarOwner
integer, save, dimension(:,:), allocatable :: IoutletFullyDevelopedScalarNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IoutletFullyDevelopedScalarNBFaces

end MODULE BoundaryConditionsScalar2
