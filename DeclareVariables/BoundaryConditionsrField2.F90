MODULE BoundaryConditionsrField2
!use User0
implicit none

integer, save, dimension(:), allocatable :: IWallDirichletrField
integer, save, dimension(:), allocatable :: IWallVonNeumannrField
integer, save, dimension(:), allocatable :: IWallRobinrField
integer, save, dimension(:), allocatable :: IinletSpecifiedValuerField
integer, save, dimension(:), allocatable :: IinletSupersonicrField
integer, save, dimension(:), allocatable :: IoutletFullyDevelopedrField
integer, save, dimension(:), allocatable :: IoutletSupersonicrField

integer, save, dimension(:,:), allocatable :: IWallDirichletrFieldOwner
integer, save, dimension(:,:), allocatable :: IWallDirichletrFieldNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IWallDirichletrFieldNBFaces

integer, save, dimension(:,:), allocatable :: IWallVonNeumannrFieldOwner
integer, save, dimension(:,:), allocatable :: IWallVonNeumannrFieldNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IWallVonNeumannrFieldNBFaces

integer, save, dimension(:,:), allocatable :: IWallRobinrFieldOwner
integer, save, dimension(:,:), allocatable :: IWallRobinrFieldNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IWallRobinrFieldNBFaces
double precision, save, dimension(:,:), allocatable :: rFieldinfinityRobin
double precision, save, dimension(:,:), allocatable :: rFieldConvectionCoefficientRobin

integer, save, dimension(:,:), allocatable :: IinletsupersonicrFieldOwner
integer, save, dimension(:,:), allocatable :: IinletsupersonicrFieldNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IinletsupersonicrFieldNBFaces

integer, save, dimension(:,:), allocatable :: IinletSpecifiedValuerFieldOwner
integer, save, dimension(:,:), allocatable :: IinletSpecifiedValuerFieldNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IinletSpecifiedValuerFieldNBFaces

integer, save, dimension(:,:), allocatable :: IoutletsupersonicrFieldOwner
integer, save, dimension(:,:), allocatable :: IoutletsupersonicrFieldNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IoutletsupersonicrFieldNBFaces

integer, save, dimension(:,:), allocatable :: IoutletFullyDevelopedrFieldOwner
integer, save, dimension(:,:), allocatable :: IoutletFullyDevelopedrFieldNumberOfBCSets
integer, save, dimension(:,:), allocatable :: IoutletFullyDevelopedrFieldNBFaces

end MODULE BoundaryConditionsrField2
