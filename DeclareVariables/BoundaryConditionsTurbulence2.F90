MODULE BoundaryConditionsTurbulence2
!use User0
implicit none

integer, save  :: IwallTurbulence
integer, save  :: IinletTurbulence
integer, save  :: IoutletTurbulence

integer, save, dimension(:), allocatable :: IWallTurbulenceOwner
integer, save, dimension(:), allocatable :: IWallTurbulenceNumberOfBCSets
integer, save, dimension(:), allocatable :: IWallTurbulenceNBFaces

integer, save, dimension(:), allocatable :: IinletTurbulenceOwner
integer, save, dimension(:), allocatable :: IinletTurbulenceNumberOfBCSets
integer, save, dimension(:), allocatable :: IinletTurbulenceNBFaces

integer, save, dimension(:), allocatable :: IoutletTurbulenceOwner
integer, save, dimension(:), allocatable :: IoutletTurbulenceNumberOfBCSets
integer, save, dimension(:), allocatable :: IoutletTurbulenceNBFaces


end MODULE BoundaryConditionsTurbulence2
