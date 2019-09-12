MODULE RoughWall1
implicit none

logical, save :: LRough=.false.
double precision, save, dimension(:), allocatable :: WallRoughnessSize
double precision, save, dimension(:,:), allocatable :: hSplus
double precision, save, dimension(:,:), allocatable :: BwallRough

end MODULE RoughWall1
