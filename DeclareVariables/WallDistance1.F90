MODULE WallDistance1

  implicit none
  double precision, save, dimension(:), allocatable :: WallDistance
  double precision, save, dimension(:,:), allocatable :: BWallDistance
  integer, save, dimension(:), allocatable :: iTau
  integer, save, dimension(:), allocatable :: jTau
  integer, save, dimension(:,:), allocatable :: BiTau
  integer, save, dimension(:,:), allocatable :: BjTau
  double precision, save :: urfWall=1.
  double precision, save :: rrfWD=1.e-7 !0.0001
  double precision, save :: WallResiduals=1.e-10
  integer, save :: IterWDMax=10000
  integer, save :: ASIterWallDistance=10
  integer, save :: MethodCalculateNormalDistance=2  !(1 or 2(solve differential equation))
  character*6,save :: ASSolverWD='sor'
  logical, save :: LMultigridWD=.true.
end MODULE WallDistance1
