MODULE Tecplot1
implicit none

double precision, dimension(:), allocatable :: SumInvd
double precision, dimension(:), allocatable :: SumPhiInvd
double precision, dimension(:), allocatable :: dummyN
double precision, dimension(:), allocatable :: dummyE
double precision, save, dimension(:,:), allocatable :: BdummyE
double precision, save, dimension(:), allocatable :: yplusPlot
double precision, save, dimension(:,:), allocatable :: ByplusPlot
double precision, save, dimension(:), allocatable :: uplusPlot
double precision, save, dimension(:,:), allocatable :: BuplusPlot


end MODULE Tecplot1
