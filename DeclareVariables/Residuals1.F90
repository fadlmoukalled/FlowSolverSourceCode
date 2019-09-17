MODULE Residuals1
implicit none
integer, save :: nVariablesPlot
double precision, save, dimension(:), allocatable :: ResorAbs
double precision, save, dimension(:), allocatable :: ResorMax
double precision, save, dimension(:), allocatable :: ResorRMS
double precision, save, dimension(:), allocatable :: ResorScaled
double precision, save, dimension(:), allocatable :: ResidualsPlot
double precision, save :: OverallImbalance
integer, save :: iMaxImbalance
end MODULE Residuals1
