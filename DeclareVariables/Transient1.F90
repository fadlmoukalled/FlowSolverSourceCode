MODULE Transient1
implicit none
 integer, save :: ndt
 double precision, save :: dtOld,dtOldOld
 double precision, save, dimension(:), allocatable :: dCpdt,dpdt
end MODULE Transient1
