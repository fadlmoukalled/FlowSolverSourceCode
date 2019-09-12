MODULE PBCG1 !Preconditioned BiConjugate Gradient method
implicit none
!
double precision, save, dimension(:), allocatable :: p,pp,r,rr,z,zz
double precision, save, dimension(:), allocatable :: sa
integer, save, dimension(:), allocatable :: ija
double precision, save :: eps=1.d-14
!double precision, save :: tol=1.d-6
!
end MODULE PBCG1