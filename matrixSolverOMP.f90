subroutine ComputeSourcesNew()
use ElementsVariables
use PhysicalVariables
use GMRESParameters
use VortexMethodParameters
use SparseMatricesVariables
use QuickSearchVariables
use GeneralSubroutine
use GeneralVariables
use CoreFunctions
use MathsTools
!$ use omp_lib
implicit none

type(SprsmatrixA) :: AC
integer :: i, l
integer :: k, k1, k2
integer :: j, jj, j1, j2
integer :: izone, jzone
integer :: CHUNK, nThreads, threadId, LastCHUNK
integer :: MaxNbOfInteraction

double precision :: Rcutoff2
double precision :: normr
double precision :: coefficient,phi
double precision :: dx, dy, dz, d, d2
double precision :: xT, yT, zT
double precision, allocatable, dimension(:) ::r, b, sol
integer*4, allocatable :: ThreadMaxNbOfInteraction(:)
integer*4, allocatable :: jjj(:, :, :)
integer, allocatable :: i1(:), i2(:)

write(*,*)
write(*,*)'You are in subroutine ComputeSourcesNew'

allocate(r(nElements))
allocate(b(nElements))
allocate(sol(nElements))

Coefficient = 1.0D0/(4.0D0*pi*sigma*sigma*sigma)

Rcutoff2 = RcutoffInversion*RcutoffInversion

!PARALLELIZE THIS BLOCK
r(1:nElements) = elements(1:nElements)%x(1)
b(1:nElements) = elements(1:nElements)%x(2)
sol(1:nElements) = elements(1:nElements)%x(3)

! Should be further considered

Zonesize = RcutoffInversion*sigma
call ReorderNew()

nThreads = OMP_get_max_Threads()
CHUNK = nElements/nThreads
LastCHUNK = nElements - (nThreads-1)*CHUNK

allocate(i1(nThreads))
allocate(i2(nThreads))

allocate(ThreadMaxNbOfInteraction(nThreads))
ThreadMaxNbOfInteraction(1:nThreads)= 0


!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nElements,elements,neighZoneRelInd,zones,r,b,sol,&
!$OMP                                  &Rcutoff2,CHUNK,nThreads,i1,i2,ThreadMaxNbOfInteraction,sigmaI)

ThreadID = omp_get_thread_num()
i1(threadID+1) = ThreadID*CHUNK+1
i2(ThreadID+1)=i1(ThreadID+1)+CHUNK-1
if((ThreadID+1).eq.nThreads) i2(ThreadID+1) = nElements


do i = i1(ThreadID+1),i2(ThreadID+1) !targets
	xT = r(i); yT = b(i); zT = sol(i)
	izone = elements(i)%indexZone
	k=0
	do l=1,27
		jzone = izone + neighZoneRelInd(l)
		if(zones(jzone)%nEls.gt.0) then
			j1 = zones(jzone)%startLoc
			j2 = j1+zones(jzone)%nEls-1
			do jj = j1,j2 !sources
				j = elements(jj)%index
				dx = (xT - r(j))*sigmaI; dy = (yT - b(j))*sigmaI; dz = (zT - sol(j))*sigmaI
				d2 = (dx*dx+dy*dy+dz*dz)
				if(d2.lt.Rcutoff2) then
					k=k+1
				end if
			end do ! jj
		end if
	end do !l
	if(ThreadMaxNbOfInteraction(threadID+1).lt.k) ThreadMaxNbOfInteraction(threadID+1)= k
end do !i
!$OMP END PARALLEL


MaxNbOfInteraction = 0
do i=1,nThreads
 if(MaxNbOfInteraction.lt.ThreadMaxNbOfInteraction(i)) MaxNbOfInteraction = ThreadMaxNbOfInteraction(i)
end do

deallocate(ThreadMaxNbOfInteraction)

allocate(jjj(nThreads,LastCHUNK,MaxNbOfInteraction))
jjj(1:nThreads,1:LastCHUNK,1:MaxNbOfInteraction)=0

allocate(AC%nentries(nElements))

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nElements,elements,neighZoneRelInd,zones,jjj,r,b,sol,&
!$OMP                                  &Rcutoff2,AC,CHUNK,nThreads,i1,i2,sigmaI)

ThreadID = omp_get_thread_num()

do i=i1(threadID+1),i2(threadID+1) !targets
	xT = r(i); yT = b(i); zT = sol(i)
	izone = elements(i)%indexZone
	k=0
	do l=1,27
		jzone = izone + neighZoneRelInd(l)
		if(zones(jzone)%nEls.gt.0) then
			j1=zones(jzone)%startLoc
			j2=j1+zones(jzone)%nEls-1
			do jj=j1,j2 !sources
				j = elements(jj)%index
				dx = (xT - r(j))*sigmaI; dy = (yT - b(j))*sigmaI; dz = (zT - sol(j))*sigmaI
				d2 = (dx*dx+dy*dy+dz*dz)
				if(d2.lt.Rcutoff2) then
					k=k+1;
					jjj(threadID+1,i-i1(threadID+1)+1,k) = j
				end if
			end do ! jj
		end if
	end do !l
	AC%nentries(i) = k
end do ! i
!$OMP END PARALLEL

do i=2,nElements
	AC%nentries(i) = AC%nentries(i) + AC%nentries(i-1)
end do


write(*,*)'total number of matrix coefficients =', AC%nentries(nElements)
write(*,*)'storage in Gigabytes is: ', 8.0D0*DBLE(AC%nentries(nElements))/DBLE(1073741824)
write(*,*)'filling matrix'

AC%length = AC%nentries(nElements)
allocate(AC%j(AC%length))

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(jjj,AC,i1,i2)

ThreadID = omp_get_thread_num()
do i=i1(ThreadID+1),i2(ThreadID+1) !targets
	k1=1
	if(i.gt.1) k1=AC%nentries(i-1)+1
	k2=AC%nentries(i)
	AC%j(k1:k2) = jjj(ThreadID+1,i-i1(ThreadID+1)+1,1:(k2-k1+1))
end do
!$OMP END PARALLEL

deallocate(jjj)
allocate(AC%value(AC%length))

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nElements,elements,AC,r,b,sol,CHUNK,nThreads,&
!$OMP                                  & i1,i2,coefOffset,sigmaI)
threadID = omp_get_thread_num()

	do i=i1(threadID+1),i2(threadID+1) !targets
		xT = r(i); yT = b(i); zT = sol(i)
		k1=1
		if(i.gt.1)k1=AC%nentries(i-1)+1
		k2=AC%nentries(i)
		do k=k1,k2 !sources
			j=AC%j(k)
			dx = (xT-r(j))*sigmaI
			dy = (yT-b(j))*sigmaI
			dz = (zT-sol(j))*sigmaI
			d2 = (dx*dx+dy*dy+dz*dz);
			d=DSQRT(d2)
			call CoreFunction(d,Phi)
			AC%value(k) = (phi-coefOffset)
		end do ! k
  end do ! i
!$OMP END PARALLEL

deallocate(i1)
deallocate(i2)

AC%value(1:AC%length) = coefficient*weight*AC%value(1:AC%length)
AC%m = nElements
AC%n = nElements
AC%rowIndexed = 1

do i = 1,nElements
  elements(i)%Gamma(1:3) = 0.0D0
end do

do i=1,3
	if(S0(i).GT.0.0D0) then
		b(1:nElements) = elements(1:nElements)%Vorticity(i)
		!call GMRES_RESTARTED_OMP_RETURNR(A,b,x,r,n,mMax,nIter,res,normr)
		call GMRES_RESTARTED_OMP_RETURNR(AC,b(1:nElements),sol(1:nElements),r(1:nElements),&
		&nElements,GMRES_ninner,GMRES_nouter,GMRES_tol,normr)
		S0r(i) = normr
		elements(1:nElements)%Gamma(i) = elements(1:nElements)%Gamma(i)+sol(1:nElements)
		write(*,*)'GMRES A',i, 1, S0(i), S0r(i), S0r(i)/S0(i)
		if(S0r(i)/S0(i)>0.01D0)then
			write(*,*)'inversion failed, stopping'
		end if
	else
		write(*,*)'GMRES A',i, 1, S0(i)
	end if
end do

deallocate(AC%value)
deallocate(AC%j)
deallocate(AC%nentries)
if(allocated(zones)) deallocate(zones)

deallocate(r)
deallocate(b)
deallocate(sol)

return
end subroutine ComputeSourcesNew

!---------------------

module SparseMatricesVariables
implicit none
save
double precision :: nmaxOvernElements = 0

type SprsmatrixA
	integer :: n             ! = nElements ( nb of rows in the ordinary matrix )
	integer :: m             ! = nElements ( nb of columns in the ordinary matrix )
	integer :: rowIndexed    ! = 1
	integer*8, allocatable :: nentries(:)  ! Ac%nentries(j)= nb of interaction of  elements 1 to j
	                                                       ! its dimension is (nElements)
	integer*8 :: length                      !  number of matrix coefficients =AC%nentries(nElements)
	integer*4, allocatable :: j(:)         ! Contains the list of interacting elements
                                                      ! its dimension is (Ac%lenght)
	double precision, allocatable :: value(:)  ! Contains the value of  interaction btw elements
                                                                 ! its dimension is (Ac%lenght)
end type SprsmatrixA

end module SparseMatricesVariables


!============================
module MathsTools
implicit none

public :: GMRES_RESTARTED_OMP_RETURNR
public :: AxOMP


contains

subroutine GMRES_RESTARTED_OMP_RETURNR(A,b,x,r,n,mMax,nIter,res,normr)
use sparseMatricesVariables
!$ use omp_lib
implicit none
! This subroutine solves the linear system Ax=b
! It gives as outpout normr=module(Ax-b)/module(b) where x is the optimal approximate solution
! CutoffFraction plays a crucial role in the convergence of the inversion scheme
! As cutoffFraction increases we will obtain more accurate solution

type(sprsmatrixA), intent(in) :: A !< coefficients matrix of size ( n x n), stored in sparse format
double precision, intent(inout) :: b(n) !< vector of size n
double precision, intent(inout) :: x(n) !< solution vector of size n
double precision, intent(out) :: r(n) !< vector of size n

integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
integer,  intent(in)  :: nIter !< maximum number of iterations
integer,  intent(in)  :: mMax !< maximum number of (inner) iterations to take.  0 < m <= m. Typically m is a small integer

double precision, intent(in) :: res !< parameter for convergence. See above.
double precision, intent(inout) :: normr !< is used as cutoff circulation when in

double precision :: w(n),x1(n),xOpt(n),rOpt(n)
double precision :: hi,hiP1,dotproduct
double precision :: d,gj,gjP1,resIter,normb,normrPrev,normmin,normb0
integer :: i,j,iter,totiter,totitermax,m,leave,mm
double precision :: cutoffCirculation,xMax,cutoffFraction
double precision, allocatable :: H(:,:),v(:,:),y(:),g(:),hj(:),hjP1(:),s(:),c(:)
double precision :: normminGLOBAL,xOptGlobal(n),rOptGlobal(n)
integer :: resetFlag=1

cutoffFraction=1.0D-3 !1.0D-3
cutoffCirculation=-1.0D0

!$OMP WORKSHARE
normb0 = dot_product(b,b)
!$OMP END WORKSHARE
normb0 = DSQRT(normb0)

b(1:n)=b(1:n)/normb0
normb = 1.0D0

mm=5
normminGLOBAL = 1.0D30

do while(mm.le.mMax)
  allocate(H(mm+1,mm))
  allocate(v(n,mm+1))
  allocate(y(mm))
  allocate(g(mm+1))
  allocate(hj(mm))
  allocate(hjP1(mm))
  allocate(s(mm))
  allocate(c(mm))

  normmin=1.0D30
  m=mm

  resIter = 1.0D6
  iter = 0; totiter = 0;

  if(resetFlag.eq.1.or.mm.eq.5)then
    x(1:n) = 0.0D0
    r=b

    !$OMP WORKSHARE
    normr = dot_product(r,r)
    !$OMP END WORKSHARE
    normr = DSQRT(normr)
  else
    x = xOptGlobal
    r = rOptGlobal
    normr=normminGLOBAL
  end if

  normrPrev = normr
  totiterMax = nIter*mm
  leave=0

  !write(*,*)'starting with',normr,normb

  do while(normmin.gt.res*normb.and.leave.eq.0)
    iter = iter+1
    v(1:n,1)=r/normr; H(1:m+1,1:m) = 0.0D0
    g(1:m+1) = 0; g(1) = normr
    do j=1,m
      totiter = totiter+1
      call AxOMP(A,v(1:n,j),w)
      do i=1,j
        !H(i,j) = dot(w,v(1:n,i),n)
        !$OMP WORKSHARE
        H(i,j) = dot_product(w,v(1:n,i))
        !$OMP END WORKSHARE
        w = w - H(i,j)*v(1:n,i)
      end do

      !H(j+1,j)=DSQRT(dot(w,w,n))
      !$OMP WORKSHARE
      dotproduct = dot_product(w,w)
      !$OMP END WORKSHARE
      H(j+1,j) = DSQRT(dotproduct)

      if(H(j+1,j).eq.0.0D0) then
        m=j; exit
      end if
      v(1:n,j+1) = w/H(j+1,j)

      ! apply j-1 Givens rotations to jth column of H
      if(j.gt.1) then
        do i=1,j-1
          hi = c(i)*H(i,j)+s(i)*H(i+1,j); hiP1 = -s(i)*H(i,j)+c(i)*H(i+1,j)
          H(i,j) = hi; H(i+1,j) = hiP1
        end do
      end if

      ! apply jth Givens rotation to H and g
      d = DSQRT(H(j,j)**2+H(j+1,j)**2);s(j) = H(j+1,j)/d;c(j) = H(j,j)/d
      !new jth and j+1 rows of H after multiplying it form the left with a Givens rotation
      hj(1:m) = c(j)*H(j,1:m)+s(j)*H(j+1,1:m); hjP1(1:m) = -s(j)*H(j,1:m)+c(j)*H(j+1,1:m)
      H(j,1:m) = hj(1:m); H(j+1,1:m) = hjP1(1:m)
      ! applying m rotations to g
      gj =  c(j)*g(j)+s(j)*g(j+1); gjP1 = -s(j)*g(j)+c(j)*g(j+1); g(j)=gj; g(j+1)=gjP1
      resIter = DABS(gjP1)
    end do

    ! getting y
    y(m) =  g(m)/H(m,m)

    do j=m-1,1,-1
      !dotproduct = dot(H(j,j+1:m),y(j+1:m),m-j)
      !$OMP WORKSHARE
      dotproduct = dot_product(H(j,j+1:m),y(j+1:m))
      !$OMP END WORKSHARE
      y(j)=(g(j)-dotproduct)/H(j,j)
    end do


    ! this is matrix multiplication
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    do i=1,n
      w(i) =dot_product(v(i,1:m),y(1:m))
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    x1 = x+w

    if(cutoffCirculation.le.0.0D0)then
      xMax=-1.0D30
      do i=1,n;
        if(DABS(x1(i)).gt.xMax)xMax=DABS(x1(i))
      end do;
      cutoffCirculation=cutoffFraction*xMax;
    end if

    do i=1,n;
      !    if(x1(i)*b(i).lt.0.0D0)x1(i)=0.0D0;
      if(x1(i)*b(i).lt.0.0D0.and.DABS(x1(i)).gt.cutoffCirculation)then
        if(b(i).lt.0)x1(i)=cutoffCirculation;
        if(b(i).gt.0)x1(i)=-cutoffCirculation;
      end if
    end do;


    call AxOMP(A,x1,r);
    r=b-r;

    !$OMP WORKSHARE
    normr = dot_product(r,r)
    !$OMP END WORKSHARE
    normr=DSQRT(normr)


    if((normr-normrPrev)/normrPrev.gt.-0.01D0.and.m.eq.1)leave=1
    if(m.eq.1.and.iter.gt.100)leave=1

    if(normr.lt.normmin) then
      normmin=normr
      xOpt=x1
      rOpt=r
    end if
    if(normr.lt.normrPrev) then
      x=x1
    else
      x = xOpt
      r = rOpt
    end if

    if((DABS((normr-normrPrev)/normr).lt.0.1D0.or.normr.gt.normrPrev).and.m.gt.1)then
      if(m.gt.10)then
        m=m-5
      else
        m=m-1
      end if
      iter=0
    end if

    normrPrev = normr

  end do

  write(*,*)'GMRES',m,mm,normr,normminGLOBAL,res,cutoffCirculation

  deallocate(H)
  deallocate(v)
  deallocate(y)
  deallocate(g)
  deallocate(hj)
  deallocate(hjP1)
  deallocate(s)
  deallocate(c)

  if(normmin.lt.normminGLOBAL)then
    normminGLOBAL=normmin
    xOptGlobal=xOpt
    rOptGlobal=rOpt
    mm=mm+5
  else
    mm=mMax+5
  end if

  if(normminGLOBAL.le.res)then
    mm=mMax+5
  end if
end do

x=xOptGlobal

call AxOMP(A,x,r)

r=b-r;
x(1:n)=x(1:n)*normb0
b(1:n)=b(1:n)*normb0
r(1:n)=r(1:n)*normb0
!$OMP WORKSHARE
normr = dot_product(r,r)
!$OMP END WORKSHARE
normr=DSQRT(normr)

return
end subroutine GMRES_RESTARTED_OMP_RETURNR


!-------------

subroutine BICGSTAB_OMP(A,b,x,n,nIter,res)
use sparseMatrices
implicit none
type(sprsmatrixA), intent(in) :: A !< coefficients matrix of size ( n x n), stored in sparse format
integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
integer,  intent(in)  :: nIter !< maximum number of iterations
double precision, intent(in) :: b(n) !< vector of size n
double precision, intent(inout) :: x(n) !< solution vector of size n
double precision, intent(in) :: res !< parameter for convergence. See above.
integer :: k
double precision :: r(n),r0(n),p(n),s(n),v(n),As(n),rho,rhoPrev,dp0,dp1,dp2
double precision :: normb,alpha,beta,omega,resnorm,regularity,regularityMin


call AxOMP(A,x,r0,0); r0 = b - r0; r = r0 !r0 = b - A x0
!$OMP WORKSHARE
rhoPrev = dot_product(r,r0)
!$OMP END WORKSHARE
resnorm = DSQRT(DABS(rhoPrev))
p = r0
!$OMP WORKSHARE
normb = dot_product(b,b)
!$OMP END WORKSHARE
normb = DSQRT(normb)

k = 1
regularityMin = 1.0D0
do while(k.lt.nIter.and.resnorm.gt.res*normb)
	k=k+1
	call AxOMP(A,p,v,0);
  !$OMP WORKSHARE
  dp0 = dot_product(v,r0)
  !$OMP END WORKSHARE

  alpha = rhoPrev/dp0
	s = r - alpha*v
	call AxOMP(A,s,As,0);
  !$OMP WORKSHARE
  dp1 = dot_product(As,s)
  !$OMP END WORKSHARE
  !$OMP WORKSHARE
  dp2 = dot_product(As,As)
  !$OMP END WORKSHARE

  omega = dp1/dp2
	x = x + alpha*p + omega*s
	r = s - omega*As
  !$OMP WORKSHARE
  rho = dot_product(r,r0)
  !$OMP END WORKSHARE

	beta = (rho/rhoPrev)*(alpha/omega); rhoPrev = rho
	p = r + beta*(p - omega*v)
	resnorm = DSQRT(DABS(rho))
  write(*,*)k,resnorm
end do
call AxOMP(A,x,r,0); r = b - r;
!$OMP WORKSHARE
dp1 = dot_product(r,r)
!$OMP END WORKSHARE
resnorm = DSQRT(dp1)
write(*,*)'BICGSTAB target residual=',k,nIter,resnorm,res*normb

return
end subroutine BICGSTAB_OMP

!----------------------------------------------------
! This subroutine takes as imput a Sparse Matrix A in Sparse Storage Format and a vector x
! and gives as an outpout the vector y=Ax
SUBROUTINE AxOMP(A,x,y)
	use sparseMatricesVariables
	!$ use omp_lib
	implicit none

	type(sprsmatrixA), intent(in) :: A         !< Sparse Matrix A in Sparse Storage Format
    double precision, intent(in)    :: x(*)    !< vector x
	double precision, intent(inout) :: y(A%m)  !< vector y

	integer(KIND=8) :: i,k
	double precision :: t
	integer(KIND=8)  :: nentries(A%m+1)

y(1:A%m) = 0.0D0

nentries(1) = 0
do i=2,A%m+1
	nentries(i) = A%nentries(i-1)
end do

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,t)
!$OMP DO SCHEDULE(STATIC)
do i=1,A%m
	t = 0.0D0
	do k=nentries(i)+1,nentries(i+1)
		t = t + A%value(k)*x(A%j(k))
	end do
	y(i) = t
end do
!$OMP END DO
!$OMP END PARALLEL

return
END SUBROUTINE AxOMP

end module MathsTools
