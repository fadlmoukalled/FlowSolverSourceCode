module GMRESParameters
    implicit none
    save
    integer,parameter :: GMRES_ninner = 10
    integer,parameter :: GMRES_nouter = 100
    double precision ::  GMRES_tol = 1.0D-9    ! Maximum Residual
    integer,parameter :: GMRES_save = 0
    integer,parameter :: GMRES_constrained = 1
    double precision,parameter :: GMRES_mu = 0.0D0
    double precision,parameter :: GMRES_threshold = 5.0E-008
end module GMRESParameters
!
module SparseMatrix
    implicit none
    save
    type SprsmatrixA
        integer :: n             ! = NumberOfElements ( nb of rows in the ordinary matrix )
        integer :: m             ! = NumberOfElements ( nb of columns in the ordinary matrix )
        integer :: rowIndexed    ! = 1
        integer*8, allocatable :: nentries(:)  ! AP%nentries(j)= nb of interaction of  elements 1 to j
    ! its dimension is (NumberOfElements)
        integer*8 :: length                      !  number of matrix coefficients =AP%nentries(NumberOfElements)
        integer*4, allocatable :: j(:)         ! Contains the list of interacting elements
    ! its dimension is (AP%lenght)
        double precision, allocatable :: value(:)  ! Contains the value of  interaction btw elements
    ! its dimension is (AP%lenght)
    end type SprsmatrixA    
end module SparseMatrix
!    
module SparseMatricesVariables
    use  SparseMatrix
    implicit none
    save
    double precision :: nmaxOverNumberOfElements = 0
    type (sprsmatrixA) :: AP 
    double precision, allocatable, dimension(:) :: b,x,r

end module SparseMatricesVariables
!    
module MathsTools
    implicit none
    public :: GMRES_RESTARTED_OMP_RETURNR
    public :: BICGSTAB_OMP
    public :: AxOMP
!
contains
!
  subroutine GMRES_RESTARTED_OMP_RETURNR(n,mMax,nIter,res,normr)
    use sparseMatricesVariables
    !$ use omp_lib
    implicit none
    ! This subroutine solves the linear system Ax=b
    ! It gives as outpout normr=module(Ax-b)/module(b) where x is the optimal approximate solution
    ! CutoffFraction plays a crucial role in the convergence of the inversion scheme
    ! As cutoffFraction increases we will obtain more accurate solution
    integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
!
    !type(sprsmatrixA), intent(in) :: A !< coefficients matrix of size ( n x n), stored in sparse format
    !double precision, intent(inout) :: b(n) !< vector of size n
    !double precision, intent(inout) :: x(n) !< solution vector of size n
    !double precision, intent(out) :: r(n) !< vector of size n
!
!    integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
    integer,  intent(in)  :: nIter !< maximum number of iterations
    integer,  intent(in)  :: mMax !< maximum number of (inner) iterations to take.  0 < m <= m. Typically m is a small integer
!
    double precision, intent(in) :: res !< parameter for convergence. See above.
    double precision, intent(inout) :: normr !< is used as cutoff circulation when in
!
    double precision, dimension (:), allocatable ::  w,x1,xOpt,rOpt
!    double precision :: w(n),x1(n),xOpt(n),rOpt(n)
    double precision :: hi,hiP1,dotproduct
    double precision :: d,gj,gjP1,resIter,normb,normrPrev,normmin,normb0
    integer :: i,j,iter,totiter,totitermax,m,leave,mm
    double precision :: cutoffCirculation,xMax,cutoffFraction
    double precision, dimension (:,:), allocatable ::  H,v
    double precision, dimension (:), allocatable ::  y,g,hj,hjP1,s,c
!    double precision, allocatable :: H(:,:),v(:,:),y(:),g(:),hj(:),hjP1(:),s(:),c(:)
    double precision, dimension (:), allocatable ::  xOptGlobal,rOptGlobal
    double precision :: normminGLOBAL !,xOptGlobal(n),rOptGlobal(n)
    integer :: resetFlag=1
!
    cutoffFraction=1.0D-3 !1.0D-3
    cutoffCirculation=-1.0D0
!    
    allocate(w(n))
    allocate(x1(n))
    allocate(xOpt(n))
    allocate(rOpt(n))
    allocate(xOptGlobal(n))
    allocate(rOptGlobal(n))
!    
    !$OMP WORKSHARE
    normb0 = dot_product(b,b)
    !$OMP END WORKSHARE
    normb0 = DSQRT(normb0)
!
    b(1:n)=b(1:n)/normb0
    normb = 1.0D0
!
    mm=5
    normminGLOBAL = 1.0D30
!
    do while(mm.le.mMax)
      allocate(H(mm+1,mm))
      allocate(v(n,mm+1))
      allocate(y(mm))
      allocate(g(mm+1))
      allocate(hj(mm))
      allocate(hjP1(mm))
      allocate(s(mm))
      allocate(c(mm))
!
      normmin=1.0D30
      m=mm
!
      resIter = 1.0D6
      iter = 0; totiter = 0;
!
      if(resetFlag.eq.1.or.mm.eq.5)then
        x(1:n) = 0.0D0
        r=b
!
        !$OMP WORKSHARE
        normr = dot_product(r,r)
        !$OMP END WORKSHARE
        normr = DSQRT(normr)
      else
        x = xOptGlobal
        r = rOptGlobal
        normr=normminGLOBAL
      end if
!
      normrPrev = normr
      totiterMax = nIter*mm
      leave=0
!
      !write(*,*)'starting with',normr,normb
!
      do while(normmin.gt.res*normb.and.leave.eq.0)
        iter = iter+1
        v(1:n,1)=r/normr; H(1:m+1,1:m) = 0.0D0
        g(1:m+1) = 0; g(1) = normr
        do j=1,m
          totiter = totiter+1
          call AxOMP(AP,v(1:n,j),w)
          do i=1,j
            !H(i,j) = dot(w,v(1:n,i),n)
            !$OMP WORKSHARE
            H(i,j) = dot_product(w,v(1:n,i))
            !$OMP END WORKSHARE
            w = w - H(i,j)*v(1:n,i)
          end do
!
          !H(j+1,j)=DSQRT(dot(w,w,n))
          !$OMP WORKSHARE
          dotproduct = dot_product(w,w)
          !$OMP END WORKSHARE
          H(j+1,j) = DSQRT(dotproduct)
!
          if(H(j+1,j).eq.0.0D0) then
            m=j; exit
          end if
          v(1:n,j+1) = w/H(j+1,j)
!
          ! apply j-1 Givens rotations to jth column of H
          if(j.gt.1) then
            do i=1,j-1
              hi = c(i)*H(i,j)+s(i)*H(i+1,j); hiP1 = -s(i)*H(i,j)+c(i)*H(i+1,j)
              H(i,j) = hi; H(i+1,j) = hiP1
            end do
          end if
!
          ! apply jth Givens rotation to H and g
          d = DSQRT(H(j,j)**2+H(j+1,j)**2);s(j) = H(j+1,j)/d;c(j) = H(j,j)/d
          !new jth and j+1 rows of H after multiplying it form the left with a Givens rotation
          hj(1:m) = c(j)*H(j,1:m)+s(j)*H(j+1,1:m); hjP1(1:m) = -s(j)*H(j,1:m)+c(j)*H(j+1,1:m)
          H(j,1:m) = hj(1:m); H(j+1,1:m) = hjP1(1:m)
          ! applying m rotations to g
          gj =  c(j)*g(j)+s(j)*g(j+1); gjP1 = -s(j)*g(j)+c(j)*g(j+1); g(j)=gj; g(j+1)=gjP1
          resIter = DABS(gjP1)
        end do
!
        ! getting y
        y(m) =  g(m)/H(m,m)
!
        do j=m-1,1,-1
          !dotproduct = dot(H(j,j+1:m),y(j+1:m),m-j)
          !$OMP WORKSHARE
          dotproduct = dot_product(H(j,j+1:m),y(j+1:m))
          !$OMP END WORKSHARE
          y(j)=(g(j)-dotproduct)/H(j,j)
        end do
!
        ! this is matrix multiplication
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
        !$OMP DO SCHEDULE(STATIC)
        do i=1,n
          w(i) =dot_product(v(i,1:m),y(1:m))
        end do
        !$OMP END DO
        !$OMP END PARALLEL
!
        x1 = x+w
!
        !if(cutoffCirculation.le.0.0D0)then
         ! xMax=-1.0D30
          !do i=1,n;
          !  if(DABS(x1(i)).gt.xMax)xMax=DABS(x1(i))
          !end do;
          !cutoffCirculation=cutoffFraction*xMax;
        !end if
!
        !do i=1,n;
        !  !    if(x1(i)*b(i).lt.0.0D0)x1(i)=0.0D0;
         ! if(x1(i)*b(i).lt.0.0D0.and.DABS(x1(i)).gt.cutoffCirculation)then
          !  if(b(i).lt.0)x1(i)=cutoffCirculation;
          !  if(b(i).gt.0)x1(i)=-cutoffCirculation;
          !end if
        !end do;
!
        call AxOMP(AP,x1,r);
        r=b-r;
!
        !$OMP WORKSHARE
        normr = dot_product(r,r)
        !$OMP END WORKSHARE
        normr=DSQRT(normr)
!
        if((normr-normrPrev)/normrPrev.gt.-0.01D0.and.m.eq.1)leave=1
        if(m.eq.1.and.iter.gt.100)leave=1
!
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
!
        if((DABS((normr-normrPrev)/normr).lt.0.1D0.or.normr.gt.normrPrev).and.m.gt.1)then
          if(m.gt.10)then
            m=m-5
          else
            m=m-1
          end if
          iter=0
        end if
!
        normrPrev = normr
!
      end do
!
      write(*,*)'GMRES',m,mm,normr,normminGLOBAL,res,cutoffCirculation
!
      deallocate(H)
      deallocate(v)
      deallocate(y)
      deallocate(g)
      deallocate(hj)
      deallocate(hjP1)
      deallocate(s)
      deallocate(c)
!
      if(normmin.lt.normminGLOBAL)then
        normminGLOBAL=normmin
        xOptGlobal=xOpt
        rOptGlobal=rOpt
        mm=mm+5
      else
        mm=mMax+5
      end if
!
      if(normminGLOBAL.le.res)then
        mm=mMax+5
      end if
    end do
!
    x=xOptGlobal
!
    call AxOMP(AP,x,r)
!
    r=b-r;
    x(1:n)=x(1:n)*normb0
    b(1:n)=b(1:n)*normb0
    r(1:n)=r(1:n)*normb0
    !$OMP WORKSHARE
    normr = dot_product(r,r)
    !$OMP END WORKSHARE
    normr=DSQRT(normr)
!
    deallocate(w)
    deallocate(x1)
    deallocate(xOpt)
    deallocate(rOpt)
    deallocate(xOptGlobal)
    deallocate(rOptGlobal)
!
    return
  end subroutine GMRES_RESTARTED_OMP_RETURNR
!
  subroutine BICGSTAB_OMP(n,nIter,res)
    use sparseMatricesVariables
    !$ use omp_lib
    implicit none
!    type(sprsmatrixA), intent(in) :: A !< coefficients matrix of size ( n x n), stored in sparse format
    integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
    integer,  intent(in)  :: nIter !< maximum number of iterations
!    double precision, intent(in) :: b(n) !< vector of size n
!    double precision, intent(inout) :: x(n) !< solution vector of size n
    double precision, intent(in) :: res !< parameter for convergence. See above.
    integer :: k
    double precision, dimension (:), allocatable ::  r0,p,s,v,As
    double precision :: rho,rhoPrev,dp0,dp1,dp2
    double precision :: normb,alpha,beta,omega,resnorm,regularity,regularityMin
!  
!    allocate(r(n))
    allocate(r0(n))
    allocate(p(n))
    allocate(s(n))
    allocate(v(n))
    allocate(As(n))
!    
    call AxOMP(AP,x,r0); r0 = b - r0; r = r0 !r0 = b - A x0
    !$OMP WORKSHARE
    rhoPrev = dot_product(r,r0)
    !$OMP END WORKSHARE
    resnorm = DSQRT(DABS(rhoPrev))
    p = r0
    !$OMP WORKSHARE
    normb = dot_product(b,b)
    !$OMP END WORKSHARE
    normb = DSQRT(normb)
!  
    k = 1
    regularityMin = 1.0D0
    do while(k.lt.nIter.and.resnorm.gt.res*normb)
      k=k+1
      call AxOMP(AP,p,v);
      !$OMP WORKSHARE
      dp0 = dot_product(v,r0)
      !$OMP END WORKSHARE
!  
      alpha = rhoPrev/dp0
      s = r - alpha*v
      call AxOMP(AP,s,As);
      !$OMP WORKSHARE
      dp1 = dot_product(As,s)
      !$OMP END WORKSHARE
      !$OMP WORKSHARE
      dp2 = dot_product(As,As)
      !$OMP END WORKSHARE
!  
      omega = dp1/dp2
      x = x + alpha*p + omega*s
      r = s - omega*As
      !$OMP WORKSHARE
      rho = dot_product(r,r0)
      !$OMP END WORKSHARE
!  
      beta = (rho/rhoPrev)*(alpha/omega); rhoPrev = rho
      p = r + beta*(p - omega*v)
      resnorm = DSQRT(DABS(rho))
      write(*,*) k,resnorm
    end do
!
    call AxOMP(AP,x,r); r = b - r;
    !$OMP WORKSHARE
    dp1 = dot_product(r,r)
    !$OMP END WORKSHARE
    resnorm = DSQRT(dp1)
    write(*,*)'BICGSTAB target residual=',k,nIter,resnorm,res*normb
!
!    deallocate(r)
    deallocate(r0)
    deallocate(p)
    deallocate(s)
    deallocate(v)
    deallocate(As)
!    
    return
  end subroutine BICGSTAB_OMP

  ! This subroutine takes as imput a Sparse Matrix A in Sparse Storage Format and a vector x
  ! and gives as an outpout the vector y=Ax
  SUBROUTINE AxOMP(A,x,y)
    use sparseMatrix
    !$ use omp_lib
    implicit none
!
    type(sprsmatrixA), intent(in) :: A         !< Sparse Matrix A in Sparse Storage Format
    double precision, intent(in)    :: x(*)    !< vector x
    double precision, intent(inout) :: y(A%m)  !< vector y
!
    integer(KIND=8) :: i,k
    double precision :: t
    integer(KIND=8)  :: nentries(A%m+1)
!
    y(1:A%m) = 0.0D0
!
    nentries(1) = 0
    do i=2,A%m+1
      nentries(i) = A%nentries(i-1)
    end do
!
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
!
    return
  END SUBROUTINE AxOMP

end module MathsTools
!
subroutine SolveAlgebraicSystemInParallel(solver)
!  
  use Geometry1, only : NumberOfElements
  use variables2, only : ac, bc, anb, dphi
  use Geometry3, only : NumberofElementNeighbors, ElementNeighbor
  use MathsTools
  use SparseMatricesVariables
  use GMRESParameters
  use Geometry1, only: MaximumNumberofElementFaces
  !$ use omp_lib
  implicit none
!
!  type(SprsmatrixA) :: AP
  integer :: i, l,kindx
  integer :: k, k1, k2
  integer :: j, jj, j1, j2
  integer :: izone, jzone, indexInGlobalList
  integer :: CHUNK, nThreads, threadId, LastCHUNK
  integer :: MaxNbOfInteraction
!  
  character*6 :: solver
!
  double precision :: Rcutoff2
  double precision :: normr
  double precision :: coefficient,phi
  double precision :: dx, dy, dz, d, d2
  double precision :: xT, yT, zT
!  double precision, allocatable, dimension(:) ::r, sol
!  integer*4, allocatable :: ThreadMaxNbOfInteraction(:)
  integer*4, allocatable :: jjj(:, :, :)
  integer*4, save, allocatable :: jjj1(:, :, :)
  integer, save, allocatable :: i1(:), i2(:)
  integer, save :: LocalIndex
  data LocalIndex/0/
!  
  if(LocalIndex.eq.0) then
  nThreads = OMP_get_max_threads();
  call OMP_SET_NUM_THREADS(nThreads)
!
  ! HERE WE ARE SOLVING AP sol = b, rthe solver yields sol and the residual r
  ! the number of unknonws is NumberOfElements
!
  CHUNK = NumberOfElements/nThreads
  LastCHUNK = NumberOfElements - (nThreads-1)*CHUNK !note that LastCHUNKshould always be >= CHUNK
!
  allocate(i1(nThreads))
  allocate(i2(nThreads))
!
!  allocate(ThreadMaxNbOfInteraction(nThreads))
!  ThreadMaxNbOfInteraction(1:nThreads)= 0
!
  !Compute the numnber of interactions per thread
  !For a grid-based solution this is the number of targets assigned to the thread (i2(ThreadID+1)-i1(ThreadID+1)+1)
  !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NumberOfElements,CHUNK,nThreads,i1,i2)
  ThreadID = omp_get_thread_num()
  i1(threadID+1) = ThreadID*CHUNK+1
  i2(ThreadID+1)=i1(ThreadID+1)+CHUNK-1
  if((ThreadID+1).eq.nThreads) i2(ThreadID+1) = NumberOfElements
  !do i = i1(ThreadID+1),i2(ThreadID+1) !targets
  !  !specify or compute ThreadMaxNbOfInteraction(threadID+1)
  !  ThreadMaxNbOfInteraction(threadID+1)=i2(ThreadID+1)-i1(ThreadID+1)+1;
  !end do !i
  !$OMP END PARALLEL
  !
  !!compute the maximum number of interactions among ther threads
  !MaxNbOfInteraction = 0
  !do i=1,nThreads
  !  if(MaxNbOfInteraction.lt.ThreadMaxNbOfInteraction(i)) MaxNbOfInteraction = ThreadMaxNbOfInteraction(i)
  !end do
  !deallocate(ThreadMaxNbOfInteraction)
!
  MaxNbOfInteraction=MaximumNumberofElementFaces+1
!
  allocate(jjj(nThreads,LastCHUNK,MaxNbOfInteraction))
  allocate(jjj1(nThreads,LastCHUNK,MaxNbOfInteraction))
!
  jjj(1:nThreads,1:LastCHUNK,1:MaxNbOfInteraction)=0
  jjj1(1:nThreads,1:LastCHUNK,1:MaxNbOfInteraction)=0
!
  allocate(AP%nentries(NumberOfElements))
!
  ! NOW we fill matrix jjj():
  ! the first index is the thread index:  threadID+1
  ! the second index is index of the target within the thread target list relative to the first one in the list: i-i1(threadID+1)+1
  ! the third one is a counter for the thread being processed which counts the nujmber of nonzero entries in AP
  ! the value of jjj() is the index of the source in the global list that interacts with target i
  !
  !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NumberOfElements,jjj,jjj1,AP,CHUNK,nThreads,i1,i2,NumberofElementNeighbors,ElementNeighbor)
  ThreadID = omp_get_thread_num()
!
  do i=i1(threadID+1),i2(threadID+1) !targets
    k=0;
    k=k+1;jjj(threadID+1,i-i1(threadID+1)+1,k) = i
    do j = 1,NumberofElementNeighbors(i)
      indexInGlobalList = ElementNeighbor(i,j)
      if(indexInGlobalList .ne.0) then
        k=k+1;
        ! assign or compute jjj(i1,i2,i3)
        jjj(threadID+1,i-i1(threadID+1)+1,k) = indexInGlobalList
        jjj1(threadID+1,i-i1(threadID+1)+1,k) = j
      end if
    end do
    AP%nentries(i) = k
  end do ! i
  !$OMP END PARALLEL
!
  do i=2,NumberOfElements
    AP%nentries(i) = AP%nentries(i) + AP%nentries(i-1)
  end do
!
  write(*,*)'total number of matrix coefficients =', AP%nentries(NumberOfElements)
  write(*,*)'storage in Gigabytes is: ', 8.0D0*DBLE(AP%nentries(NumberOfElements))/DBLE(1073741824)
  write(*,*)'filling matrix'
!
  AP%length = AP%nentries(NumberOfElements)
  allocate(AP%j(AP%length))
!
  !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(jjj,AP,i1,i2)
  ThreadID = omp_get_thread_num()
  do i=i1(ThreadID+1),i2(ThreadID+1) !targets
    k1=1
    if(i.gt.1) k1=AP%nentries(i-1)+1
    k2=AP%nentries(i)
!    AP%j(k1:k2) = jjj(ThreadID+1,i-i1(ThreadID+1)+1,1:(k2-k1+1))
    do kindx=k1,k2
        AP%j(kindx) = jjj(ThreadID+1,i-i1(ThreadID+1)+1,kindx-k1+1)
    end do
  end do
  !$OMP END PARALLEL
!
  deallocate(jjj)
  allocate(AP%value(AP%length))
!
  allocate(x(NumberOfElements))
  allocate(r(NumberOfElements))
  allocate(b(NumberOfElements))
!  
  LocalIndex=1
!  
  endif  !LocalIndex
!  
  !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NumberOfElements,AP,CHUNK,nThreads,&
  !$OMP                                   i1,i2, ac, anb, jjj1)
  threadID = omp_get_thread_num()
  do i=i1(threadID+1),i2(threadID+1) !targets
    k1=1
    if(i.gt.1)k1=AP%nentries(i-1)+1
    k2=AP%nentries(i)
    j=AP%j(k1) ! j is the index of source k in the global list
    AP%value(k1) = ac(j);
!
    do k=k1+1,k2 !sources
      !j=AP%j(k) ! j is the index of source k in the global list
      j = jjj1(threadID+1,i-i1(threadID+1)+1,k-k1+1)
      AP%value(k) = anb(i,j)
    end do ! k
  end do ! i
  !$OMP END PARALLEL
!
!  deallocate(i1)
!  deallocate(i2)
!  deallocate(jjj1)
!  
  AP%m = NumberOfElements
  AP%n = NumberOfElements
  AP%rowIndexed = 1
  write(*,*) NumberOfElements,GMRES_ninner,GMRES_nouter,GMRES_tol,normr
  !call GMRES_RESTARTED_OMP_RETURNR(A,b,x,r,n,mMax,nIter,res,normr)
  !call GMRES_RESTARTED_OMP_RETURNR(AP,bc(1:NumberOfElements),sol(1:NumberOfElements),r(1:NumberOfElements),&
  !               NumberOfElements,GMRES_ninner,GMRES_nouter,GMRES_tol,normr)
  !allocate(x(NumberOfElements))
  !allocate(r(NumberOfElements))
  !allocate(b(NumberOfElements))
!  
  b=bc;
!
  if(solver.eq.'gmres') then
    call GMRES_RESTARTED_OMP_RETURNR(NumberOfElements,GMRES_ninner,GMRES_nouter,GMRES_tol,normr)
    write(*,*) normr
  elseif(solver.eq.'bicgs') then
    call BICGSTAB_OMP(NumberOfElements,100,1.d-5)
  endif
!  
dphi=x
!
!deallocate(AP%value)
!deallocate(AP%j)
!deallocate(AP%nentries)
!
!deallocate(r)
!deallocate(x)
!deallocate(b)
  return
end subroutine SolveAlgebraicSystemInParallel