c
C#############################################################################################
c
      SUBROUTINE allocateForPBCG
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use PBCG1
c
C#############################################################################################
c
      implicit none
c********************************************************************************************
      integer, save :: ArraySize
      integer,save :: index
      data index/0/
c********************************************************************************************
c
      if(index.eq.1) return
c      
      allocate(p(NumberOfElements))
      allocate(pp(NumberOfElements))
      allocate(r(NumberOfElements))
      allocate(rr(NumberOfElements))
      allocate(z(NumberOfElements))
      allocate(zz(NumberOfElements))
c
      if(index.eq.0) call GetSizeOfPBCGarrays(ArraySize)
      index=1
c
      allocate(sa(ArraySize))
      allocate(ija(ArraySize))
c
c---- Initialize variables
c
      p=0.
      pp=0.
      r=0.
      rr=0.
      z=0.
      zz=0.
      sa=0.
      ija=0
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE GetSizeOfPBCGarrays(ArraySize)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry3, only: ElementNeighbor,NumberofElementNeighbors
      use Variables2, only: anb
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,k1,ArraySize
      double precision thresh
      data thresh/0.d0/
c********************************************************************************************
c
      k=NumberOfElements+1
c
      do i=1,NumberOfElements
        do j=1,NumberofElementNeighbors(i)
c
          k1 = ElementNeighbor(i,j)
c
          if(k1.ne.0) then   
c
            if(dabs(anb(i,j)).ge.thresh)then
c
              if(i.ne.k1)then
c
                k=k+1
c
              endif
c
            endif
c
          endif
c
         enddo
       enddo
       ArraySize=k
      return
      end
c
C#############################################################################################
c
      SUBROUTINE SolveEquationsUsingPBCG(iter)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Variables2, only: bc,dphi
      use PBCG1
      use constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer iter,j
      double precision ak,akden,bk,bkden,bknum
c
C********************************************************************************************
C     USES atimes,asolve,snrm
C***********************************************************************************************************
!     Solves A · x = b for x(1:n), given b(1:n), by the iterative biconjugate gradient method.             !
!     On input x(1:n) should be set to an initial guess of the solution (or all zeros);                    !
!     itol is 1,2,3, or 4, specifying which convergence test is applied (see text);                        !
!     itmax is the maximum number of allowed iterations; and tol is the desired convergence tolerance.     !
!     On output, x(1:n) is reset to the improved solution, iter is the number of iterations actually taken,!
!     and err is the estimated error. The matrix A is referenced only through the user-supplied routines   !
!     atimes, which computes the product of either A or its transpose on a vector; and                     !
!     asolve, which solves A  · x = b or A T · x = b for some preconditioner matrix A                      !
!     (possibly the trivial diagonal part of A).                                                           !
C***********************************************************************************************************
c
c--- Interfaces
c
      interface
c---------------------------------------------------------------------------------------------
        SUBROUTINE atimes(itrnsp,x1,b1)
c---------------------------------------------------------------------------------------------
          integer itrnsp
          double precision, dimension(:) :: x1
          double precision, dimension(:) :: b1
c---------------------------------------------------------------------------------------------
        end SUBROUTINE atimes
c---------------------------------------------------------------------------------------------
        SUBROUTINE asolve(itrnsp,x2,b2)
c---------------------------------------------------------------------------------------------
          integer itrnsp
          double precision, dimension (:) :: x2
          double precision, dimension (:) :: b2
c---------------------------------------------------------------------------------------------
        end SUBROUTINE asolve
c---------------------------------------------------------------------------------------------
      end interface
c---------------------------------------------------------------------------------------------
c
      if(iter.eq.0) then
c
        call atimes(0,r,dphi)
c
        do j=1,NumberOfElements
c
          r(j)=bc(j)-r(j)
          rr(j)=r(j)
c
        enddo
c
        call atimes(0,rr,r)  ! By uncommentting this line we get the 
                             ! “minimum residual” variant of the algorithm.
c
      endif
c
      call asolve(0,z,r)
c
      call asolve(1,zz,rr)
      bknum=0.d0
c
      do j=1,NumberOfElements
        bknum=bknum+z(j)*rr(j)
      enddo
c
      if(iter.eq.0) then
c
        do j=1,NumberOfElements
c
          p(j)=z(j)
          pp(j)=zz(j)
c
        enddo
c
      else
c
        bk=bknum/(bkden+1.d-10)
c
        do j=1,NumberOfElements
c
          p(j)=bk*p(j)+z(j)
          pp(j)=bk*pp(j)+zz(j)
c
        enddo
c
      endif
c
      bkden=bknum
      call atimes(0,z,p)
      akden=0.d0
c
      do j=1,NumberOfElements
c
        akden=akden+z(j)*pp(j)
c
      enddo
c
      ak=bknum/(akden+1.d-10)
      call atimes(1,zz,pp)
c
      do j=1,NumberOfElements
c
        dphi(j)=dphi(j)+ak*p(j)
        r(j)=r(j)-ak*z(j)
        rr(j)=rr(j)-ak*zz(j)
c
      enddo
c
      return
      end
c
C  Adopted and modified from Copr. 1986-92 Numerical Recipes Software.
c
C#############################################################################################
c
      SUBROUTINE atimes(itrnsp,x1,b1)
c
C#############################################################################################
c
      implicit none
c********************************************************************************************
      integer itrnsp
      double precision, dimension(:) :: x1
      double precision, dimension(:) :: b1
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE dsprsax(x3,b3)
c---------------------------------------------------------------------------------------------
          double precision, dimension(:) :: x3
          double precision, dimension(:) :: b3
c---------------------------------------------------------------------------------------------
        end SUBROUTINE dsprsax
c---------------------------------------------------------------------------------------------
        SUBROUTINE dsprstx(x4,b4)
c---------------------------------------------------------------------------------------------
          double precision, dimension(:) :: x4
          double precision, dimension(:) :: b4
c---------------------------------------------------------------------------------------------
        end SUBROUTINE dsprstx
c---------------------------------------------------------------------------------------------
      end interface
c
C********************************************************************************************
C     USES dsprsax,dsprstx
C********************************************************************************************
c
      if (itrnsp.eq.0) then
c
        call dsprsax(x1,b1)
c
      else
c
        call dsprstx(x1,b1)
c
      endif
c
      return
      end
c
C  Adopted and modified from Copr. 1986-92 Numerical Recipes Software.
c
C#############################################################################################
c
      SUBROUTINE asolve(itrnsp,x2,b2)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Variables2, only: dc
      use PBCG1, only: sa
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,itrnsp
      double precision, dimension (:) :: x2
      double precision, dimension (:) :: b2
c********************************************************************************************
c
c   use dc(i) for using the diagonal obtained from the ILU factorization as the preconditioning
c   matrix, and use sa(i) for using the diagonal of A as simply the preconditioning matrix. A
c   parameter can be added to make it optional for the user.
c
      do i=1,NumberOfElements
c
        x2(i)=b2(i)*dc(i) !/sa(i)
c
      enddo
c
      return
      end
c
C  Adopted and modified from Copr. 1986-92 Numerical Recipes Software.
c
C#############################################################################################
c
      SUBROUTINE dsprsax(x3,b3)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use PBCG1, only: sa,ija
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,k
      double precision, dimension(:) :: x3
      double precision, dimension(:) :: b3
c********************************************************************************************
c
      if (ija(1).ne.NumberOfElements+2) then 
c
        print*, 'mismatched vector and matrix in sprstx'
        pause
c
      endif
c
      do i=1,NumberOfElements
        x3(i)=sa(i)*b3(i)
        do k=ija(i),ija(i+1)-1
          x3(i)=x3(i)+sa(k)*b3(ija(k))
        enddo
      enddo
c
      return
      end
c
C  Adopted and modified from Copr. 1986-92 Numerical Recipes Software.
c
C#############################################################################################
c
      SUBROUTINE dsprstx(x4,b4)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use PBCG1, only: sa,ija
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k
      double precision, dimension(:) :: x4
      double precision, dimension(:) :: b4
c********************************************************************************************
c
      if (ija(1).ne.NumberOfElements+2) then 
c
        print*, 'mismatched vector and matrix in sprstx'
        pause
c
      endif
c
      do i=1,NumberOfElements
        x4(i)=sa(i)*b4(i)
      enddo
c
      do i=1,NumberOfElements
        do k=ija(i),ija(i+1)-1
          j=ija(k)
          x4(j)=x4(j)+sa(k)*b4(i)
        enddo
      enddo
c
      return
      end
c
C  Adopted and modified from Copr. 1986-92 Numerical Recipes Software.
c
C#############################################################################################
c
      SUBROUTINE sprsin
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry3, only: ElementNeighbor,NumberofElementNeighbors
      use PBCG1
      use Variables2, only: ac,anb
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer index,i,j,k,k1,ArraySize
      double precision thresh
      data thresh/0.d0/
c********************************************************************************************
c
      ArraySize=size(sa)
c
      do j=1,NumberOfElements
c
        sa(j)=ac(j)
c
      enddo
c
      ija(1)=NumberOfElements+2
      k=NumberOfElements+1
c
      do i=1,NumberOfElements
        do j=1,NumberofElementNeighbors(i)
c
          k1 = ElementNeighbor(i,j)
c
          if(k1.ne.0) then   
c
            if(dabs(anb(i,j)).ge.thresh) then
c
              if(i.ne.k1) then
c
                k=k+1
                if(k.gt.ArraySize) then
c
                  print*, 'nmax too small in sprsin'
                  print*, i,k,ArraySize
                  pause
c
                endif
c
                sa(k)=anb(i,j)
                ija(k)= k1 
c
              endif
c
            endif
c
          endif
c
         enddo
         ija(i+1)=k+1
       enddo
c
      return
      end