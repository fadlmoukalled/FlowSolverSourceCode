c
c  Adopted and modified from Copr. 1986-92 Numerical Recipes Software.
c
c*********************************************************************************************
c
c--- A direct algebraic equation solver
c
C#############################################################################################
c
      SUBROUTINE SetupMatrixforDirectSolver(NumberOfElements)
c
C#############################################################################################
c
      use DirectSolver1
      use Variables2, only:ac,anb,bc
      use Geometry3, only:ElementNeighbor,NumberofElementNeighbors
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      integer :: NumberOfElements
c********************************************************************************************
c      
      allocate(a(NumberOfElements,NumberOfElements))
      allocate(b(NumberOfElements))
      allocate(indx(NumberOfElements))
c
c--- Initialize
c
      a=0.
      b=0.
      indx=0
c
c-- Start computations from first to the last element
c
      do i=1,NumberOfElements
c
        b(i)=bc(i)
        a(i,i)=ac(i)
      enddo
c
      do i=1,NumberOfElements
        do j=1,NumberofElementNeighbors(i)
c
          k = ElementNeighbor(i,j)
c
          if(k.ne.0) then
c
            a(i,k)=anb(i,j)
c
          endif
        enddo
      enddo
c               
      return
      end
C#############################################################################################
c
      SUBROUTINE deallocateMatrixStorage
c
C#############################################################################################
c
      use DirectSolver1
c
c*********************************************************************************************
      implicit none      
c*********************************************************************************************
c
      deallocate(a)
      deallocate(b)
      deallocate(indx)
c
      return
      end
c
c*********************************************************************************************
c
c--- A direct algebraic equation solver (Adopted from Numerical Recipes Fortran 77)
c
C#############################################################################################
c
      SUBROUTINE ludcmpF77
c
C#############################################################################################
c
      use DirectSolver1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision, dimension(size(a,1)) :: vv   !vv stores the implicit scaling of each row.
      double precision, parameter :: tiny=1.0e-20    !A small number.
      integer :: i,imax,j,k,n
      double precision :: aamax,dum,sum,d
      
c*********************************************************************************************
      interface
c*********************************************************************************************
        FUNCTION asserteq(n1,n2,n3,string)
          character(len=*), intent(in) :: string
          integer, intent(in) :: n1,n2,n3
          integer :: asserteq
        end FUNCTION asserteq
      end interface
c
      interface
        SUBROUTINE nrerror(string)
          character(LEN=*), intent(in) :: string
        end SUBROUTINE nrerror
      end interface
c
c*********************************************************************************************
c
      n=asserteq(size(a,1),size(a,2),size(indx),"ludcmp")
      d=1.
c
      do i=1,n
c
        aamax=0.
c
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0.) call nrerror("singular matrix in ludcmp")
        vv(i)=1./aamax
      enddo
c
      do j=1,n
c
        do i=1,j-1
          sum=a(i,j)
          do  k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
c
        aamax=0.
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
c
        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.) a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
c
      enddo
c
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software
c
C#############################################################################################
      SUBROUTINE lubksbF77
C#############################################################################################
c
      use DirectSolver1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,n,ii,ll,j
      double precision :: sum
c*********************************************************************************************
      interface
c*********************************************************************************************
        FUNCTION asserteq(n1,n2,n3,string)
c*********************************************************************************************
          character(len=*), intent(in) :: string
          integer, intent(in) :: n1,n2,n3
          integer :: asserteq
c*********************************************************************************************
        end FUNCTION asserteq
c*********************************************************************************************
      end interface
c*********************************************************************************************
c      
      n=asserteq(size(a,1),size(a,2),size(indx),"lubksb")
c
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
      enddo
c
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
c
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software P$,-15+<Z$.
c
c*********************************************************************************************
c
c--- A direct algebraic equation solver (Adopted from Numerical Recipes Fortran 90)
c
C#############################################################################################
      SUBROUTINE ludcmp
C#############################################################################################
c
      use DirectSolver1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision :: d
      double precision, dimension(size(a,1)) :: vv   !vv stores the implicit scaling of each row.
      double precision, parameter :: tiny=1.0e-20    !A small number.
      integer :: j,n,imax
c*********************************************************************************************
      interface
c*********************************************************************************************
        SUBROUTINE nrerror(string)
          character(LEN=*), intent(in) :: string
        end SUBROUTINE nrerror
      end interface
c
      interface
        FUNCTION asserteq(n1,n2,n3,string)
          character(len=*), intent(in) :: string
          integer, intent(in) :: n1,n2,n3
          integer :: asserteq
        end FUNCTION asserteq
      end interface
c
      interface
        FUNCTION outerprod(e,f)
          double precision, dimension(:), intent(in) :: e,f
          double precision, dimension(size(e),size(f)) :: outerprod
        end FUNCTION outerprod
      end interface
c
      interface
        FUNCTION imaxloc(arr)
          double precision, dimension(:), intent(in) :: arr
          integer :: imaxloc
          integer, dimension(1) :: imax
        end FUNCTION imaxloc
      end interface
c
      interface
        SUBROUTINE swap1dDParrays(a,b)
          double precision, dimension(:), intent(INOUT) :: a,b
          double precision, dimension(size(a)) :: dum
        end SUBROUTINE swap1dDParrays
      end interface
c
c*********************************************************************************************
c     Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
c     rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
c     output vector of length N that records the row permutation effected by the partial pivoting;
c     d is output as ±1 depending on whether the number of row interchanges was even or odd,
c     respectively. This routine is used in combination with lubksb to solve linear equations or
c     invert a matrix.
c
c*********************************************************************************************
c
      n=asserteq(size(a,1),size(a,2),size(indx),"ludcmp")
      d=1.0                                            !No row interchanges yet.
      vv=maxval(dabs(a),dim=2)                         !Loop over rows to get the implicit scaling
      if (any(vv == 0.0)) call nrerror("singular matrix in ludcmp")      !information.
                                                       !There is a row of zeros.
      vv=1.0/vv                                        !Save the scaling.
      do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*dabs(a(j:n,j)))     !Find the pivot row.
        if (j /= imax) then                            !Do we need to interchange rows?
          call swap1dDParrays(a(imax,:),a(j,:))                  !Yes, do so...
          d=-d                                         !...and change the parity of d.
          vv(imax)=vv(j)                               !Also interchange the scale factor.
        endif
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
c
c      If the pivot element is zero the matrix is singular (at least to 
c      the precision of the algorithm).For some applications on singular 
c      matrices, it is desirable to substitute TINY for zero.
c
        a(j+1:n,j)=a(j+1:n,j)/a(j,j)                   !Divide by the pivot element.
c        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
                                                        !Reduce remaining submatrix.
      enddo
      return
      end 
c
C#############################################################################################
      SUBROUTINE lubksb
C#############################################################################################
c
      use DirectSolver1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,n,ii,ll
      double precision :: summ
c*********************************************************************************************
      interface
c*********************************************************************************************
        FUNCTION asserteq(n1,n2,n3,string)
c*********************************************************************************************
          character(len=*), intent(in) :: string
          integer, intent(in) :: n1,n2,n3
          integer :: asserteq
c*********************************************************************************************
        end FUNCTION asserteq
c*********************************************************************************************
      end interface
c*********************************************************************************************
c
c     Solves the set of N linear equations A · X = B. Here the N × N 
c     matrix a is input, not as the original matrix A, but rather as its
c     LU decomposition, determined by the routine ludcmp. indx is input 
c     as the permutation vector of length N returned by ludcmp. b is 
c     input as the right-hand-side vector B, also of length N, and 
c     returns with the solution vectorX. a and indx are not modified by 
c     this routine and can be left in place for successive calls with 
c     different right-hand sides b. This routine takes into account the 
c     possibility that b will begin with many zero elements, 
c     so it is efficient for use in matrix inversion.
c
c*********************************************************************************************
c
      n=asserteq(size(a,1),size(a,2),size(indx),"lubksb")
      ii=0                                    
c
c*********************************************************************************************
c
c     When ii is set to a positive value, it will become the index
c     of the first nonvanishing element of b. We now do the forward 
c     substitution. The only new  wrinkle is to unscramble the 
c     permutation as we go.
c
c*********************************************************************************************
c
      do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
      if (ii /= 0) then
        summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
      elseif (summ /= 0.0) then
        ii=i                  !A nonzero element was encountered, 
      endif                   !so from now on we will  
        b(i)=summ             !have to do the dot product above.
      enddo
      do i=n,1,-1             !Now we do the backsubstitution
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
      enddo
      return
      end 
C#############################################################################################
      FUNCTION imaxloc(arr)
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision, dimension(:), intent(in) :: arr
      integer :: imaxloc
      integer, dimension(1) :: imax
c*********************************************************************************************
      imax=maxloc(arr(:))
      imaxloc=imax(1)
c
      end FUNCTION imaxloc
c
C#############################################################################################
      FUNCTION outerprod(e,f)
C#############################################################################################
c
      implicit none
c*********************************************************************************************
      double precision, dimension(:), intent(in) :: e,f
      double precision, dimension(size(e),size(f)) :: outerprod
c*********************************************************************************************
      outerprod = spread(e,dim=2,ncopies=size(f))*
     *            spread(f,dim=1,ncopies=size(e))
      
      end FUNCTION outerprod
c
C#############################################################################################
      FUNCTION asserteq(n1,n2,n3,string)
C#############################################################################################
c
      implicit none
c*********************************************************************************************
      character(len=*), intent(in) :: string
      integer, intent(in) :: n1,n2,n3
      integer :: asserteq
c*********************************************************************************************
      if (n1 == n2 .and. n2 == n3) then
        asserteq=n1
      else
        write (*,*) "nrerror: an assert_eq failed with this tag:", 
     *               "program terminated by asserteq"
         stop
      endif
      end FUNCTION asserteq

