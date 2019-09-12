c
C#############################################################################################
      FUNCTION reallocate2dDPpointer(p,n,m)
C#############################################################################################
c
  ! Reallocate a 2-dimensional double precision pointer to a new size, preserving its previous contents.
c
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision, dimension(:,:), 
     *                         pointer :: p, reallocate2dDPpointer
      integer, intent(in) :: n,m
      integer :: nold,mold,ierr
c*********************************************************************************************
      interface
c*********************************************************************************************
        SUBROUTINE nrerror(string)
          CHARACTER(LEN=*), INTENT(IN) :: string
        end SUBROUTINE nrerror
c*********************************************************************************************
      end interface
c*********************************************************************************************
c
      allocate(reallocate2dDPpointer(n,m),stat=ierr)
      if (ierr /= 0) call nrerror
     *  ("reallocate2dDPpointer: problem in attempt to allocate memory")
      if (.not. associated(p)) return
      nold=size(p,1)
      mold=size(p,2)
      reallocate2dDPpointer(1:min(nold,n),1:min(mold,m))=
     *                              p(1:min(nold,n),1:min(mold,m))
      deallocate(p)
c
      return
      end
