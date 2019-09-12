c
C#############################################################################################
      FUNCTION reallocate2dINpointer(p,n,m)
C#############################################################################################
c
c  Reallocate a 2-dimensional integer pointer to a new size, preserving its previous contents.
c
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer, dimension(:,:), pointer :: p, reallocate2dINpointer
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
      allocate(reallocate2dINpointer(n,m),stat=ierr)
      if (ierr /= 0) call nrerror
     *  ("reallocate2dINpointer: problem in attempt to allocate memory")
      if (.not. associated(p)) return
      nold=size(p,1)
      mold=size(p,2)
      reallocate2dINpointer(1:min(nold,n),1:min(mold,m))=
     *     p(1:min(nold,n),1:min(mold,m))
      deallocate(p)
c
      return
      end
