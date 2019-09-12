c
C#############################################################################################
      FUNCTION reallocate1dINpointer(p,n)
C#############################################################################################
c
c  Reallocate a 1-dimensional intger pointer to a new size, preserving its previous contents.
c
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer, dimension(:), pointer :: p, reallocate1dINpointer
      integer, intent(in) :: n
      integer :: nold,ierr
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
      allocate(reallocate1dINpointer(n),stat=ierr)
      if (ierr /= 0) call nrerror
     *  ("reallocate1dINpointer: problem in attempt toallocate memory")
      if (.not. associated(p)) return
      nold=size(p)
      reallocate1dINpointer(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
c
      return
      end
