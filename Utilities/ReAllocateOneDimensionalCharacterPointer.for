c
C#############################################################################################
      FUNCTION reallocate1dCHpointer(p,n)
C#############################################################################################
c
c  Reallocate a 1-dimensional character pointer to a new size, preserving its previous contents.
c
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      character(1), dimension(:), pointer :: p, reallocate1dCHpointer
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
      allocate(reallocate1dCHpointer(n),stat=ierr)
      if (ierr /= 0) call nrerror
     *  ("reallocate1dCHpointer: problem in attempt to allocate memory")
      if (.not. associated(p)) return
      nold=size(p)
      reallocate1dCHpointer(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
c
      return
      end
