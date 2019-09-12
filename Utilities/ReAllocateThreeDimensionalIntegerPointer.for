c
C#############################################################################################
      FUNCTION reallocate3dINpointer(p,n,m,k)
C#############################################################################################
c
c  Reallocate a 3-dimensional integer pointer to a new size, preserving its previous contents.
c
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer, dimension(:,:,:), pointer :: p,reallocate3dINpointer
      integer, intent(in) :: n,m,k
      integer :: nold,mold,kold,ierr
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
      allocate(reallocate3dINpointer(n,m,k),stat=ierr)
      if (ierr /= 0) call nrerror
     *  ("reallocate3dINpointer: problem in attempt to allocate memory")
      if (.not. associated(p)) return
      nold=size(p,1)
      mold=size(p,2)
      kold=size(p,3)
      reallocate3dINpointer(1:min(nold,n),1:min(mold,m),1:min(kold,k))=
     *     p(1:min(nold,n),1:min(mold,m),1:min(kold,k))
      deallocate(p)
c
      return
      end
