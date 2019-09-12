MODULE reallocateData
implicit none
  interface
    FUNCTION reallocate1dDPpointer(p,n)
!
  ! Reallocate a 1-dimensional double precision pointer to a new size, preserving its previous contents.
!
      double precision, dimension(:), pointer :: p, reallocate1dDPointer
      integer, intent(in) :: n
      integer :: nold,ierr
    end FUNCTION reallocate1dDPpointer


    FUNCTION reallocate1dINpointer(p,n)
!
  ! Reallocate a 1-dimensional intger pointer to a new size, preserving its previous contents.
!
      integer, dimension(:), pointer :: p, reallocate1dINpointer
      integer, intent(IN) :: n
      integer :: nold,ierr
    end FUNCTION reallocate1dINpointer


    FUNCTION reallocate1dCHpointer(p,n)
!
  ! Reallocate a 1-dimensional character pointer to a new size, preserving its previous contents.
!
      character(1), dimension(:), pointer :: p, reallocate1dCHpointer
      integer, intent(IN) :: n
      integer :: nold,ierr
    end FUNCTION reallocate1dCHpointer


    FUNCTION reallocate2dDPpointer(p,n,m)
!
  ! Reallocate a 2-dimensional double precision pointer to a new size, preserving its previous contents.
!
      double precision, dimension(:,:), pointer :: p, reallocate2dDPpointer
      integer, intent(IN) :: n,m
      integer :: nold,mold,ierr
    end FUNCTION reallocate2dDPpointer

    FUNCTION reallocate2dINpointer(p,n,m)
!
  ! Reallocate a 2-dimensional integer pointer to a new size, preserving its previous contents.
!
      integer, dimension(:,:), pointer :: p, reallocate2dINpointer
      integer, intent(IN) :: n,m
      integer :: nold,mold,ierr
    end FUNCTION reallocate2dINpointer
!
    FUNCTION reallocate3dDPpointer(p,n,m,k)
!
  ! Reallocate a 3-dimensional double precision pointer to a new size, preserving its previous contents.
!
      implicit none
      double precision, dimension(:,:,:),pointer :: p, reallocate3dDPpointer
      integer, intent(in) :: n,m,k
      integer :: nold,mold,kold,ierr
      end FUNCTION reallocate3dDPpointer
!      
    FUNCTION reallocate3dINpointer(p,n,m,k)
!
!  Reallocate a 3-dimensional integer pointer to a new size, preserving its previous contents.
!
      implicit none
      integer, dimension(:,:,:), pointer :: p, reallocate3dINpointer
      integer, intent(in) :: n,m,k
      integer :: nold,mold,kold,ierr
    end FUNCTION reallocate3dINpointer
  end interface
end MODULE reallocateData 