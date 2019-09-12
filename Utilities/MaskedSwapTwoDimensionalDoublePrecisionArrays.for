c
C#############################################################################################
      SUBROUTINE maskedSwap2dDParrays(a,b,mask)
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision, dimension(:,:), intent(INOUT) :: a,b
      logical, dimension(:,:), intent(IN) :: mask
      double precision, dimension(size(a,1),size(a,2)) :: swp
c*********************************************************************************************
c
      where (mask)
c
        swp=a
        a=b
        b=swp
c
      end where
c
      return
      end
