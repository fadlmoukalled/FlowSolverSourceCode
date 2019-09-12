c
C#############################################################################################
      SUBROUTINE maskedSwap1dDParrays(a,b,mask)
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision, dimension(:), intent(INOUT) :: a,b
      logical, dimension(:), intent(IN) :: mask
      double precision, dimension(size(a)) :: swp
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
