MODULE swapData

implicit none

  interface
    SUBROUTINE swapINvariables(a,b)
      integer, intent(INOUT) :: a,b
      integer :: dum
    end SUBROUTINE swapINvariables
!
    SUBROUTINE swapDPvariables(a,b)
      double precision, intent(INOUT) :: a,b
      double precision :: dum
    end SUBROUTINE swapDPvariables
!
    SUBROUTINE swap1dDParrays(a,b)
      double precision, dimension(:), intent(INOUT) :: a,b
      double precision, dimension(SIZE(a)) :: dum
    end SUBROUTINE swap1dDParrays
!
    SUBROUTINE maskedSwapDPvariables(a,b,mask)
      double precision, intent(INOUT) :: a,b
      logical, intent(IN) :: mask
      double precision :: swp
    end SUBROUTINE maskedSwapDPvariables
!
    SUBROUTINE maskedSwap1dDParrays(a,b,mask)
      double precision, dimension(:), intent(INOUT) :: a,b
      logical, dimension(:), intent(IN) :: mask
      double precision, dimension(size(a)) :: swp
    end SUBROUTINE maskedSwap1dDParrays
!
    SUBROUTINE maskedSwap2dDParrays(a,b,mask)
      double precision, dimension(:,:), intent(INOUT) :: a,b
      logical, dimension(:,:), intent(IN) :: mask
      double precision, dimension(size(a,1),size(a,2)) :: swp
    end SUBROUTINE maskedSwap2dDParrays
  end interface

end MODULE swapData