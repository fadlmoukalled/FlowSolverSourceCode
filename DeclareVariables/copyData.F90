MODULE copyData
implicit none
!
  interface
    SUBROUTINE copy1dDParrayUS(src,dest,ncopied,nNOTcopied)
      double precision, dimension(:), intent(IN) :: src
      double precision, dimension(:), intent(OUT) :: dest
      integer, intent(OUT) :: ncopied, nNOTcopied
    END SUBROUTINE copy1dDParrayUS
!
    SUBROUTINE copy1dINarrayUS(src,dest,ncopied,nNOTcopied)
      integer, dimension(:), intent(IN) :: src
      integer, dimension(:), intent(OUT) :: dest
      integer, intent(OUT) :: ncopied, nNOTcopied
    END SUBROUTINE copy1dINarrayUS
  end interface

end MODULE copyData
