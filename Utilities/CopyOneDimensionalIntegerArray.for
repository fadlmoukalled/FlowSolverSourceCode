c
C#############################################################################################
      SUBROUTINE copy1dINarrayUS(src,dest,ncopied,nNOTcopied)
C#############################################################################################
      implicit none
c*********************************************************************************************
      integer, dimension(:), intent(IN) :: src
      integer, dimension(:), intent(OUT) :: dest
      integer, intent(OUT) :: ncopied, nNOTcopied
c*********************************************************************************************
c
      ncopied=min(size(src),size(dest))
      nNOTcopied=size(src)-ncopied
      dest(1:ncopied)=src(1:ncopied)
c
      return
      end
