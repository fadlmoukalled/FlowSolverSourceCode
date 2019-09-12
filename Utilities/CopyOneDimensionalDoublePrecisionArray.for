c
C#############################################################################################
      SUBROUTINE copy1dDParrayUS(src,dest,ncopied,nNOTcopied)
C#############################################################################################
c---  Copy array where size of source not known in advance.
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision, dimension(:), intent(IN) :: src
      double precision, dimension(:), intent(OUT) :: dest
      integer, intent(OUT) :: ncopied, nNOTcopied
c*********************************************************************************************
c
      ncopied=min(size(src),size(dest))
      nNOTcopied=size(src)-ncopied
      dest(1:ncopied)=src(1:ncopied)
c
      return
      end
