c
C#############################################################################################
      SUBROUTINE swapDPvariables(a,b)
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision, intent(INOUT) :: a,b
      double precision :: dum
c*********************************************************************************************
c
      dum=a
      a=b
      b=dum
c
      return
      end
