c
C#############################################################################################
c
      SUBROUTINE InitializeCoefficients
c
C#############################################################################################
c
      use User0, only: LUnsteady
      use Variables2
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      ac=0.
      anb=0.
      bc=0.
c
      if(LUnsteady) then
        acold=0.
        acoldold=0.
      endif
c
      return
      end
