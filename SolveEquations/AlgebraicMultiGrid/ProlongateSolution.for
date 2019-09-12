c
C#############################################################################################
c
      SUBROUTINE Prolongate(NF,iLevel,rrF,itmax,solver)
c
C#############################################################################################
c
      implicit none
c********************************************************************************************
      character*6 solver
      integer NF,itmax,iLevel 
      double precision rrF
c********************************************************************************************
c
      call CorrectFinerLevelSolution(iLevel)
      call SolveEquationMG(NF,iLevel-1,rrF,itmax,solver)
c
      return
      end
