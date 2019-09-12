c
C#############################################################################################
c
      SUBROUTINE Restrict(NF,iLevel,rrF,itmax,solver)
c
C#############################################################################################
c
      use User0
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*6 solver
      integer NF,itmax,iLevel
      double precision rrF
c********************************************************************************************
c
      call SolveEquationMG(NF,iLevel,rrF,itmax,solver)
c
      call CalculateResidualsMG(iLevel)
      call AssembleAgglomeratedRHS(iLevel+1)
c
      return
      end
