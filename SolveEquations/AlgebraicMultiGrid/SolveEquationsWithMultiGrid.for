c
C#############################################################################################
c
      SUBROUTINE SolveEquationMG(NF,iLevel,rrF,itmax,solver)
c
C#############################################################################################
c
      use Residuals1
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*6 solver
      integer NF,itmax,iter,iLevel 
      double precision rrF
      double precision IResiduals,FResiduals
c********************************************************************************************
c
      if(NF.eq.-1) then     !wall distance
        call CheckConvergence(NF,IResiduals)
      else
        IResiduals=ResorAbs(NF)
      endif
c      
      FResiduals=IResiduals
c
      if(solver.eq.'sor') then      
c
        iter=0
        do while(FResiduals.GT.rrF*IResiduals.and.iter.lt.itmax)
c
          call SolveEquationsUsingSORMG(iLevel)
          call CheckConvergenceMG(NF,iLevel,FResiduals)
c          
          iter=iter+1
        enddo
c
      elseif(solver.eq.'ilu') then      
c
        iter=0
        call ILUFactorizationMG(iLevel)
        do while(FResiduals.GT.rrF*IResiduals.and.iter.lt.itmax)
c
          call SolveEquationsUsingILUMG(iLevel)
          call CheckConvergenceMG(NF,iLevel,FResiduals)
c
          iter=iter+1
        enddo
c
      endif
c
      return
      end
