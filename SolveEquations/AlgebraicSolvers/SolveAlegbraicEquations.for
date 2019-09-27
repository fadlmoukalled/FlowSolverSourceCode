c
C#############################################################################################
c
      SUBROUTINE SolveEquation(Variable,FiT,rrF,itmax,solver,LMultigrid)
c
C#############################################################################################
c
      use user0, only: nIterStartApplyingMG,MGVariable,MGType
      use MultiGrid2, only: nIter
      use Residuals1
      use Variables2, only: dphi
      use DirectSolver1, only: b
      use Geometry1, only: NumberOfElements
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      character*35 Variable1
      character*6 solver
      logical LMultigrid
      integer NF,itmax,iter
      double precision rrF,InitialResorAbs,InitialResorMax
      double precision InitialResorRMS,InitialResorScaled
      double precision IResiduals,FResiduals
      double precision, dimension(:)  :: FiT
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE CalculateResiduals(NF,FiT)
C--------------------------------------------------------------------------------------
          integer NF
          double precision, dimension (:) :: FiT
C--------------------------------------------------------------------------------------
        end SUBROUTINE CalculateResiduals
C--------------------------------------------------------------------------------------
        SUBROUTINE UpdateSourceResistance(iter,FiT)
C--------------------------------------------------------------------------------------
          integer iter
          double precision, dimension (:) :: FiT
C--------------------------------------------------------------------------------------
        end SUBROUTINE UpdateSourceResistance
C--------------------------------------------------------------------------------------
      end interface
C********************************************************************************************
c
      call GetNFfromVariable(Variable,NF,Variable1)
      dphi=0.
      call CalculateResiduals(NF,FiT)
c
      if(nIter.eq.nIterStartApplyingMG.and.
     *                   LMultigrid.and.MGType.eq.'algebraic') then
c
        if(Variable1.eq.MGVariable) then 
c
          IResiduals=ResorAbs(NF)
          FResiduals=IResiduals
c
          call AlgebraicMultigrid
     *          (NF,rrF,itmax,solver,IResiduals,FResiduals,Variable1)
c
c---- Update solution
c
          FiT=FiT+dphi
c
        endif
c
      elseif(nIter.gt.nIterStartApplyingMG.and.
     *                   LMultigrid.and.MGType.eq.'algebraic') then
c
        IResiduals=ResorAbs(NF)
        FResiduals=IResiduals
c
        call AlgebraicMultigrid
     *          (NF,rrF,itmax,solver,IResiduals,FResiduals,Variable1)
c
c---- Update solution
c
        FiT=FiT+dphi
c
      elseif(nIter.ge.nIterStartApplyingMG.and.
     *                LMultigrid.and.MGType.eq.'geometricelement') then
c
        IResiduals=ResorAbs(NF)
        FResiduals=IResiduals
c
        call AlgebraicMultigrid
     *          (NF,rrF,itmax,solver,IResiduals,FResiduals,Variable1)
c
c---- Update solution
c
        FiT=FiT+dphi
c
      elseif(nIter.ge.nIterStartApplyingMG.and.
     *                LMultigrid.and.MGType.eq.'geometricnode') then
c
        IResiduals=ResorAbs(NF)
        FResiduals=IResiduals
c
        call AlgebraicMultigrid
     *          (NF,rrF,itmax,solver,IResiduals,FResiduals,Variable1)
c
c---- Update solution
c
        FiT=FiT+dphi
c
      else
c
        IResiduals=ResorAbs(NF)
        FResiduals=IResiduals
c
        if(solver.eq.'sor') then      
c
          iter=0
          do while(FResiduals.GT.rrF*IResiduals.and.iter.lt.itmax)
c
            if(NF.eq.4.and.iter.ne.0) 
     *                call UpdateSourceResistance(iter,FiT)   !should be implemented for multigrid
            call SolveEquationsUsingSOR                       ! and direct solver
            call CheckConvergence(NF,FResiduals)
c          
            iter=iter+1
          enddo
c
        elseif(solver.eq.'ilu') then      
c
          iter=0
          call ILUFactorization
          do while(FResiduals.GT.rrF*IResiduals.and.iter.lt.itmax)
c
            if(NF.eq.4.and.iter.ne.0) 
     *                call UpdateSourceResistance(iter,FiT)
            call SolveEquationsUsingILU
            call CheckConvergence(NF,FResiduals)
c
            iter=iter+1
          enddo
c
        elseif(solver.eq.'pbcg') then      
c
          iter=0
          call allocateForPBCG
          call sprsin
          call ILUFactorization
c
          do while(FResiduals.GT.rrF*IResiduals.and.iter.lt.itmax)
c
            if(NF.eq.4.and.iter.ne.0) 
     *                call UpdateSourceResistance(iter,FiT)
            call SolveEquationsUsingPBCG(iter)
            call CheckConvergence(NF,FResiduals)
c
            iter=iter+1
          enddo
c
        elseif(solver.eq.'direct') then 
c
          call SetupMatrixforDirectSolver(NumberOfElements)
          call ludcmpF77
          call lubksbF77
          dphi=b
          call deallocateMatrixStorage
c
        endif
c
c---- Update solution
c
c        print*,dphi
        FiT=FiT+dphi
c
      endif
c
      return
      end
