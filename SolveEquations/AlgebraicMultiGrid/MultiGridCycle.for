c
C#############################################################################################
c
      SUBROUTINE vMultiGridCycle(NF,rrF,itmax,solver)
c
C###################################################################################################
C  Start V cycle
C
C             1     @               @
C                    \             /
C             2       @           @
C                      \         /
C             3         @       @
C                        \     /
C             4           @   @
C                          \ /
C             5             @
C###################################################################################################
C
      use MultiGrid2, only: iLevelMax
      use User0, only: MultiGridpreSweep,MultiGridpostSweep
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*6 solver
      integer NF,itmax,iLevel 
      double precision rrF
C--------------------------------------------------------------------------------------
c
c.....Restriction
c
      do iLevel=1,iLevelMax-1
c
        call Restrict(NF,iLevel,rrF,MultiGridpreSweep,solver)
c
      enddo
c
      call SolveEquationMG(NF,iLevelMax,rrF,MultiGridpostSweep,solver)
      call CalculateResidualsMG(iLevelMax)
c
c.....Prolongation
c
      do iLevel=iLevelMax,2,-1
c
        call Prolongate(NF,iLevel,rrF,MultiGridpostSweep,solver)
c
      enddo


      return
      end
c
C#############################################################################################
c
      SUBROUTINE fMultiGridCycle(NF,rrF,itmax,solver)
c
C
C###################################################################################################
C  Start F cycle
C
C             1     @                                       @
C                    \                                     /
C             2       @                       @           @
C                      \                     / \         /
C             3         @           @       @   @       @
C                        \         / \     /     \     /
C             4           @   @   @   @   @       @   @
C                          \ / \ /     \ /         \ /
C             5             @   @       @           @
C###################################################################################################
c
      use MultiGrid2, only: iLevelMax
      use User0, only: MultiGridpreSweep,MultiGridpostSweep
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*6 solver
      integer NF,itmax,iLevel,iLevel1
      double precision rrF
c********************************************************************************************
c
c.....Restriction
c
      do iLevel=1,iLevelMax-1
c
        call Restrict(NF,iLevel,rrF,MultiGridpreSweep,solver)
c
      enddo
c
      call SolveEquationMG(NF,iLevelMax,rrF,MultiGridpostSweep,solver)
      call CalculateResidualsMG(iLevelMax)
c
      iLevel1=iLevelMax-1
c
      do while(iLevel1.ge.2)
c
c.....Prolongation
c
        do iLevel=iLevelMax,iLevel1+1,-1
C
          call Prolongate(NF,iLevel,rrF,MultiGridpostSweep,solver)
C
        enddo
c
c.....Restriction
c
        do iLevel=iLevel1,iLevelMax-1
c
          call Restrict(NF,iLevel,rrF,MultiGridpreSweep,solver)
c
        enddo
c
        call SolveEquationMG(NF,iLevelMax,rrF,MultiGridpostSweep,solver)
        call CalculateResidualsMG(iLevelMax)
c
        iLevel1=iLevel1-1
c
      enddo 
c
c.....Prolongation
c
      do iLevel=iLevelMax,2,-1
c
        call Prolongate(NF,iLevel,rrF,MultiGridpostSweep,solver)
c
      enddo
c      
      return
      end
c
C#############################################################################################
c
      SUBROUTINE wMultiGridCycle(NF,rrF,itmax,solver)
c
C
C###################################################################################################
C  Start W cycle
C
C             1     @                                                   @
C                    \                                                 /
C             2       @                       @                       @
C                      \                     / \                     /
C             3         @           @       @   @       @           @
C                        \         / \     /     \     / \         /
C             4           @   @   @   @   @       @   @   @   @   @
C                          \ / \ /     \ /         \ /     \ / \ /
C             5             @   @       @           @       @   @
C###################################################################################################
c
      use MultiGrid2, only: iLevelMax
      use User0, only: MultiGridpreSweep,MultiGridpostSweep
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*6 solver
      integer NF,itmax,iLevel,iLevel1
      double precision rrF
c********************************************************************************************
c
c.....Restriction
c
      do iLevel=1,iLevelMax-1
c
        call Restrict(NF,iLevel,rrF,MultiGridpreSweep,solver)
c
      enddo
c
      call SolveEquationMG(NF,iLevelMax,rrF,MultiGridpostSweep,solver)
      call CalculateResidualsMG(iLevelMax)
c
      iLevel1=iLevelMax-1
c
      do while(iLevel1.ge.2)
c
c.....Prolongation
c
        do iLevel=iLevelMax,iLevel1+1,-1
C
          call Prolongate(NF,iLevel,rrF,MultiGridpostSweep,solver)
C
        enddo
c
c.....Restriction
c
        do iLevel=iLevel1,iLevelMax-1
c
          call Restrict(NF,iLevel,rrF,MultiGridpreSweep,solver)
c
        enddo
c
        call SolveEquationMG(NF,iLevelMax,rrF,MultiGridpostSweep,solver)
        call CalculateResidualsMG(iLevelMax)
c
        iLevel1=iLevel1-1
c
      enddo
c
      iLevel1=3
c
      Do while(iLevel1.le.iLevelMax-1)
c
c.....Prolongation
c
        do iLevel=iLevelMax,iLevel1+1,-1
c
          call Prolongate(NF,iLevel,rrF,MultiGridpostSweep,solver)
c
        enddo
c
c.....Restriction
c
        do iLevel=iLevel1,iLevelMax-1
c
          call Restrict(NF,iLevel,rrF,MultiGridpreSweep,solver)
c
        enddo
c
        call SolveEquationMG(NF,iLevelMax,rrF,MultiGridpostSweep,solver)
        call CalculateResidualsMG(iLevelMax)
c
        iLevel1=iLevel1+1
c
      enddo
c
c.....Prolongation
c
      do iLevel=iLevelMax,2,-1
c
        call Prolongate(NF,iLevel,rrF,MultiGridpostSweep,solver)
c
      enddo
c      
      return
      end