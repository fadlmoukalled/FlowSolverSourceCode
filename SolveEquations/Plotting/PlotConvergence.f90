!-------------------------------------------------------------------------------
!    GnuPlot Interface
!-------------------------------------------------------------------------------
subroutine PlotConvergenceHistory
        use Residuals1
        use MultiGrid2, only: nIter
        use User0, only: LSolveMomentum,LSolveContinuity,LSolveEnergy,LSolveTurbulenceKineticEnergy,&
                         LSolveTurbulenceDissipationRate,LSolveTurbulenceSpecificDissipationRate,&
                         LSolveTurbulentKL,LSolveModifiedED,LSolveTurbulenceGammaEquation,&
                         LSolveTurbulenceReynoldsThetaEquation,LSolveTurbulencefRelaxationEquation,&
                         LSolveTurbulenceV2Equation,LSolveTurbulenceZetaEquation,LSolveLambdaELEEquation,&
                         LSolverField,LSolveScalar,NumberOfrFieldsToSolve,NumberOfScalarsToSolve,NstopType,&
                         rFieldName,ScalarName,LReadOldSolution
  implicit none
  integer :: i,j,k
  integer, save :: nIterOld
  integer :: IOstatus
  character*30, dimension(:), allocatable :: titleOld
  double precision, dimension(:), allocatable :: xOld
  character(len = 200) :: PlotLine
  character (len = 30) :: nVariablesChar
!
1       format (i7,<nVariablesPlot>e20.12)
5       format(1x,' "iteration"')
10      format(1x,' "x-momentum"')
11      format(1x,' "y-momentum"')
12      format(1x,' "z-momentum"')
20      format(1x,' "continuity"')
30      format(1x,' "energy"')
40      format(1x,' "k"')
50      format(1x,' "epsilon"')
60      format(1x,' "omega"')
70      format(1x,' "KL"')
80      format(1x,' "EddyDiffusivity"')
90      format(1x,' "Gamma"')
100     format(1x,' "ReTheta"')
110     format(1x,' "fRelaxation"')
120     format(1x,' "v2"')
130     format(1x,' "zeta"')
140     format(1x,' "lambda"')
150     format(1x,'"',A10,'" ')
!  
   if(nIter.eq.1) then 
        k=nVariablesPlot+1
        write(nVariablesChar , *) k
!    
        PlotLine = 'plot for [i=2:'//trim(adjustl(nVariablesChar))//'] "data.txt" using 1:i with lines title columnheader(i)'
!--------------------------------------------------------------------------!
        nIterOld=0
        if(LReadOldSolution) then
          open ( unit = 24, file = "data.txt")
!          
          write(24,5,advance='no') 
!
          if(LSolveMomentum) then
            write(24,10,advance='no') 
            write(24,11,advance='no') 
            write(24,12,advance='no') 
          endif
!      
          if(LSolveContinuity) then
            write(24,20,advance='no') 
          endif
!        
          if(LSolveEnergy) then
            write(24,30,advance='no') 
          endif
!        
          if(LSolveTurbulenceKineticEnergy) then
            i=i+1
            write(24,40,advance='no') 
          endif
!        
          if(LSolveTurbulenceDissipationRate) then
            write(24,50,advance='no') 
          endif
!        
          if(LSolveTurbulenceSpecificDissipationRate) then
            write(24,60,advance='no') 
          endif
!        
          if(LSolveTurbulentKL) then
            write(24,70,advance='no') 
          endif
!        
          if(LSolveModifiedED) then
            write(24,80,advance='no') 
          endif
!        
          if(LSolveTurbulenceGammaEquation) then
            write(24,90,advance='no') 
          endif
!        
          if(LSolveTurbulenceReynoldsThetaEquation) then
            write(24,100,advance='no') 
          endif
!        
          if(LSolveTurbulencefRelaxationEquation) then
            write(24,110,advance='no') 
          endif
!        
          if(LSolveTurbulenceV2Equation) then
            write(24,120,advance='no') 
          endif
!        
          if(LSolveTurbulenceZetaEquation) then
            write(24,130,advance='no') 
          endif
!        
          if(LSolveLambdaELEEquation) then
            write(24,140,advance='no') 
          endif
!        
          do j=1,NumberOfrFieldsToSolve
            if(LSolverField(j)) then
              write(24,150,advance='no') rFieldName(j) 
            endif
          enddo
!
          do j=1,NumberOfScalarsToSolve
            if(LSolveScalar(j)) then
              write(24,150,advance='no') ScalarName(j) 
            endif
          enddo
!   
          write(24,*) ""
          
          allocate(xOld(nVariablesPlot))
          allocate(titleOld(k))
          read(25,*,iostat=IOstatus) (titleOld(i),i=1,k)
          do while(IOstatus==0)
            read(25,1,iostat=IOstatus) nIterOld,(xOld(i),i=1,nVariablesPlot)
            if(nIterOld.ne.0.and.IOstatus==0) write(24,1) nIterOld,(xOld(i),i=1,nVariablesPlot)
          enddo
          deallocate(xOld)
          deallocate(titleOld)
        else
          open ( unit = 24, file = "data.txt")
!write titles of the graphs here
!
          write(24,5,advance='no') 
!
          if(LSolveMomentum) then
            write(24,10,advance='no') 
            write(24,11,advance='no') 
            write(24,12,advance='no') 
          endif
!      
          if(LSolveContinuity) then
            write(24,20,advance='no') 
          endif
!        
          if(LSolveEnergy) then
            write(24,30,advance='no') 
          endif
!        
          if(LSolveTurbulenceKineticEnergy) then
            i=i+1
            write(24,40,advance='no') 
          endif
!        
          if(LSolveTurbulenceDissipationRate) then
            write(24,50,advance='no') 
          endif
!        
          if(LSolveTurbulenceSpecificDissipationRate) then
            write(24,60,advance='no') 
          endif
!        
          if(LSolveTurbulentKL) then
            write(24,70,advance='no') 
          endif
!        
          if(LSolveModifiedED) then
            write(24,80,advance='no') 
          endif
!        
          if(LSolveTurbulenceGammaEquation) then
            write(24,90,advance='no') 
          endif
!        
          if(LSolveTurbulenceReynoldsThetaEquation) then
            write(24,100,advance='no') 
          endif
!        
          if(LSolveTurbulencefRelaxationEquation) then
            write(24,110,advance='no') 
          endif
!        
          if(LSolveTurbulenceV2Equation) then
            write(24,120,advance='no') 
          endif
!        
          if(LSolveTurbulenceZetaEquation) then
            write(24,130,advance='no') 
          endif
!        
          if(LSolveLambdaELEEquation) then
            write(24,140,advance='no') 
          endif
!        
          do j=1,NumberOfrFieldsToSolve
            if(LSolverField(j)) then
              write(24,150,advance='no') rFieldName(j) 
            endif
          enddo
!
          do j=1,NumberOfScalarsToSolve
            if(LSolveScalar(j)) then
              write(24,150,advance='no') ScalarName(j) 
            endif
          enddo
!   
          write(24,*) ""
!
        endif
!    
!--------------------------------------------------------------------------!
        open ( unit = 22, file = "commands2.txt")
        write(22,*) 'a = 0'
        write(22,*) 'trim(PlotLine)'
        close(22) 
!--------------------------------------------------------------------------!
!Write the options here
        open ( unit = 21, file = "commands1.txt")
        write(21,*) 'set title "Convergence History" '
        write(21,*) 'set xlabel "Iteration" '
        if(NstopType.eq.1) then
          write(21,*) 'set ylabel "Absolute Residuals" '
        elseif(NstopType.eq.2) then
          write(21,*) 'set ylabel "Maximum Residuals" '
        elseif(NstopType.eq.3) then
          write(21,*) 'set ylabel "RMS Residuals" '
        elseif(NstopType.eq.4) then
          write(21,*) 'set ylabel "Normalized Residuals" '
        endif
        k=nVariablesPlot
        write(nVariablesChar , *) k
        write(21,*) 'set for [i=1:'//trim(adjustl(nVariablesChar))//'] linetype i lw 3'
        write(21,*) 'set grid'
        write(21,*) 'set log y'
        write(21,*) 'load "commands2.txt" ' 
        write(21,*) 'if(a==1) '//trim(PlotLine)
        write(21,*) 'reread'
        close(21)
!--------------------------------------------------------------------------!
        call system('start "" gnuplot commands1.txt ')
!--------------------------------------------------------------------------!
    endif
!
    if(NstopType.eq.1) then
!
        i=0
        if(LSolveMomentum) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(1)
            i=i+1
            ResidualsPlot(i)=ResorAbs(2)
            i=i+1
            ResidualsPlot(i)=ResorAbs(3)
        endif
!        
        if(LSolveContinuity) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(5)
        endif
!        
        if(LSolveEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(6)
        endif
!        
        if(LSolveTurbulenceKineticEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(8)
        endif
!        
        if(LSolveTurbulenceDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(9)
        endif
!        
        if(LSolveTurbulenceSpecificDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(10)
        endif
!        
        if(LSolveTurbulentKL) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(12)
        endif
!        
        if(LSolveModifiedED) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(11)
        endif
!        
        if(LSolveTurbulenceGammaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(13)
        endif
!        
        if(LSolveTurbulenceReynoldsThetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(14)
        endif
!        
        if(LSolveTurbulencefRelaxationEquation) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(17)
        endif
!        
        if(LSolveTurbulenceV2Equation) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(15)
        endif
!        
        if(LSolveTurbulenceZetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(16)
        endif
!        
        if(LSolveLambdaELEEquation) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(18)
        endif
!        
        do j=1,NumberOfrFieldsToSolve
          if(LSolverField(j)) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(18+j)
          endif
        enddo
!
        do j=1,NumberOfScalarsToSolve
          if(LSolveScalar(j)) then
            i=i+1
            ResidualsPlot(i)=ResorAbs(18+NumberOfrFieldsToSolve+j)
          endif
        enddo
!
    elseif(NstopType.eq.2) then
!
        i=0
        if(LSolveMomentum) then
            i=i+1
            ResidualsPlot(i)=ResorMax(1)
            i=i+1
            ResidualsPlot(i)=ResorMax(2)
            i=i+1
            ResidualsPlot(i)=ResorMax(3)
        endif
!        
        if(LSolveContinuity) then
            i=i+1
            ResidualsPlot(i)=ResorMax(5)
        endif
!        
        if(LSolveEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorMax(6)
        endif
!        
        if(LSolveTurbulenceKineticEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorMax(8)
        endif
!        
        if(LSolveTurbulenceDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorMax(9)
        endif
!        
        if(LSolveTurbulenceSpecificDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorMax(10)
        endif
!        
        if(LSolveTurbulentKL) then
            i=i+1
            ResidualsPlot(i)=ResorMax(12)
        endif
!        
        if(LSolveModifiedED) then
            i=i+1
            ResidualsPlot(i)=ResorMax(11)
        endif
!        
        if(LSolveTurbulenceGammaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorMax(13)
        endif
!        
        if(LSolveTurbulenceReynoldsThetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorMax(14)
        endif
!        
        if(LSolveTurbulencefRelaxationEquation) then
            i=i+1
            ResidualsPlot(i)=ResorMax(17)
        endif
!        
        if(LSolveTurbulenceV2Equation) then
            i=i+1
            ResidualsPlot(i)=ResorMax(15)
        endif
!        
        if(LSolveTurbulenceZetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorMax(16)
        endif
!        
        if(LSolveLambdaELEEquation) then
            i=i+1
            ResidualsPlot(i)=ResorMax(18)
        endif
!        
        do j=1,NumberOfrFieldsToSolve
          if(LSolverField(j)) then
            i=i+1
            ResidualsPlot(i)=ResorMax(18+j)
          endif
        enddo
!
        do j=1,NumberOfScalarsToSolve
          if(LSolveScalar(j)) then
            i=i+1
            ResidualsPlot(i)=ResorMax(18+NumberOfrFieldsToSolve+j)
          endif
        enddo
!
    elseif(NstopType.eq.3) then
!
        i=0
        if(LSolveMomentum) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(1)
            i=i+1
            ResidualsPlot(i)=ResorRMS(2)
            i=i+1
            ResidualsPlot(i)=ResorRMS(3)
        endif
!        
        if(LSolveContinuity) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(5)
        endif
!        
        if(LSolveEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(6)
        endif
!        
        if(LSolveTurbulenceKineticEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(8)
        endif
!        
        if(LSolveTurbulenceDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(9)
        endif
!        
        if(LSolveTurbulenceSpecificDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(10)
        endif
!        
        if(LSolveTurbulentKL) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(12)
        endif
!        
        if(LSolveModifiedED) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(11)
        endif
!        
        if(LSolveTurbulenceGammaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(13)
        endif
!        
        if(LSolveTurbulenceReynoldsThetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(14)
        endif
!        
        if(LSolveTurbulencefRelaxationEquation) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(17)
        endif
!        
        if(LSolveTurbulenceV2Equation) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(15)
        endif
!        
        if(LSolveTurbulenceZetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(16)
        endif
!        
        if(LSolveLambdaELEEquation) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(18)
        endif
!        
        do j=1,NumberOfrFieldsToSolve
          if(LSolverField(j)) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(18+j)
          endif
        enddo
!
        do j=1,NumberOfScalarsToSolve
          if(LSolveScalar(j)) then
            i=i+1
            ResidualsPlot(i)=ResorRMS(18+NumberOfrFieldsToSolve+j)
          endif
        enddo
!
    elseif(NstopType.eq.4) then
!
        i=0
        if(LSolveMomentum) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(1)
            i=i+1
            ResidualsPlot(i)=ResorScaled(2)
            i=i+1
            ResidualsPlot(i)=ResorScaled(3)
        endif
!        
        if(LSolveContinuity) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(5)
        endif
!        
        if(LSolveEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(6)
        endif
!        
        if(LSolveTurbulenceKineticEnergy) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(8)
        endif
!        
        if(LSolveTurbulenceDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(9)
        endif
!        
        if(LSolveTurbulenceSpecificDissipationRate) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(10)
        endif
!        
        if(LSolveTurbulentKL) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(12)
        endif
!        
        if(LSolveModifiedED) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(11)
        endif
!        
        if(LSolveTurbulenceGammaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(13)
        endif
!        
        if(LSolveTurbulenceReynoldsThetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(14)
        endif
!        
        if(LSolveTurbulencefRelaxationEquation) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(17)
        endif
!        
        if(LSolveTurbulenceV2Equation) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(15)
        endif
!        
        if(LSolveTurbulenceZetaEquation) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(16)
        endif
!        
        if(LSolveLambdaELEEquation) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(18)
        endif
!        
        do j=1,NumberOfrFieldsToSolve
          if(LSolverField(j)) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(18+j)
          endif
        enddo
!
        do j=1,NumberOfScalarsToSolve
          if(LSolveScalar(j)) then
            i=i+1
            ResidualsPlot(i)=ResorScaled(18+NumberOfrFieldsToSolve+j)
          endif
        enddo
!        
    endif
!
    write ( 24,1) nIter+nIterOld,(ResidualsPlot(i),i=1,nVariablesPlot)
!--------------------------------------------------------------!
!Should be added directly after writing residuals to the files\
!
        open ( unit = 22, file = "commands2.txt")
        write(22,*) 'a=1'
        close(22)
        open ( unit = 22, file = "commands2.txt")
        write(22,*) 'a=0'
        close(22)
!--------------------------------------------------------------!
        return   
!    call sleep(1/100)
 end subroutine PlotConvergenceHistory
!    
 Subroutine CalculateNumberOfVariablesToPlot
!
        use User0, only: LSolveMomentum,LSolveContinuity,LSolveEnergy,LSolveTurbulenceKineticEnergy,&
                         LSolveTurbulenceDissipationRate,LSolveTurbulenceSpecificDissipationRate,&
                         LSolveTurbulentKL,LSolveModifiedED,LSolveTurbulenceGammaEquation,&
                         LSolveTurbulenceReynoldsThetaEquation,LSolveTurbulencefRelaxationEquation,&
                         LSolveTurbulenceV2Equation,LSolveTurbulenceZetaEquation,LSolveLambdaELEEquation,&
                         LSolverField,LSolveScalar,NumberOfrFieldsToSolve,NumberOfScalarsToSolve
        use Residuals1
!
    implicit none   
!    
        integer :: i
!            
        nVariablesPlot=0
!      
        if(LSolveMomentum) nVariablesPlot=nVariablesPlot+3
        if(LSolveContinuity) nVariablesPlot=nVariablesPlot+1
        if(LSolveEnergy) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulenceKineticEnergy) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulenceDissipationRate) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulenceSpecificDissipationRate) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulentKL) nVariablesPlot=nVariablesPlot+1
        if(LSolveModifiedED) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulenceGammaEquation) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulenceReynoldsThetaEquation) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulencefRelaxationEquation) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulenceV2Equation) nVariablesPlot=nVariablesPlot+1
        if(LSolveTurbulenceZetaEquation) nVariablesPlot=nVariablesPlot+1
        if(LSolveLambdaELEEquation) nVariablesPlot=nVariablesPlot+1
!      
        do i=1,NumberOfrFieldsToSolve
            if(LSolverField(i)) nVariablesPlot=nVariablesPlot+1
        enddo
        do i=1,NumberOfScalarsToSolve
            if(LSolveScalar(i)) nVariablesPlot=nVariablesPlot+1
        enddo
        allocate(ResidualsPlot(nVariablesPlot))
 end Subroutine CalculateNumberOfVariablesToPlot   
!
subroutine CommitTognufile
        use Residuals1, only: nVariablesPlot
        use User0, only: LSolveMomentum,LSolveContinuity,LSolveEnergy,LSolveTurbulenceKineticEnergy,&
                         LSolveTurbulenceDissipationRate,LSolveTurbulenceSpecificDissipationRate,&
                         LSolveTurbulentKL,LSolveModifiedED,LSolveTurbulenceGammaEquation,&
                         LSolveTurbulenceReynoldsThetaEquation,LSolveTurbulencefRelaxationEquation,&
                         LSolveTurbulenceV2Equation,LSolveTurbulenceZetaEquation,LSolveLambdaELEEquation,&
                         LSolverField,LSolveScalar,NumberOfrFieldsToSolve,NumberOfScalarsToSolve,&
                         rFieldName,ScalarName
  implicit none
  integer :: i,j,k
  integer, save :: nIterOld
  integer :: IOstatus
  character*30, dimension(:), allocatable :: titleOld
  double precision, dimension(:), allocatable :: xOld
  character (len = 30) :: nVariablesChar
!
1       format (i7,<nVariablesPlot>e20.12)
5       format(1x,' "iteration"')
10      format(1x,' "x-momentum"')
11      format(1x,' "y-momentum"')
12      format(1x,' "z-momentum"')
20      format(1x,' "continuity"')
30      format(1x,' "energy"')
40      format(1x,' "k"')
50      format(1x,' "epsilon"')
60      format(1x,' "omega"')
70      format(1x,' "KL"')
80      format(1x,' "EddyDiffusivity"')
90      format(1x,' "Gamma"')
100     format(1x,' "ReTheta"')
110     format(1x,' "fRelaxation"')
120     format(1x,' "v2"')
130     format(1x,' "zeta"')
140     format(1x,' "lambda"')
150     format(1x,'"',A10,'" ')
!  
        k=nVariablesPlot+1
        write(nVariablesChar , *) k
        rewind(25)          
        write(25,5,advance='no') 
!
          if(LSolveMomentum) then
            write(25,10,advance='no') 
            write(25,11,advance='no') 
            write(25,12,advance='no') 
          endif
!      
          if(LSolveContinuity) then
            write(25,20,advance='no') 
          endif
!        
          if(LSolveEnergy) then
            write(25,30,advance='no') 
          endif
!        
          if(LSolveTurbulenceKineticEnergy) then
            i=i+1
            write(25,40,advance='no') 
          endif
!        
          if(LSolveTurbulenceDissipationRate) then
            write(25,50,advance='no') 
          endif
!        
          if(LSolveTurbulenceSpecificDissipationRate) then
            write(25,60,advance='no') 
          endif
!        
          if(LSolveTurbulentKL) then
            write(25,70,advance='no') 
          endif
!        
          if(LSolveModifiedED) then
            write(25,80,advance='no') 
          endif
!        
          if(LSolveTurbulenceGammaEquation) then
            write(25,90,advance='no') 
          endif
!        
          if(LSolveTurbulenceReynoldsThetaEquation) then
            write(25,100,advance='no') 
          endif
!        
          if(LSolveTurbulencefRelaxationEquation) then
            write(25,110,advance='no') 
          endif
!        
          if(LSolveTurbulenceV2Equation) then
            write(25,120,advance='no') 
          endif
!        
          if(LSolveTurbulenceZetaEquation) then
            write(25,130,advance='no') 
          endif
!        
          if(LSolveLambdaELEEquation) then
            write(25,140,advance='no') 
          endif
!        
          do j=1,NumberOfrFieldsToSolve
            if(LSolverField(j)) then
              write(25,150,advance='no') rFieldName(j) 
            endif
          enddo
!
          do j=1,NumberOfScalarsToSolve
            if(LSolveScalar(j)) then
              write(25,150,advance='no') ScalarName(j) 
            endif
          enddo
!   
          write(25,*) ""
          allocate(xOld(nVariablesPlot))
          allocate(titleOld(k))
          rewind(24)
          read(24,*,iostat=IOstatus) (titleOld(i),i=1,k)
          do while(IOstatus==0)
            read(24,1,iostat=IOstatus) nIterOld,(xOld(i),i=1,nVariablesPlot)
            if(nIterOld.ne.0.and.IOstatus==0) write(25,1) nIterOld,(xOld(i),i=1,nVariablesPlot)
          enddo
          deallocate(xOld)
          deallocate(titleOld)
          return
end subroutine CommitTognufile