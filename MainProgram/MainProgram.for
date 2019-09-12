C#############################################################################################
c
c    A collocated three-dimensional, Finite Volume-based, flow solver on unstructured grids. The
c    code deals with both compressible and incompressible flow and transfer phenomena problems. 
c    Convection schemes are implemented following the TVD and NVF formulations. The algebraic
c    system of equations may be solved using the SOR, ILU, PBCG (Preconditioned biconjugate 
c    gradient method), or a direct solver. In addition algebraic systems may be solved using   
c    an algebraic multigrid solver with any of the above iterative methods as a smoother. 
c    A variety of boundary conditions are implemented for both incompressible and compressible 
c    flows. Both steady and unsteady flow and transfer phenomena flow problems can be handeled
c    using a first or a second order transient schemes with fixed or adaptive time steps. 
c    The SIMPLE and SIMPLEC algorithms are implemented. High Reynolds and low Reynolds number
c    turbulence models are implemented. Incompressible free surface flows based on the volume 
c    of fluid method.
c
C#############################################################################################
c
      PROGRAM Unstructured
c
C#############################################################################################
c
      use User0
      use Geometry1, only: NumberOfElements
      use Variables1
      use MultiGrid2, only: nIter,nIterMG,nIterTotal
      use PhysicalProperties1
      use Residuals1
      use Transient1
      use Geometry4, only:xc,yc,zc
      use ReferencePressure1, only: LSetReferencePressure
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*14, save :: TransientSchemeTemp,TransientSchemerFieldT
      integer :: nTemp,i,index
      integer, save :: iScalar,irField
      double precision source
      double precision t1
c
c********************************************************************************************
c
      print *, ' Enter Problem Name (up to seventy Characters):  '
      read(*,'(A70)') name
c
      call OpenFiles
      print*,'Setting variables to solve for...Please wait'
      call SetVariablesToSolve
c
      call OpenInternalFiles
c
      print*,'Computations started: '
c
      write(11,*),'Computations started: '
      call timestamp
      call cpu_time(t1)
c
      call Grid
      print*,'Allocate storage...Please wait' 
      call allocateBCStorage
      call allocateStorage
      print*,'Reading initial and boundary conditions...Please wait'
      call SetGlobalVariables
      call setBoundaryConditions
      call SetBCPatches
      call RotationMatrixCoefficients
      call allocateCentrifugalForce
      print*,'checking for point sources...Please wait'
      call LocatePointSources
      call SetReferencePressure
      if(LSetReferencePressure) call LocateReferencePressureElement
      if(LfixPressure) call LocateFixedPressureElement
      call LocateMonitoringElement
      if(LTurbulentFlow) then 
        print*,'Setting turbulence model coefficients...Please wait'
        call SetTurbulenceModelCoefficients
        call CalculateNormalDistanceToWall
        call CheckPositivityOfNormalDistanceToWall
      endif
c
      print*,'......................................................'
      print*,'..............Starting computations...................'
      print*,'......................................................'
c
c-----------------------------------------------------------------------------------
c--- Make sure boundary conditions and properties are assigned the proper values 
c-----------------------------------------------------------------------------------
      nIter=0
      nIterTotal=0
      call MakeChangesDuringComputations
      call updateBoundaryConditions
      call updateDensity
      call updateViscosity
      call updateConductivity
      call updateSpecificHeat
      call updateSources
      call updateScalars
      call SetTimeStepValue
c
c--- Calculate the far field velocities from given data
c
      call calculateFarFieldVariables
c
      call CalculatePhysicalProperties
c
c-----------------------------------------------------------------------------------
c
      if(LSolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                call InitializeHTotalFromTemperature
c
      if(LsolveMomentum) then
c
        if(Lcompressible) then
c
          call CalculateMachNumber  
          call Calculatedrhodp  
          call CalculateDensity  
          If(LsolveEnergy) call UpdateBoundaryTemperature
          call CalculateFaceDensity  
c
        endif
c
        call saveDensityForFalseTransient
        call UpdateBoundaryPressure
c
c---- Update values along symmetry boundaries
c
        call UpdateSymmetry
c
c---- Update velocity field along slip walls
c
        call UpdateSlipWallBoundary
c
      endif
c
      if(LSolveLambdaELEEquation) then
        call CalculateLambdaDiffusionTensor
        call IntializeLambdaVelocityField
        call CalculateLambdaELESources
      endif
c
      if(.not.Lcompressible) call CalculateFaceDensity
      if(LsolveMomentum.or.LConvectScalar) call CalculateInitialmdot
      if(LsolveMomentum.and.LBuoyancy) call reDistributeTheBuoyancyTerm
c
c--- check if transient
c
      if(.not.LUnsteady) then
c
        timemax=0.
c
      endif
c
      time=0.
      ndt=0
c
      index=1
c
      nIterMG=0
      TransientSchemeTemp=TransientScheme
      TransientSchemerFieldT=TransientSchemerField
c
      if(LUnsteady) then
c
        dtOld=dt
        dtOldOld=dtOld
        call InitializeOldVariables
c
      endif
c
      if(LReadOldSolution) then
c
        print*,'Reading saved solution...Please wait'
        call readSolution
        close(4)
        print*,'Finished reading saved solution..........'
c
      endif
c
      call WriteInputDataToFile
c
      do while(time.le.timemax.and.index.eq.1)
c
        index=0
        if(LUnsteady) then
c
          index=1
c
          if(TransientScheme.eq.'cranknicolson1') then
c
            if(ndt.ne.0) call CrankNicolson
c
          endif
c
          if(TransientSchemerField.eq.'cranknicolson1') then
c
            if(ndt.ne.0) call CrankNicolsonrField
c
          endif
c
          if(ndt.eq.0) then
c
            If(TransientSchemeTemp.eq.'cranknicolson2'.or.
     *            TransientSchemeTemp.eq.'adamsmoulton') then
c
              TransientScheme='euler'
c
            endif
c
            If(TransientSchemerFieldT.eq.'cranknicolson2'.or.
     *            TransientSchemerFieldT.eq.'adamsmoulton'.or.
     *            TransientSchemerFieldT.eq.'tics1.75'.or.
     *            TransientSchemerFieldT.eq.'tics2.5') then
c
              TransientSchemerField='euler'
c
            endif
c
          else
c
            TransientScheme=TransientSchemeTemp   
            TransientSchemerField=TransientSchemerFieldT   
c          
          endif
c
          call updateVariablesInTime
c
          dtOldOld=dtOld
          dtOld=dt
c
          call PrintMaximumTauWall
          call PrintBoundaryMassFlowRatesInTime
          call PrintVariablesInTime
c
          if(LadaptiveTimeStep.and.ndt.ne.0) then
c
            call AdaptiveTimeStep
c
          else
c
            call SetTimeStepValue
c
          endif
c
          time=time+dt
          ndt=ndt+1
          call UpdateUnsteadyBoundaryConditions
          call UpdateArteryResistance
c 
          write(*,*) 'solution at time= ', time
          write(13,*) 'solution at time= ', time
c 
        endif
c
c---  Start outer iterations
c
        do nIter=1,IterMax
c
          nIterMG=nIterMG+1
          nIterTotal=nIterTotal+1
c
          call saveDensityForFalseTransient
          call CalculatePhysicalProperties
c
          if(LsolveMomentum) then
c
            call SolveMomentumX
            call SolveMomentumY
            call SolveMomentumZ
c
          endif
c
          if(LsolveContinuity) call SolveContinuity
c
          if(LTurbulentFlow) then
c
            call CalculateEffectsOfTurbulence
c
          endif
c
          if(LsolveEnergy) call SolveEnergy 
c
          do iScalar=1,NumberOfScalarsToSolve
c
            if(LsolveScalar(iScalar)) call SolveScalar(iScalar)
c
          enddo
c
          if(LFreeSurfaceFlow) then
c
            do irField=1,NumberOfrFieldsToSolve
c
              if(LsolverField(irField)) Call SolveVolumeOfFluid(irField)
c
            enddo
c
            call CalculateLastrField
c
          endif
c
          if(LSolveLambdaELEEquation) then
            call SolveLambdaEulerLagrangeEquation
          endif
c
          if(LsolveContinuity) then
c
            print*,'Overall mass conservation = ', OverallImbalance
            print*,'Location of maximum mass imbalance = ',iMaxImbalance
            write(*,14) xc(iMaxImbalance),
     *                     yc(iMaxImbalance),zc(iMaxImbalance)
 14         format(2x,'xLocation = '1P1E10.3,2x,
     *        'yLocation = '1P1E10.3,2x,'zLocation = '1P1E10.3)
c
            write(13,*) "Overall mass conservation = ", OverallImbalance
            write(13,*) 
     *          "Location of maximum mass imbalance = ",iMaxImbalance
            write(13,14) xc(iMaxImbalance),
     *                     yc(iMaxImbalance),zc(iMaxImbalance)
c
          endif
c
          call PrintConvergenceHistory(nIter,nIterTotal,source)
c          call PlotConvergenceHistory
c
          if(source.lt.maximumResidual) then
c
            exit
c
          else
c
            if(LsolveContinuity) then 
              call CalculateBoundaryMassFlowRates
              call PrintBoundaryMassFlowRates
              call CalculatePressureAtOutlets
              if(LUnsteady.and..not.LCompressible) call WindKesselModel
            endif
            call MakeChangesDuringComputations
            call updateBoundaryConditions
            call updateDensity
            call updateViscosity
            call updateConductivity
            call updateSpecificHeat
            call updateSources
            call updateScalars
c
          endif
c      
        enddo
c
        if(LUnsteady) then
c
          if(MeshType.eq.'neutral') then
c
              print*,time,timemax
            if(mod(ndt,nTimePrint).eq.0.or.time.ge.timemax) then
c
              call GetUnitName(ndt,'.dat',Filprint,12) 
c
              call tecplot
c
              close(12)
c
            endif
c
		elseif(MeshType.eq.'polymesh') then
c
            if(mod(ndt,nTimePrint).eq.0.or.time.ge.timemax) then
c
              call GetUnitName(ndt,'.vtu',Filprint,12) 
c
              call WriteParaviewFile
c
              close(12)
c
            endif
c
		endif
c
          if(mod(ndt,nTimeSave).eq.0) then
c
            call GetUnitName(ndt,'.old',Filold,4) 
c
            call writeSolution
c
            close(4)
c
          endif
c
        else
c 
          if(LSolveLambdaELEEquation) call UpdateEulerLagrangeVelocity
c
          if(MeshType.eq.'neutral') then
            call tecplot
            close(12)
          elseif(MeshType.eq.'polymesh') then
            if(LprintParaviewFile) then
              call WriteParaviewFile
              close(12)
            endif
          endif
c
          call GetUnitName(ndt,'.old',Filold,4) 
          call writeSolution
          close(4)
c
        endif
c
      enddo
c
      call printresults(t1)
c
      pause
      end
c
C#############################################################################################
      SUBROUTINE OpenInternalFiles
C#############################################################################################
c
      use User0, only: MeshType,name,FilNeu,Filprint,Filin,Filres,
     *                 Filout,Filwdist,Filold,FilMG,LTestMultiGrid,
     *                 NeutralMeshDirectory,SolutionDirectory,
     *                 LReadOldSolution,LUnsteady
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      if(MeshType.eq.'neutral') then
        FilNeu=trim(name)//'.neu'
        Filprint=trim(name)//'.dat'
      elseif(MeshType.eq.'polymesh') then
        Filprint=trim(name)//'.vtu'
      endif
c
      Filin=trim(name)//'.cin'
      Filres=trim(name)//'.res'
      Filout=trim(name)//'.out'
      Filwdist=trim(name)//'.norm'
      Filold=trim(name)//'.old'
C
      if(MeshType.eq.'neutral') 
     *  open (unit=1,file=trim(NeutralMeshDirectory)//'/'//trim(FilNeu))
      if(LReadOldSolution)
     *  open (unit=4,file=trim(SolutionDirectory)//'/'//trim(Filold),
     *                                               form='unformatted')
      open (unit=10,file=trim(SolutionDirectory)//'/'//trim(Filin))
      open (unit=11,file=trim(SolutionDirectory)//'/'//trim(Filout))
      open (unit=12,file=trim(SolutionDirectory)//'/'//trim(Filprint))
      open (unit=13,file=trim(SolutionDirectory)//'/'//trim(Filres))
      open (unit=15,file=trim(SolutionDirectory)//'/'//trim(Filwdist))
c
      if(MeshType.eq.'neutral') rewind 1
      if(LReadOldSolution) rewind 4
      rewind 10
      rewind 11
      rewind 12
      rewind 13
      rewind 15
c
      if(LTestMultiGrid) then
        FilMG=trim(name)//'.MG'
        open (unit=14,file=trim(SolutionDirectory)//'/'//trim(FilMG))
        rewind 14
      endif
c
      if(LUnsteady) then
        open(unit=16,status='unknown',
     *    file=trim(SolutionDirectory)//"/MassFlowRate",
     *                access='sequential',form='formatted')
        open(unit=18,status='unknown',
     *    file=trim(SolutionDirectory)//"/Shear",
     *                access='sequential',form='formatted')
      endif
c      
      return
      end
c
C#############################################################################################
c
      SUBROUTINE GetUnitName(ndt,extension,FileTemp,i)
c
C#############################################################################################
      use User0, only: name,SolutionDirectory
C#############################################################################################
      implicit none
C#############################################################################################
      integer :: ndt,i
      character*80 FileTemp
      character*4 extension
      character(len=8) :: fmt ! format descriptor
      character(len=8) :: x1
C#############################################################################################
c
      fmt = '(I5.5)' ! an integer of width 5 with zeros at the left
      write (x1,fmt) ndt ! converting integer to string using a 'internal file'
      FileTemp=trim(name)//trim(x1)//trim(extension)
      if(i.eq.4) then
        if(ndt.eq.0) FileTemp=trim(name)//trim(extension)
        open (unit=i,file=trim(SolutionDirectory)//'/'//trim(FileTemp),
     *       form='unformatted')
      else
        open (unit=i,file=trim(SolutionDirectory)//'/'//trim(FileTemp))
      endif
      rewind i
c
      return
      end