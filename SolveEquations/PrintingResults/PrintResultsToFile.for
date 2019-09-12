c
C#############################################################################################
c
      SUBROUTINE printresults(t1)
c
C#############################################################################################
c
      use User0, only: LsolveMomentum,LsolveContinuity,LsolveEnergy,
     *                 LCompressible,ScalarName,LsolveScalar,
     *                 NumberOfScalarsToSolve,LSolveModifiedED,
     *                 LSolveTurbulenceKineticEnergy,
     *                 LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LSolveTurbulentKL,LTurbulentFlow,EnergyEquation,
     *                 rFieldName,LsolverField,NumberOfrFieldsToSolve,
     *                 LSolveTurbulenceGammaEquation,
     *                 LSolveTurbulenceReynoldsThetaEquation,
     *                 LSolveTurbulencefRelaxationEquation,
     *                 LSolveTurbulenceV2Equation,
     *                 LSolveTurbulenceZetaEquation,
     *                 LSolveLambdaELEEquation,IAbsDivergence,
     *                 IMaxDivergence,IrmsDivergence,FAbsDivergence,
     *                 FMaxDivergence,FrmsDivergence,LUnsteady
      use Geometry1, only: NumberOfElements,NumberOfNodes,x,y,z,
     *                     NumberOfBCSets
      use Geometry3, only: NIFaces,NBFacesTotal,NBFaces
      use Geometry4, only: xc,yc,zc,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      Pressure,temperature,
     *                      MachNumber,BuVelocity,BvVelocity,
     *                      BwVelocity,BPressure,
     *                      Btemperature,BMachNumber,mdot,Bmdot,
     *                      TurbulentKE,TurbulentED,TurbulentOmega,
     *                      BTurbulentKE,BTurbulentED,BTurbulentOmega,
     *                      ModifiedED,BModifiedED,TurbulentKL,
     *                      BTurbulentKL,Htotal,BHtotal,
     *                      TGamma,TReTheta,BTGamma,BTReTheta,
     *                      TurbulentV2,TurbulentZeta,TfRelaxation,
     *                      BTurbulentV2,BTurbulentZeta,BTfRelaxation,
     *                      LambdaELE,BLambdaELE,InitialVelDivergence,
     *                      FinalVelDivergence
      use Scalar1, only: Scalar,BScalar
      use VolumeOfFluid1, only: rField,BrField
      use PhysicalProperties1, only: Density,BDensity,SpecificHeat,
     *                               BSpecificHeat,ReferenceTemperature
      use WallDistance1, only: WallDistance,BWallDistance
      use FlowInOut1, only: LPrintMassFlowRate,MassFlowRate
      use WallStress1
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,jbeg,jend,iskip,k,iScalar,irField
      double precision cpu
      double precision t1,t2,Tref
      data iskip/7/
c********************************************************************************************
c
   7  Format(1x,'z-centroids of faces for boundary set # ',I4)
   8  Format(1x,'y-centroids of faces for boundary set # ',I4)
   9  Format(1x,'x-centroids of faces for boundary set # ',I4)
  10  format(1P7E9.2)
  11  format(1x,'number of nodes= ',I7)
  12  format(1x,'number of elements= ',I7)
  13  format(1x,'number of internal faces= ',I7)
  14  format(1x,'number of boundary faces= ',I7)
  17  format(1x,'CPU time taken = ',F9.2,' Seconds')
  30  format(1x,/)
  40  Format(1x,'u-velocity at face centroids for boundary set # ',I4)
  50  Format(1x,'v-velocity at face centroids for boundary set # ',I4)
  55  Format(1x,'w-velocity at face centroids for boundary set # ',I4)
  51  Format(1x,'vel. magnitude at centroids for boundary set # ',I4)
  60  Format(1x,'Pressure at face centroids for boundary set # ',I4)
  61  Format(1x,'mass flow rate at boundary centroids for set # ',I4)
  65  Format(1x,'Turbulent kinetic energy at face centroids for 
     *boundary set # ',I4)
  41  Format(1x,'v2 at face centroids for boundary set # ',I4)
  42  Format(1x,'Zeta at face centroids for boundary set # ',I4)
  43  Format(1x,'f at face centroids for boundary set # ',I4)
  66  Format(1x,'Turbulent eddy diffusivity at face centroids for 
     *boundary set # ',I4)
  67  Format(1x,'Turbulent specific dissipation rate at face centroids 
     *for boundary set # ',I4)
  81  Format(1x,'Gamma (Turbulence intermittency) at face centroids 
     *for boundary set # ',I4)
  82  Format(1x,'Turbulent Reynolds Theta at face centroids 
     *for boundary set # ',I4)
  68  Format(1x,'Modified eddy diffusivity at face centroids 
     *for boundary set # ',I4)
  69  Format(1x,'Turbulent KL at face centroids for boundary set # ',I4)
  71  Format(1x,'Normal distance to wall at face centroids 
     *for boundary set # ',I4)
  70  Format(1x,'T-static at face centroids for boundary set # ',I4)
  75  Format(1x,'T-total at face centroids for boundary set # ',I4)
  72  Format(1x,'h-static at face centroids for boundary set # ',I4)
  74  Format(1x,'H-total at face centroids for boundary set # ',I4)
  80  Format(1x,'Density at face centroids for boundary set # ',I4)
  90  Format(1x,'Mach Number at face centroids for boundary set # ',I4)
  91  Format(1x,'Lambda at face centroids for boundary set # ',I4)
c********************************************************************************************
c
      Tref=ReferenceTemperature
c
      if(LSolveLambdaELEEquation) then  
        print*,'*******************************************************'
        print*, '*     Initial velocity divergence parameters         *'
        print*,'*******************************************************'
        print*, 'Sum of Initial absolute value=',IAbsDivergence
        print*, 'Maximum Initial absolute value=',IMaxDivergence
        print*, 'rms Initial value=',IrmsDivergence
        print*,'*******************************************************'
        print*, '*     Final velocity divergence parameters           *'
        print*,'*******************************************************'
        print*, 'Sum of Final absolute value=',FAbsDivergence
        print*, 'Maximum Final absolute value=',FMaxDivergence
        print*, 'rms Final value=',FrmsDivergence
        print*,'*******************************************************'
      endif
c      
      print*,''
      print*,'Computations finished: '
      write(11,*),'Computations finished: '
      call timestamp
c
      call cpu_time(t2)
C
      cpu=t2-t1
      print*,'CPU time taken =',cpu,' Seconds'
      write(11,17) ,cpu
c
      write(11,11) NumberOfNodes
      write(11,12) NumberOfElements
      write(11,13) NIFaces
      write(11,14) NBFacesTotal
      write(11,30)
c
      if(LSolveLambdaELEEquation) then
        write(11,*)'***************************************************'
        write(11,*)'*     Initial velocity divergence parameters      *'
        write(11,*)'***************************************************'
        write(11,*) 'Initial velocity divergence parameters'
        write(11,*) 'Sum of Initial absolute value=',IAbsDivergence
        write(11,*) 'Maximum Initial absolute value=',IMaxDivergence
        write(11,*) 'rms Initial value=',IrmsDivergence
        write(11,*)'***************************************************'
        write(11,*)'*     Final velocity divergence parameters        *'
        write(11,*)'***************************************************'
        write(11,*) 'Sum of Final absolute value=',FAbsDivergence
        write(11,*) 'Maximum Final absolute value=',FMaxDivergence
        write(11,*) 'rms Final value=',FrmsDivergence
        write(11,*)'***************************************************'
      endif
c      
      write(11,30)
      write(11,*) 'x-values of nodes in sequence'
      write(11,30)
c
      do i=1,NumberOfNodes,iskip
c      
        jbeg=i
        jend=min(jbeg+iskip-1,NumberOfNodes)
        write(11,10) (x(j),j=jbeg,jend)
c
      enddo
c
      write(11,30)
      write(11,*) 'y-values of nodes in sequence'
      write(11,30)
c
      do i=1,NumberOfNodes,iskip
c      
        jbeg=i
        jend=min(jbeg+iskip-1,NumberOfNodes)
        write(11,10) (y(j),j=jbeg,jend)
c
      enddo
c
      write(11,30)
      write(11,*) 'z-values of nodes in sequence'
      write(11,30)
c
      do i=1,NumberOfNodes,iskip
c      
        jbeg=i
        jend=min(jbeg+iskip-1,NumberOfNodes)
        write(11,10) (z(j),j=jbeg,jend)
c
      enddo
c
      write(11,30)
      write(11,*) 'x-centroids of elements in sequence'
      write(11,30)
c
      do i=1,NumberOfElements,iskip
c      
        jbeg=i
        jend=min(jbeg+iskip-1,NumberOfElements)
        write(11,10) (xc(j),j=jbeg,jend)
c
      enddo
c
      write(11,30)
      write(11,*) 'x-centroids of boundary faces in sequence'
      write(11,30)
c
      do i=1,NumberOfBCSets
c
        write(11,9) i
c        
        do j=1,NBFaces(i),iskip      
c      
          jbeg=j
          jend=min(jbeg+iskip-1,NBFaces(i))
          write(11,10) (BFaceCentroidx(i,k),k=jbeg,jend)
c
        enddo
c
      enddo
c
      write(11,30)
      write(11,*) 'y-centroids of elements in sequence'
      write(11,30)
c
      do i=1,NumberOfElements,iskip
c      
        jbeg=i
        jend=min(jbeg+iskip-1,NumberOfElements)
        write(11,10) (yc(j),j=jbeg,jend)
c
      enddo
c
      do i=1,NumberOfBCSets
c
        write(11,8) i
c        
        do j=1,NBFaces(i),iskip      
c      
          jbeg=j
          jend=min(jbeg+iskip-1,NBFaces(i))
          write(11,10) (BFaceCentroidy(i,k),k=jbeg,jend)
c
        enddo
c
      enddo
c
      write(11,30)
      write(11,*) 'z-centroids of elements in sequence'
      write(11,30)
c
      do i=1,NumberOfElements,iskip
c      
        jbeg=i
        jend=min(jbeg+iskip-1,NumberOfElements)
        write(11,10) (zc(j),j=jbeg,jend)
c
      enddo
c
      do i=1,NumberOfBCSets
c
        write(11,7) i
c        
        do j=1,NBFaces(i),iskip      
c      
          jbeg=j
          jend=min(jbeg+iskip-1,NBFaces(i))
          write(11,10) (BFaceCentroidz(i,k),k=jbeg,jend)
c
        enddo
c
      enddo
c
      if(LsolveMomentum.or.LSolveLambdaELEEquation) then
c
        write(11,30)
        write(11,*) 'u-Velocity at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (uVelocity(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,40) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BuVelocity(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
        write(11,30)
        write(11,*) 'v-Velocity at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (vVelocity(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,50) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BvVelocity(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
        write(11,30)
        write(11,*) 'w-Velocity at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (wVelocity(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,55) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BwVelocity(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
        write(11,30)
        write(11,*)'Velocity Magnitude at element centroids in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) 
     *        (dsqrt(uVelocity(j)**2+vVelocity(j)**2+
     *                          wVelocity(j)**2),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,51) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) 
     *        (dsqrt(BuVelocity(i,k)**2+BvVelocity(i,k)**2+
     *                              BwVelocity(i,k)**2),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LsolveContinuity) then
c
        write(11,30)
        write(11,*) 'Pressure at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (Pressure(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,60) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BPressure(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
        write(11,30)
        write(11,*) 'Mass flow rate at centroids of internal faces' 
        write(11,30)
c
        do i=1,NIFaces,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NIFaces)
          write(11,10) (mdot(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,61) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (Bmdot(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulenceKineticEnergy) then
c
        write(11,30)
        write(11,*) 'Turbulent kinetic energy at centroids of elements 
     *in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TurbulentKE(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,65) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTurbulentKE(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulenceDissipationRate) then
c
        write(11,30)
        write(11,*) 'Turbulent dissipation rate at centroids of elements
     * in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TurbulentED(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,66) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTurbulentED(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulenceV2Equation) then
c
        write(11,30)
        write(11,*) 'v2 values at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TurbulentV2(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,41) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTurbulentV2(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulenceZetaEquation) then
c
        write(11,30)
        write(11,*) 'Zeta values at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TurbulentZeta(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,42) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTurbulentZeta(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulencefRelaxationEquation) then
c
        write(11,30)
        write(11,*) 'f values at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TfRelaxation(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,43) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTfRelaxation(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulenceSpecificDissipationRate) then
c
        write(11,30)
        write(11,*) 'Turbulent specific dissipation rate at centroids
     * of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TurbulentOmega(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,67) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTurbulentOmega(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulenceGammaEquation) then
c
        write(11,30)
        write(11,*) 'Turbulent Gamma at centroids
     * of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TGamma(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,81) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTGamma(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulenceReynoldsThetaEquation) then
c
        write(11,30)
        write(11,*) 'Turbulent Reynolds Theta at centroids
     * of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TReTheta(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,82) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTReTheta(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveTurbulentKL) then
c
        write(11,30)
        write(11,*) 'Turbulent KL at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (TurbulentKL(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,69) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTurbulentKL(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveModifiedED) then
c
        write(11,30)
        write(11,*) 'Modified eddy diffusivity at centroids
     * of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (ModifiedED(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,68) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BModifiedED(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LSolveLambdaELEEquation) then
c
        write(11,30)
        write(11,*) 'Lambda at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (LambdaELE(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,91) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BLambdaELE(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
        write(11,30)
        write(11,*) 'Initial velocity divergence at centroids
     * of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (InitialVelDivergence(j),j=jbeg,jend)
c
        enddo
c
        write(11,30)
        write(11,*) 'Final velocity divergence at centroids
     * of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (FinalVelDivergence(j),j=jbeg,jend)
c
        enddo
c
      endif
c
      if(LTurbulentFlow) then
c
        write(11,30)
        write(11,*) 'Normal distance to wall at centroids
     * of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (WallDistance(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,71) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BWallDistance(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LsolveEnergy) then
c
        write(11,30)
        write(11,*) 'T-static at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (Temperature(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,70) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTemperature(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LsolveEnergy) then
c
        write(11,30)
        write(11,*) 'T-total at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (Temperature(j)+0.5*(uVelocity(j)**2+
     *                        vVelocity(j)**2+wVelocity(j)**2)/
     *                                  SpecificHeat(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,75) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BTemperature(i,k)+0.5*(BuVelocity(i,k)**2+
     *                          BvVelocity(i,k)**2+BwVelocity(i,k)**2)/
     *                                   BSpecificHeat(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LsolveEnergy) then
c
        write(11,30)
        write(11,*) 'h-static at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (SpecificHeat(j)*(Temperature(j)-
     *                            Tref),j=jbeg,jend)
     
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,72) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BSpecificHeat(i,k)*(BTemperature(i,k)-
     *                               Tref),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') then
c
        write(11,30)
        write(11,*) 'H-total at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (Htotal(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,74) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
            jbeg=j
            jend=min(jbeg+iskip-1,NBFaces(i))
            write(11,10) (BHtotal(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        do irField=1,NumberOfrFieldsToSolve
c
          if(LSolverField(irField)) then
c
            write(11,30)
            write(11,*) rFieldName(irField),
     *          ' at centroids of elements in sequence'
            write(11,30)
c
            do i=1,NumberOfElements,iskip
c      
              jbeg=i
              jend=min(jbeg+iskip-1,NumberOfElements)
              write(11,10) (rField(j,irField),j=jbeg,jend)
c
            enddo
c
            do i=1,NumberOfBCSets
c  70  Format(1x,'T-static at face centroids for boundary set # ',I4)
              write(11,30)
              write(11,*) rFieldName(irField),
     *          'at face centroids for boundary set # ',i
              write(11,30)
c        
              do j=1,NBFaces(i),iskip      
c      
                jbeg=j
                jend=min(jbeg+iskip-1,NBFaces(i))
                write(11,10) (BrField(i,k,irField),k=jbeg,jend)
c
              enddo
c
            enddo
c
          endif
c
        enddo
c
      endif
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        do iScalar=1,NumberOfScalarsToSolve
c
          if(LSolveScalar(iScalar)) then
c
            write(11,30)
            write(11,*) ScalarName(iScalar),
     *          ' at centroids of elements in sequence'
            write(11,30)
c
            do i=1,NumberOfElements,iskip
c      
              jbeg=i
              jend=min(jbeg+iskip-1,NumberOfElements)
              write(11,10) (Scalar(j,iScalar),j=jbeg,jend)
c
            enddo
c
            do i=1,NumberOfBCSets
c
              write(11,30)
              write(11,*) ScalarName(iScalar),
     *          'at face centroids for boundary set # ',i
              write(11,30)
c        
              do j=1,NBFaces(i),iskip      
c      
                jbeg=j
                jend=min(jbeg+iskip-1,NBFaces(i))
                write(11,10) (BScalar(i,k,iScalar),k=jbeg,jend)
c
              enddo
c
            enddo
c
          endif
c
        enddo
c
      endif
c
      if(LCompressible) then
c
        write(11,30)
        write(11,*) 'Density at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (Density(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,80) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
          jbeg=j
          jend=min(jbeg+iskip-1,NBFaces(i))
          write(11,10) (BDensity(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
        write(11,30)
        write(11,*) 'Mach number at centroids of elements in sequence'
        write(11,30)
c
        do i=1,NumberOfElements,iskip
c      
          jbeg=i
          jend=min(jbeg+iskip-1,NumberOfElements)
          write(11,10) (MachNumber(j),j=jbeg,jend)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          write(11,30)
          write(11,90) i
          write(11,30)
c        
          do j=1,NBFaces(i),iskip      
c      
          jbeg=j
          jend=min(jbeg+iskip-1,NBFaces(i))
          write(11,10) (BMachNumber(i,k),k=jbeg,jend)
c
          enddo
c
        enddo
c
      endif
c
      do i=1,NumberOfBCSets
        If(LPrintMassFlowRate(i)) then
c              
         write(11,*) 'Mass flow rate through boundary ',i,' = ',
     *       MassFlowRate(i),'Kg/s'
c
        endif
      enddo
c
      if(.not.LUnsteady) then
        call PrintMaximumTauWall
        do i=1,NumberOfBCSets
          write(11,*) 
     *     'x,y,z, and shear stress value along wall ',i,' = ',
     *      'x=',xLocationMaxWallShearStress(i),
     *      'y=',yLocationMaxWallShearStress(i),
     *      'z=',zLocationMaxWallShearStress(i),
     *      'Maximu Stress=',MaximumShearStress(i)
        enddo
      endif
c
      call CalculateFluxes
c
      return
      end