
module WriteFoamModule
contains
!************************************************************************************************
     SUBROUTINE PrintFoam
!************************************************************************************************
    use user0, only:  time
    use user0, only:  FoamCasedirectory
    
     use User0, only: LSolveMomentum,LSolveContinuity,LSolveTurbulenceKineticEnergy,&
                     LSolveModifiedED,LSolveTurbulenceDissipationRate,LSolveTurbulenceV2Equation,&
                     LSolveTurbulenceZetaEquation,LSolveTurbulencefRelaxationEquation,&
                     LSolveTurbulenceSpecificDissipationRate,LSolveTurbulenceGammaEquation,&              
                     LSolveTurbulenceReynoldsThetaEquation,LSolveTurbulentKL,LSolveEnergy,&
                     TurbulenceModel,LTurbulentFlow,LBuoyancy,LtestTurbulenceModel,EnergyEquation,&
                     NumberOfrFieldsToSolve,LSolverField,NumberOfScalarsToSolve,LSolveScalar,&
                     LCompressible,LFreeSurfaceFlow,rFieldName,ScalarName,LUnsteady,LPrintGradients,&
                     LSolveLambdaELEEquation
     use Geometry1, only: NumberOfElements, NumberOfBCSets
     use Geometry3, only: NBFaces,NBFacesmax,NBFaceOwner
     !----------------------------------------------------------------------------------------------------!
     use Variables1, only: uVelocity,vVelocity,WVelocity,Pressure,TurbulentKE,ModifiedED,TurbulentED,&
                          TurbulentV2,TurbulentZeta,TfRelaxation,TurbulentOmega,TGamma,TGammaEff,&
                          TReTheta,TurbulentKL,TurbulenceProduction,TurbulenceProductionB,Temperature,&
                          Htotal,MachNumber,uVelGradx,uVelGrady,uVelGradz,vVelGradx,vVelGrady,vVelGradz,&
                          wVelGradx,wVelGrady,wVelGradz,PressGradx,PressGrady,PressGradz,&
                          TKEGradx,TKEGrady,TKEGradz,TEDGradx,TEDGrady,TEDGradz,TurbulentV2Gradx,&
                          TurbulentV2Grady,TurbulentV2Gradz,TurbulentZetaGradx,TurbulentZetaGrady,&
                          TurbulentZetaGradz,TfRelaxationGradx,TfRelaxationGrady,TfRelaxationGradz,&
                          TOmegaGradx,TOmegaGrady,TOmegaGradz,TGammaGradx,TGammaGrady,TGammaGradz,&
                          TReThetaGradx,TReThetaGrady,TReThetaGradz,TurbulentKLGradx,TurbulentKLGrady,&
                          TurbulentKLGradz,ModifiedEDGradx,ModifiedEDGrady,ModifiedEDGradz,TempGradx,TempGrady,&
                          TempGradz,HtotalGradx,HtotalGrady,HtotalGradz,LambdaELE,LambdaELEGradx,&
                          LambdaELEGrady,LambdaELEGradz,InitialVelDivergence,FinalVelDivergence, &
                          BInitialVelDivergence,BFinalVelDivergence,BTGammaEff,&
                          BuVelocity,BvVelocity,BWVelocity,BPressure,BTurbulentKE,BModifiedED,BTurbulentED,&
                          BTurbulentV2,BTurbulentZeta,BTfRelaxation,BTurbulentOmega,BTGamma,&
                          BTReTheta,BTurbulentKL,BTurbulenceProduction,BTurbulenceProductionB,BTemperature,&
                          BHtotal,BMachNumber,BuVelGradx,BuVelGrady,BuVelGradz,BvVelGradx,BvVelGrady,BvVelGradz,&
                          BwVelGradx,BwVelGrady,BwVelGradz,BPressGradx,BPressGrady,BPressGradz,&
                          BTKEGradx,BTKEGrady,BTKEGradz,BTEDGradx,BTEDGrady,BTEDGradz,BTurbulentV2Gradx,&
                          BTurbulentV2Grady,BTurbulentV2Gradz,BTurbulentZetaGradx,BTurbulentZetaGrady,&
                          BTurbulentZetaGradz,BTfRelaxationGradx,BTfRelaxationGrady,BTfRelaxationGradz,&
                          BTOmegaGradx,BTOmegaGrady,BTOmegaGradz,BTGammaGradx,BTGammaGrady,BTGammaGradz,&
                          BTReThetaGradx,BTReThetaGrady,BTReThetaGradz,BTurbulentKLGradx,BTurbulentKLGrady,&
                          BTurbulentKLGradz,BModifiedEDGradx,BModifiedEDGrady,BModifiedEDGradz,BTempGradx,BTempGrady,&
                          BTempGradz,BHtotalGradx,BHtotalGrady,BHtotalGradz,BLambdaELE,BLambdaELEGradx,&
                          BLambdaELEGrady,BLambdaELEGradz
                          
     use PhysicalProperties1, only: Density,Viscosity,TurbulentViscosity,SpecificHeat,ReferenceTemperature,&
                                   eDiffCoefficient, &
                                   
                                   BDensity,BViscosity,BTurbulentViscosity,BSpecificHeat,&
                                   BeDiffCoefficient
                                   
     use Turbulence1, only: TurbulentViscosityTl,TurbulentViscosityTs,ProductionKT,ProductionKL,ReT,&
                            StrainRate,Vorticity,BProductionKT,BProductionKL,  &
                            BTurbulentViscosityTl,BTurbulentViscosityTs,BReT,&
                            BStrainRate,BVorticity
     !----------------------------------------------------------------------------------------------------!
                            
     use ReferenceValues1, only: UInfinity
     use constants1, only: twothird,tiny
     use WallDistance1, only: WallDistance, BWallDistance
     use VolumeOfFluid1, only: rField,rFieldGradx,rFieldGrady,rFieldGradz, &
                               brField,brFieldGradx,brFieldGrady,brFieldGradz
     use Scalar1, only: Scalar,ScalarGradx,ScalarGrady,ScalarGradz, &
                        bScalar,bScalarGradx,bScalarGrady,bScalarGradz
     use Tecplot1, only: yplusPlot,uplusPlot,ByplusPlot,BuplusPlot
!************************************************************************************************
     implicit none 
!************************************************************************************************
     character*10 :: Variable1
     integer :: iScalar,irField,i,j,k
     integer :: VTK_unit
     character*10 :: part1
     character*5 :: part2,part3,part4
     character*15 :: name1,name2,name3
     integer ::   namesize
     double precision, allocatable :: Momentum(:,:), BMomentum(:,:,:)
     double precision, allocatable :: ScalarToPrint(:), BScalarToPrint(:,:)
     character(len=15) :: timeChar
     
1    format(1x,',"',A15,'" ')
!********************************************************************************************
     
     write(timeChar , '(f10.5)') time
     
     call system('mkdir -p '//trim(FoamCasedirectory)//'\'//adjustl(trim(timeChar)))
!
     if(LUnsteady) LtestTurbulenceModel=.false.   !in order not to modify normal distance to wall
!    
     allocate(ScalarToPrint(NumberOfElements))
     allocate(BScalarToPrint(NumberOfBCSets,NBFacesmax))
!
     if(LsolveMomentum.or.LSolveLambdaELEEquation) then    
       allocate(Momentum(NumberOfElements,3))
       allocate(BMomentum(NumberOfBCSets,NBFacesmax,3))
       
       Momentum(:,1)=uVelocity(:)
       Momentum(:,2)=vVelocity(:)
       Momentum(:,3)=WVelocity(:)
       
	   BMomentum(:,:,1) = BuVelocity(:,:) 
	   BMomentum(:,:,2) = BvVelocity(:,:) 
	   BMomentum(:,:,3) = BwVelocity(:,:) 
       
       call printFoamVector(Momentum,BMomentum , 'Velocity', timeChar)
       deallocate(Momentum)
       deallocate(BMomentum)
     endif
     
     if(LSolveLambdaELEEquation) then    
       call printFoamScalar(LambdaELE,BLambdaELE,'Lambda-Eulerian-Lagrangian-Equation', timeChar)
       call printFoamScalar(InitialVelDivergence,BInitialVelDivergence,'Initial-velocity-divergence', timeChar)
       call printFoamScalar(FinalVelDivergence,BFinalVelDivergence,'Final-velocity-divergence', timeChar)
     endif
     if(LSolveContinuity) then
       call printFoamScalar(Pressure,BPressure,'Pressure', timeChar)
     endif
     if(LSolveTurbulenceKineticEnergy) then
       call printFoamScalar(TurbulentKE,BTurbulentKE,'Turbulent-kinetic-energy', timeChar)
     endif
     if(LSolveModifiedED) then
       call printFoamScalar(ModifiedED,BModifiedED,'Modified-eddy-diffusivity', timeChar)
     endif
     if(LSolveTurbulenceDissipationRate) then
       call printFoamScalar(TurbulentED,BTurbulentED,'Turbulent-dissipation-rate', timeChar)
     endif
     if(LSolveTurbulenceV2Equation) then
       call printFoamScalar(TurbulentV2,BTurbulentV2,'Turbulent-v2', timeChar)
     endif
     if(LSolveTurbulenceZetaEquation) then
       call printFoamScalar(TurbulentZeta,BTurbulentZeta,'Turbulent-zeta', timeChar)
     endif
     if(LSolveTurbulencefRelaxationEquation) then
       call printFoamScalar(TfRelaxation,BTfRelaxation,'f-relaxation', timeChar)
     endif
     if(LSolveTurbulenceSpecificDissipationRate) then
       call printFoamScalar(TurbulentOmega,BTurbulentOmega,'Turbulent-specific-dissipation-rate', timeChar)
     endif
     if(LSolveTurbulenceGammaEquation) then
       call printFoamScalar(TGamma,BTGamma,'Intermittency(Turbulent-gama)', timeChar)
     endif
     if(TurbulenceModel.eq.'sstgamaretheta') then
       do i=1,NumberOfBCSets
         do j=1,NBFaces(i)
           k=NBFaceOwner(i,j)
           BTGammaEff(i,j)=TGammaEff(k)
         enddo
       enddo
       call printFoamScalar(TGammaEff,BTGammaEff,'Intermittency-effective', timeChar)
     endif
     if(LSolveTurbulenceReynoldsThetaEquation) then
       call printFoamScalar(TReTheta,BTReTheta,'Momentum-thickness', timeChar)
     endif
     if(LSolveTurbulentKL) then
       call printFoamScalar(TurbulentKL,BTurbulentKL,'Laminar-kinetic-energy', timeChar)
     endif
     if(LTurbulentFlow.and.TurbulenceModel.eq.'kklomega') then
       ScalarToPrint=TurbulentKE+TurbulentKL
       BScalarToPrint = BTurbulentKE+BTurbulentKL
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Total-fluctuation-energy', timeChar)  
     endif
     if(LSolveMomentum) then
       call printFoamScalar(Viscosity,BViscosity,'Laminar-viscosity', timeChar)
     endif
     if(LTurbulentFlow) then
       call printFoamScalar(TurbulentViscosity,BTurbulentViscosity,'Turbulent-viscosity', timeChar)
     endif
     if(LTurbulentFlow.and.TurbulenceModel.eq.'kklomega') then
       ScalarToPrint=Density*TurbulentViscosityTl
       BScalarToPrint = BDensity*BTurbulentViscosityTl
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Turbulent-viscosity-(large-scale)', timeChar)
       
       ScalarToPrint=Density*TurbulentViscosityTs
       BScalarToPrint=BDensity*BTurbulentViscosityTs
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Turbulent-viscosity(small-scale)', timeChar)
     endif
     if(LTurbulentFlow) then
       ScalarToPrint=Viscosity+TurbulentViscosity
       BScalarToPrint=BViscosity+BTurbulentViscosity
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Effective-viscosity', timeChar)
       
       ScalarToPrint=TurbulentViscosity/Viscosity
       BScalarToPrint=BTurbulentViscosity/BViscosity
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Turbulent-viscosity-ratio', timeChar)
     endif
     if(LSolveTurbulenceKineticEnergy) then
       if(TurbulenceModel.eq.'kklomega') then         
         ScalarToPrint=100.*dsqrt(twothird*(TurbulentKE+TurbulentKL))/dmax1(Uinfinity,tiny)
         BScalarToPrint=100.*dsqrt(twothird*(BTurbulentKE+BTurbulentKL))/dmax1(Uinfinity,tiny)
       else
         ScalarToPrint=100.*dsqrt(twothird*TurbulentKE)/dmax1(Uinfinity,tiny)
         BScalarToPrint=100.*dsqrt(twothird*BTurbulentKE)/dmax1(Uinfinity,tiny) 
       endif
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Turbulent-intensity(%)', timeChar)
     endif
     if(LTurbulentFlow) then
       if(TurbulenceModel.ne.'spalartallmaras'.and.TurbulenceModel.ne.'wrayagarwal') then
         if(TurbulenceModel.eq.'kklomega') then
           ScalarToPrint=Density*ProductionKT   ! Production of k
           do i=1,NumberOfBCSets
             do j=1,NBFaces(i)
               k=NBFaceOwner(i,j)
               BScalarToPrint(i,j)=Density(k)*ProductionKT(k)
             enddo
           enddo
           call printFoamScalar(ScalarToPrint,BScalarToPrint,'Production-of-k)', timeChar)
           
           ScalarToPrint=Density*ProductionKL   ! Production of laminar k
           do i=1,NumberOfBCSets
             do j=1,NBFaces(i)
               k=NBFaceOwner(i,j)
               BScalarToPrint(i,j)=Density(k)*ProductionKL(k)
             enddo
           enddo
           call printFoamScalar(ScalarToPrint,BScalarToPrint,'Production-of-laminar-k)', timeChar)
         elseif(TurbulenceModel.eq.'sstgamaretheta') then
           call printFoamScalar(TurbulenceProduction,BTurbulenceProduction,'Production-of-k)', timeChar)
         else       
           call printFoamScalar(TurbulenceProduction,BTurbulenceProduction,'Production-of-k)', timeChar)
         endif
       endif
       if(LBuoyancy) then 
           call printFoamScalar(TurbulenceProductionB,BTurbulenceProductionB,'Production-of-k-by-buoyancy', timeChar)
       endif
     endif
     if(LTurbulentFlow.and.TurbulenceModel.ne.'wrayagarwal'.and.TurbulenceModel.ne.'spalartallmaras') then 
       call CalculateTurbulentReynoldsNumber
       call printFoamScalar(ReT,BReT,'Turbulent-Reynolds-number', timeChar)
     endif
     if(LTurbulentFlow.and.LtestTurbulenceModel) then
       call CalculateYplusPlotting
       call printFoamScalar(yplusPlot,ByplusPlot,'yplus', timeChar) 
       call printFoamScalar(uplusPlot,BuplusPlot,'uplus', timeChar)
       deallocate(yplusPlot)
       deallocate(uplusPlot)
       deallocate(ByplusPlot)
       deallocate(BuplusPlot)

     endif
     if(LTurbulentFlow) then
       call printFoamScalar(WallDistance, BWallDistance,'Normal-distance-to-wall(m)', timeChar)
     endif
     if(LSolveEnergy) then
       call printFoamScalar(Temperature,BTemperature,'Temperature(K)', timeChar)
       ScalarToPrint=Temperature+0.5*(uVelocity*uVelocity+vVelocity*vVelocity+wVelocity*wVelocity)/SpecificHeat   
       BScalarToPrint=BTemperature+0.5*(BuVelocity*BuVelocity+BvVelocity*BvVelocity+BwVelocity*BwVelocity)/BSpecificHeat   
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Total-Temperature(K)', timeChar) 
       ScalarToPrint=SpecificHeat*(Temperature-ReferenceTemperature)
       BScalarToPrint=BSpecificHeat*(BTemperature-ReferenceTemperature)
       call printFoamScalar(ScalarToPrint,BScalarToPrint,'Static-enthalpy', timeChar)
     endif
     if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
       call printFoamScalar(Htotal,BHtotal,'Total-enthalpy', timeChar)
     endif
     if(LSolveEnergy.and.LTurbulentFlow) then
       Variable1='temp'
       call CalculateEffectiveDiffusionCoefficient(Variable1)
       call printFoamScalar(eDiffCoefficient,BeDiffCoefficient,'Effective-thermal-conductivity', timeChar)
     endif
!     
      do irField=1,NumberOfrFieldsToSolve
        if(LSolverField(irField)) then
          ScalarToPrint(:)=rField(:,irField)  !rfield
          bScalarToPrint(:,:)=brField(:,:,irField)  !rfield 
          call printFoamScalar(ScalarToPrint,BScalarToPrint,rFieldName(irField), timeChar)
        endif
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LSolveScalar(iScalar)) then
          ScalarToPrint(:)=Scalar(:,iScalar)  !scalar
          bScalarToPrint(:,:)=bScalar(:,:,iScalar)  !scalar
          call printFoamScalar(ScalarToPrint,BScalarToPrint,ScalarName(iScalar), timeChar)
        endif
      enddo
      if(LCompressible) then
        call printFoamScalar(Density,BDensity,'Density', timeChar)
      endif
      if(.not.LCompressible) then
        if(LFreeSurfaceFlow.and.LSolveMomentum) then
          call printFoamScalar(Density,BDensity,'Density', timeChar)
        endif
      endif
      if(LCompressible) then
        call printFoamScalar(MachNumber,BMachNumber,'Mach-Number', timeChar) 
      endif
!
      if(LPrintGradients) then
!
        if(LSolveMomentum) then
          call printFoamScalar(uVelGradx,buVelGradx,'uVelGradx', timeChar)
          call printFoamScalar(uVelGrady,buVelGrady,'uVelGrady', timeChar)
          call printFoamScalar(uVelGradz,buVelGradz,'uVelGradz', timeChar)
          call printFoamScalar(vVelGradx,bvVelGradx,'vVelGradx', timeChar)
          call printFoamScalar(vVelGrady,bvVelGrady,'vVelGrady', timeChar)
          call printFoamScalar(vVelGradz,bvVelGradz,'vVelGradz', timeChar)
          call printFoamScalar(wVelGradx,bwVelGradx,'wVelGradx', timeChar)
          call printFoamScalar(wVelGrady,bwVelGrady,'wVelGrady', timeChar)
          call printFoamScalar(wVelGradz,bwVelGradz,'wVelGradz', timeChar)
          if(.not.LTurbulentFlow) then
            call AllocateStrainRateTensor
            call CalculateStrainRateTensor
          endif
          call printFoamScalar(StrainRate,bStrainRate,'Strain-rate', timeChar)
          if(.not.LTurbulentFlow) call deAllocateStrainRateTensor
          if(.not.LTurbulentFlow) then 
            call AllocateVorticityTensor
            call CalculateVorticityTensor
          endif
          call printFoamScalar(Vorticity,bVorticity,'Vorticity', timeChar)
          if(.not.LTurbulentFlow) call deAllocateVorticityTensor
        endif
        if(LSolveLambdaELEEquation) then
          call printFoamScalar(uVelGradx,buVelGradx,'uVelGradx', timeChar)
          call printFoamScalar(uVelGrady,buVelGrady,'uVelGrady', timeChar)
          call printFoamScalar(uVelGradz,buVelGradz,'uVelGradz', timeChar)
          call printFoamScalar(vVelGradx,bvVelGradx,'vVelGradx', timeChar)
          call printFoamScalar(vVelGrady,bvVelGrady,'vVelGrady', timeChar)
          call printFoamScalar(vVelGradz,bvVelGradz,'vVelGradz', timeChar)
          call printFoamScalar(wVelGradx,bwVelGradx,'wVelGradx', timeChar)
          call printFoamScalar(wVelGrady,bwVelGrady,'wVelGrady', timeChar)
          call printFoamScalar(wVelGradz,bwVelGradz,'wVelGradz', timeChar)
          call printFoamScalar(LambdaELEGradx,bLambdaELEGradx,'LambdaGradx', timeChar)
          call printFoamScalar(LambdaELEGrady,bLambdaELEGrady,'LambdaGrady', timeChar)
          call printFoamScalar(LambdaELEGradz,bLambdaELEGradz,'LambdaGradz', timeChar)
        endif
        if(LSolveContinuity) then
          call printFoamScalar(PressGradx,bPressGradx,'PressureGradx', timeChar)
          call printFoamScalar(PressGrady,bPressGrady,'PressureGrady', timeChar)
          call printFoamScalar(PressGradz,bPressGradz,'PressureGradz', timeChar)
        endif
        if(LSolveTurbulenceKineticEnergy) then
          call printFoamScalar(TKEGradx,bTKEGradx,'TKEGradx', timeChar)
          call printFoamScalar(TKEGrady,bTKEGrady,'TKEGrady', timeChar)
          call printFoamScalar(TKEGradz,bTKEGrady,'TKEGradz', timeChar)
        endif
        if(LSolveTurbulenceDissipationRate) then
          call printFoamScalar(TEDGradx,bTEDGradx,'TEDGradx', timeChar)
          call printFoamScalar(TEDGrady,bTEDGrady,'TEDGrady', timeChar)
          call printFoamScalar(TEDGradz,bTEDGradz,'TEDGradz', timeChar)
        endif
        if(LSolveTurbulenceV2Equation) then
          call printFoamScalar(TurbulentV2Gradx,bTurbulentV2Gradx,'v2Gradx', timeChar)
          call printFoamScalar(TurbulentV2Grady,bTurbulentV2Grady,'v2Grady', timeChar)
          call printFoamScalar(TurbulentV2Gradz,bTurbulentV2Gradz,'v2Gradz', timeChar)
        endif
        if(LSolveTurbulenceZetaEquation) then
          call printFoamScalar(TurbulentZetaGradx,bTurbulentZetaGradx,'ZetaGradx', timeChar)
          call printFoamScalar(TurbulentZetaGrady,bTurbulentZetaGrady,'ZetaGrady', timeChar)
          call printFoamScalar(TurbulentZetaGrady,bTurbulentZetaGrady,'ZetaGradz', timeChar)
        endif
        if(LSolveTurbulencefRelaxationEquation) then
          call printFoamScalar(TfRelaxationGradx,bTfRelaxationGradx,'fGradx', timeChar)
          call printFoamScalar(TfRelaxationGrady,bTfRelaxationGrady,'fGrady', timeChar)
          call printFoamScalar(TfRelaxationGradz,bTfRelaxationGradz,'fGradz', timeChar)
        endif
        if(LSolveTurbulenceSpecificDissipationRate) then
          call printFoamScalar(TOmegaGradx,bTOmegaGradx,'TOmegaGradx', timeChar)
          call printFoamScalar(TOmegaGrady,bTOmegaGrady,'TOmegaGrady', timeChar)
          call printFoamScalar(TOmegaGradz,bTOmegaGradz,'TOmegaGradz', timeChar)
        endif
        if(LSolveTurbulenceGammaEquation) then
          call printFoamScalar(TGammaGradx,bTGammaGradx,'TGammaGradx', timeChar)
          call printFoamScalar(TGammaGrady,bTGammaGrady,'TGammaGrady', timeChar)
          call printFoamScalar(TGammaGradz,bTGammaGradz,'TGammaGradz', timeChar)
        endif
        if(LSolveTurbulenceReynoldsThetaEquation) then
          call printFoamScalar(TReThetaGradx,bTReThetaGradx,'TReThetaGradx', timeChar)
          call printFoamScalar(TReThetaGrady,bTReThetaGrady,'TReThetaGrady', timeChar)
          call printFoamScalar(TReThetaGradz,bTReThetaGradz,'TReThetaGradz', timeChar)
        endif
        if(LSolveTurbulentKL) then
          call printFoamScalar(TurbulentKLGradx,bTurbulentKLGradx,'TKLGradx', timeChar)
          call printFoamScalar(TurbulentKLGrady,bTurbulentKLGrady,'TKLGrady', timeChar)
          call printFoamScalar(TurbulentKLGradz,bTurbulentKLGradz,'TKLGradz', timeChar)
        endif
        if(LSolveModifiedED) then
          call printFoamScalar(ModifiedEDGradx,bModifiedEDGradx,'MEDGradx', timeChar)
          call printFoamScalar(ModifiedEDGrady,bModifiedEDGrady,'MEDGrady', timeChar)
          call printFoamScalar(ModifiedEDGradz,bModifiedEDGradz,'MEDGradz', timeChar)
        endif
        if(LSolveEnergy) then
          call printFoamScalar(TempGradx,bTempGradx,'TempGradx', timeChar)
          call printFoamScalar(TempGrady,bTempGrady,'TempGrady', timeChar)
          call printFoamScalar(TempGradz,bTempGradz,'TempGradz', timeChar)
        endif
        if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
          call printFoamScalar(HtotalGradx,bHtotalGradx,'HTotalGradx', timeChar)
          call printFoamScalar(HtotalGrady,bHtotalGrady,'HTotalGrady', timeChar)
          call printFoamScalar(HtotalGradz,bHtotalGradz,'HTotalGradz', timeChar)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!problem here!!!!!!!!!!!!!!!!!!!!!!!!!
        do irField=1,NumberOfrFieldsToSolve
          if(LSolverField(irField)) then
            part1=trim(rFieldName(irField))
            part2="Gradx"
            part3="Grady"
            part4="Gradz"
            namesize=len(part1)+5
            name1=trim(part1)//part2
            name2=trim(part1)//part3
            name3=trim(part1)//part4
            
            ScalarToPrint=rFieldGradx(:,irField) 
            bScalarToPrint=BrFieldGradx(:,:,irField)
            call printFoamScalar(ScalarToPrint,BScalarToPrint,name1, timeChar)
            
            ScalarToPrint=rFieldGrady(:,irField)
            bScalarToPrint=BrFieldGrady(:,:,irField)
            call printFoamScalar(ScalarToPrint,BScalarToPrint,name2, timeChar)
            
            ScalarToPrint=rFieldGradz(:,irField)
            bScalarToPrint=BrFieldGradz(:,:,irField)
            call printFoamScalar(ScalarToPrint,BScalarToPrint,name3, timeChar)
          endif
        enddo
        
        do iScalar=1,NumberOfScalarsToSolve
          if(LSolveScalar(iScalar)) then
            part1=trim(ScalarName(iScalar))
            part2="Gradx"
            part3="Grady"
            part4="Gradz"
            namesize=len(part1)+5
            name1=trim(part1)//part2
            name2=trim(part1)//part3
            name3=trim(part1)//part4
            
            ScalarToPrint=ScalarGradx(:,iScalar) 
            bScalarToPrint=bScalarGradx(:,:,iScalar) 
            call printFoamScalar(ScalarToPrint,BScalarToPrint,name1, timeChar)
            
            ScalarToPrint=ScalarGrady(:,iScalar) 
            bScalarToPrint=bScalarGrady(:,:,iScalar) 
            call printFoamScalar(ScalarToPrint,BScalarToPrint,name2, timeChar)
            
            ScalarToPrint=ScalarGradz(:,iScalar) 
            bScalarToPrint=bScalarGradz(:,:,iScalar) 
            call printFoamScalar(ScalarToPrint,BScalarToPrint,name3, timeChar)
          endif
        enddo
      endif
!      
      deallocate(ScalarToPrint)
      deallocate(bScalarToPrint)
      return
    end SUBROUTINE PrintFoam
!     
!************************************************************************************************
subroutine printFoamScalar(Variable,BFVariable , NameV, timeChar)
use Geometry1, only: NumberOfBCSets, NumberOfElements
use Geometry3, only: NBFaces
use Geometry2, only:BoundaryName
use user0, only:FoamCasedirectory

implicit none
double precision, allocatable :: Variable(:)
double precision, allocatable :: BFVariable(:,:)
CHARACTER(LEN=*), INTENT(IN) :: NameV, timeChar
integer :: i, j 

open(unit = 30, file = trim(FoamCasedirectory)//'\'//trim(adjustl(timeChar))//'\'//trim(adjustl(NameV)))

write(30,'(a)') "/*--------------------------------*- C++ -*----------------------------------*\"
write(30,'(a)') "| =========                 |                                                 |"
write(30,'(a)') "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |"
write(30,'(a)') "|  \\    /   O peration     | Version:  2.0.1                                 |"
write(30,'(a)') "|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |"
write(30,'(a)') "|    \\/     M anipulation  |                                                 |"
write(30,'(a)') "\*---------------------------------------------------------------------------*/"
write(30,'(a)') "FoamFile"
write(30,'(a)') "{"
write(30,'(a)') "    version     2.0;"
write(30,'(a)') "    format      ascii;"
write(30,'(a)') "    class       volScalarField;"
write(30,*) "    object      ",NameV,";"
write(30,'(a)') "}"
write(30,'(a)') "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
write(30,'(a)') "internalField   nonuniform List<scalar> "
write(30,*) NumberOfElements
write(30,'(a)') "("
do i=1,NumberOfElements
   write(30,*) Variable(i)
enddo
write(30,'(a)') ");"

write(30,'(a)') "boundaryField"
write(30,'(a)') "{"

do i=1,NumberOfBCSets
   write(30,*) trim(adjustl(BoundaryName(i)))
   write(30,'(a)') "{"
   write(30,'(a)') 'value  nonuniform List<scalar>'
   write(30,*) NBFaces(i)
   write(30,'(a)') "("
   do j=1,NBFaces(i)
      write(30,*) BFVariable(i,j)
   enddo
   write(30,'(a)') ");"
   write(30,'(a)') "}"
enddo
write(30,'(a)') "}"
close(30)
end subroutine printFoamScalar
!************************************************************************************************
subroutine printFoamVector(Variable,BFVariable , NameV, timeChar)
use Geometry1, only: NumberOfBCSets, NumberOfElements
use Geometry3, only: NBFaces
use Geometry2, only:BoundaryName
use user0, only:FoamCasedirectory

implicit none
double precision, allocatable :: Variable(:,:) !Variable(5,1) cell 5 1>xcomponent
double precision, allocatable :: BFVariable(:,:,:) !Variable(5,6,2) patch 5 face 6 variabe 2>y-component 
CHARACTER(LEN=*), INTENT(IN) :: NameV, timeChar
integer :: i, j 

open(unit = 30, file = trim(FoamCasedirectory)//'\'//trim(adjustl(timeChar))//'\'//trim(adjustl(NameV)))

write(30,'(a)') "/*--------------------------------*- C++ -*----------------------------------*\"
write(30,'(a)') "| =========                 |                                                 |"
write(30,'(a)') "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |"
write(30,'(a)') "|  \\    /   O peration     | Version:  2.0.1                                 |"
write(30,'(a)') "|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |"
write(30,'(a)') "|    \\/     M anipulation  |                                                 |"
write(30,'(a)') "\*---------------------------------------------------------------------------*/"
write(30,'(a)') "FoamFile"
write(30,'(a)') "{"
write(30,'(a)') "    version     2.0;"
write(30,'(a)') "    format      ascii;"
write(30,'(a)') "    class       volVectorField;"
write(30,*) "    object      ",NameV,";"
write(30,'(a)') "}"
write(30,'(a)') "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
write(30,'(a)') "internalField   nonuniform List<vector> "
write(30,*) size(Variable,1)
write(30,'(a)') "("

do i=1,NumberOfElements
   write(30,*) '(',Variable(i,1),Variable(i,2),Variable(i,3),')'
enddo
write(30,'(a)') ");"

write(30,'(a)') "boundaryField"
write(30,'(a)') "{"

do i=1,NumberOfBCSets
   write(30,*) trim(adjustl(BoundaryName(i)))
   write(30,'(a)') "{"
   write(30,'(a)') 'value  nonuniform List<vector>'
   write(30,*) NBFaces(i)
   write(30,'(a)') "("
   do j=1,NBFaces(i)
      write(30,*) '(',BFVariable(i,j,1),BFVariable(i,j,2),BFVariable(i,j,3),')'
   enddo
   write(30,'(a)') ");"
   write(30,'(a)') "}"
enddo
write(30,'(a)') "}"
close(30)
end subroutine printFoamVector

end module WriteFoamModule
