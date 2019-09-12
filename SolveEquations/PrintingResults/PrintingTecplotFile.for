c
C#############################################################################################
c
      SUBROUTINE Tecplot
c
C#############################################################################################
      use Geometry1
      use Geometry3, only:NBFacesMax
      use User0, only: LSolveMomentum,LSolveContinuity,LSolveEnergy,
     *                 LCompressible,NumberOfScalarsToSolve,
     *                 LSolveScalar,ScalarName,LTurbulentFlow,
     *                 LSolveTurbulenceKineticEnergy,LBuoyancy,
     *                 LSolveTurbulenceDissipationRate,LTransitionalSA,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LsolveModifiedED,LtestTurbulenceModel,
     *                 LSolveTurbulentKL,EnergyEquation,
     *                 NumberOfrFieldsToSolve,LSolverField,rFieldName,
     *                 LFreeSurfaceFlow,LSolveTurbulenceGammaEquation,
     *                 LSolveTurbulenceReynoldsThetaEquation,
     *                 LSolveTurbulencefRelaxationEquation,
     *                 LSolveTurbulenceV2Equation,
     *                 LSolveTurbulenceZetaEquation,
     *                 LPrintGradients,TurbulenceModel,LUnsteady,
     *                 LSolveLambdaELEEquation
      use Variables1, only: TurbulentKE,BTurbulentKE
      use tecplot1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      logical :: LPrintGradientsTemp
      character*10 :: part1
      character*5 :: part2,part3,part4
      character*15 :: name1,name2,name3
      character*15, save :: zonetype,zonetype1
      integer n,iScalar,irField,namesize
      integer i,j,k,j1
      integer j2,j3,j4
c*********************************************************************************************
  10  format(1x,'Variables= "x (m)", "y (m)", "z (m)" ')
  11  format(1x,',"Lambda Eulerian Lagrangian Equation"') 
  12  format(1x,',"Initial velocity divergence (1/s)"') 
  13  format(1x,',"Final velocity divergence (1/s)"') 
  20  format(1x,',"u-Velocity (m/s)", "v-Velocity (m/s)",
     * "w-Velocity (m/s)" ')
  21  format(1x,',"Velocity-Magnitude (m/s)"')
  30  format(1x,',"Pressure (Pa)"')
  35  format(1x,',"Turbulent kinetic energy (m2/s2)" ')
  36  format(1x,',"Turbulent dissipation rate (m2/s3)" ')
  96  format(1x,',"Turbulent v2 (m2/s2)" ')
  97  format(1x,',"Turbulent zeta" ')
  98  format(1x,',"f relaxation (m2/s)" ')
  37  format(1x,',"Turbulent specific dissipation rate (1/s)" ')
  61  format(1x,',"Intermittency (Turbulent gama)" ')
  62  format(1x,',"Momentum thickness (turbulent Re Theta)" ')
  15  format(1x,',"Laminar kinetic energy (m2/s2)" ')
  17  format(1x,',"Total fluctuation energy (m2/s2)" ')
  39  format(1x,',"Modified eddy diffusivity (m2/s)" ')
  22  format(1x,',"Laminar viscosity (Pa.s)" ')
  38  format(1x,',"Turbulent viscosity (Pa.s)" ')
  42  format(1x,',"Production of k (Kg/(m.s3))" ')
  51  format(1x,',"Production of laminar k (Kg/(m.s3))" ')
  52  format(1x,',"Turbulent viscosity (large-scale) (Pa.s)" ')
  53  format(1x,',"Turbulent viscosity (small-scale) (Pa.s)" ')
  54  format(1x,',"Effective viscosity (Pa.s)" ')
  55  format(1x,',"Turbulent viscosity ratio" ')
  56  format(1x,',"Effective thermal conductivity (W/(m.K))" ')
  57  format(1x,',"Turbulent intensity (%)" ')
  58  format(1x,',"Intermittency effective" ')
  43  format(1x,',"Production of k by buoyancy (Kg/(m.s3))" ')
  44  format(1x,',"Turbulent Reynolds number" ')
  45  format(1x,',"yplus" ')
  46  format(1x,',"uplus" ')
  47  format(1x,',"Normal distance to wall (m)" ')
  40  format(1x,',"Temperature (K)" ')
  16  format(1x,',"Total Temperature (K)" ')
  48  format(1x,',"Static enthalpy (J/Kg)" ')
  49  format(1x,',"Total enthalpy (J/Kg)" ')
  41  format(1x,',"',A10,'" ')
  50  format(1x,',"Density (Kg/m3)" ')
  60  format(1x,',"Mach" ')
  70  format(1x,\)      
  80  format(1x,',"uVelGradx (1/s)", "uVelGrady (1/s)",
     * "uVelGradz (1/s)" ')
  81  format(1x,',"vVelGradx (1/s)", "vVelGrady (1/s)",
     * "vVelGradz (1/s)" ')
  82  format(1x,',"wVelGradx (1/s)", "wVelGrady (1/s)",
     * "wVelGradz (1/s)" ')
  83  format(1x,',"PressureGradx (Pa/m)", "PressureGrady (Pa/m)",
     * "PressureGradz (Pa/m)" ')
  84  format(1x,',"TKEGradx (m/s2)", "TKEGrady (m/s2)",
     * "TKEGradz (m/s2)" ')
  85  format(1x,',"TEDGradx (m/s3)", "TEDGrady (m/s3)",
     * "TEDGradz (m/s3)" ')
  86  format(1x,',"TOmegaGradx (1/(m.s))", "TOmegaGrady (1/(m.s))",
     * "TOmegaGradz (1/(m.s))" ')
  87  format(1x,',"TGammaGradx (1/m)", "TGammaGrady (1/m)",
     * "TGammaGradz (1/m)" ')
  88  format(1x,',"TReThetaGradx (1/m)", "TReThetaGrady (1/m)",
     * "TReThetaGradz (1/m)" ')
  89  format(1x,',"TKLGradx (m/s2)", "TKLGrady (m/s2)",
     * "TKLGradz (m/s2)" ')
  90  format(1x,',"MEDGradx (m/s)", "MEDGrady (m/s)",
     * "MEDGradz (m/s)" ')
  91  format(1x,',"TempGradx (K/m)", "TempGrady (K/m)",
     * "TempGradz (K/m)" ')
  92  format(1x,',"HTotalGradx (J/(Kg.m))", "HTotalGrady (J/(Kg.m))",
     * "HTotalGradz (J/(Kg.m))" ')
  93  format(1x,',"',A15,'" ')
  94  format(1x,',"Strain rate (1/s)" ')
  95  format(1x,',"Vorticity (1/s)" ')
 101  format(1x,',"v2Gradx (m/s2)", "v2Grady (m/s2)",
     * "v2Gradz (m/s2)" ')
 102  format(1x,',"ZetaGradx (1/m)", "ZetaGrady (1/m)",
     * "ZetaGradz (1/m)" ')
 103  format(1x,',"fGradx (m/s)", "fGrady (m/s)", "fGradz (m/s)" ')
 104  format(1x,',"LambdaGradx", "LambdaGrady", "LambdaGradz" ')
c
c*********************************************************************************************
c
      if(LUnsteady) LtestTurbulenceModel=.false. !in order not to modify
                                                 !normal distance to wall
c
      LPrintGradientsTemp=LPrintGradients
c
      allocate(dummyE(NumberOfElements))
      allocate(BdummyE(NumberOfBCSets,NBFacesMax))
c
      write(12,*) 'Title = "Example" '
      write(12,10,advance='no') 
      if(LSolveMomentum.or.LSolveLambdaELEEquation) 
     *                                write(12,20,advance='no') 
      if(LSolveMomentum.or.LSolveLambdaELEEquation) 
     *                                write(12,21,advance='no') 
      if(LSolveLambdaELEEquation) write(12,11,advance='no') 
      if(LSolveLambdaELEEquation) write(12,12,advance='no') 
      if(LSolveLambdaELEEquation) write(12,13,advance='no') 
      if(LSolveContinuity) write(12,30,advance='no') 
      if(LSolveTurbulenceKineticEnergy) write(12,35,advance='no')
      if(LTurbulentFlow.and.TurbulenceModel.eq.'spalartallmaras')
     *                                   write(12,35,advance='no')
      if(LTurbulentFlow.and.TurbulenceModel.eq.'wrayagarwal')
     *                                   write(12,35,advance='no')
      if(LTurbulentFlow.and.TurbulenceModel.eq.'nut92')
     *                                   write(12,35,advance='no')
      if(LSolveTurbulenceDissipationRate) write(12,36,advance='no') 
      if(LSolveTurbulenceV2Equation) write(12,96,advance='no')
      if(LSolveTurbulenceZetaEquation) write(12,97,advance='no')
      if(LSolveTurbulencefRelaxationEquation) write(12,98,advance='no')
      if(LSolveTurbulenceSpecificDissipationRate) 
     *                                   write(12,37,advance='no') 
      if(LTurbulentFlow.and.TurbulenceModel.eq.
     *            'spalartallmaras'.and.LTransitionalSA) then
        write(12,61,advance='no')
      endif
      if(LSolveTurbulenceGammaEquation) then
        write(12,61,advance='no') 
        if(TurbulenceModel.eq.'sstgamaretheta') 
     *                        write(12,58,advance='no')
      endif
      if(LSolveTurbulenceReynoldsThetaEquation)write(12,62,advance='no')
      if(LSolveTurbulentKL) write(12,15,advance='no') 
      if(LTurbulentFlow.and.TurbulenceModel.eq.'kklomega')
     *                                         write(12,17,advance='no')
      if(LSolveModifiedED) write(12,39,advance='no') 
      if(LSolveMomentum) write(12,22,advance='no') 
      if(LTurbulentFlow) write(12,38,advance='no')
      if(LTurbulentFlow.and.TurbulenceModel.eq.
     *                    'kklomega') then
        write(12,52,advance='no')
        write(12,53,advance='no')
      endif  
      if(LTurbulentFlow) then
        write(12,54,advance='no')
        write(12,55,advance='no')
      endif
      if(LSolveTurbulenceKineticEnergy) write(12,57,advance='no')
      if(LTurbulentFlow) then
        if(TurbulenceModel.ne.'spalartallmaras'.and.
     *                   TurbulenceModel.ne.'wrayagarwal') then
          if(TurbulenceModel.eq.'kklomega') then
            write(12,42,advance='no')
            write(12,51,advance='no')
          elseif(TurbulenceModel.eq.'sstgamaretheta') then
            write(12,42,advance='no')
          else       
            write(12,42,advance='no')
          endif
        endif
        if(LBuoyancy) write(12,43,advance='no')
      endif
      if(LTurbulentFlow.and.TurbulenceModel.ne.'wrayagarwal'.and.
     *  TurbulenceModel.ne.'spalartallmaras') write(12,44,advance='no')
      if(LTurbulentFlow.and.LtestTurbulenceModel) then
        write(12,45,advance='no')
        write(12,46,advance='no')
      endif
      if(LTurbulentFlow) write(12,47,advance='no')
      if(LSolveEnergy) write(12,40,advance='no')
      if(LSolveEnergy) write(12,16,advance='no')
      if(LSolveEnergy) write(12,48,advance='no')
      if(LSolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                         write(12,49,advance='no')
      if(LSolveEnergy.and.LTurbulentFlow) write(12,56,advance='no')
      do irField=1,NumberOfrFieldsToSolve
        if(LSolverField(irField)) write(12,41,advance='no') 
     *                                      rFieldName(irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LSolveScalar(iScalar)) write(12,41,advance='no') 
     *                                      ScalarName(iScalar)
      enddo
      if(LCompressible) write(12,50,advance='no') 
      if(.not.LCompressible) then
        if(LFreeSurfaceFlow.and.LSolveMomentum) then
          write(12,50,advance='no')
        endif
      endif
      if(LCompressible) write(12,60,advance='no') 
c
      if(LPrintGradientsTemp) then
c
        if(LSolveMomentum) then
c
          write(12,80,advance='no') 
          write(12,81,advance='no') 
          write(12,82,advance='no') 
          write(12,94,advance='no') 
          write(12,95,advance='no') 
c
        endif
c
        if(LSolveLambdaELEEquation) then
c
          write(12,80,advance='no') 
          write(12,81,advance='no') 
          write(12,82,advance='no') 
c
        endif
c
        if(LSolveLambdaELEEquation) then
          write(12,104,advance='no') 
        endif
c
        if(LSolveContinuity) then
          write(12,83,advance='no') 
        endif
        if(LSolveTurbulenceKineticEnergy) then
          write(12,84,advance='no') 
        endif
        if(LSolveTurbulenceDissipationRate) then
          write(12,85,advance='no') 
        endif
        if(LSolveTurbulenceV2Equation) then
          write(12,101,advance='no')
        endif
        if(LSolveTurbulenceZetaEquation) then
          write(12,102,advance='no')
        endif
        if(LSolveTurbulencefRelaxationEquation) then
          write(12,103,advance='no')
        endif
        if(LSolveTurbulenceSpecificDissipationRate) then
          write(12,86,advance='no') 
        endif
        if(LSolveTurbulenceGammaEquation) then
          write(12,87,advance='no') 
        endif
        if(LSolveTurbulenceReynoldsThetaEquation) then
          write(12,88,advance='no') 
        endif
        if(LSolveTurbulentKL) then
          write(12,89,advance='no') 
        endif
        if(LSolveModifiedED) then
          write(12,90,advance='no') 
        endif
        if(LSolveEnergy) then
          write(12,91,advance='no') 
        endif
        if(LSolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                         write(12,92,advance='no')
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
            write(12,93,advance='no') name1(1:namesize)
            write(12,93,advance='no') name2(1:namesize)
            write(12,93,advance='no') name3(1:namesize)
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
            write(12,93,advance='no') name1(1:namesize)
            write(12,93,advance='no') name2(1:namesize)
            write(12,93,advance='no') name3(1:namesize)
          endif
        enddo
      endif
c
      write(12,70)
c
      do k=1,NumberOfElementGroups
c
        LPrintGradientsTemp=.false.
c
        if(MaxNodesPerGroupElement(k).ne.
     *              MinNodesPerGroupElement(k)) then
          zonetype='BRICK'
          zonetype1='MIXED'
        elseif(MaxNodesPerGroupElement(k).eq.4) then
          zonetype='TETRAHEDRON'
          zonetype1='TETRAHEDRON'
        elseif(MaxNodesPerGroupElement(k).eq.8) then
          zonetype='BRICK'
          zonetype1='BRICK'
        endif
c
        write(12,*) 'Zone ','N=', NumberOfNodesInGroup(k),'  ,E=', 
     *  NumberOfGroupElements(k),' ,F=FEBLOCK,  ET=', zonetype
c
        do i=1,NumberOfNodesInGroup(k)
          j=ListOfNodesInGroup(i,k)
          write(12,*) x(j)
        enddo
c
        do i=1,NumberOfNodesInGroup(k)
          j=ListOfNodesInGroup(i,k)
          write(12,*) y(j)
        enddo
c
        do i=1,NumberOfNodesInGroup(k)
          j=ListOfNodesInGroup(i,k)
          write(12,*) z(j)
        enddo
c
        If(LSolveMomentum.or.LSolveLambdaELEEquation) then
c
          call InterpolateFromElementsToNodes(1,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
c
          deallocate(dummyN)
          call InterpolateFromElementsToNodes(2,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
c
          deallocate(dummyN)
          call InterpolateFromElementsToNodes(3,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
c
          deallocate(dummyN)
          call InterpolateFromElementsToNodes(14,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        If(LSolveLambdaELEEquation) then
c
          call InterpolateFromElementsToNodes(40,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
c
          deallocate(dummyN)
c
          call InterpolateFromElementsToNodes(41,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
c
          deallocate(dummyN)
c
          call InterpolateFromElementsToNodes(42,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
c
          deallocate(dummyN)
c          
        endif
c      
        If(LSolveContinuity) then
c
          call InterpolateFromElementsToNodes(4,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        If(LSolveTurbulenceKineticEnergy) then
c
          call InterpolateFromElementsToNodes(8,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        if(LTurbulentFlow.and.TurbulenceModel.eq.'spalartallmaras') then
c
          call CalculateOneEquationModelTurbulentKE
          call InterpolateFromElementsToNodes(8,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        if(LTurbulentFlow.and.TurbulenceModel.eq.'wrayagarwal') then
c
          call CalculateOneEquationModelTurbulentKE
          call InterpolateFromElementsToNodes(8,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        if(LTurbulentFlow.and.TurbulenceModel.eq.'nut92') then
c
          call CalculateOneEquationModelTurbulentKE
          call InterpolateFromElementsToNodes(8,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        If(LSolveTurbulenceDissipationRate) then
c
          call InterpolateFromElementsToNodes(9,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        if(LSolveTurbulenceV2Equation) then
c
          call InterpolateFromElementsToNodes(37,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        if(LSolveTurbulenceZetaEquation) then
c
          call InterpolateFromElementsToNodes(38,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        if(LSolveTurbulencefRelaxationEquation) then
c
          call InterpolateFromElementsToNodes(39,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        If(LSolveTurbulenceSpecificDissipationRate) then
c
          call InterpolateFromElementsToNodes(10,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
      if(LTurbulentFlow.and.TurbulenceModel.eq.
     *            'spalartallmaras'.and.LTransitionalSA) then
c
          call InterpolateFromElementsToNodes(25,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
      endif
c
        If(LSolveTurbulenceGammaEquation) then
c
          call InterpolateFromElementsToNodes(25,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
          if(TurbulenceModel.eq.'sstgamaretheta') then
            call InterpolateFromElementsToNodes
     *                           (36,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
          endif      
c
        endif      
c
        If(LSolveTurbulenceReynoldsThetaEquation) then
c
          call InterpolateFromElementsToNodes(26,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        If(LSolveTurbulentKL) then
c
          call InterpolateFromElementsToNodes(20,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        if(LTurbulentFlow.and.TurbulenceModel.eq.
     *                            'kklomega') then
c
          call InterpolateFromElementsToNodes(27,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        If(LSolveModifiedED) then
c
          call InterpolateFromElementsToNodes(18,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        If(LSolveMomentum) then
c
          call InterpolateFromElementsToNodes(19,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        If(LTurbulentFlow) then
c
          call InterpolateFromElementsToNodes(11,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        if(LTurbulentFlow.and.TurbulenceModel.eq.
     *                              'kklomega') then
c
          call InterpolateFromElementsToNodes(28,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
          call InterpolateFromElementsToNodes(29,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        if(LTurbulentFlow) then
c
          call InterpolateFromElementsToNodes(30,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
          call InterpolateFromElementsToNodes(31,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          call InterpolateFromElementsToNodes(32,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        if(LTurbulentFlow) then
c
          if(TurbulenceModel.ne.'spalartallmaras'.and.
     *                   TurbulenceModel.ne.'wrayagarwal') then
c
            if(TurbulenceModel.eq.'kklomega') then
c
              call InterpolateFromElementsToNodes
     *                              (33,1,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
              call InterpolateFromElementsToNodes
     *                               (34,1,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
            elseif(TurbulenceModel.eq.'sstgamaretheta') then
c
              call InterpolateFromElementsToNodes
     *                              (12,1,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
            else
c
              call InterpolateFromElementsToNodes
     *                              (12,1,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
            endif
c
          endif
c
          if(LBuoyancy) then
c
            call InterpolateFromElementsToNodes
     *                               (13,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
           endif      
c
          if(TurbulenceModel.ne.'wrayagarwal'.and.
     *                TurbulenceModel.ne.'spalartallmaras') then
            call CalculateTurbulentReynoldsNumber
c
            call InterpolateFromElementsToNodes
     *                                 (15,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
          endif
c
          if(LtestTurbulenceModel) then
c
            call CalculateYplusPlotting
c
            call InterpolateFromElementsToNodes
     *                                    (16,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
            call InterpolateFromElementsToNodes
     *                                    (17,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
          endif
c
          call InterpolateFromElementsToNodes(21,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        If(LSolveEnergy) then
c
          call InterpolateFromElementsToNodes(5,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
          call InterpolateFromElementsToNodes(24,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
          call InterpolateFromElementsToNodes(22,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
          if(EnergyEquation.eq.'htotal') then
c
            call InterpolateFromElementsToNodes
     *                                      (23,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
          endif
c
          if(LTurbulentFlow) then
c
            call InterpolateFromElementsToNodes
     *                                      (35,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
c
          endif
c
        endif      
c
        do irField=1,NumberOfrFieldsToSolve
c
          if(LSolverField(irField)) then
c
            call InterpolateFromElementsToNodes
     *                                (42+irField,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
          endif
c
        enddo
c
        do iScalar=1,NumberOfScalarsToSolve
c
          if(LSolveScalar(iScalar)) then
c
            call InterpolateFromElementsToNodes
     *       (42+NumberOfrFieldsToSolve+iScalar,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
          endif
c
        enddo
c
        If(LCompressible) then
c
          call InterpolateFromElementsToNodes(6,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif
c
        if(.not.LCompressible) then
          if(LFreeSurfaceFlow.and.LSolveMomentum) then
c
            call InterpolateFromElementsToNodes(6,1,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
            deallocate(dummyN)
c
          endif
        endif
c
        If(LCompressible) then
c
          call InterpolateFromElementsToNodes(7,1,LPrintGradientsTemp)
c
          do i=1,NumberOfNodesInGroup(k)
            j=ListOfNodesInGroup(i,k)
            write(12,*) dummyN(j)
          enddo
          deallocate(dummyN)
c
        endif      
c
        LPrintGradientsTemp=LPrintGradients
c
        if(LPrintGradientsTemp) then
c
          If(LSolveMomentum) then
c
            call InterpolateFromElementsToNodes(1,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(2,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(3,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
            call InterpolateFromElementsToNodes(4,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(5,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(6,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
            call InterpolateFromElementsToNodes(7,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(8,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(9,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
            if(.not.LTurbulentFlow) call AllocateStrainRateTensor
c
            call CalculateStrainRateTensor
            call InterpolateFromElementsToNodes
     *                     (40,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            if(.not.LTurbulentFlow) call deAllocateStrainRateTensor
c
            if(.not.LTurbulentFlow) call AllocateVorticityTensor
c
            call CalculateVorticityTensor
            call InterpolateFromElementsToNodes
     *                     (41,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            if(.not.LTurbulentFlow) call deAllocateVorticityTensor
c
          endif
c
          If(LSolveLambdaELEEquation) then
c
            call InterpolateFromElementsToNodes(1,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(2,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(3,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
            call InterpolateFromElementsToNodes(4,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(5,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(6,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
            call InterpolateFromElementsToNodes(7,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(8,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes(9,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveLambdaELEEquation) then
c
            call InterpolateFromElementsToNodes
     *                                (51,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
            call InterpolateFromElementsToNodes
     *                               (52,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
            call InterpolateFromElementsToNodes
     *                               (53,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
          endif
c
          if(LSolveContinuity) then
c
            call InterpolateFromElementsToNodes
     *                                     (10,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                    (11,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (12,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulenceKineticEnergy) then
c
            call InterpolateFromElementsToNodes
     *                                   (13,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (14,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (15,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulenceDissipationRate) then
c
            call InterpolateFromElementsToNodes
     *                                   (16,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (17,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (18,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulenceV2Equation) then
c
            call InterpolateFromElementsToNodes
     *                                   (42,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (43,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (44,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulenceZetaEquation) then
c
            call InterpolateFromElementsToNodes
     *                                   (45,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (46,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (47,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulencefRelaxationEquation) then
c
            call InterpolateFromElementsToNodes
     *                                   (48,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (49,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (50,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulenceSpecificDissipationRate) then
c
            call InterpolateFromElementsToNodes
     *                                    (19,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                     (20,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                    (21,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulenceGammaEquation) then
c
            call InterpolateFromElementsToNodes
     *                                   (22,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (23,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (24,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulenceReynoldsThetaEquation) then
c
            call InterpolateFromElementsToNodes
     *                                   (25,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (26,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                   (27,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveTurbulentKL) then
c
            call InterpolateFromElementsToNodes
     *                                  (28,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                 (29,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                (30,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveModifiedED) then
c
            call InterpolateFromElementsToNodes
     *                                  (31,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                 (32,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                 (33,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveEnergy) then
c
            call InterpolateFromElementsToNodes
     *                                  (34,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                  (35,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                  (36,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
c
            call InterpolateFromElementsToNodes
     *                                   (37,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                  (38,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
            call InterpolateFromElementsToNodes
     *                                  (39,2,LPrintGradientsTemp)
c
            do i=1,NumberOfNodesInGroup(k)
              j=ListOfNodesInGroup(i,k)
              write(12,*) dummyN(j)
            enddo
c
            deallocate(dummyN)
c
          endif
c
          do irField=1,NumberOfrFieldsToSolve
c
            if(LSolverField(irField)) then
c
              call InterpolateFromElementsToNodes
     *                        (53+3*irField-2,2,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
              call InterpolateFromElementsToNodes
     *                        (53+3*irField-1,2,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
              call InterpolateFromElementsToNodes
     *                        (53+3*irField,2,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
            endif
c
          enddo
c
          do iScalar=1,NumberOfScalarsToSolve
c
            if(LSolveScalar(iScalar)) then
c
              call InterpolateFromElementsToNodes
     *                (53+3*irField+3*iScalar-2,2,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
              call InterpolateFromElementsToNodes
     *                (53+3*irField+3*iScalar-1,2,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
              call InterpolateFromElementsToNodes
     *                 (53+3*irField+3*iScalar,2,LPrintGradientsTemp)
c
              do i=1,NumberOfNodesInGroup(k)
                j=ListOfNodesInGroup(i,k)
                write(12,*) dummyN(j)
              enddo
              deallocate(dummyN)
c
            endif
c
          enddo
c
        endif
c
        do i=1,NumberOfGroupElements(k)
          j=NElementsInGroup(k,i)  
c
          if(zonetype1.eq.'TETRAHEDRON') then
c
            write(12,*) ListOfElementNodeslocal(i,4,k),
     *                  ListOfElementNodeslocal(i,1,k),
     *                  ListOfElementNodeslocal(i,3,k),
     *                  ListOfElementNodeslocal(i,2,k)
c
          elseif(zonetype1.eq.'BRICK') then
c
            write(12,*) ListOfElementNodeslocal(i,5,k),
     *                  ListOfElementNodeslocal(i,6,k),
     *                  ListOfElementNodeslocal(i,2,k),
     *                  ListOfElementNodeslocal(i,1,k),
     *                  ListOfElementNodeslocal(i,7,k),
     *                  ListOfElementNodeslocal(i,8,k),
     *                  ListOfElementNodeslocal(i,4,k),
     *                  ListOfElementNodeslocal(i,3,k)
c
          elseif(zonetype1.eq.'MIXED') then
c          
            if(NumbOfElementNodes(j).eq.4) then         
c
            write(12,*) ListOfElementNodeslocal(i,4,k),
     *                  ListOfElementNodeslocal(i,1,k),
     *                  ListOfElementNodeslocal(i,3,k),
     *                  ListOfElementNodeslocal(i,2,k)
c
            elseif(NumbOfElementNodes(j).eq.8) then         
c
            write(12,*) ListOfElementNodeslocal(i,5,k),
     *                  ListOfElementNodeslocal(i,6,k),
     *                  ListOfElementNodeslocal(i,2,k),
     *                  ListOfElementNodeslocal(i,1,k),
     *                  ListOfElementNodeslocal(i,7,k),
     *                  ListOfElementNodeslocal(i,8,k),
     *                  ListOfElementNodeslocal(i,4,k),
     *                  ListOfElementNodeslocal(i,3,k)
c
            endif
c
          endif
c
        enddo
c
      enddo
c
      deallocate(dummyE)
      deallocate(BdummyE)
c
      return
      end