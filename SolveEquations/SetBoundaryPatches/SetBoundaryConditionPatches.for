c
C#############################################################################################
      SUBROUTINE SetBCPatches
C#############################################################################################
      use Geometry1, only : NumberOfBCSets
      use Geometry3, only : NBFaces,NBFaceOwner 
      use BoundaryConditions1
      use BoundaryConditions2
      use BoundaryConditionsScalar1
      use BoundaryConditionsScalar2
      use BoundaryConditionsrField1
      use BoundaryConditionsrField2
      use BoundaryConditionsTurbulence1
      use BoundaryConditionsTurbulence2
      use BoundaryFluxes
      use Turbulence1, only: yplus,yplusT,WallViscosity,
     *                       uplus,ustar,duplusdyplus,
     *                       WallViscosityT,TauWall,ystar,
     *                       WallViscosityS,yplusS,uTau,
     *                       WallViscosityR,yplusR,KsPlus,
     *                       RoughccM,RoughccE,RoughccS
      use User0
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k,i1,i2,i3
c*********************************************************************************************
c
      Iwall=0
      IWallSlip=0
      IWallnoSlip=0
      IWallDirichlet=0
      IWallVonNeumann=0
      IWallRobin=0
      IWallTurbulence=0
c      
      Iinlet=0
      Iinletsupersonic=0
      IinletSpecifiedVelocity=0
      IinletSpecifiedMassFlowRate=0
      IinletSpecifiedStaticPressure=0
      IinletSpecifiedStagnationPressure=0
      IinletSpecifiedStaticTemperature=0
      IinletSpecifiedStagnationTemperature=0
      IinletTurbulence=0
c
      Ioutlet=0
      Ioutletsupersonic=0
      IoutletspecifiedVelocity=0
      IoutletSpecifiedStaticPressure=0
      IoutletSpecifiedAverageStaticPressure=0
      IoutletSpecifiedResistance=0
      IoutletSpecifiedMassFlowRate=0
      IoutletFullyDeveloped=0
      IoutletFullyDevelopedEnergy=0
      IoutletTurbulence=0
c
      IpressureFarField=0
      IoutletTransmissive=0
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        i=NumberOfScalarsToSolve
c
        allocate(IWallDirichletScalar(i))
        allocate(IWallVonNeumannScalar(i))
        allocate(IWallRobinScalar(i))
        allocate(IinletsupersonicScalar(i))
        allocate(IinletSpecifiedValueScalar(i))
        allocate(IoutletsupersonicScalar(i))
        allocate(IoutletFullyDevelopedScalar(i))

        IWallDirichletScalar=0
        IWallVonNeumannScalar=0
        IWallRobinScalar=0
        IinletsupersonicScalar=0
        IinletSpecifiedValueScalar=0
        IoutletsupersonicScalar=0
        IoutletFullyDevelopedScalar=0
c
      endif
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        i=NumberOfrFieldsToSolve
c
        allocate(IWallDirichletrField(i))
        allocate(IWallVonNeumannrField(i))
        allocate(IWallRobinrField(i))
        allocate(IinletsupersonicrField(i))
        allocate(IinletSpecifiedValuerField(i))
        allocate(IoutletsupersonicrField(i))
        allocate(IoutletFullyDevelopedrField(i))

        IWallDirichletrField=0
        IWallVonNeumannrField=0
        IWallRobinrField=0
        IinletsupersonicrField=0
        IinletSpecifiedValuerField=0
        IoutletsupersonicrField=0
        IoutletFullyDevelopedrField=0
c
      endif
c
      Isymmetry=0
      Iperiodic=0
      Iaxis=0
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          if(BCType(i,j).eq.'wall') then
c
            Iwall=Iwall+1
            if(wallTypeMomentum(i,j).eq.'slip') then
              IWallSlip=IWallSlip+1
            elseif(wallTypeMomentum(i,j).eq.'noslip') then
              IWallnoSlip=IWallnoSlip+1
              if(LTurbulentFlow) IWallTurbulence=IWallTurbulence+1
            endif
            if(wallTypeEnergy(i,j).eq.'dirichlet') then
              IWallDirichlet=IWallDirichlet+1
            elseif(wallTypeEnergy(i,j).eq.'vonneumann') then
              IWallVonNeumann=IWallVonNeumann+1
            elseif(wallTypeEnergy(i,j).eq.'robin') then
              IWallRobin=IWallRobin+1
            endif
            if(wallTypeLambda(i,j).eq.'dirichlet') then
              IWallDirichlet=IWallDirichlet+1
            elseif(wallTypeLambda(i,j).eq.'vonneumann') then
              IWallVonNeumann=IWallVonNeumann+1
            endif
            do k=1,NumberOfScalarsToSolve
              if(LSolveScalar(k)) then
                if(wallTypeScalar(i,j,k).eq.'dirichletscalar') then
                  IWallDirichletScalar(k)=IWallDirichletScalar(k)+1
                elseif(wallTypeScalar(i,j,k).eq.'vonneumannscalar') then
                  IWallVonNeumannScalar(k)=IWallVonNeumannScalar(k)+1
                elseif(wallTypeScalar(i,j,k).eq.'robinscalar') then
                  IWallRobinScalar(k)=IWallRobinScalar(k)+1
                endif
              endif
            enddo
            do k=1,NumberOfrFieldsToSolve
              if(LSolverField(k)) then
                if(wallTyperField(i,j,k).eq.'dirichletrfield') then
                  IWallDirichletrField(k)=IWallDirichletrField(k)+1
                elseif(wallTyperField(i,j,k).eq.'vonneumannrfield') then
                  IWallVonNeumannrField(k)=IWallVonNeumannrField(k)+1
                elseif(wallTyperField(i,j,k).eq.'robinrfield') then
                  IWallRobinrField(k)=IWallRobinrField(k)+1
                endif
              endif
            enddo
c            
          elseif(BCType(i,j).eq.'inlet') then
c
            Iinlet=Iinlet+1
            if(LTurbulentFlow) IinletTurbulence=IinletTurbulence+1
            if(inletTypeMomentum(i,j).EQ.'supersonic') then
              Iinletsupersonic=Iinletsupersonic+1
            elseif(inletTypeMomentum(i,j).EQ.'specifiedvelocity') then
              IinletSpecifiedVelocity=IinletSpecifiedVelocity+1
            elseif(inletTypeMomentum(i,j).EQ.
     *                            'specifiedmassflowrate')then
              IinletSpecifiedMassFlowRate=IinletSpecifiedMassFlowRate+1
            elseif(inletTypeMomentum(i,j).EQ.
     *                   'specifiedstaticpressure')then
              IinletSpecifiedStaticPressure=
     *                              IinletSpecifiedStaticPressure+1
            elseif(inletTypeMomentum(i,j).EQ.
     *                   'specifiedstagnationpressure')then
              IinletSpecifiedStagnationPressure=
     *                              IinletSpecifiedStagnationPressure+1
            endif
c    
            if(inletTypeEnergy(i,j).EQ.'supersonic')then
               if(.not.LsolveMomentum) then
                 Iinletsupersonic=Iinletsupersonic+1
               endif
            elseif(inletTypeEnergy(i,j).EQ.
     *                              'specifiedstatictemperature')then
              IinletSpecifiedStaticTemperature=
     *                              IinletSpecifiedStaticTemperature+1
            elseif(inletTypeEnergy(i,j).EQ.
     *                  'specifiedstagnationtemperature')then
              IinletSpecifiedStagnationTemperature=
     *                        IinletSpecifiedStagnationTemperature+1
            endif
c
            do k=1,NumberOfScalarsToSolve
              if(LSolveScalar(k)) then
                if(inletTypeScalar(i,j,k).eq.'supersonic') then
                  IinletSupersonicScalar(k)=
     *                            IinletSupersonicScalar(k)+1
                elseif(inletTypeScalar(i,j,k).eq.
     *                               'specifiedvaluescalar') then
                  IinletSpecifiedValueScalar(k)=
     *                            IinletSpecifiedValueScalar(k)+1
                endif
              endif
            enddo
c
            do k=1,NumberOfrFieldsToSolve
              if(LSolverField(k)) then
                if(inletTyperField(i,j,k).eq.'supersonic') then
                  IinletSupersonicrField(k)=
     *                            IinletSupersonicrField(k)+1
                elseif(inletTyperField(i,j,k).eq.
     *                               'specifiedvaluerfield') then
                  IinletSpecifiedValuerField(k)=
     *                            IinletSpecifiedValuerField(k)+1
                endif
              endif
            enddo
c
          elseif(BCType(i,j).eq.'outlet') then
c
            Ioutlet=Ioutlet+1
            if(LTurbulentFlow) IoutletTurbulence=IoutletTurbulence+1
            if(outletTypeMomentum(i,j).EQ.'supersonic') then
              Ioutletsupersonic=Ioutletsupersonic+1
            elseif(outletTypeMomentum(i,j).EQ.'specifiedvelocity') then
              IoutletspecifiedVelocity=IoutletspecifiedVelocity+1
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedstaticpressure') then
              IoutletSpecifiedStaticPressure=
     *                      IoutletSpecifiedStaticPressure+1
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedaveragestaticpressure') then
              IoutletSpecifiedAverageStaticPressure=
     *                      IoutletSpecifiedAverageStaticPressure+1
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedresistance') then
              IoutletSpecifiedResistance=
     *                       IoutletSpecifiedResistance+1
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedmassflowrate') then
              IoutletSpecifiedMassFlowRate=
     *                      IoutletSpecifiedMassFlowRate+1
            elseif(outletTypeMomentum(i,j).EQ.'fullydeveloped') then
              IoutletFullyDeveloped=IoutletFullyDeveloped+1
            elseif(outletTypeMomentum(i,j).EQ.'transmissive') then
              IoutletTransmissive=IoutletTransmissive+1
            endif
            if(outletTypeEnergy(i,j).EQ.'supersonic') then
               if(.not.LsolveMomentum) then
                 Ioutletsupersonic=Ioutletsupersonic+1
               endif
            elseif(outletTypeEnergy(i,j).EQ.'fullydevelopedenergy') then
              IoutletFullyDevelopedEnergy=IoutletFullyDevelopedEnergy+1
            elseif(outletTypeEnergy(i,j).EQ.'Transmissive') then
               if(.not.LsolveMomentum) then
                 IoutletTransmissive=IoutletTransmissive+1
               endif
            endif
c
            do k=1,NumberOfScalarsToSolve
              if(LSolveScalar(k)) then
                if(outletTypeScalar(i,j,k).eq.'supersonic') then
                  IoutletSupersonicScalar(k)=
     *                            IoutletSupersonicScalar(k)+1
                elseif(outletTypeScalar(i,j,k).eq.
     *                               'fullydevelopedscalar') then
                  IoutletFullyDevelopedScalar(k)=
     *                            IoutletFullyDevelopedScalar(k)+1
                endif
              endif
            enddo
c
            do k=1,NumberOfrFieldsToSolve
              if(LSolverField(k)) then
                if(outletTyperField(i,j,k).eq.'supersonic') then
                  IoutletSupersonicrField(k)=
     *                            IoutletSupersonicrField(k)+1
                elseif(outletTyperField(i,j,k).eq.
     *                               'fullydevelopedrfield') then
                  IoutletFullyDevelopedrField(k)=
     *                            IoutletFullyDevelopedrField(k)+1
                endif
              endif
            enddo
c 
          elseif(BCType(i,j).eq.'pressurefarfield') then
c
            IpressureFarField=IpressureFarField+1
c
          elseif(BCType(i,j).eq.'symmetry') then
c
            Isymmetry=Isymmetry+1
c
          elseif(BCType(i,j).eq.'periodic') then
c
            Iperiodic=Iperiodic+1
c 
          elseif(BCType(i,j).eq.'axis') then
c
            Iaxis=Iaxis+1
c
          endif
c
        enddo
      enddo
c
c--- Allocate storage for walls
c
      allocate(IWallSlipOwner(IWallSlip))
      allocate(IWallSlipNumberOfBCSets(IWallSlip))
      allocate(IWallSlipNBFaces(IWallSlip))
c
      allocate(IWallnoSlipOwner(IWallnoSlip))
      allocate(IWallnoSlipNumberOfBCSets(IWallnoSlip))
      allocate(IWallnoSlipNBFaces(IWallnoSlip))
c
      allocate(IWallDirichletOwner(IWallDirichlet))
      allocate(IWallDirichletNumberOfBCSets(IWallDirichlet))
      allocate(IWallDirichletNBFaces(IWallDirichlet))
c
      allocate(IWallVonNeumannOwner(IWallVonNeumann))
      allocate(IWallVonNeumannNumberOfBCSets(IWallVonNeumann))
      allocate(IWallVonNeumannNBFaces(IWallVonNeumann))
c
      allocate(IWallRobinOwner(IWallRobin))
      allocate(IWallRobinNumberOfBCSets(IWallRobin))
      allocate(IWallRobinNBFaces(IWallRobin))
      allocate(TinfinityRobin(IWallRobin))
      allocate(HinfinityRobin(IWallRobin))
c
      allocate(IWallTurbulenceOwner(IWallTurbulence))
      allocate(IWallTurbulenceNumberOfBCSets(IWallTurbulence))
      allocate(IWallTurbulenceNBFaces(IWallTurbulence))
      allocate(yplus(IWallTurbulence))
      allocate(yplusT(IWallTurbulence))
      allocate(WallViscosity(IWallTurbulence))
      allocate(WallViscosityT(IWallTurbulence))
      allocate(uplus(IWallTurbulence))
      allocate(ustar(IWallTurbulence))
      allocate(duplusdyplus(IWallTurbulence))
      allocate(TauWall(IWallTurbulence))
      allocate(KsPlus(IWallTurbulence))
      allocate(uTau(IWallTurbulence))
      allocate(ystar(IWallTurbulence))
      allocate(RoughccM(IWallTurbulence))
      allocate(RoughccE(IWallTurbulence))
c
      if(NumberOfScalarsToSolve.gt.0) then
        i=NumberOfScalarsToSolve
        allocate(Roughccs(IWallTurbulence,i))
      endif
c
c---- Initialize wall values
c
      RoughccM=0.
      RoughccE=0.
      KsPlus=0.
      uTau=0.
      ystar=0.
      yplus=1.
      yplusT=1.
      uplus=1.
      ustar=1.
      duplusdyplus=1.
      TauWall=0.
      WallViscosity=0.
      WallViscosityT=0.
c
      if(NumberOfScalarsToSolve.gt.0) then
        Roughccs=0.
      endif
c
      i1=0
      i2=0
      i3=0
c
      do i=1,NumberOfScalarsToSolve
c
        i1=max(i1,IWallDirichletScalar(i))
        i2=max(i2,IWallVonNeumannScalar(i))
        i3=max(i3,IWallRobinScalar(i))
c
      enddo
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        i=NumberOfScalarsToSolve
c
        allocate(IWallDirichletScalarOwner(i1,i))
        allocate(IWallDirichletScalarNumberOfBCSets(i1,i))
        allocate(IWallDirichletScalarNBFaces(i1,i))
c
        allocate(IWallVonNeumannScalarOwner(i2,i))
        allocate(IWallVonNeumannScalarNumberOfBCSets(i2,i))
        allocate(IWallVonNeumannScalarNBFaces(i2,i))
c
        allocate(IWallRobinScalarOwner(i3,i))
        allocate(IWallRobinScalarNumberOfBCSets(i3,i))
        allocate(IWallRobinScalarNBFaces(i3,i))
        allocate(PhiinfinityRobin(i3,i))
        allocate(ConvectionCoefficientRobin(i3,i))
c
        allocate(WallViscosityS(IWallTurbulence,i))
        allocate(yplusS(IWallTurbulence,i))
c
      endif
c
      i1=0
      i2=0
      i3=0
c
      do i=1,NumberOfrFieldsToSolve
c
        i1=max(i1,IWallDirichletrField(i))
        i2=max(i2,IWallVonNeumannrField(i))
        i3=max(i3,IWallRobinrField(i))
c
      enddo
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        i=NumberOfrFieldsToSolve
c
        allocate(IWallDirichletrFieldOwner(i1,i))
        allocate(IWallDirichletrFieldNumberOfBCSets(i1,i))
        allocate(IWallDirichletrFieldNBFaces(i1,i))
c
        allocate(IWallVonNeumannrFieldOwner(i2,i))
        allocate(IWallVonNeumannrFieldNumberOfBCSets(i2,i))
        allocate(IWallVonNeumannrFieldNBFaces(i2,i))
c
        allocate(IWallRobinrFieldOwner(i3,i))
        allocate(IWallRobinrFieldNumberOfBCSets(i3,i))
        allocate(IWallRobinrFieldNBFaces(i3,i))
        allocate(rFieldinfinityRobin(i3,i))
        allocate(rFieldConvectionCoefficientRobin(i3,i))
c
        allocate(WallViscosityR(IWallTurbulence,i))
        allocate(yplusR(IWallTurbulence,i))
c
      endif
c
c--- Allocate storage for inlets
c
      allocate(IinletsupersonicOwner(Iinletsupersonic))
      allocate(IinletsupersonicNumberOfBCSets(Iinletsupersonic))
      allocate(IinletsupersonicNBFaces(Iinletsupersonic))
c
      allocate(IinletSpecifiedVelocityOwner(IinletSpecifiedVelocity))
      allocate(IinletSpecifiedVelocityNumberOfBCSets
     *                               (IinletSpecifiedVelocity))
      allocate(IinletSpecifiedVelocityNBFaces(IinletSpecifiedVelocity))
c
      allocate(IinletSpecifiedMassFlowRateOwner
     *                         (IinletSpecifiedMassFlowRate))
      allocate(IinletSpecifiedMassFlowRateNumberOfBCSets
     *                         (IinletSpecifiedMassFlowRate))
      allocate(IinletSpecifiedMassFlowRateNBFaces
     *                         (IinletSpecifiedMassFlowRate))
c
      allocate(IinletSpecifiedStaticPressureOwner
     *                         (IinletSpecifiedStaticPressure))
      allocate(IinletSpecifiedStaticPressureNumberOfBCSets
     *                         (IinletSpecifiedStaticPressure))
      allocate(IinletSpecifiedStaticPressureNBFaces
     *                         (IinletSpecifiedStaticPressure))
c
      allocate(IinletSpecifiedStagnationPressureOwner
     *                      (IinletSpecifiedStagnationPressure))
      allocate(IinletSpecifiedStagnationPressureNumberOfBCSets
     *                      (IinletSpecifiedStagnationPressure))
      allocate(IinletSpecifiedStagnationPressureNBFaces
     *                      (IinletSpecifiedStagnationPressure))
c
      allocate(IinletSpecifiedStaticTemperatureOwner
     *                         (IinletSpecifiedStaticTemperature))
      allocate(IinletSpecifiedStaticTemperatureNumberOfBCSets
     *                         (IinletSpecifiedStaticTemperature))
      allocate(IinletSpecifiedStaticTemperatureNBFaces
     *                         (IinletSpecifiedStaticTemperature))
c
      allocate(IinletSpecifiedStagnationTemperatureOwner
     *                         (IinletSpecifiedStagnationTemperature))
      allocate(IinletSpecifiedStagnationTemperatureNumberOfBCSets
     *                         (IinletSpecifiedStagnationTemperature))
      allocate(IinletSpecifiedStagnationTemperatureNBFaces
     *                         (IinletSpecifiedStagnationTemperature))
c
      allocate(IinletTurbulenceOwner(IinletTurbulence))
      allocate(IinletTurbulenceNumberOfBCSets(IinletTurbulence))
      allocate(IinletTurbulenceNBFaces(IinletTurbulence))
c
      i1=0
      i2=0
c
      do i=1,NumberOfScalarsToSolve
c
        i1=max(i1,IinletsupersonicScalar(i))
        i2=max(i2,IinletSpecifiedValueScalar(i))
c
      enddo
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        i=NumberOfScalarsToSolve
c
        allocate(IinletsupersonicScalarOwner(i1,i))
        allocate(IinletsupersonicScalarNumberOfBCSets(i1,i))
        allocate(IinletsupersonicScalarNBFaces(i1,i))
c
        allocate(IinletSpecifiedValueScalarOwner(i2,i))
        allocate(IinletSpecifiedValueScalarNumberOfBCSets(i2,i))
        allocate(IinletSpecifiedValueScalarNBFaces(i2,i))
c
      endif
c
      i1=0
      i2=0
c
      do i=1,NumberOfrFieldsToSolve
c
        i1=max(i1,IinletsupersonicrField(i))
        i2=max(i2,IinletSpecifiedValuerField(i))
c
      enddo
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        i=NumberOfrFieldsToSolve
c
        allocate(IinletsupersonicrFieldOwner(i1,i))
        allocate(IinletsupersonicrFieldNumberOfBCSets(i1,i))
        allocate(IinletsupersonicrFieldNBFaces(i1,i))
c
        allocate(IinletSpecifiedValuerFieldOwner(i2,i))
        allocate(IinletSpecifiedValuerFieldNumberOfBCSets(i2,i))
        allocate(IinletSpecifiedValuerFieldNBFaces(i2,i))
c
      endif
cc--- Allocate storage for outlets
c
      allocate(IoutletsupersonicOwner(Ioutletsupersonic))
      allocate(IoutletsupersonicNumberOfBCSets(Ioutletsupersonic))
      allocate(IoutletsupersonicNBFaces(Ioutletsupersonic))
c
      allocate(IoutletspecifiedVelocityOwner(IoutletspecifiedVelocity))
      allocate(IoutletspecifiedVelocityNumberOfBCSets
     *                                   (IoutletspecifiedVelocity))
      allocate(IoutletspecifiedVelocityNBFaces
     *                                   (IoutletspecifiedVelocity))
c
      allocate(IoutletSpecifiedStaticPressureOwner
     *                             (IoutletSpecifiedStaticPressure))
      allocate(IoutletSpecifiedStaticPressureNumberOfBCSets
     *                             (IoutletSpecifiedStaticPressure))
      allocate(IoutletSpecifiedStaticPressureNBFaces
     *                             (IoutletSpecifiedStaticPressure))
c
      allocate(IoutletSpecifiedAverageStaticPressureOwner
     *                         (IoutletSpecifiedAverageStaticPressure))
      allocate(IoutletSpecifiedAverageStaticPressureNumberOfBCSets
     *                         (IoutletSpecifiedAverageStaticPressure))
      allocate(IoutletSpecifiedAverageStaticPressureNBFaces
     *                         (IoutletSpecifiedAverageStaticPressure))
c
      allocate(IoutletSpecifiedResistanceOwner
     *                                    (IoutletSpecifiedResistance))
      allocate(IoutletSpecifiedResistanceNumberOfBCSets
     *                                    (IoutletSpecifiedResistance))
      allocate(IoutletSpecifiedResistanceNBFaces
     *                                    (IoutletSpecifiedResistance))
c
      allocate(IoutletSpecifiedMassFlowRateOwner
     *                             (IoutletSpecifiedMassFlowRate))
      allocate(IoutletSpecifiedMassFlowRateNumberOfBCSets
     *                             (IoutletSpecifiedMassFlowRate))
      allocate(IoutletSpecifiedMassFlowRateNBFaces
     *                             (IoutletSpecifiedMassFlowRate))
c
      allocate(IoutletFullyDevelopedOwner(IoutletFullyDeveloped))
      allocate(IoutletFullyDevelopedNumberOfBCSets
     *                             (IoutletFullyDeveloped))
      allocate(IoutletFullyDevelopedNBFaces(IoutletFullyDeveloped))
c
      allocate(IoutletFullyDevelopedEnergyOwner
     *                             (IoutletFullyDevelopedEnergy))
      allocate(IoutletFullyDevelopedEnergyNumberOfBCSets
     *                             (IoutletFullyDevelopedEnergy))
      allocate(IoutletFullyDevelopedEnergyNBFaces
     *                             (IoutletFullyDevelopedEnergy))
c
      allocate(IoutletTransmissiveOwner(IoutletTransmissive))
      allocate(IoutletTransmissiveNumberOfBCSets(IoutletTransmissive))
      allocate(IoutletTransmissiveNBFaces(IoutletTransmissive))
c
      allocate(IoutletTurbulenceOwner(IoutletTurbulence))
      allocate(IoutletTurbulenceNumberOfBCSets(IoutletTurbulence))
      allocate(IoutletTurbulenceNBFaces(IoutletTurbulence))
c
      i1=0
      i2=0
c
      do i=1,NumberOfScalarsToSolve
c
        i1=max(i1,IoutletsupersonicScalar(i))
        i2=max(i2,IoutletFullyDevelopedScalar(i))
c
      enddo
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        i=NumberOfScalarsToSolve
c
        allocate(IoutletsupersonicScalarOwner(i1,i))
        allocate(IoutletsupersonicScalarNumberOfBCSets(i1,i))
        allocate(IoutletsupersonicScalarNBFaces(i1,i))
c
        allocate(IoutletFullyDevelopedScalarOwner(i2,i))
        allocate(IoutletFullyDevelopedScalarNumberOfBCSets(i2,i))
        allocate(IoutletFullyDevelopedScalarNBFaces(i2,i))
c
      endif
c
      i1=0
      i2=0
c
      do i=1,NumberOfrFieldsToSolve
c
        i1=max(i1,IoutletsupersonicrField(i))
        i2=max(i2,IoutletFullyDevelopedrField(i))
c
      enddo
c
      if(NumberOfrFieldsToSolve.gt.0) then
c
        i=NumberOfrFieldsToSolve
c
        allocate(IoutletsupersonicrFieldOwner(i1,i))
        allocate(IoutletsupersonicrFieldNumberOfBCSets(i1,i))
        allocate(IoutletsupersonicrFieldNBFaces(i1,i))
c
        allocate(IoutletFullyDevelopedrFieldOwner(i2,i))
        allocate(IoutletFullyDevelopedrFieldNumberOfBCSets(i2,i))
        allocate(IoutletFullyDevelopedrFieldNBFaces(i2,i))
c
      endif
c
c--- Allocate storage for pressure far field
c
      allocate(IpressurefarfieldOwner(Ipressurefarfield))
      allocate(IpressurefarfieldNumberOfBCSets(Ipressurefarfield))
      allocate(IpressurefarfieldNBFaces(Ipressurefarfield))
c
c--- Allocate storage for symmetry planes
c
      allocate(IsymmetryOwner(Isymmetry))
      allocate(IsymmetryNumberOfBCSets(Isymmetry))
      allocate(IsymmetryNBFaces(Isymmetry))
c
c--- Allocate storage for periodic boundaries
c
      allocate(IperiodicOwner(Iperiodic))
      allocate(IperiodicNumberOfBCSets(Iperiodic))
      allocate(IperiodicNBFaces(Iperiodic))
c
c--- Allocate storage for Axis boundaries (Axisymmetric)
c
      allocate(IaxisOwner(Iaxis))
      allocate(IaxisNumberOfBCSets(Iaxis))
      allocate(IaxisNBFaces(Iaxis))
c
c--- Assign correct parameter values along different boundaries
c
      IWallSlip=0
      IWallnoSlip=0
      IWallDirichlet=0
      IWallVonNeumann=0
      IWallRobin=0
      IWallDirichletScalar=0
      IWallVonNeumannScalar=0
      IWallRobinScalar=0
      IWallDirichletrField=0
      IWallVonNeumannrField=0
      IWallRobinrField=0
      IwallTurbulence=0
c      
      Iinletsupersonic=0
      IinletSpecifiedVelocity=0
      IinletSpecifiedMassFlowRate=0
      IinletSpecifiedStaticPressure=0
      IinletSpecifiedStagnationPressure=0
      IinletSpecifiedStaticTemperature=0
      IinletSpecifiedStagnationTemperature=0
      IinletsupersonicScalar=0
      IinletSpecifiedValueScalar=0
      IinletsupersonicrField=0
      IinletSpecifiedValuerField=0
      IinletTurbulence=0
c
      Ioutletsupersonic=0
      IoutletspecifiedVelocity=0
      IoutletSpecifiedStaticPressure=0
      IoutletSpecifiedAverageStaticPressure=0
      IoutletSpecifiedResistance=0
      IoutletSpecifiedMassFlowRate=0
      IoutletFullyDeveloped=0
      IoutletFullyDevelopedEnergy=0
      IoutletsupersonicScalar=0
      IoutletFullyDevelopedScalar=0
      IoutletsupersonicrField=0
      IoutletFullyDevelopedrField=0
      IoutletTurbulence=0
c
      IoutletTransmissive=0
      IpressureFarField=0
c
      Isymmetry=0
      Iperiodic=0
      Iaxis=0
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          if(BCType(i,j).eq.'wall') then
c
            if(wallTypeMomentum(i,j).eq.'slip') then
c
              IWallSlip=IWallSlip+1
c
              IWallSlipOwner(IWallSlip)=NBFaceOwner(i,j)
              IWallSlipNumberOfBCSets(IWallSlip)=i
              IWallSlipNBFaces(IWallSlip)=j
c
            elseif(wallTypeMomentum(i,j).eq.'noslip') then
c
              IWallnoSlip=IWallnoSlip+1
c
              IWallnoSlipOwner(IWallnoSlip)=NBFaceOwner(i,j)
              IWallnoSlipNumberOfBCSets(IWallnoSlip)=i
              IWallnoSlipNBFaces(IWallnoSlip)=j
c           
              if(LTurbulentFlow) then
c
                IwallTurbulence=IwallTurbulence+1
c
                IwallTurbulenceOwner(IwallTurbulence)=NBFaceOwner(i,j)
                IwallTurbulenceNumberOfBCSets(IwallTurbulence)=i
                IwallTurbulenceNBFaces(IwallTurbulence)=j
c
              endif
c 
            endif
c
            if(wallTypeEnergy(i,j).eq.'dirichlet') then
c
              IWallDirichlet=IWallDirichlet+1
c
              IWallDirichletOwner(IWallDirichlet)=NBFaceOwner(i,j)
              IWallDirichletNumberOfBCSets(IWallDirichlet)=i
              IWallDirichletNBFaces(IWallDirichlet)=j
c
           elseif(wallTypeEnergy(i,j).eq.'vonneumann') then
c
              IWallVonNeumann=IWallVonNeumann+1
c
              IWallVonNeumannOwner(IWallVonNeumann)=NBFaceOwner(i,j)
              IWallVonNeumannNumberOfBCSets(IWallVonNeumann)=i
              IWallVonNeumannNBFaces(IWallVonNeumann)=j
c
            elseif(wallTypeEnergy(i,j).eq.'robin') then
c
              IWallRobin=IWallRobin+1
c
              IWallRobinOwner(IWallRobin)=NBFaceOwner(i,j)
              IWallRobinNumberOfBCSets(IWallRobin)=i
              IWallRobinNBFaces(IWallRobin)=j
              HinfinityRobin(IWallRobin)=Hinfinity(i,j)
              TinfinityRobin(IwallRobin)=Tinfinity(i,j)
c
            endif
c
            if(wallTypeLambda(i,j).eq.'dirichlet') then
c
              IWallDirichlet=IWallDirichlet+1
c
              IWallDirichletOwner(IWallDirichlet)=NBFaceOwner(i,j)
              IWallDirichletNumberOfBCSets(IWallDirichlet)=i
              IWallDirichletNBFaces(IWallDirichlet)=j
c
           elseif(wallTypeLambda(i,j).eq.'vonneumann') then
c
              IWallVonNeumann=IWallVonNeumann+1
c
              IWallVonNeumannOwner(IWallVonNeumann)=NBFaceOwner(i,j)
              IWallVonNeumannNumberOfBCSets(IWallVonNeumann)=i
              IWallVonNeumannNBFaces(IWallVonNeumann)=j
c
            endif
c
            do k=1,NumberOfScalarsToSolve
c
              if(LSolveScalar(k)) then
c
                if(wallTypeScalar(i,j,k).eq.'dirichletscalar') then
c
                  IWallDirichletScalar(k)=IWallDirichletScalar(k)+1
c
                  IWallDirichletScalarOwner
     *                   (IWallDirichletScalar(k),k)=NBFaceOwner(i,j)
                  IWallDirichletScalarNumberOfBCSets
     *                                 (IWallDirichletScalar(k),k)=i
                  IWallDirichletScalarNBFaces
     *                                 (IWallDirichletScalar(k),k)=j
c
                elseif(wallTypeScalar(i,j,k).eq.'vonneumannscalar') then
c
                  IWallVonNeumannScalar(k)=IWallVonNeumannScalar(k)+1
c
                  IWallVonNeumannScalarOwner
     *                   (IWallVonNeumannScalar(k),k)=NBFaceOwner(i,j)
                  IWallVonNeumannScalarNumberOfBCSets
     *                                 (IWallVonNeumannScalar(k),k)=i
                  IWallVonNeumannScalarNBFaces
     *                                 (IWallVonNeumannScalar(k),k)=j
c
                elseif(wallTypeScalar(i,j,k).eq.'robinscalar') then
c
                  IWallRobinScalar(k)=IWallRobinScalar(k)+1
c
                  IWallRobinScalarOwner
     *                   (IWallRobinScalar(k),k)=NBFaceOwner(i,j)
                  IWallRobinScalarNumberOfBCSets
     *                                 (IWallRobinScalar(k),k)=i
                  IWallRobinScalarNBFaces(IWallRobinScalar(k),k)=j
c
                endif
c
              endif
c
            enddo
c
            do k=1,NumberOfrFieldsToSolve
c
              if(LSolverField(k)) then
c
                if(wallTyperField(i,j,k).eq.'dirichletrfield') then
c
                  IWallDirichletrField(k)=IWallDirichletrField(k)+1
c
                  IWallDirichletrFieldOwner
     *                   (IWallDirichletrField(k),k)=NBFaceOwner(i,j)
                  IWallDirichletrFieldNumberOfBCSets
     *                                 (IWallDirichletrField(k),k)=i
                  IWallDirichletrFieldNBFaces
     *                                 (IWallDirichletrField(k),k)=j
c
                elseif(wallTyperField(i,j,k).eq.'vonneumannrfield') then
c
                  IWallVonNeumannrField(k)=IWallVonNeumannrField(k)+1
c
                  IWallVonNeumannrFieldOwner
     *                   (IWallVonNeumannrField(k),k)=NBFaceOwner(i,j)
                  IWallVonNeumannrFieldNumberOfBCSets
     *                                 (IWallVonNeumannrField(k),k)=i
                  IWallVonNeumannrFieldNBFaces
     *                                 (IWallVonNeumannrField(k),k)=j
c
                elseif(wallTyperField(i,j,k).eq.'robinrfield') then
c
                  IWallRobinrField(k)=IWallRobinrField(k)+1
c
                  IWallRobinrFieldOwner
     *                   (IWallRobinrField(k),k)=NBFaceOwner(i,j)
                  IWallRobinrFieldNumberOfBCSets
     *                                 (IWallRobinrField(k),k)=i
                  IWallRobinrFieldNBFaces(IWallRobinrField(k),k)=j
c
                endif
c
              endif
c
            enddo
c
          elseif(BCType(i,j).eq.'inlet') then
c
            if(inletTypeMomentum(i,j).EQ.'supersonic') then
c
              Iinletsupersonic=Iinletsupersonic+1
c
              IinletsupersonicOwner(Iinletsupersonic)=NBFaceOwner(i,j)
              IinletsupersonicNumberOfBCSets(Iinletsupersonic)=i
              IinletsupersonicNBFaces(Iinletsupersonic)=j
c
            elseif(inletTypeMomentum(i,j).EQ.'specifiedvelocity') then
c
              IinletSpecifiedVelocity=IinletSpecifiedVelocity+1
c
              IinletSpecifiedVelocityOwner(IinletSpecifiedVelocity)=
     *                                                  NBFaceOwner(i,j)
              IinletSpecifiedVelocityNumberOfBCSets
     *                                 (IinletSpecifiedVelocity)=i
              IinletSpecifiedVelocityNBFaces(IinletSpecifiedVelocity)=j
c
            elseif(inletTypeMomentum(i,j).EQ.
     *                            'specifiedmassflowrate')then
c
              IinletSpecifiedMassFlowRate=IinletSpecifiedMassFlowRate+1
c
              IinletSpecifiedMassFlowRateOwner
     *                  (IinletSpecifiedMassFlowRate)=NBFaceOwner(i,j)
              IinletSpecifiedMassFlowRateNumberOfBCSets
     *                                 (IinletSpecifiedMassFlowRate)=i
              IinletSpecifiedMassFlowRateNBFaces
     *                                 (IinletSpecifiedMassFlowRate)=j
c
            elseif(inletTypeMomentum(i,j).EQ.
     *                   'specifiedstaticpressure')then
c
              IinletSpecifiedStaticPressure=
     *                              IinletSpecifiedStaticPressure+1
c
              IinletSpecifiedStaticPressureOwner
     *                  (IinletSpecifiedStaticPressure)=NBFaceOwner(i,j)
              IinletSpecifiedStaticPressureNumberOfBCSets
     *                                 (IinletSpecifiedStaticPressure)=i
              IinletSpecifiedStaticPressureNBFaces
     *                                 (IinletSpecifiedStaticPressure)=j
c
            elseif(inletTypeMomentum(i,j).EQ.
     *                   'specifiedstagnationpressure')then
c
              IinletSpecifiedStagnationPressure=
     *                              IinletSpecifiedStagnationPressure+1
c
              IinletSpecifiedStagnationPressureOwner
     *            (IinletSpecifiedStagnationPressure)=NBFaceOwner(i,j)
              IinletSpecifiedStagnationPressureNumberOfBCSets
     *                          (IinletSpecifiedStagnationPressure)=i
              IinletSpecifiedStagnationPressureNBFaces
     *                          (IinletSpecifiedStagnationPressure)=j
c
            endif
c           
            if(LTurbulentFlow) then
c
              IinletTurbulence=IinletTurbulence+1
c
              IinletTurbulenceOwner(IinletTurbulence)=NBFaceOwner(i,j)
              IinletTurbulenceNumberOfBCSets(IinletTurbulence)=i
              IinletTurbulenceNBFaces(IinletTurbulence)=j
c
            endif
c 
            if(inletTypeEnergy(i,j).EQ.'supersonic')then
              if(.not.LsolveMomentum) then
c
                Iinletsupersonic=Iinletsupersonic+1
c
                IinletsupersonicOwner(Iinletsupersonic)=NBFaceOwner(i,j)
                IinletsupersonicNumberOfBCSets(Iinletsupersonic)=i
                IinletsupersonicNBFaces(Iinletsupersonic)=j
c
              endif            
            elseif(inletTypeEnergy(i,j).EQ.
     *                                 'specifiedstatictemperature')then
c
              IinletSpecifiedStaticTemperature=
     *                              IinletSpecifiedStaticTemperature+1
c
              IinletSpecifiedStaticTemperatureOwner
     *               (IinletSpecifiedStaticTemperature)=NBFaceOwner(i,j)
              IinletSpecifiedStaticTemperatureNumberOfBCSets
     *                              (IinletSpecifiedStaticTemperature)=i
              IinletSpecifiedStaticTemperatureNBFaces
     *                              (IinletSpecifiedStaticTemperature)=j
c
            elseif(inletTypeEnergy(i,j).EQ.
     *                  'specifiedstagnationtemperature')then
c
              IinletSpecifiedStagnationTemperature=
     *                        IinletSpecifiedStagnationTemperature+1
c
              IinletSpecifiedStagnationTemperatureOwner
     *           (IinletSpecifiedStagnationTemperature)=NBFaceOwner(i,j)
              IinletSpecifiedStagnationTemperatureNumberOfBCSets
     *                          (IinletSpecifiedStagnationTemperature)=i
              IinletSpecifiedStagnationTemperatureNBFaces
     *                          (IinletSpecifiedStagnationTemperature)=j
c
            endif
c            
            do k=1,NumberOfScalarsToSolve
c
              if(LSolveScalar(k)) then
c
                if(inletTypeScalar(i,j,k).eq.'specifiedvaluescalar')then
c
                  IinletSpecifiedValueScalar(k)=
     *                             IinletSpecifiedValueScalar(k)+1
c
                  IinletSpecifiedValueScalarOwner
     *               (IinletSpecifiedValueScalar(k),k)=NBFaceOwner(i,j)
                  IinletSpecifiedValueScalarNumberOfBCSets
     *                              (IinletSpecifiedValueScalar(k),k)=i
                  IinletSpecifiedValueScalarNBFaces
     *                              (IinletSpecifiedValueScalar(k),k)=j
c
                elseif(inletTypeScalar(i,j,k).eq.'supersonic') then
c
                  IinletSupersonicScalar(k)=IinletSupersonicScalar(k)+1
c
                  IinletSupersonicScalarOwner
     *                   (IinletSupersonicScalar(k),k)=NBFaceOwner(i,j)
                  IinletSupersonicScalarNumberOfBCSets
     *                                 (IinletSupersonicScalar(k),k)=i
                  IinletSupersonicScalarNBFaces
     *                                 (IinletSupersonicScalar(k),k)=j
c
                endif
c
              endif
c
            enddo
c            
            do k=1,NumberOfrFieldsToSolve
c
              if(LSolverField(k)) then
c
                if(inletTyperField(i,j,k).eq.'specifiedvaluerfield')then
c
                  IinletSpecifiedValuerField(k)=
     *                             IinletSpecifiedValuerField(k)+1
c
                  IinletSpecifiedValuerFieldOwner
     *               (IinletSpecifiedValuerField(k),k)=NBFaceOwner(i,j)
                  IinletSpecifiedValuerFieldNumberOfBCSets
     *                              (IinletSpecifiedValuerField(k),k)=i
                  IinletSpecifiedValuerFieldNBFaces
     *                              (IinletSpecifiedValuerField(k),k)=j
c
                elseif(inletTyperField(i,j,k).eq.'supersonic') then
c
                  IinletSupersonicrField(k)=IinletSupersonicrField(k)+1
c
                  IinletSupersonicrFieldOwner
     *                   (IinletSupersonicrField(k),k)=NBFaceOwner(i,j)
                  IinletSupersonicrFieldNumberOfBCSets
     *                                 (IinletSupersonicrField(k),k)=i
                  IinletSupersonicrFieldNBFaces
     *                                 (IinletSupersonicrField(k),k)=j
c
                endif
c
              endif
c
            enddo
c
          elseif(BCType(i,j).eq.'outlet') then
c
            if(outletTypeMomentum(i,j).EQ.'supersonic') then
c
              Ioutletsupersonic=Ioutletsupersonic+1
c
              IoutletsupersonicOwner(Ioutletsupersonic)=NBFaceOwner(i,j)
              IoutletsupersonicNumberOfBCSets(Ioutletsupersonic)=i
              IoutletsupersonicNBFaces(Ioutletsupersonic)=j
c
            elseif(outletTypeMomentum(i,j).EQ.'specifiedvelocity') then
c
              IoutletspecifiedVelocity=IoutletspecifiedVelocity+1
c
              IoutletspecifiedVelocityOwner
     *                       (IoutletspecifiedVelocity)=NBFaceOwner(i,j)
              IoutletspecifiedVelocityNumberOfBCSets
     *                                 (IoutletspecifiedVelocity)=i
              IoutletspecifiedVelocityNBFaces
     *                                  (IoutletspecifiedVelocity)=j
c
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedstaticpressure') then
c
              IoutletSpecifiedStaticPressure=
     *                      IoutletSpecifiedStaticPressure+1
c
              IoutletSpecifiedStaticPressureOwner
     *                (IoutletSpecifiedStaticPressure)=NBFaceOwner(i,j)
              IoutletSpecifiedStaticPressureNumberOfBCSets
     *                               (IoutletSpecifiedStaticPressure)=i
              IoutletSpecifiedStaticPressureNBFaces
     *                               (IoutletSpecifiedStaticPressure)=j
c
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedaveragestaticpressure') then
c
              IoutletSpecifiedAverageStaticPressure=
     *                      IoutletSpecifiedAverageStaticPressure+1
c
              IoutletSpecifiedAverageStaticPressureOwner
     *          (IoutletSpecifiedAverageStaticPressure)=NBFaceOwner(i,j)
              IoutletSpecifiedAverageStaticPressureNumberOfBCSets
     *                         (IoutletSpecifiedAverageStaticPressure)=i
              IoutletSpecifiedAverageStaticPressureNBFaces
     *                         (IoutletSpecifiedAverageStaticPressure)=j
c
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedresistance') then
c
              IoutletSpecifiedResistance=
     *                      IoutletSpecifiedResistance+1
c
              IoutletSpecifiedResistanceOwner
     *                (IoutletSpecifiedResistance)=NBFaceOwner(i,j)
              IoutletSpecifiedResistanceNumberOfBCSets
     *                               (IoutletSpecifiedResistance)=i
              IoutletSpecifiedResistanceNBFaces
     *                               (IoutletSpecifiedResistance)=j
c
            elseif(outletTypeMomentum(i,j).EQ.
     *                      'specifiedmassflowrate') then
c
              IoutletSpecifiedMassFlowRate=
     *                      IoutletSpecifiedMassFlowRate+1
c
              IoutletSpecifiedMassFlowRateOwner
     *                (IoutletSpecifiedMassFlowRate)=NBFaceOwner(i,j)
              IoutletSpecifiedMassFlowRateNumberOfBCSets
     *                               (IoutletSpecifiedMassFlowRate)=i
              IoutletSpecifiedMassFlowRateNBFaces
     *                               (IoutletSpecifiedMassFlowRate)=j
c
            elseif(outletTypeMomentum(i,j).EQ.'fullydeveloped') then
c
              IoutletFullyDeveloped=IoutletFullyDeveloped+1
c
              IoutletFullyDevelopedOwner
     *                (IoutletFullyDeveloped)=NBFaceOwner(i,j)
              IoutletFullyDevelopedNumberOfBCSets
     *                               (IoutletFullyDeveloped)=i
              IoutletFullyDevelopedNBFaces
     *                               (IoutletFullyDeveloped)=j
c
            elseif(outletTypeMomentum(i,j).EQ.'transmissive') then
c
              IoutletTransmissive=IoutletTransmissive+1
c
              IoutletTransmissiveOwner(IoutletTransmissive)=
     *                                                NBFaceOwner(i,j)
              IoutletTransmissiveNumberOfBCSets(IoutletTransmissive)=i
              IoutletTransmissiveNBFaces(IoutletTransmissive)=j
c
           endif
c           
            if(LTurbulentFlow) then
c
              IoutletTurbulence=IoutletTurbulence+1
c
              IoutletTurbulenceOwner(IoutletTurbulence)=NBFaceOwner(i,j)
              IoutletTurbulenceNumberOfBCSets(IoutletTurbulence)=i
              IoutletTurbulenceNBFaces(IoutletTurbulence)=j
c
            endif
c 
            if(outletTypeEnergy(i,j).EQ.'supersonic') then
c
              if(.not.LSolveMomentum) then
                Ioutletsupersonic=Ioutletsupersonic+1
c
                IoutletsupersonicOwner(Ioutletsupersonic)=
     *                                             NBFaceOwner(i,j)
                IoutletsupersonicNumberOfBCSets(Ioutletsupersonic)=i
                IoutletsupersonicNBFaces(Ioutletsupersonic)=j
              endif
c
            elseif(outletTypeEnergy(i,j).EQ.'fullydevelopedenergy') then
c
              IoutletFullyDevelopedEnergy=IoutletFullyDevelopedEnergy+1
c
              IoutletFullyDevelopedEnergyOwner
     *                (IoutletFullyDevelopedEnergy)=NBFaceOwner(i,j)
              IoutletFullyDevelopedEnergyNumberOfBCSets
     *                               (IoutletFullyDevelopedEnergy)=i
              IoutletFullyDevelopedEnergyNBFaces
     *                               (IoutletFullyDevelopedEnergy)=j
c
            elseif(outletTypeEnergy(i,j).EQ.'transmissive') then
c
              if(.not.LsolveMomentum) then
                IoutletTransmissive=IoutletTransmissive+1
c
                IoutletTransmissiveOwner(IoutletTransmissive)=
     *                                                NBFaceOwner(i,j)
                IoutletTransmissiveNumberOfBCSets(IoutletTransmissive)=i
                IoutletTransmissiveNBFaces(IoutletTransmissive)=j
              endif
c
            endif
c            
            do k=1,NumberOfScalarsToSolve
c
              if(LSolveScalar(k)) then
c
                if(outletTypeScalar(i,j,k).eq.
     *                                  'fullydevelopedscalar') then
c
                  IoutletFullyDevelopedScalar(k)=
     *                             IoutletFullyDevelopedScalar(k)+1
c
                  IoutletFullyDevelopedScalarOwner
     *               (IoutletFullyDevelopedScalar(k),k)=NBFaceOwner(i,j)
                  IoutletFullyDevelopedScalarNumberOfBCSets
     *                              (IoutletFullyDevelopedScalar(k),k)=i
                  IoutletFullyDevelopedScalarNBFaces
     *                              (IoutletFullyDevelopedScalar(k),k)=j
c
                elseif(outletTypeScalar(i,j,k).eq.'supersonic') then
c
                  IoutletSupersonicScalar(k)=
     *                                 IoutletSupersonicScalar(k)+1
c
                  IoutletSupersonicScalarOwner
     *                   (IoutletSupersonicScalar(k),k)=NBFaceOwner(i,j)
                  IoutletSupersonicScalarNumberOfBCSets
     *                                 (IoutletSupersonicScalar(k),k)=i
                  IoutletSupersonicScalarNBFaces
     *                                 (IoutletSupersonicScalar(k),k)=j
c
                endif
c
              endif
c
            enddo
c            
            do k=1,NumberOfrFieldsToSolve
c
              if(LSolverField(k)) then
c
                if(outletTyperField(i,j,k).eq.
     *                                  'fullydevelopedrfield') then
c
                  IoutletFullyDevelopedrField(k)=
     *                             IoutletFullyDevelopedrField(k)+1
c
                  IoutletFullyDevelopedrFieldOwner
     *               (IoutletFullyDevelopedrField(k),k)=NBFaceOwner(i,j)
                  IoutletFullyDevelopedrFieldNumberOfBCSets
     *                              (IoutletFullyDevelopedrField(k),k)=i
                  IoutletFullyDevelopedrFieldNBFaces
     *                              (IoutletFullyDevelopedrField(k),k)=j
c
                elseif(outletTyperField(i,j,k).eq.'supersonic') then
c
                  IoutletSupersonicrField(k)=
     *                                 IoutletSupersonicrField(k)+1
c
                  IoutletSupersonicrFieldOwner
     *                   (IoutletSupersonicrField(k),k)=NBFaceOwner(i,j)
                  IoutletSupersonicrFieldNumberOfBCSets
     *                                 (IoutletSupersonicrField(k),k)=i
                  IoutletSupersonicrFieldNBFaces
     *                                 (IoutletSupersonicrField(k),k)=j
c
                endif
c
              endif
c
            enddo
c
          elseif(BCType(i,j).eq.'pressurefarfield') then
c
            IpressureFarField=IpressureFarField+1
            IpressureFarFieldOwner(IpressureFarField)=NBFaceOwner(i,j)
            IpressureFarFieldNumberOfBCSets(IpressureFarField)=i
            IpressureFarFieldNBFaces(IpressureFarField)=j
c
          elseif(BCType(i,j).eq.'symmetry') then
c
            Isymmetry=Isymmetry+1
            IsymmetryOwner(Isymmetry)=NBFaceOwner(i,j)
            IsymmetryNumberOfBCSets(Isymmetry)=i
            IsymmetryNBFaces(Isymmetry)=j
c
          elseif(BCType(i,j).eq.'periodic') then
c
            Iperiodic=Iperiodic+1
            IperiodicOwner(Iperiodic)=NBFaceOwner(i,j)
            IperiodicNumberOfBCSets(Iperiodic)=i
            IperiodicNBFaces(Iperiodic)=j
c
          elseif(BCType(i,j).eq.'axis') then
c
            Iaxis=Iaxis+1
            IaxisOwner(Iaxis)=NBFaceOwner(i,j)
            IaxisNumberOfBCSets(Iaxis)=i
            IaxisNBFaces(Iaxis)=j
c
          endif
c
        enddo
      enddo
c
      return
      end
