c
c#############################################################################################
c
      SUBROUTINE CalculateEffectsOfTurbulence
c
C#############################################################################################
      use User0, only: LCompressibilityTemperatureCorrection,
     *                 WallTreatment
      use Turbulence1, only: ModelNumber
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      call CalculateStrainRateTensor
      call CalculateVorticityTensor
      if(WallTreatment.eq.'lowreynoldsnumber') 
     *                          call CalculateWallTurbulence
c
      SolveEquations:select case (ModelNumber)
c-----------------------------------------------------------------------------------------------      
        case(1) SolveEquations           !kepsilon
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(2) SolveEquations           !kepsilonchien
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateKEChienCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(3) SolveEquations           !kepsilonsharma
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateKELaunderSharmaCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(4) SolveEquations           !kepsilonchc
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateKEchcCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(5) SolveEquations           !kepsilonkasagi
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateKEkasagiCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(6) SolveEquations           !kepsilontagawa
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateKEtagawaCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(7) SolveEquations           !kepsilonhishida
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateKEhishidaCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(8) SolveEquations           !kelambremhorst
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculatekelambremhorstCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(9) SolveEquations           !kelambremhorstm
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculatekelambremhorstmCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(10) SolveEquations           !realizable
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateRealizableCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(11) SolveEquations           !komega
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(12) SolveEquations           !komegaepsilon
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(13) SolveEquations           !komegabsl
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(14) SolveEquations           !komegasst
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(15) SolveEquations           !sstgamaretheta
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          call CalculateTurbulentViscosity
          call SolveTurbulenceGammaEquation             ! Intermittency equation
          call SolveTurbulenceReynoldsThetaEquation
c
c-----------------------------------------------------------------------------------------------      
        case(16) SolveEquations           !komega2006
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(17) SolveEquations           !komega2006lrn
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          call CalculateKOmega2006Coefficients
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(18) SolveEquations           !kklmodel
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceKLequation
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(19) SolveEquations           !spalartallmaras
c-----------------------------------------------------------------------------------------------      
c
          call SolveSpalartEddyViscosityModel
          call CalculateSpalartAllmarasCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(20) SolveEquations           !wrayagarwal
c-----------------------------------------------------------------------------------------------      
c
          call SolveWrayAgarwalEddyViscosityModel
          call CalculateWrayAgarwalCoefficients
          if(.not.LCompressibilityTemperatureCorrection) then
            call CalculateTurbulentViscosity
          else
            call ModifiyTViscosityByTemperatureCorrection
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(21) SolveEquations           !kklomega
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceKLequation
          call SolveTurbulenceSpecificDissipationRate
          call CalculateKKLWCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(22) SolveEquations           !kepsilonrt
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call SolvePseudoEddyViscosityEquation
          call CalculatePseudoEddyViscosityCoefficients
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(23) SolveEquations           !sstgama
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceSpecificDissipationRate
          call SolveTurbulenceGammaEquation             ! Intermittency equation
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
        case(24) SolveEquations           !nut92
c-----------------------------------------------------------------------------------------------      
c
          call Solvenut92EddyViscosityModel
          call CalculateTurbulentViscosity
c-----------------------------------------------------------------------------------------------      
        case(25) SolveEquations           !kepsilonrng
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateRNGCoefficients
          call CalculateTurbulentViscosity
c-----------------------------------------------------------------------------------------------      
        case(26) SolveEquations           !kepsilonv2f
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateV2fCoefficients
          call SolveTurbulenceV2fEquation
          call SolveTurbulencefRelaxationEquation
          call CalculateTurbulentViscosity
c-----------------------------------------------------------------------------------------------      
        case(27) SolveEquations           !kepsilonzetaf
c-----------------------------------------------------------------------------------------------      
c
          call SolveTurbulenceKineticEnergy
          call SolveTurbulenceDissipationRate
          call CalculateZetafCoefficients
          call SolveTurbulenceZetaEquation
          call SolveTurbulencefRelaxationEquation
          call CalculateTurbulentViscosity
c
c-----------------------------------------------------------------------------------------------      
      end select SolveEquations 
c-----------------------------------------------------------------------------------------------      
c
      return
      end