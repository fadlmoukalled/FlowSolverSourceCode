c
c#############################################################################################
c
      SUBROUTINE SolveEnergy
c
C#############################################################################################
      use User0, only: EnergyEquation
c********************************************************************************************
      implicit none
c********************************************************************************************
c
      if(EnergyEquation.eq.'temperature') then
c
        call SolveTemperature
c
      elseif(EnergyEquation.eq.'htotal') then
c
        call SolveTotalEnthalpy
c
      endif
c
      return
      end
