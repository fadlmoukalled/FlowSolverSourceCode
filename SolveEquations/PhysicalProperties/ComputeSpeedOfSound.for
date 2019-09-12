c
c#############################################################################################
c
      FUNCTION SpeedOfSound(i)
c
C#############################################################################################
c
      use User0, only: BulkModulus
      use PhysicalProperties1, only: GammaGas,RGas,Density,
     *                               EquationOfState

      use Variables1, only: Temperature
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: SpeedOfSound
c********************************************************************************************
c
        if(EquationOfState.eq.'idealgas') then
c
          SpeedOfSound=dsqrt(GammaGas*RGas*Temperature(i))
c
        elseif(EquationOfState.eq.'constant') then
c
          SpeedOfSound=dsqrt(BulkModulus/Density(i))
c
        elseif(EquationOfState.eq.'tait') then
c
          SpeedOfSound=dsqrt(BulkModulus/Density(i))
c
        endif
c
      end FUNCTION SpeedOfSound
