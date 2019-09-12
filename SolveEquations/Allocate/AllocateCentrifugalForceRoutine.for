c
C#############################################################################################
      SUBROUTINE allocateCentrifugalForce
C#############################################################################################
      use User0, only: LCoriolis
      use Geometry1, only:NumberOfElements
      use Centrifugal1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
c
      if(LCoriolis) then
c
        allocate(xCentrifugal(NumberOfElements))
        allocate(yCentrifugal(NumberOfElements))
        allocate(zCentrifugal(NumberOfElements))
c
        call CalculateCentrifugalForce
c
      endif
c
      return
      end
