c
C********************************************************************************************
      SUBROUTINE CalculateLengthScale
C********************************************************************************************
c       
      use User0
      use Turbulence1, only: LengthScale
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use constants1, only: oneThird
C********************************************************************************************
      implicit none
C********************************************************************************************
      integer :: i
C********************************************************************************************
      LengthScale=0.
      do i=1,NumberOfElements
        LengthScale=LengthScale+Volume(i)
      enddo
      LengthScale=(LengthScale)**oneThird
c
      return
      end