      MODULE ArteryResistance1
      Implicit none
!
      logical, save, dimension(:), allocatable :: LArteryExplicit
      double precision, save, dimension(:), allocatable :: urfPressureResistance
      double precision, save, dimension(:), allocatable :: ArteryResistance  
      double precision, save, dimension(:), allocatable :: geoDiffSum
      double precision, save, dimension(:,:), allocatable :: geoDiffB
      double precision, save, dimension(:), allocatable :: PressureCResistance
!
      end MODULE ArteryResistance1