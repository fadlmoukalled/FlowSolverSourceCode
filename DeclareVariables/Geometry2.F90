      MODULE Geometry2
      Implicit none

      character*20 Control_info
      character*20 Version
      character*22 Gambit_file
      character*80, save :: Problem_title
      character*80 Program_version
      character*80 Date
      character*40 NodalCoordinates
      character*12 EndOfSection
      character*40 ElementsCells
      character*40 ElementGroup
      character*32, dimension(:), save, allocatable :: GroupName
      character*7, save :: NGroup
      character*11, save :: NElementRead
      character*10, save :: Material
      character*8, save ::  NFlagsRead
      character*40 BoundaryConditions
      character*5  NUMNP,NELEM,NGRPS,NDFCD,NDFVL
      character*6  NBSETS
      character*32, dimension(:), allocatable :: BoundaryName

      integer :: n,ntemp
      integer,  dimension(:,:), allocatable  :: ListOfElementNodesTemp
      integer, save :: NumberOfGroup
      integer, dimension(:), save, allocatable  :: NumberOfGroupFlags
      integer, dimension(:), save, allocatable  :: NGroupFlags
      integer, dimension(:,:), save, allocatable  :: NElementsInGroupTemp
      integer :: NBCValuesPerRecord,NBCCode

      integer, dimension(:,:), allocatable  :: NodeBCTemp
      integer, dimension(:,:), allocatable  :: NElementBCTemp
      integer, dimension(:,:), allocatable  :: NElementBCTypeTemp
      integer, dimension(:,:), allocatable  :: NElementBCFaceTemp

      end MODULE Geometry2