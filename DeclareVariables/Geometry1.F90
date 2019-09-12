      MODULE Geometry1
      Implicit none
      integer, save :: NumberOfNodes
      integer, save :: NumberOfElements
      integer, save :: NumberOfElementGroups
      integer, save :: NumberOfBCSets
      integer, save :: NumberOfCoordDir
      integer, save :: NumberOfVelComp
      integer, save :: MaximumNumberofElementNodes
      integer, save :: MaximumNumberofElementFaces

      integer, save, dimension(:), allocatable :: NTypeGeometry
      integer, save, dimension(:), allocatable  :: NumbOfElementNodes
      integer, save, dimension(:), allocatable  :: NumbOfElementFaces
!      integer, save, dimension(:), allocatable  :: NumbOfElementEdges
      integer, save, dimension(:,:), allocatable  :: ListOfElementNodes
      integer, save, dimension(:), allocatable  :: NumberOfGroupElements
      integer, save, dimension(:), allocatable  :: MaterialGroupType
      integer, save, dimension(:,:), allocatable  :: NElementsInGroup
      integer, save :: NBCSets
      integer, save, dimension(:), allocatable  :: NBCDataType
      integer, save, dimension(:), allocatable  :: NBCDataRecords
      integer, save, dimension(:,:), allocatable  :: NodeBC
      integer, save, dimension(:,:), allocatable  :: NElementBC
      integer, save, dimension(:,:), allocatable  :: NElementBCType
      integer, save, dimension(:,:), allocatable  :: NElementBCFace
      integer, save, dimension(:), allocatable  :: NodeFlag

      integer, save, dimension(:,:), allocatable :: ListOfNodesInGroup
      integer, save, dimension(:), allocatable :: NumberOfNodesInGroup
      integer, save, dimension(:), allocatable :: MaxNodesPerGroupElement
      integer, save, dimension(:), allocatable :: MinNodesPerGroupElement
      integer, save, dimension(:,:,:), allocatable :: ListOfElementNodesLocal
      double precision, save, dimension(:), allocatable :: x,y,z    

      end MODULE Geometry1