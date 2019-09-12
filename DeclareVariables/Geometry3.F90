      MODULE Geometry3
      Implicit none
!      integer, save, dimension(:,:,:), allocatable :: LocalElementEdgeNodes
      integer, save, dimension(:,:,:), allocatable :: LocalElementFaceNodes
!      integer, save, dimension(:,:,:), allocatable :: LocalElementFaceEdges
!      integer, save, dimension(:,:,:), allocatable :: ElementFaceEdges
!      integer, save :: NElementEdges     
      integer, save :: NElementFaces      

!      integer, save, dimension(:), allocatable :: EdgeNode1
!      integer, save, dimension(:), allocatable :: EdgeNode2
!      integer, save, dimension(:,:), allocatable :: ElementEdges
!      integer, save :: NIEdges
      integer, save :: MaximumNumberOfFaceNodes

      integer, save, dimension(:,:), allocatable :: NumberOfElementFaceNodes
      integer, save, dimension(:), allocatable :: GlobalFaceNumberOfNodesT
      integer, save, dimension(:), allocatable :: GlobalFaceNumberOfNodes

      
      integer, save :: NIFaces,NFacesTotal
      integer, save, dimension(:,:),allocatable :: NIFaceNodes
      integer, save, dimension(:),allocatable :: NIFaceOwner
      integer, save, dimension(:),allocatable :: NIFaceNeighbor
      integer, save :: NBFacesMax,NBFacesTotal
      integer, save, dimension(:), Allocatable :: NBFaces 
      integer, save, dimension(:,:,:),allocatable :: NBFaceNodes
      integer, save, dimension(:,:),allocatable :: NBFaceOwner
      integer, save, dimension(:,:),allocatable :: NBFaceNeighbor
      integer, save, dimension(:,:),allocatable :: NGlobalEFaces

      integer, save, dimension(:),allocatable :: iOwnerNeighbor
      integer, save, dimension(:),allocatable :: iNeighborOwner
      integer, save, dimension(:,:),allocatable :: ElementNeighbor
      integer, save, dimension(:),allocatable :: NumberofElementNeighbors
      integer, save, dimension(:),allocatable :: BCSet
      integer, save, dimension(:),allocatable :: BCRecord

      end MODULE Geometry3