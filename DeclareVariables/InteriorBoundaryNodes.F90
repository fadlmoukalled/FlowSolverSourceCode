      MODULE InteriorBoundaryNodes
      Implicit none
      integer, save, dimension(:), allocatable :: InteriorNodeTemp
      integer, save, dimension(:), allocatable :: BoundaryNodeTemp
      integer, save, dimension(:), allocatable :: InteriorNode
      integer, save, dimension(:), allocatable :: BoundaryNode
!      
      integer, save :: NumberOfInteriorNodes,NumberOfBoundaryNodes     
!
      end MODULE InteriorBoundaryNodes