c
C#############################################################################################
      SUBROUTINE SaveGrid
C#############################################################################################
c
      use Geometry1
      use Geometry2
      use Geometry3
      use InteriorBoundaryNodes
      use User0
      use ReadPolymesh
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k,i1,j1
c*********************************************************************************************
c
      FilGrd=trim(name)//'.Grd'
      open (unit=3,file=trim(SolutionDirectory)//'/'//trim(FilGrd),
     *                                             form='unformatted')
      rewind 3
c
C--- Write the neutral file
c
      write(3) NumberOfNodes,NumberOfElements,NumberOfElementGroups,
     *           NumberOfBCSets,NumberOfCoordDir,NumberOfVelComp,
     *           MaximumNumberofElementNodes
c
c---- Nodal Coordinates
c
      write(3) x
      write(3) y
      write(3) z
c
c---- Element/Cell Connectivity
c
      do i=1,NumberOfElements
c
!        write(3) NTypeGeometry(i),NumbOfElementNodes(i),
!     *               NumbOfElementEdges(i),NumbOfElementFaces(i)
        write(3) NTypeGeometry(i),NumbOfElementNodes(i),
     *               NumbOfElementFaces(i)
c
      enddo
c
      do i=1,NumberOfElements
c
        write(3)(ListOfElementNodes(i,j),j=1,NumbOfElementNodes(i))
c
      enddo
c
c---- Element Group Information
c
      do i=1,NumberOfElementGroups
c
        write(3) ElementGroup    
        write(3) NGroup,i,
     *             NElementRead,NumberOfGroupElements(i),
     *             Material,MaterialGroupType(i),
     *             NFlagsRead,NumberOfGroupFlags(i)
c      
        write(3) GroupName(i)      
        write(3)(NGroupFlags(j),j=1,NumberOfGroupFlags(i))
        write(3)(NElementsInGroup(i,j),j=1,
     *             NumberOfGroupElements(i))
c        
      enddo
c
c---- Boundary Conditions
c
      do i=1,NumberOfBCSets
c
        write(3) BoundaryConditions      
        write(3) BoundaryName(i),NBCDataType(i),NBCDataRecords(i),
     *             NBCValuesPerRecord,NBCCode
c
        if(NBCDataType(i).eq.0) then
c
          write(3) (NodeBC(i,j),j=1,NBCDataRecords(i))        
c
        elseif(NBCDataType(i).eq.1) then
c
          do j=1,NBCDataRecords(i)
c
            write(3) NElementBC(i,j),NElementBCType(i,j),
     *                  NElementBCFace(i,j)
c        
          enddo
c
        endif
c        
      enddo
c
c--- Write the grid connectivity
c
      !write(3) NElementEdges
      write(3) NElementFaces
      if(MeshType.eq.'polymesh') write(3) MaximumNumberOfFaceNodes
c 
      do i=1,NumberOfElements
        do j=1,NumbOfElementFaces(i)
c        
        write(3) NumberOfElementFaceNodes(i,j)
c
        enddo
      enddo
!c
!      do i=1,NumberOfElements
!        do j=1,NumbOfElementedges(i)
!        
!        write(3) LocalElementEdgeNodes(i,j,1),
!     *                 LocalElementEdgeNodes(i,j,2)
!        enddo
!      enddo
!c 
      do i=1,NumberOfElements
        do j=1,NumbOfElementFaces(i)
          do k=1,NumberOfElementFaceNodes(i,j)
            write(3) LocalElementFaceNodes(i,j,k)
          enddo
        enddo
      enddo
c
c---- write Global Edge Numbering
c
!      write(3) NIEdges
!c
!      do i=1,NIEdges
!c
!        write(3) EdgeNode1(i),EdgeNode2(i)
!c
!      enddo
!c
!      do i=1,NumberOfElements
!        do j=1,NumbOfElementEdges(i)
!c
!          write(3) ElementEdges(i,j)
!c
!        enddo
!      enddo






c
c---  Develop Global Face Numbering and Neighbors
c
      write(3) NIFaces,NFacesTotal,NBFacesMax,NBFacesTotal
c
      if(MeshType.eq.'neutral') MaximumNumberOfFaceNodes=4
      do i=1,NIFaces
c      
        write(3) (NIFaceNodes(i,j),j=1,MaximumNumberOfFaceNodes), 
     *           NIFaceOwner(i),NIFaceNeighbor(i)
c
      enddo
c
      do i=1,NFacesTotal
        write(3) GlobalFaceNumberOfNodes(i)
      enddo
c      
      do i=1,NumberOfElements
        do j=1,NumbOfElementFaces(i)
c
          write(3) NGlobalEFaces(i,j)
c
        enddo
      enddo
c
c---- Process boundary faces (Assumes BC type 1, i.e. element based)
c
      do i=1,NumberOfBCSets
c
        write(3) NBFaces(i)
c
      enddo
c
      write(3) NBFacesTotal
c
      do i=1,NBFacesTotal
c
        write(3) BCSet(i),BCRecord(i)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
          write(3) (NBFaceNodes(i,j,k),k=1,MaximumNumberOfFaceNodes),
     *    NBFaceOwner(i,j),NBFaceNeighbor(i,j)
c               
        enddo
      enddo
c
      write(3) MaximumNumberofElementFaces
c
      do i=1,NIFaces
c
        write(3) iOwnerNeighbor(i),iNeighborOwner(i) 
c
      enddo
c
      do i=1,NumberOfElements
c
        do j=1,NumbOfElementFaces(i)
c
            write(3) ElementNeighbor(i,j)
c
        enddo
c
      enddo
c
      do i=1,NumberOfElements
c
        write(3) NumberofElementNeighbors(i)
c
      enddo
c
      do i=1,NumberOfNodes     
c      
        write(3) NodeFlag(i)
c
      enddo
c
      write(3) NumberOfInteriorNodes,NumberOfBoundaryNodes
c
      do i=1,NumberOfInteriorNodes     
c      
        write(3) InteriorNode(i)
c
      enddo
c
      do i=1,NumberOfBoundaryNodes     
c      
        write(3) BoundaryNode(i)
c
      enddo
c
c---- Save connectivity of nodes in groups
c
      write(3) (NumberOfNodesInGroup(i),i=1,NumberOfElementGroups)
c
      i1=-1
      do i=1,NumberOfElementGroups
        i1=max(NumberOfNodesInGroup(i),i1)
      enddo
c
      do i=1,i1
c
        write(3) (ListOfNodesInGroup(i,j),j=1,NumberOfElementGroups)
c
      enddo
c
      write(3) (MaxNodesPerGroupElement(i),i=1,NumberOfElementGroups)
      write(3) (MinNodesPerGroupElement(i),i=1,NumberOfElementGroups)
c
      do i=1,NumberOfElementGroups
        do j=1,NumberOfGroupElements(i)
          j1=NElementsInGroup(i,j)
          do k=1,NumbOfElementNodes(j1)
c
             write(3) ListOfElementNodesLocal(j,k,i)
c
          enddo
        enddo
      enddo
c
      if(MeshType.eq.'polymesh') then
c       
        write(3) nPoints
        write(3) points
        write(3) nFaces
        do i=1,nFaces
          write(3) faces(i)%nPoints
        enddo
        do i=1,nFaces
          do j=1,faces(i)%nPoints
            write(3) faces(i)%PointsID(j)
          enddo
        enddo
        do i=1,nFaces
          write(3) faces(i)%owner
        enddo
        do i=1,nFaces
          write(3) faces(i)%neighbour
        enddo
        write(3) nBoundaries
        do i=1,nBoundaries
          write(3) Boundaries(i)%Boundary_name
        enddo
        do i=1,nBoundaries
          write(3) Boundaries(i)%StartFace
        enddo
        do i=1,nBoundaries
          write(3) Boundaries(i)%nFaces
        enddo
        write(3) nCells
        do i=1,nCells
          write(3) Cells(i)%nFaces
        enddo
        do i=1,nCells
          do j=1,Cells(i)%nFaces
            write(3) Cells(i)%FacesID(j)
          enddo
        enddo
        do i=1,nCells
          write(3) Cells(i)%nPoints
        enddo
        do i=1,nCells
          do j=1,Cells(i)%nPoints
            write(3) Cells(i)%PointsID(j)
          enddo
        enddo
c
      endif
c
      close(3)
      return
c
      end