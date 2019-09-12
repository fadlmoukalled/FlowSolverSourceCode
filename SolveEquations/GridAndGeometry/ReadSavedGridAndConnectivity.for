c
C#############################################################################################
      SUBROUTINE ReadSavedGrid
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
      integer i,j,k,i1,j1,ElementMax,NodeMax
c*********************************************************************************************
c
      FilGrd=trim(name)//'.Grd'
      open (unit=3,file=trim(SolutionDirectory)//'/'//trim(FilGrd),
     *                                             form='unformatted')
      rewind 3
c
C--- Write the neutral file
c
      read(3) NumberOfNodes,NumberOfElements,NumberOfElementGroups,
     *           NumberOfBCSets,NumberOfCoordDir,NumberOfVelComp,
     *            MaximumNumberofElementNodes
c
c---- Nodal Coordinates
c
c
      allocate(x(NumberOfNodes))
      allocate(y(NumberOfNodes))
      allocate(z(NumberOfNodes))
c
      read(3) x
      read(3) y
      read(3) z
      x=x*GridScalex
      y=y*GridScaley
      z=z*GridScalez
c
c---- Element/Cell Connectivity
c
      allocate(NTypeGeometry(NumberOfElements))
      allocate(NumbOfElementNodes(NumberOfElements))
c      allocate(NumbOfElementEdges(NumberOfElements))
      allocate(NumbOfElementFaces(NumberOfElements))
c
      do i=1,NumberOfElements
c
c        read(3) NTypeGeometry(i),NumbOfElementNodes(i),
c     *                NumbOfElementEdges(i),NumbOfElementFaces(i)
        read(3) NTypeGeometry(i),NumbOfElementNodes(i),
     *                NumbOfElementFaces(i)
c
      enddo
c
      ntemp=-1
      do i=1,NumberOfElements
c
        ntemp=max(ntemp,NumbOfElementNodes(i))
c
      enddo
c
      allocate(ListOfElementNodes(NumberOfElements,ntemp))
c
      do i=1,NumberOfElements
c
        read(3)(ListOfElementNodes(i,j),j=1,NumbOfElementNodes(i))
c
      enddo
c
c---- Element Group Information
c
      allocate(NumberOfGroupElements(NumberOfElementGroups))
      allocate(MaterialGroupType(NumberOfElementGroups))
      allocate(NumberOfGroupFlags(NumberOfElementGroups))
      allocate(GroupName(NumberOfElementGroups))
      allocate(NElementsInGroupTemp
     *          (NumberOfElementGroups,NumberOfElements))
      allocate(NGroupFlags(1))
c
      do i=1,NumberOfElementGroups
c
        read(3) ElementGroup    
        read(3) NGroup,NumberOfGroup,
     *             NElementRead,NumberOfGroupElements(NumberOfGroup),
     *             Material,MaterialGroupType(NumberOfGroup),
     *             NFlagsRead,NumberOfGroupFlags(NumberOfGroup)
c      
        read(3) GroupName(NumberOfGroup)      
        read(3)(NGroupFlags(j),j=1,NumberOfGroupFlags(NumberOfGroup))
        read(3)(NElementsInGroupTemp(NumberOfGroup,j),j=1,
     *             NumberOfGroupElements(NumberOfGroup))
c        
      enddo
c
      ntemp=-1
      do i=1,NumberOfElementGroups
c
        ntemp=max(ntemp,NumberOfGroupElements(i))
c
      enddo
c
      allocate(NElementsInGroup(NumberOfElementGroups,ntemp))
c
      do i=1,NumberOfElementGroups
        do j=1,NumberOfGroupElements(i)
c
          NElementsInGroup(i,j)=NElementsInGroupTemp(i,j)
c
        enddo
      enddo
c
      deallocate(NElementsInGroupTemp)
c
c---- Boundary Conditions
c
      allocate(BoundaryName(NumberOfBCSets))
      allocate(NBCDataType(NumberOfBCSets))
      allocate(NBCDataRecords(NumberOfBCSets))
c
      ntemp=NumberOfElements*NumberOfBCSets
      allocate(NodeBCTemp(NumberOfBCSets,ntemp))
      allocate(NElementBCTemp(NumberOfBCSets,ntemp))
      allocate(NElementBCTypeTemp(NumberOfBCSets,ntemp))
      allocate(NElementBCFaceTemp(NumberOfBCSets,ntemp))
c
      do i=1,NumberOfBCSets
c
        read(3) BoundaryConditions      
        read(3) BoundaryName(i),NBCDataType(i),NBCDataRecords(i),
     *             NBCValuesPerRecord,NBCCode
c
        if(NBCDataType(i).eq.0) then
c
          read(3) (NodeBCTemp(i,j),j=1,NBCDataRecords(i))        
c
        elseif(NBCDataType(i).eq.1) then
c
          do j=1,NBCDataRecords(i)
c
            read(3) NElementBCTemp(i,j),NElementBCTypeTemp(i,j),
     *                  NElementBCFaceTemp(i,j)
c        
          enddo
c
        endif
c        
      enddo
c
      ntemp=-1
      do i=1,NumberOfBCSets
c
        ntemp=max(ntemp,NBCDataRecords(i))
c
      enddo
c
      allocate(NodeBC(NumberOfBCSets,ntemp))
      allocate(NElementBC(NumberOfBCSets,ntemp))
      allocate(NElementBCType(NumberOfBCSets,ntemp))
      allocate(NElementBCFace(NumberOfBCSets,ntemp))
c
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
          NodeBC(i,j)=NodeBCTemp(i,j)
          NElementBC(i,j)=NElementBCTemp(i,j)
          NElementBCType(i,j)=NElementBCTypeTemp(i,j)
          NElementBCFace(i,j)=NElementBCFaceTemp(i,j)
c
        enddo
      enddo
c
      deallocate(NodeBCTemp)
      deallocate(NElementBCTemp)
      deallocate(NElementBCTypeTemp)
      deallocate(NElementBCFaceTemp)
c
c--- read the grid connectivity
!c
!      read(3) NElementEdges
      read(3) NElementFaces
!c
      ntemp=4
      if(MeshType.eq.'polymesh') read(3) ntemp
      MaximumNumberOfFaceNodes=ntemp
          
      allocate(NumberOfElementFaceNodes(NumberOfElements,NElementFaces))
      allocate(LocalElementFaceNodes
     *                        (NumberOfElements,NElementFaces,ntemp))
!      allocate(LocalElementEdgeNodes(NumberOfElements,NElementEdges,2))
c 
      do i=1,NumberOfElements
        do j=1,NumbOfElementFaces(i)
c        
        read(3) NumberOfElementFaceNodes(i,j)
c
        enddo
      enddo
!c
!      do i=1,NumberOfElements
!        do j=1,NumbOfElementedges(i)
!        
!        read(3) LocalElementEdgeNodes(i,j,1),
!     *                 LocalElementEdgeNodes(i,j,2)
!        enddo
!      enddo
!c 
      do i=1,NumberOfElements
        do j=1,NumbOfElementFaces(i)
          do k=1,NumberOfElementFaceNodes(i,j)
            read(3) LocalElementFaceNodes(i,j,k)
          enddo
        enddo
      enddo
c
c---- read Global Edge Numbering
c
!      read(3) NIEdges
!      allocate(EdgeNode1(NIEdges))
!      allocate(EdgeNode2(NIEdges))
!      allocate(ElementEdges(NumberOfElements,NElementEdges))
!c
!      do i=1,NIEdges
!c
!        read(3) EdgeNode1(i),EdgeNode2(i)
!c
!      enddo
!c
!      do i=1,NumberOfElements
!        do j=1,NumbOfElementEdges(i)
!c
!          read(3) ElementEdges(i,j)
!c
!        enddo
!      enddo




c
c---  Develop Global Face Numbering and Neighbors
c
      read(3) NIFaces,NFacesTotal,NBFacesMax,NBFacesTotal
c
      allocate(NIFaceNodes(NIFaces,MaximumNumberOfFaceNodes))
      allocate(NIFaceOwner(NIFaces))
      allocate(NIFaceNeighbor(NIFaces))
      allocate(GlobalFaceNumberOfNodes(NFacesTotal))
c
      do i=1,NIFaces
c 
        read(3)  (NIFaceNodes(i,j),j=1,MaximumNumberOfFaceNodes),
     *           NIFaceOwner(i),NIFaceNeighbor(i)
c
      enddo
c
      do i=1,NFacesTotal
        read(3) GlobalFaceNumberOfNodes(i)
      enddo
c      
      allocate(NGlobalEFaces(NumberOfElements,NElementFaces))
c
      do i=1,NumberOfElements
        do j=1,NumbOfElementFaces(i)
c
          read(3) NGlobalEFaces(i,j)
c
        enddo
      enddo
c
c---- Process boundary faces (Assumes BC type 1, i.e. element based)
c
      ntemp=-1
      do i=1,NumberOfBCSets
        ntemp=max(ntemp,NBCDataRecords(i))
      enddo
c
      allocate(NBFaces(NumberOfBCSets))
      allocate(NBFaceNodes
     *              (NumberOfBCSets,ntemp,MaximumNumberOfFaceNodes))
      allocate(NBFaceOwner(NumberOfBCSets,ntemp))
      allocate(NBFaceNeighbor(NumberOfBCSets,ntemp))
c
      do i=1,NumberOfBCSets
c
        read(3) NBFaces(i)
c
      enddo
c
      read(3) NBFacesTotal
c
      allocate(BCSet(NBFacesTotal))
      allocate(BCRecord(NBFacesTotal))
c
      do i=1,NBFacesTotal
c
        read(3) BCSet(i),BCRecord(i)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
          read(3) (NBFaceNodes(i,j,k),k=1,MaximumNumberOfFaceNodes),
     *               NBFaceOwner(i,j),NBFaceNeighbor(i,j)
c               
        enddo
      enddo
c
      read(3) MaximumNumberofElementFaces
c
      allocate(iOwnerNeighbor(NIFaces))
      allocate(iNeighborOwner(NIFaces))
      allocate(ElementNeighbor(NumberOfElements,NElementFaces))
c
      do i=1,NIFaces
c
        read(3) iOwnerNeighbor(i),iNeighborOwner(i) 
c
      enddo
c
      do i=1,NumberOfElements
c
        do j=1,NumbOfElementFaces(i)
c
            read(3) ElementNeighbor(i,j)
c
        enddo
c
      enddo
c
      allocate(NumberofElementNeighbors(NumberOfElements))
c
      do i=1,NumberOfElements
c
        read(3) NumberofElementNeighbors(i)
c
      enddo
c
      allocate(NodeFlag(NumberOfNodes))
c      
      do i=1,NumberOfNodes     
c      
        read(3) NodeFlag(i)
c
      enddo
c
      read(3) NumberOfInteriorNodes,NumberOfBoundaryNodes
c
      allocate(InteriorNode(NumberOfInteriorNodes))
      allocate(BoundaryNode(NumberOfBoundaryNodes))
c
      do i=1,NumberOfInteriorNodes     
c      
        read(3) InteriorNode(i)
c
      enddo
c
      do i=1,NumberOfBoundaryNodes     
c      
        read(3) BoundaryNode(i)
c
      enddo
c
c---- Read connectivity of nodes in groups
c
      allocate(NumberOfNodesInGroup(NumberOfElementGroups))
      read(3) (NumberOfNodesInGroup(i),i=1,NumberOfElementGroups)
c
      i1=-1
      do i=1,NumberOfElementGroups
        i1=max(NumberOfNodesInGroup(i),i1)
      enddo
c
      allocate(ListOfNodesInGroup(i1,NumberOfElementGroups))
      do i=1,i1
c
        read(3) (ListOfNodesInGroup(i,j),j=1,NumberOfElementGroups)
c
      enddo
c
      allocate(MaxNodesPerGroupElement(NumberOfElementGroups))
      allocate(MinNodesPerGroupElement(NumberOfElementGroups))
c
      read(3) (MaxNodesPerGroupElement(i),i=1,NumberOfElementGroups)
      read(3) (MinNodesPerGroupElement(i),i=1,NumberOfElementGroups)
c
      ElementMax=1
      NodeMax=1
      do i=1,NumberOfElementGroups
        ElementMax=max(ElementMax,NumberOfGroupElements(i))
        NodeMax=max(NodeMax,MaxNodesPerGroupElement(i))
      enddo

      allocate(ListOfElementNodesLocal
     *        (ElementMax,NodeMax,NumberOfElementGroups))
c
      do i=1,NumberOfElementGroups
        do j=1,NumberOfGroupElements(i)
          j1=NElementsInGroup(i,j)
          do k=1,NumbOfElementNodes(j1)
c
             read(3) ListOfElementNodesLocal(j,k,i)
c
          enddo
        enddo
      enddo
c
      if(MeshType.eq.'polymesh') then
c       
        read(3) nPoints
c
        allocate(points(nPoints,3))
c
        read(3) points
        read(3) nFaces
c        
        allocate(Faces(nFaces))
c        
        do i=1,nFaces
          read(3) faces(i)%nPoints
        enddo
        do i=1,nFaces
          allocate(faces(i)%PointsID(faces(i)%nPoints))
          do j=1,faces(i)%nPoints
            read(3) faces(i)%PointsID(j)
          enddo
        enddo
        do i=1,nFaces
          read(3) faces(i)%owner
        enddo
        do i=1,nFaces
          read(3) faces(i)%neighbour
        enddo
        read(3) nBoundaries
c        
        allocate(Boundaries(nBoundaries))
c        
        do i=1,nBoundaries
          read(3) Boundaries(i)%Boundary_name
        enddo
        do i=1,nBoundaries
          read(3) Boundaries(i)%StartFace
        enddo
        do i=1,nBoundaries
          read(3) Boundaries(i)%nFaces
        enddo
        read(3) nCells
        allocate(cells(nCells))
        do i=1,nCells
          read(3) Cells(i)%nFaces
        enddo
        do i=1,nCells
          allocate(Cells(i)%FacesID(Cells(i)%nFaces))
          do j=1,Cells(i)%nFaces
            read(3) Cells(i)%FacesID(j)
          enddo
        enddo
        do i=1,nCells
          read(3) Cells(i)%nPoints
        enddo
        do i=1,nCells
          allocate(Cells(i)%PointsID(Cells(i)%nPoints))
          do j=1,Cells(i)%nPoints
            read(3) Cells(i)%PointsID(j)
          enddo
        enddo
c
      endif
c    
      close(3)
      return
c
      end
