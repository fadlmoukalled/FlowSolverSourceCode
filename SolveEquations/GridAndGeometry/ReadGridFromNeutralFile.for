c
C#############################################################################################
      SUBROUTINE ReadNeutralGrid
C#############################################################################################
      use User0, only: GridScalex,GridScaley,GridScalez
      use Geometry1
      use Geometry2
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k,m
c*********************************************************************************************
c
  1   Format(A20,1x,A20)
  10  Format(5X,A5,5X,A5,5X,A5,4X,A6,5X,A5,5X,A5)
  20  Format(6(1X,I9))
  30  Format(I10,3E20.11)
  40  Format(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))
  50  Format(A7,I10,A11,I10,A11,I10,A8,I10)      
  60  Format(10I8)
  70  Format(10I8)
  80  Format(A32,4I8)
  90  Format(I10)
  100 Format(I10,I5,I5)
c
c********************************************************************************************
c
      read(1,1) Control_info,Version
      read(1,'(A22)') Gambit_file
      read(1,'(A80)') Problem_title
      read(1,'(A80)') Program_version
      read(1,'(A80)') Date
      read(1,10) NUMNP,NELEM,NGRPS,NBSETS,NDFCD,NDFVL
      read(1,20) NumberOfNodes,NumberOfElements,NumberOfElementGroups,
     *           NumberOfBCSets,NumberOfCoordDir,NumberOfVelComp
      read(1,'(A12)') EndOfSection    
c
c---- Nodal Coordinates
c
      read(1,'(A40)') NodalCoordinates
c
      allocate(x(NumberOfNodes))
      allocate(y(NumberOfNodes))
      allocate(z(NumberOfNodes))
c
      do i=1,NumberOfNodes
c
        read(1,30) j,x(j),y(j),z(j)
c
      enddo
c      
      x=x*GridScalex
      y=y*GridScaley
      z=z*GridScalez
c
      read(1,'(A12)') EndOfSection    
c
c---- Element/Cell Connectivity
c
      read(1,'(A40)') ElementsCells    
c
c---- Set the maximum possible number of nodes
c
      n=8
c
c---- Allocate storage based on the maximum possible number of nodes
c
      allocate(NTypeGeometry(NumberOfElements))
      allocate(NumbOfElementNodes(NumberOfElements))
      allocate(ListOfElementNodesTemp(NumberOfElements,n))
c      allocate(NumbOfElementEdges(NumberOfElements))
      allocate(NumbOfElementFaces(NumberOfElements))
c
      do i=1,NumberOfElements
c
        read(1,40) j,NTypeGeometry(j),NumbOfElementNodes(j),
     *       (ListOfElementNodesTemp(j,k),k=1,NumbOfElementNodes(j))
c
c---- Number of element edges and faces: Brick
c
        if(NTypeGeometry(j).eq.4) then
c
c          NumbOfElementedges(j)=12
          NumbOfElementFaces(j)=6
c
c---- Number of element edges and faces: Wedge (Prism)
c
        elseif(NTypeGeometry(j).eq.5) then
c
c          NumbOfElementedges(j)=9
          NumbOfElementFaces(j)=5
c
c---- Number of element edges and faces: Tetrahedron
c
        elseif(NTypeGeometry(j).eq.6) then
c
c          NumbOfElementedges(j)=6
          NumbOfElementFaces(j)=4
c
c---- Number of element edges and faces: Pyramid
c
        elseif(NTypeGeometry(j).eq.7) then
c
c          NumbOfElementedges(j)=8
          NumbOfElementFaces(j)=5
c
        endif
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
        do j=1,NumbOfElementNodes(i)
c
          ListOfElementNodes(i,j)=ListOfElementNodesTemp(i,j)
c
        enddo
      enddo
c
      deallocate(ListOfElementNodesTemp)
c
      read(1,'(A12)') EndOfSection    
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
        read(1,'(A40)') ElementGroup    
        read(1,50) NGroup,NumberOfGroup,
     *             NElementRead,NumberOfGroupElements(NumberOfGroup),
     *             Material,MaterialGroupType(NumberOfGroup),
     *             NFlagsRead,NumberOfGroupFlags(NumberOfGroup)
c      
        read(1,'(A32)') GroupName(NumberOfGroup)      
        read(1,60)(NGroupFlags(j),j=1,NumberOfGroupFlags(NumberOfGroup))
        read(1,70)(NElementsInGroupTemp(NumberOfGroup,j),j=1,
     *             NumberOfGroupElements(NumberOfGroup))
c        
        read(1,'(A12)') EndOfSection    
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
        read(1,'(A40)') BoundaryConditions      
        read(1,80) BoundaryName(i),NBCDataType(i),NBCDataRecords(i),
     *             NBCValuesPerRecord,NBCCode
c
        if(NBCDataType(i).eq.0) then
c
          read(1,90) (NodeBCTemp(i,j),j=1,NBCDataRecords(i))        
c
        elseif(NBCDataType(i).eq.1) then
c
          do j=1,NBCDataRecords(i)
c
            read(1,100) NElementBCTemp(i,j),NElementBCTypeTemp(i,j),
     *                  NElementBCFaceTemp(i,j)
c        
          enddo
c
        endif
c
        read(1,'(A12)') EndOfSection    
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
      close(1)
c
      return
c
      end