c
C#############################################################################################
      SUBROUTINE CalculateNumberOfNodesInGroups
C#############################################################################################
      use Geometry1, only:NumberOfNodes,NumberOfElementGroups,
     *              NumberOfGroupElements,NElementsInGroup,
     *              NumbOfElementNodes,ListOfElementNodes,
     *              ListOfNodesInGroup,NumberOfNodesInGroup,
     *              MaxNodesPerGroupElement,ListOfElementNodesLocal,
     *              MinNodesPerGroupElement
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer index,i,i1,j,j1,j2,k,k1,k2,k3,k4
      integer ElementMax,NodeMax,NodeMin
      integer, save, dimension(:,:), allocatable :: ListOfNodesInGroupT
c*********************************************************************************************
c
      allocate(ListOfNodesInGroupT(NumberOfNodes,NumberOfElementGroups))
      allocate(NumberOfNodesInGroup(NumberOfElementGroups))
      allocate(MaxNodesPerGroupElement(NumberOfElementGroups))
      allocate(MinNodesPerGroupElement(NumberOfElementGroups))
c
      do i=1,NumberOfElementGroups
c
        index=1
        k4=0
        do j=1,NumberOfGroupElements(i)
          j1=NElementsInGroup(i,j)
          do k=1,NumbOfElementNodes(j1)
            k1=ListOfElementNodes(j1,k)
c
            k3=k4
            index=1
            k2=1
            do while(k2.le.k3.and.index.eq.1)
              if(k1.eq.ListOfNodesInGroupT(k2,i)) then
                index=0
              endif
              k2=k2+1
            enddo
            if(index.eq.1) then
              k4=k4+1
              ListOfNodesInGroupT(k4,i)=k1
            endif
          enddo
        enddo
c
        NumberOfNodesInGroup(i)=k4
c
      enddo
c
      i1=-1
      do i=1,NumberOfElementGroups
        i1=max(NumberOfNodesInGroup(i),i1)
      enddo
c
      allocate(ListOfNodesInGroup(i1,NumberOfElementGroups))
c
      do i=1,i1
        do j=1,NumberOfElementGroups
          ListOfNodesInGroup(i,j)=ListOfNodesInGroupT(i,j)
        enddo
      enddo
c
      do i=1,NumberOfElementGroups
c
        k=-1
        j2=1000000
        do j=1,NumberOfGroupElements(i)
          j1=NElementsInGroup(i,j)
          k=max(k,NumbOfElementNodes(j1))
          j2=min(j2,NumbOfElementNodes(j1))
        enddo  
        MaxNodesPerGroupElement(i)=k
        MinNodesPerGroupElement(i)=j2
      enddo 
c
      deallocate(ListOfNodesInGroupT)
c
      ElementMax=1
      NodeMax=1
      NodeMin=100000000
      do i=1,NumberOfElementGroups
        ElementMax=max(ElementMax,NumberOfGroupElements(i))
        NodeMax=max(NodeMax,MaxNodesPerGroupElement(i))
        NodeMin=min(NodeMin,MinNodesPerGroupElement(i))
      enddo

      allocate(ListOfElementNodesLocal
     *        (ElementMax,NodeMax,NumberOfElementGroups))
c
      do i=1,NumberOfElementGroups
        do j=1,NumberOfGroupElements(i)
          j1=NElementsInGroup(i,j)
          do k=1,NumbOfElementNodes(j1)
            k1=ListOfElementNodes(j1,k)
            do i1=1,NumberOfNodesInGroup(i)
c            
               if(k1.eq.ListOfNodesInGroup(i1,i)) then
                 ListOfElementNodesLocal(j,k,i)=i1
               endif
c
            enddo
          enddo
        enddo
      enddo
c
      return
      end
c
C#############################################################################################
      SUBROUTINE CalculateFoamNumberOfNodesInGroups
C#############################################################################################
      use User0, only: MeshType
      use Geometry1, only:NumberOfNodes,NumberOfElementGroups,
     *              NumberOfGroupElements,NElementsInGroup,
     *              NumbOfElementNodes,ListOfElementNodes,
     *              ListOfNodesInGroup,NumberOfNodesInGroup,
     *              MaxNodesPerGroupElement,ListOfElementNodesLocal,
     *              MinNodesPerGroupElement
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer index,i,i1,j,j1,j2,k,k1,k2,k3,k4
      integer ElementMax,NodeMax,NodeMin
      integer, save, dimension(:,:), allocatable :: ListOfNodesInGroupT
c*********************************************************************************************
c
      allocate(ListOfNodesInGroupT(NumberOfNodes,NumberOfElementGroups))
      allocate(NumberOfNodesInGroup(NumberOfElementGroups))
      allocate(MaxNodesPerGroupElement(NumberOfElementGroups))
      allocate(MinNodesPerGroupElement(NumberOfElementGroups))
c
      do i=1,NumberOfElementGroups
         do j=1,NumberOfNodes
         ListOfNodesInGroupT(j,i)=j
         enddo
c
         NumberOfNodesInGroup(i)=NumberOfNodes
c
      enddo
c
!      do i=1,NumberOfElementGroups
!c
!        index=1
!        k4=0
!        do j=1,NumberOfGroupElements(i)
!          j1=NElementsInGroup(i,j)
!          do k=1,NumbOfElementNodes(j1)
!            k1=ListOfElementNodes(j1,k)
!c
!            k3=k4
!            index=1
!            k2=1
!            do while(k2.le.k3.and.index.eq.1)
!              if(k1.eq.ListOfNodesInGroupT(k2,i)) then
!                index=0
!              endif
!              k2=k2+1
!            enddo
!            if(index.eq.1) then
!              k4=k4+1
!              ListOfNodesInGroupT(k4,i)=k1
!            endif
!          enddo
!        enddo
!c
!        NumberOfNodesInGroup(i)=k4
!c
!      enddo
!c
      i1=-1
      do i=1,NumberOfElementGroups
        i1=max(NumberOfNodesInGroup(i),i1)
      enddo
c
      allocate(ListOfNodesInGroup(i1,NumberOfElementGroups))
c
      do i=1,i1
        do j=1,NumberOfElementGroups
          ListOfNodesInGroup(i,j)=ListOfNodesInGroupT(i,j)
        enddo
      enddo
c
      do i=1,NumberOfElementGroups
c
        k=-1
        j2=100000000
        do j=1,NumberOfGroupElements(i)
          j1=NElementsInGroup(i,j)
          k=max(k,NumbOfElementNodes(j1))
          j2=min(j2,NumbOfElementNodes(j1))
        enddo  
        MaxNodesPerGroupElement(i)=k
        MinNodesPerGroupElement(i)=j2
      enddo 
c
      deallocate(ListOfNodesInGroupT)
c
      ElementMax=1
      NodeMax=1
      NodeMin=100000000
      do i=1,NumberOfElementGroups
        ElementMax=max(ElementMax,NumberOfGroupElements(i))
        NodeMax=max(NodeMax,MaxNodesPerGroupElement(i))
        NodeMin=min(NodeMin,MinNodesPerGroupElement(i))
      enddo

      allocate(ListOfElementNodesLocal
     *        (ElementMax,NodeMax,NumberOfElementGroups))
c
      if(MeshType.eq.'polymesh') then
        ListOfElementNodesLocal=1
        return
      endif
c
      do i=1,NumberOfElementGroups
        do j=1,NumberOfGroupElements(i)
          j1=NElementsInGroup(i,j)
          do k=1,NumbOfElementNodes(j1)
            k1=ListOfElementNodes(j1,k)
            do i1=1,NumberOfNodesInGroup(i)
c            
               if(k1.eq.ListOfNodesInGroup(i1,i)) then
                 ListOfElementNodesLocal(j,k,i)=i1
               endif
c
            enddo
          enddo
        enddo
      enddo
c
      return
      end