C
C#############################################################################################
      SUBROUTINE GlobalConnectivity
C#############################################################################################
c
      use Geometry1
      use Geometry3
      use InteriorBoundaryNodes
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer, dimension(:,:), allocatable :: NFaceNodesT
      integer, dimension(:), allocatable :: NFaceOwnerT
      integer, dimension(:), allocatable :: NFaceNeighborT
c
      integer :: i,j,k,m,ntemp,index1,index2,i1
      integer :: iElement,iNeighbor,iEdge,indexEdge
      integer :: k1,k2,k3
c
c---      NTypeGeometry=4  Brick
c---      NTypeGeometry=5  Wedge (Prism)
c---      NTypeGeometry=6  Tetrahedron
c---      NTypeGeometry=7  Pyramid
c
c*********************************************************************************************
c
c      NElementEdges=-1
c      
c      do i=1,NumberOfElements
c
c        NElementEdges=max(NElementEdges,NumbOfElementEdges(i))
c
c      enddo
c                      
      NElementFaces=-1
c      
      do i=1,NumberOfElements
c
        NElementFaces=max(NElementFaces,NumbOfElementFaces(i))
c
      enddo
c                      
c      allocate(LocalElementEdgeNodes(NumberOfElements,NElementEdges,2))
      allocate(LocalElementFaceNodes(NumberOfElements,NElementFaces,4))
      allocate(NumberOfElementFaceNodes(NumberOfElements,NElementFaces))
c      allocate(LocalElementFaceEdges(NumberOfElements,NElementFaces,4))
c      allocate(ElementFaceEdges(NumberOfElements,NElementFaces,4))
c 
      LocalElementFaceNodes=-1
c
      do i=1,NumberOfElements
c
        if(NTypeGeometry(i).eq.4) then
!c
!          LocalElementEdgeNodes(i,1,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,1,2)=ListOfElementNodes(i,5)  
!          LocalElementEdgeNodes(i,2,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,2,2)=ListOfElementNodes(i,2)  
!          LocalElementEdgeNodes(i,3,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,3,2)=ListOfElementNodes(i,6)  
!          LocalElementEdgeNodes(i,4,1)=ListOfElementNodes(i,5) 
!          LocalElementEdgeNodes(i,4,2)=ListOfElementNodes(i,6)  
!          LocalElementEdgeNodes(i,5,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,5,2)=ListOfElementNodes(i,4)  
!          LocalElementEdgeNodes(i,6,1)=ListOfElementNodes(i,4) 
!          LocalElementEdgeNodes(i,6,2)=ListOfElementNodes(i,8)  
!          LocalElementEdgeNodes(i,7,1)=ListOfElementNodes(i,6) 
!          LocalElementEdgeNodes(i,7,2)=ListOfElementNodes(i,8)  
!          LocalElementEdgeNodes(i,8,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,8,2)=ListOfElementNodes(i,4)  
!          LocalElementEdgeNodes(i,9,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,9,2)=ListOfElementNodes(i,7)  
!          LocalElementEdgeNodes(i,10,1)=ListOfElementNodes(i,7) 
!          LocalElementEdgeNodes(i,10,2)=ListOfElementNodes(i,8)  
!          LocalElementEdgeNodes(i,11,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,11,2)=ListOfElementNodes(i,3)  
!          LocalElementEdgeNodes(i,12,1)=ListOfElementNodes(i,5) 
!          LocalElementEdgeNodes(i,12,2)=ListOfElementNodes(i,7)  
!c
          NumberOfElementFaceNodes(i,1)=4
          NumberOfElementFaceNodes(i,2)=4
          NumberOfElementFaceNodes(i,3)=4
          NumberOfElementFaceNodes(i,4)=4
          NumberOfElementFaceNodes(i,5)=4
          NumberOfElementFaceNodes(i,6)=4
c
          LocalElementFaceNodes(i,1,1)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,1,2)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,1,3)=ListOfElementNodes(i,6) 
          LocalElementFaceNodes(i,1,4)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,2,1)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,2,2)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,2,3)=ListOfElementNodes(i,8) 
          LocalElementFaceNodes(i,2,4)=ListOfElementNodes(i,6) 
          LocalElementFaceNodes(i,3,1)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,3,2)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,3,3)=ListOfElementNodes(i,7) 
          LocalElementFaceNodes(i,3,4)=ListOfElementNodes(i,8) 
          LocalElementFaceNodes(i,4,1)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,4,2)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,4,3)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,4,4)=ListOfElementNodes(i,7) 
          LocalElementFaceNodes(i,5,1)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,5,2)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,5,3)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,5,4)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,6,1)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,6,2)=ListOfElementNodes(i,6) 
          LocalElementFaceNodes(i,6,3)=ListOfElementNodes(i,8) 
          LocalElementFaceNodes(i,6,4)=ListOfElementNodes(i,7) 
c
        elseif(NTypeGeometry(i).eq.5) then
!c
!          LocalElementEdgeNodes(i,1,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,1,2)=ListOfElementNodes(i,2)  
!          LocalElementEdgeNodes(i,2,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,2,2)=ListOfElementNodes(i,3)  
!          LocalElementEdgeNodes(i,3,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,3,2)=ListOfElementNodes(i,1)  
!          LocalElementEdgeNodes(i,4,1)=ListOfElementNodes(i,4) 
!          LocalElementEdgeNodes(i,4,2)=ListOfElementNodes(i,5)  
!          LocalElementEdgeNodes(i,5,1)=ListOfElementNodes(i,5) 
!          LocalElementEdgeNodes(i,5,2)=ListOfElementNodes(i,6)  
!          LocalElementEdgeNodes(i,6,1)=ListOfElementNodes(i,6) 
!          LocalElementEdgeNodes(i,6,2)=ListOfElementNodes(i,4)  
!          LocalElementEdgeNodes(i,7,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,7,2)=ListOfElementNodes(i,4)  
!          LocalElementEdgeNodes(i,8,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,8,2)=ListOfElementNodes(i,5)  
!          LocalElementEdgeNodes(i,9,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,9,2)=ListOfElementNodes(i,6)  
!c
          NumberOfElementFaceNodes(i,1)=4
          NumberOfElementFaceNodes(i,2)=4
          NumberOfElementFaceNodes(i,3)=4
          NumberOfElementFaceNodes(i,4)=3
          NumberOfElementFaceNodes(i,5)=3
          LocalElementFaceNodes(i,1,1)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,1,2)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,1,3)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,1,4)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,2,1)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,2,2)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,2,3)=ListOfElementNodes(i,6) 
          LocalElementFaceNodes(i,2,4)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,3,1)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,3,2)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,3,3)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,3,4)=ListOfElementNodes(i,6) 
          LocalElementFaceNodes(i,4,1)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,4,2)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,4,3)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,5,1)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,5,2)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,5,3)=ListOfElementNodes(i,6) 
c
        elseif(NTypeGeometry(i).eq.6) then
!c
!          LocalElementEdgeNodes(i,1,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,1,2)=ListOfElementNodes(i,2)  
!          LocalElementEdgeNodes(i,2,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,2,2)=ListOfElementNodes(i,3)  
!          LocalElementEdgeNodes(i,3,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,3,2)=ListOfElementNodes(i,1)  
!          LocalElementEdgeNodes(i,4,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,4,2)=ListOfElementNodes(i,4)  
!          LocalElementEdgeNodes(i,5,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,5,2)=ListOfElementNodes(i,4)  
!          LocalElementEdgeNodes(i,6,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,6,2)=ListOfElementNodes(i,4)  
!c
          NumberOfElementFaceNodes(i,1)=3
          NumberOfElementFaceNodes(i,2)=3
          NumberOfElementFaceNodes(i,3)=3
          NumberOfElementFaceNodes(i,4)=3
          LocalElementFaceNodes(i,1,1)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,1,2)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,1,3)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,2,1)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,2,2)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,2,3)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,3,1)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,3,2)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,3,3)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,4,1)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,4,2)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,4,3)=ListOfElementNodes(i,4) 
c
        elseif(NTypeGeometry(i).eq.7) then
!c
!          LocalElementEdgeNodes(i,1,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,1,2)=ListOfElementNodes(i,2)  
!          LocalElementEdgeNodes(i,2,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,2,2)=ListOfElementNodes(i,4)  
!          LocalElementEdgeNodes(i,3,1)=ListOfElementNodes(i,4) 
!          LocalElementEdgeNodes(i,3,2)=ListOfElementNodes(i,3)  
!          LocalElementEdgeNodes(i,4,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,4,2)=ListOfElementNodes(i,1)  
!          LocalElementEdgeNodes(i,5,1)=ListOfElementNodes(i,1) 
!          LocalElementEdgeNodes(i,5,2)=ListOfElementNodes(i,5)  
!          LocalElementEdgeNodes(i,6,1)=ListOfElementNodes(i,2) 
!          LocalElementEdgeNodes(i,6,2)=ListOfElementNodes(i,5)  
!          LocalElementEdgeNodes(i,7,1)=ListOfElementNodes(i,4) 
!          LocalElementEdgeNodes(i,7,2)=ListOfElementNodes(i,5)  
!          LocalElementEdgeNodes(i,8,1)=ListOfElementNodes(i,3) 
!          LocalElementEdgeNodes(i,8,2)=ListOfElementNodes(i,5)  
!c
          NumberOfElementFaceNodes(i,1)=4
          NumberOfElementFaceNodes(i,2)=3
          NumberOfElementFaceNodes(i,3)=3
          NumberOfElementFaceNodes(i,4)=3
          NumberOfElementFaceNodes(i,5)=3
          LocalElementFaceNodes(i,1,1)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,1,2)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,1,3)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,1,4)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,2,1)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,2,2)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,2,3)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,3,1)=ListOfElementNodes(i,2) 
          LocalElementFaceNodes(i,3,2)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,3,3)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,4,1)=ListOfElementNodes(i,4) 
          LocalElementFaceNodes(i,4,2)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,4,3)=ListOfElementNodes(i,5) 
          LocalElementFaceNodes(i,5,1)=ListOfElementNodes(i,3) 
          LocalElementFaceNodes(i,5,2)=ListOfElementNodes(i,1) 
          LocalElementFaceNodes(i,5,3)=ListOfElementNodes(i,5) 
c
        endif
      enddo
c
c---- Develop Global Edge Numbering
c
      !ntemp=NElementEdges*NumberOfElements
      !allocate(EdgeNode1(ntemp))
      !allocate(EdgeNode2(ntemp))
      !allocate(ElementEdges(NumberOfElements,NElementEdges))
c      
!      iEdge=0
!c
!      do j=1,NumbOfElementEdges(1)
!c
!        iEdge=iEdge+1      
!        EdgeNode1(iEdge)=LocalElementEdgeNodes(1,j,1)
!        EdgeNode2(iEdge)=LocalElementEdgeNodes(1,j,2)
!c
!        ElementEdges(1,j)=iEdge
!c
!      enddo
!c
!      do i=2,NumberOfElements
!c
!        if(mod(i,1000).eq.0) 
!     *       print*,'Global edge connectivity for element ',i,
!     *                   ' out of ', NumberOfElements,' elements'
!c
!        do j=1,NumbOfElementEdges(i)
!c
!          k=1
!          indexEdge=1
!          do while(k.le.iEdge.and.indexEdge.eq.1)
!c
!            if(EdgeNode1(k).eq.LocalElementEdgeNodes(i,j,1).or.
!     *             EdgeNode1(k).eq.LocalElementEdgeNodes(i,j,2)) then
!              if(EdgeNode2(k).eq.LocalElementEdgeNodes(i,j,1).or.
!     *             EdgeNode2(k).eq.LocalElementEdgeNodes(i,j,2)) then
!c
!                ElementEdges(i,j)=k
!                indexEdge=0
!c
!              endif
!c
!            endif
!c
!            k=k+1
!c
!          enddo 
!c
!          if(indexEdge.eq.1) then
!c
!            iEdge=iEdge+1      
!            EdgeNode1(iEdge)=LocalElementEdgeNodes(i,j,1)
!            EdgeNode2(iEdge)=LocalElementEdgeNodes(i,j,2)
!c
!            ElementEdges(i,j)=iEdge
!c
!          endif        
!c
!        enddo
!      enddo
!c
!      NIEdges=iEdge          








c
c---  Develop Global Face Numbering and Neighbors
c
      ntemp=6*NumberOfElements
      allocate(NFaceNodesT(ntemp,4))
      allocate(NFaceOwnerT(ntemp))
      allocate(NFaceNeighborT(ntemp))
      allocate(GlobalFaceNumberOfNodesT(ntemp))
      allocate(NGlobalEFaces(NumberOfElements,NElementFaces))
c
      NIFaces=0
c
      do i=1,NumberOfElements
c
        if(mod(i,1000).eq.0) 
     *  print*,'Global face connectivity for element ',i,
     *                   ' out of ', NumberOfElements,' elements'
c
        do j=1,NumbOfElementFaces(i)
c
          index1=1
          index2=1
          NIFaces=NIFaces+1
c
          k1=NumberOfElementFaceNodes(i,j)
          do k=1,k1
c
            NFaceNodesT(NIFaces,k)=LocalElementFaceNodes(i,j,k)
c
          enddo
c
          NFaceOwnerT(NIFaces)=i
          GlobalFaceNumberOfNodesT(NIFaces)=k1
c
c--- Check if this is a boundary face. If it is do not accept
c
          do k=1,NumberOfBCSets
            do m=1,NBCDataRecords(k)
c
              if(NElementBC(k,m).eq.i) then
                if(NElementBCFace(k,m).eq.j) then
c
                  NIFaces=NIFaces-1
                  index1=0
                  exit
c
                endif
              endif
c               
            enddo
          enddo
c
c--- Process internal boundary faces
c
          if(index1.eq.1) then
c
c--- Exclude faces that have been visited
c
            do k=1,NIFaces-1
c
              if(k1.eq.GlobalFaceNumberOfNodesT(k)) then
c
                if(k1.eq.3) then
c
                  if(NFaceNodesT(NIFaces,1).eq.NFaceNodesT(k,1).or.
     *                NFaceNodesT(NIFaces,1).eq.NFaceNodesT(k,2).or.
     *                NFaceNodesT(NIFaces,1).eq.NFaceNodesT(k,3)) then
c
                    if(NFaceNodesT(NIFaces,2).eq.NFaceNodesT(k,1).or.
     *                NFaceNodesT(NIFaces,2).eq.NFaceNodesT(k,2).or.
     *                NFaceNodesT(NIFaces,2).eq.NFaceNodesT(k,3)) then
c
                      if(NFaceNodesT(NIFaces,3).eq.NFaceNodesT(k,1).or.
     *                  NFaceNodesT(NIFaces,3).eq.NFaceNodesT(k,2).or.
     *                  NFaceNodesT(NIFaces,3).eq.NFaceNodesT(k,3)) then
c
                        NIFaces=NIFaces-1
                        index2=0
                        exit
c
                      endif
c
                    endif
c
                  endif
c
                elseif(k1.eq.4) then
c
                  if(NFaceNodesT(NIFaces,1).eq.NFaceNodesT(k,1).or.
     *                NFaceNodesT(NIFaces,1).eq.NFaceNodesT(k,2).or.
     *                NFaceNodesT(NIFaces,1).eq.NFaceNodesT(k,3).or.
     *                NFaceNodesT(NIFaces,1).eq.NFaceNodesT(k,4)) then
c
                    if(NFaceNodesT(NIFaces,2).eq.NFaceNodesT(k,1).or.
     *                NFaceNodesT(NIFaces,2).eq.NFaceNodesT(k,2).or.
     *                NFaceNodesT(NIFaces,2).eq.NFaceNodesT(k,3).or.
     *                NFaceNodesT(NIFaces,2).eq.NFaceNodesT(k,4)) then
c
                      if(NFaceNodesT(NIFaces,3).eq.NFaceNodesT(k,1).or.
     *                  NFaceNodesT(NIFaces,3).eq.NFaceNodesT(k,2).or.
     *                  NFaceNodesT(NIFaces,3).eq.NFaceNodesT(k,3).or.
     *                  NFaceNodesT(NIFaces,3).eq.NFaceNodesT(k,4)) then
c
                        if(NFaceNodesT(NIFaces,4).eq.NFaceNodesT(k,1)
     *                      .or.NFaceNodesT(NIFaces,4).eq.
     *                      NFaceNodesT(k,2).or.NFaceNodesT(NIFaces,4)
     *                    .eq.NFaceNodesT(k,3).or.NFaceNodesT(NIFaces,4)
     *                                        .eq.NFaceNodesT(k,4)) then
c
                          NIFaces=NIFaces-1
                          index2=0
                          exit
c
                        endif
c
                      endif
c
                    endif
c
                  endif
c
                endif
c
              endif
c
            enddo
c
c--- Process faces that have not been visited
c
            if(index2.eq.1) then
c
              do k=1,NumberOfElements
                do m=1,NumbOfElementFaces(k)
c
                  if(k.eq.i.and.m.eq.j) then
                  else
                    if(k1.eq.NumberOfElementFaceNodes(k,m)) then
                      if(k1.eq.3) then
c
                        if(NFaceNodesT(NIFaces,1).eq.
     *                         LocalElementFaceNodes(k,m,1).or.
     *                         NFaceNodesT(NIFaces,1).eq.
     *                         LocalElementFaceNodes(k,m,2).or.
     *                         NFaceNodesT(NIFaces,1).eq.
     *                         LocalElementFaceNodes(k,m,3)) then
                          if(NFaceNodesT(NIFaces,2).eq.
     *                         LocalElementFaceNodes(k,m,1).or.
     *                         NFaceNodesT(NIFaces,2).eq.
     *                         LocalElementFaceNodes(k,m,2).or.
     *                         NFaceNodesT(NIFaces,2).eq.
     *                         LocalElementFaceNodes(k,m,3)) then
                            if(NFaceNodesT(NIFaces,3).eq.
     *                         LocalElementFaceNodes(k,m,1).or.
     *                         NFaceNodesT(NIFaces,3).eq.
     *                         LocalElementFaceNodes(k,m,2).or.
     *                         NFaceNodesT(NIFaces,3).eq.
     *                         LocalElementFaceNodes(k,m,3)) then
c
                              NFaceNeighborT(NIFaces)=k
                              NGlobalEFaces(k,m)=NIFaces    
                              exit
c
                            endif 
                          endif 
                        endif
                      elseif(k1.eq.4) then
c
                        if(NFaceNodesT(NIFaces,1).eq.
     *                         LocalElementFaceNodes(k,m,1).or.
     *                         NFaceNodesT(NIFaces,1).eq.
     *                         LocalElementFaceNodes(k,m,2).or.
     *                         NFaceNodesT(NIFaces,1).eq.
     *                         LocalElementFaceNodes(k,m,3).or.
     *                         NFaceNodesT(NIFaces,1).eq.
     *                         LocalElementFaceNodes(k,m,4)) then
                          if(NFaceNodesT(NIFaces,2).eq.
     *                         LocalElementFaceNodes(k,m,1).or.
     *                         NFaceNodesT(NIFaces,2).eq.
     *                         LocalElementFaceNodes(k,m,2).or.
     *                         NFaceNodesT(NIFaces,2).eq.
     *                         LocalElementFaceNodes(k,m,3).or.
     *                         NFaceNodesT(NIFaces,2).eq.
     *                         LocalElementFaceNodes(k,m,4)) then
                            if(NFaceNodesT(NIFaces,3).eq.
     *                         LocalElementFaceNodes(k,m,1).or.
     *                         NFaceNodesT(NIFaces,3).eq.
     *                         LocalElementFaceNodes(k,m,2).or.
     *                         NFaceNodesT(NIFaces,3).eq.
     *                         LocalElementFaceNodes(k,m,3).or.
     *                         NFaceNodesT(NIFaces,3).eq.
     *                         LocalElementFaceNodes(k,m,4)) then
                              if(NFaceNodesT(NIFaces,4).eq.
     *                           LocalElementFaceNodes(k,m,1).or.
     *                           NFaceNodesT(NIFaces,4).eq.
     *                           LocalElementFaceNodes(k,m,2).or.
     *                           NFaceNodesT(NIFaces,4).eq.
     *                           LocalElementFaceNodes(k,m,3).or.
     *                           NFaceNodesT(NIFaces,4).eq.
     *                           LocalElementFaceNodes(k,m,4)) then
c
                                NFaceNeighborT(NIFaces)=k
                                NGlobalEFaces(k,m)=NIFaces    
                                exit
c
                              endif 
                            endif 
                          endif 
                        endif 
                      endif 
                    endif 
                  endif 
c
                enddo
              enddo
c
            endif
c
          endif
c
          if(index1.eq.1.and.index2.eq.1) NGlobalEFaces(i,j)=NIFaces
        enddo
              enddo
c
c--- Define correct size of faces     
c
      NFacesTotal=NIFaces
      do i=1,NumberOfBCSets
        NFacesTotal=NFacesTotal+NBCDataRecords(i)
      enddo
      allocate(NIFaceNodes(NIFaces,4))
      allocate(NIFaceOwner(NIFaces))
      allocate(NIFaceNeighbor(NIFaces))
      allocate(GlobalFaceNumberOfNodes(NFacesTotal))
c
      GlobalFaceNumberOfNodes(1:NIFaces)=
     *            GlobalFaceNumberOfNodesT(1:NIFaces)
c
      deallocate(GlobalFaceNumberOfNodesT)
c
      NIFaceNodes=-1
      do i=1,NIFaces
c
        do j=1,GlobalFaceNumberOfNodes(i)
c
          NIFaceNodes(i,j)=NFaceNodesT(i,j) 
c
        enddo
c
        if(NFaceOwnerT(i).lt.NFaceNeighborT(i)) then
c
          NIFaceOwner(i)=NFaceOwnerT(i) 
          NIFaceNeighbor(i)=NFaceNeighborT(i)
c
        else
c
          NIFaceOwner(i)=NFaceNeighborT(i)
          NIFaceNeighbor(i)=NFaceOwnerT(i) 
c
        endif
c
      enddo
c
      deallocate(NFaceNodesT)
      deallocate(NFaceOwnerT)
      deallocate(NFaceNeighborT)
c
c---- Process boundary faces (Assumes BC type 1, i.e. element based)
c
      index1=-1
      i1=0
      do i=1,NumberOfBCSets
        index1=max(index1,NBCDataRecords(i))
        i1=i1+NBCDataRecords(i)
      enddo
      NBFacesMax=index1
c
      allocate(NBFaces(NumberOfBCSets))
      allocate(NBFaceNodes(NumberOfBCSets,index1,4))
      allocate(NBFaceOwner(NumberOfBCSets,index1))
      allocate(NBFaceNeighbor(NumberOfBCSets,index1))
      allocate(BCSet(i1))
      allocate(BCRecord(i1))
c
      NBFaces=0
      NBFacesTotal=0
      NFacesTotal=NIFaces
      NBFaceNodes=-1
c
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
          NBFacesTotal=NBFacesTotal+1
          NBFaces(i)=NBFaces(i)+1
          k=NElementBC(i,j)
          m=NElementBCFace(i,j)
          NBFaceNodes(i,j,1)=LocalElementFaceNodes(k,m,1)
          NBFaceNodes(i,j,2)=LocalElementFaceNodes(k,m,2)
          NBFaceNodes(i,j,3)=LocalElementFaceNodes(k,m,3)
          NBFaceNodes(i,j,4)=LocalElementFaceNodes(k,m,4)
          NBFaceOwner(i,j)=NElementBC(i,j)
          NBFaceNeighbor(i,j)=-1
c
          NFacesTotal=NFacesTotal+1
          GlobalFaceNumberOfNodes(NFacesTotal)=0
          do k1=1,4
            if(NBFaceNodes(i,j,k1).ne.-1) then
                 GlobalFaceNumberOfNodes(NFacesTotal)=
     *                      GlobalFaceNumberOfNodes(NFacesTotal)+1
            endif
          enddo
          NGlobalEFaces(k,m)=NFacesTotal
c
          BCSet(NBFacesTotal)=i
          BCRecord(NBFacesTotal)=j
c
        enddo
      enddo
c
c--- Locate the element neighbors, OwnerNeighbor, and NeighborOwner connectivity
c
      MaximumNumberofElementFaces=-1
      do i=1,NumberOfElements     
c
        MaximumNumberofElementFaces=
     *     max(MaximumNumberofElementFaces,NumbOfElementfaces(i))
c
        enddo
c
      allocate(iOwnerNeighbor(NIFaces))
      allocate(iNeighborOwner(NIFaces))
      allocate(ElementNeighbor(NumberOfElements,
     *                   MaximumNumberofElementFaces))
c
      iOwnerNeighbor=0
      iNeighborOwner=0
      ElementNeighbor=0
      
c
      do i=1,NIFaces
c
        iElement=NIFaceOwner(i)
        iNeighbor=NIFaceNeighbor(i)
c
        do j=1,NumbOfElementFaces(iElement)
c
          k=NGlobalEFaces(iElement,j)
c
          if(k.eq.i) then
c
            iOwnerNeighbor(i)=j
            ElementNeighbor(iElement,j)=iNeighbor
            exit
c
          endif
c
        enddo
c
        do j=1,NumbOfElementFaces(iNeighbor)
c
          k=NGlobalEFaces(iNeighbor,j)
c
          if(k.eq.i) then
c
            iNeighborOwner(i)=j 
            ElementNeighbor(iNeighbor,j)=iElement
            exit
c
          endif
c
        enddo
c
      enddo
c
      allocate(NumberofElementNeighbors(NumberOfElements))
c
      do i=1,NumberOfElements
c
        NumberofElementNeighbors(i)=NumbOfElementFaces(i)
c
      enddo
c
c--- Create flags for nodes (interior node ---> flag=1,boundary node ---> flag=0) 
c
      allocate(NodeFlag(NumberOfNodes))
      NodeFlag=1
c      
      k1=NIFaces
      do i=1,NumberOfBCSets     
        do j=1,NBCDataRecords(i)
c          k1=NElementBC(i,j)
c          k2=NElementBCFace(i,j)
          k1=k1+1
          do k=1,GlobalFaceNumberOfNodes(k1)
c          do k=1,NumberOfElementFaceNodes(k1,k2)
c 
            k2=NBFaceNodes(i,j,k)
            if(k2.gt.0) NodeFlag(k2)=0
c
          enddo
        enddo
      enddo
c
C********************************************************************************************
c
      allocate(InteriorNodeTemp(NumberOfNodes))
      allocate(BoundaryNodeTemp(NumberOfNodes))
      j=0
      k=0   
      do i=1,NumberOfNodes
c
        if(NodeFlag(i).eq.1) then
c
          j=j+1
          InteriorNodeTemp(j)=i     
c
        elseif(NodeFlag(i).eq.0) then
c
          k=k+1
          BoundaryNodeTemp(k)=i
c
        endif
c
      enddo
c
      NumberOfInteriorNodes=j
      NumberOfBoundaryNodes=k
c
      allocate(InteriorNode(NumberOfInteriorNodes))
      allocate(BoundaryNode(NumberOfBoundaryNodes))
c
      InteriorNode(1:NumberOfInteriorNodes)=
     *                     InteriorNodeTemp(1:NumberOfInteriorNodes)
      BoundaryNode(1:NumberOfBoundaryNodes)=
     *                     BoundaryNodeTemp(1:NumberOfBoundaryNodes)
c
      deallocate(InteriorNodeTemp)
      deallocate(BoundaryNodeTemp)
c
!      do i=1,NumberOfNodes
!        write(80,*) i,x(i),y(i),z(i)
!      enddo
!
!      do i=1,NumberOfElements
!        do j=1,NumbOfElementNodes(i)
!c
!          write(80,*) i,NumbOfElementNodes(i),ListOfElementNodes(i,j)
!        enddo
!      enddo
!      
!      write(80,*) 'owner neighbor'
!      do i=1,NIFaces
!c
!        write(80,*) i,NIFaceOwner(i),NIFaceNeighbor(i)
!c
!      enddo
!c
!      do i=1,NumberOfBCSets
!        do j=1,NBCDataRecords(i)
!          write(80,*) i,j,NBFaceOwner(i,j),NbFaceNeighbor(i,j) 
!        enddo
!      enddo   
!        
!        
!        
!      write(80,*) 'ElementNeighbor'
!      do i=1,NumberOfElements
!        do j=1,NumbOfElementfaces(i)
!         write(80,*) i,ElementNeighbor(i,j)
!        enddo
!      enddo
!
!      pause
!      stop
!
!
!
      return
c
      end