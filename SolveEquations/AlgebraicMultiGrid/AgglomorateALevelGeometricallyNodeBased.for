c
C#############################################################################################
c
      SUBROUTINE AgglomorateLevelGeometricallyNode(iLevel)
c
C#############################################################################################
c
      use User0
      use Constants1, only:tiny,big
      use Geometry1, only:NumberOfElements
      use MultiGrid1
      use reallocateData
      use Variables2, only:anb
      use Geometry3, only:NElementFaces,ElementNeighbor,
     *                    NumberofElementNeighbors
      use GeometricMG1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      logical DoIt,agglomorated
c
      integer iLevel,iLevel1,iElement,itemp,index
      integer iSeed,iSeed1,iSeed2
      integer iParent,iParentT,iChildT,iChildNBT
      integer NumberOfFineElements,theNumberOfParents,NumberOfParents
      integer NumberOfChildren0,iNBchildren
      integer iNBlocal,iNBlocal1,iNB,iNBchildren2,iChildNBlocal
      integer iChildlocal,iChild,iChildNB,iOrphan
      integer maxNeighbors,minNeighbors,NChildren,iNBuse
      integer maxChildren,minChildren
      integer OldNumberOfParents,DroppedParents
      integer minChildrenNew,minChildrenOld
      integer iChildren,INB1,INBParent
      integer k,k1,k2,k3,k4,k5,k6,i,j,i1,i2,iMax
      integer minimumNumberOfChildren,maximumNumberOfChildren
      integer maxNeighborsOld,iOfiNB,iOfiOrphan
      integer ij,ij1,ij2,ij3,ij4,ij5,ij6,ij7,jn,ij2n,ij2n1
      integer indexn
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9
      integer ijn1,ijn2,ijn3,ijn4,ijn5,ijn6,ijn7
c
      double precision ratio,strength,ratiomin
c--------------------------------------------------------------------------------------------
      integer, dimension(:), allocatable :: TransitionArray
      double precision, dimension(:), allocatable :: TransitionArrayR
c
      data ratiomin/0.5/
      
      integer NumberOfFineNodes,m,nmax,mmax,n,index1
      integer Localmmax,indexnode,iternode
      double precision Volmin
c*********************************************************************************************
c
      iLevel1=iLevel-1
c
      NumberOfFineElements=NumberOfElementsMG(iLevel1)
      NumberOfFineNodes=NumberOfNodesMG(iLevel1)
c
      do i=1,NumberOfFineElements
c
        ij=ijbeginE(iLevel1)+i
        Parents(ij)=0
        NumberOfChildren(ij)=0
        ElementParent(ij)=0
c
        do j=1,NBChildrenMaxMG(iLevel1)
c
          ij1=ijbeginN(iLevel1)+i+(j-1)*NumberOfFineElements
          children(ij1)=0
c
        enddo
c
      enddo
c
c---- Find number of elements connected to a node
c
      if(size(NumberOfElementsConnectedToNodeMG).eq.1) then
c
        deallocate(NumberOfElementsConnectedToNodeMG)
        allocate(NumberOfElementsConnectedToNodeMG(NumberOfFineNodes))
        NumberOfElementsConnectedToNodeMG=0
c
      else
c
        allocate(TransitionArray
     *           (size(NumberOfElementsConnectedToNodeMG)))
        TransitionArray=0
        TransitionArray=NumberOfElementsConnectedToNodeMG  
        deallocate(NumberOfElementsConnectedToNodeMG)
        ij1=size(TransitionArray)
        ij2=ij1+NumberOfFineNodes
        allocate(NumberOfElementsConnectedToNodeMG(ij2))
        NumberOfElementsConnectedToNodeMG=0
        NumberOfElementsConnectedToNodeMG(1:ij1)=TransitionArray(1:ij1)
        deallocate(TransitionArray)
c
      endif
c
      do i=1,NumberOfFineNodes
c
        ij=ijbeginNode(iLevel1)+i
c
        if(FlagNodesMG(ij).eq.iLevel1) then
c
          m=0
c
          do j=1,NumberOfFineElements
c
            ij1=ijbeginE(iLevel1)+j
c
            do k=1,NumbOfElementNodesMG(ij1)
c
              ij2=ijbeginENode(iLevel1)+j+(k-1)*NumberOfFineElements
c
              if(ListOfElementNodesMG(ij2).eq.NodesMG(ij)) then
c
                m=m+1
c
              endif
c
            enddo
c
          enddo
c
          NumberOfElementsConnectedToNodeMG(ij)=m
c
        endif
c
      enddo
c
c---- Find the maximum number of elements connected to a node
c
      mmax=-1
      do i=1,NumberOfFineNodes
c
        ij=ijbeginNode(iLevel1)+i
        mmax=max(mmax,NumberOfElementsConnectedToNodeMG(ij))
c
      enddo
c
      if(size(ListOfElementsConnectedToNodeMG).eq.1) then
c
        deallocate(ListOfElementsConnectedToNodeMG)
        allocate(ListOfElementsConnectedToNodeMG
     *                           (NumberOfFineNodes*mmax))
c
        allocate(TransitionArray(iLevel))
        TransitionArray(1:iLevel1)=ijbeginNElement(1:iLevel1)
        deallocate(ijbeginNElement)
        allocate(ijbeginNElement(iLevel))
        ijbeginNElement(1:iLevel1)=TransitionArray(1:iLevel1)
        ijbeginNElement(iLevel)=NumberOfFineNodes*mmax
        deallocate(TransitionArray)
c
      else
c
        ij1=size(ListOfElementsConnectedToNodeMG)
        allocate(TransitionArray(ij1))
        TransitionArray=ListOfElementsConnectedToNodeMG
        deallocate(ListOfElementsConnectedToNodeMG)
        allocate(ListOfElementsConnectedToNodeMG
     *                          (ij1+NumberOfFineNodes*mmax))
        ListOfElementsConnectedToNodeMG(1:ij1)=TransitionArray(1:ij1)
        deallocate(TransitionArray)
c
        allocate(TransitionArray(iLevel))
        TransitionArray(1:iLevel1)=ijbeginNElement(1:iLevel1)
        deallocate(ijbeginNElement)
        allocate(ijbeginNElement(iLevel))
        ijbeginNElement(1:iLevel1)=TransitionArray(1:iLevel1)
        ijbeginNElement(iLevel)=
     *       ijbeginNElement(iLevel1)+NumberOfFineNodes*mmax
        deallocate(TransitionArray)
c
      endif
c
c---- Find elements connected to a node
c
      do i=1,NumberOfFineNodes
c
        ij=ijbeginNode(iLevel1)+i
c
        if(FlagNodesMG(ij).eq.iLevel1) then
c
          m=0
c
          do j=1,NumberOfFineElements
c
            ij1=ijbeginE(iLevel1)+j
c
            do k=1,NumbOfElementNodesMG(ij1)
c
              ij2=ijbeginENode(iLevel1)+j+(k-1)*NumberOfFineElements
c
              if(ListOfElementNodesMG(ij2).eq.NodesMG(ij)) then
c
                m=m+1
                ij3=ijbeginNElement(iLevel1)+i+(m-1)*NumberOfFineNodes
                ListOfElementsConnectedToNodeMG(ij3)=ij1
c
              endif
c
            enddo
c
          enddo
c
        endif
c
      enddo
c       
c---- Agglomorate all elements attached to a node
c
      allocate(VolumeMGTemp(NumberOfFineElements))
c
      Localmmax=mmax
c
      iParent=0
      NumberOfChildren0=0
c
      do while(Localmmax.ge.1)
c
        do i=1,NumberOfFineNodes
c
          ij=ijbeginNode(iLevel1)+i
c
          if(NumberOfElementsConnectedToNodeMG(ij).eq.Localmmax) then
c
            if(FlagNodesMG(ij).eq.iLevel1) then
c
              iParent=iParent+1
              ij1=ijbeginE(iLevel1)+iParent
              NumberOfChildren(ij1)=
     *                  NumberOfElementsConnectedToNodeMG(ij)
c
              VolumeMGTemp(iParent)=0.
c
              do j=1,NumberOfElementsConnectedToNodeMG(ij)
c
                ij2=ijbeginNElement(iLevel1)+i+(j-1)*NumberOfFineNodes
                k=ListOfElementsConnectedToNodeMG(ij2)
                ichildren=ijbeginN(iLevel1)+iParent+
     *                                (j-1)*NumberOfFineElements
                children(ichildren)=k-ijbeginE(iLevel1)
                Parents(k)=iParent
                VolumeMGTemp(iParent)=VolumeMGTemp(iParent)+VolumeMG(k)
c
                do m=1,NumbOfElementNodesMG(k)
c
                  ij3=ijbeginENode(iLevel1)+k+
     *                   (m-1)*NumberOfFineElements-ijbeginE(iLevel1)
                  ij4=ListOfElementNodesMG(ij3)
                  k2=0
                  k3=0
c
                  do while(k2.le.NumberOfFineNodes.and.k3.eq.0)
c
                    k2=k2+1
c
                    if(NodesMG(ijbeginNode(iLevel1)+k2).eq.ij4) then
c
                      ij4=k2+ijbeginNode(iLevel1)
                      k3=1
c
                    endif
c
                  enddo
c
                  if(FlagNodesMG(ij4).eq.iLevel1) 
     *                                  FlagNodesMG(ij4)=iLevel
c
                enddo
c
              enddo
c
              FlagNodesMG(ij)=-iLevel1
              NumberOfChildren0=NumberOfChildren0+NumberOfChildren(ij1)
c
            endif

          endif
c
        enddo
c
        Localmmax=Localmmax-1
c
      enddo
c
      NumberOfParents = iParent
c
c---- Check for orphan elements
c
      indexnode=1
      iternode=0
c
      do while(indexnode.eq.1)
c      
        indexnode=0
        iternode=iternode+1
c
        do i=1,NumberOfFineElements
c              
          ij=ijbeginE(iLevel1)+i 
c
          if(Parents(ij).eq.0) then
c
            Volmin=big
            k=0
c
            do j=1,NumberofElementNeighborsMG(ij)
c
              ij1=ijbeginN(iLevel1)+i+(j-1)*NumberOfFineElements
              iNB=ElementNeighborMG(ij1)
c
              if(iNB.ne.0) then
c
                ij2=ijbeginE(iLevel1)+iNB
c
                if(Parents(ij2).ne.0) then
c
                  if(VolumeMGTemp(Parents(ij2)).lt.Volmin) then
c
                    Volmin=VolumeMGTemp(Parents(ij2))
                    k=Parents(ij2)
c
                  endif
c
                endif
c
              endif
c
            enddo
c               
            if(k.ne.0) then
c
              VolumeMGTemp(k)=VolumeMGTemp(k)+VolumeMG(ij)
              Parents(ij)=k
              ij3=ijbeginE(iLevel1)+k
              itemp=NumberOfChildren(ij3)
              ij4=ijbeginN(iLevel1)+k+
     *                        itemp*NumberOfFineElements
              children(ij4)=i
              NumberOfChildren(ij3)= itemp+1
              NumberOfChildren0=NumberOfChildren0+1
c
            else
c
              indexnode=1
c
            endif
c
          endif    
c
        enddo                     
c
      enddo                     
c
c---- Find coarse element nodes
c
c          
c---- Allocate storage 
c 
      ij1=NumberOfParents*AllocateMGGeometricStorage*iLevel*mmax
      ij2=NumberOfParents
      allocate(ListOfElementNodesMGTemp(ij1))
      allocate(NumbOfElementNodesMGTemp(ij2))
c
      nmax=-1
c
      do i=1,NumberOfParents
c
        ij1=ijbeginE(iLevel1)+i
c
        n=0
        do j=1,NumberOfChildren(ij1)
c
          ij3=ijbeginN(iLevel1)+i+(j-1)*NumberOfFineElements
          ij2=children(ij3)
c
          if(ij2.ne.0) then
c          
            do k1=1,NumbOfElementNodesMG(ij2+ijbeginE(iLevel1))
c
              ij4=ijbeginENode(iLevel1)+ij2+
     *                        (k1-1)*NumberOfFineElements 
              k= ListOfElementNodesMG(ij4)
              k2=0
              k3=0
c
              do while(k2.le.NumberOfFineNodes.and.k3.eq.0)
c
                k2=k2+1
c
                if(NodesMG(ijbeginNode(iLevel1)+k2).eq.k) then
c
                  k4=ijbeginNode(iLevel1)+k2
                  k3=1
c
                endif
c
              enddo
c                     
c              if(FlagNodesMG(k4).eq.iLevel) then
              if(FlagNodesMG(k4).ne.-iLevel1) then
c
                index1=0
c
                do i1=1,n          
c
                  ij5=i+(i1-1)*NumberOfParents
c
                  if(k.eq.ListOfElementNodesMGTemp(ij5)) then
c
                    index1=1    
c
                  endif         
c                  
                enddo
c          
                if(index1.eq.0) then
c          
                  n=n+1
                  ij6=i+(n-1)*NumberOfParents
c
                  ijn1=NumberOfParents*
     *                  AllocateMGGeometricStorage*iLevel*mmax
c
                  if(ij6.gt.ijn1) then
c 
                    allocate(TransitionArray(ijn1))
                    TransitionArray=ListOfElementNodesMGTemp
                    deallocate(ListOfElementNodesMGTemp)
                    AllocateMGGeometricStorage=
     *                          AllocateMGGeometricStorage+50
                    ijn2=NumberOfParents*
     *                  AllocateMGGeometricStorage*iLevel*mmax
                     allocate(ListOfElementNodesMGTemp(ijn2))
                     ListOfElementNodesMGTemp(1:ijn1)=
     *                             TransitionArray(1:ijn1)
                     deallocate(TransitionArray)
c                      
                  endif                  
c
                  ListOfElementNodesMGTemp(ij6)=k
c          
                endif
c
              endif
c
            enddo
c
          endif
c
        enddo
c
        NumbOfElementNodesMGTemp(i)=n
        nmax=max(nmax,n)
c
      enddo
c
c---- Calculate number of interior nodes on coarse grid
c
c      write(86,*) Parents
      ij1=size(NodesMG)
      allocate(NodesMGTemp(NumberOfFineNodes))
      allocate(FlagNodesMGTemp(NumberOfFineNodes))
c
      j=0
c
      do i=1,NumberOfFineNodes
c 
        ij=ijbeginNode(iLevel1)+i
c
        if(FlagNodesMG(ij).eq.iLevel.or.FlagNodesMG(ij).eq.0) then
c
          j=j+1
          NodesMGTemp(j)=NodesMG(ij)
          FlagNodesMGTemp(j)=FlagNodesMG(ij)
c
        endif
c        
      enddo
c
      allocate(TransitionArray(ij1))
      TransitionArray=NodesMG
      deallocate(NodesMG)
      allocate(NodesMG(ij1+j))
      NodesMG(1:ij1)=TransitionArray(1:ij1)
      NodesMG(ij1+1:ij1+j)=NodesMGTemp(1:j)
      deallocate(NodesMGTemp)
      deallocate(TransitionArray)
c
      allocate(TransitionArray(iLevel1))
      TransitionArray=NumberOfNodesMG
      deallocate(NumberOfNodesMG)
      allocate(NumberOfNodesMG(iLevel))
      NumberOfNodesMG(1:iLevel1)=TransitionArray(1:iLevel1)
      NumberOfNodesMG(iLevel)=j
      deallocate(TransitionArray)
c
      allocate(TransitionArray(ij1))
      TransitionArray=FlagNodesMG
      deallocate(FlagNodesMG)
      allocate(FlagNodesMG(ij1+j))
      FlagNodesMG(1:ij1)=TransitionArray(1:ij1)
      FlagNodesMG(ij1+1:ij1+j)=FlagNodesMGTemp(1:j)
      deallocate(FlagNodesMGTemp)
      deallocate(TransitionArray)
c
      allocate(TransitionArray(iLevel1))
      TransitionArray=ijbeginNode
      deallocate(ijbeginNode)
      allocate(ijbeginNode(iLevel))
      ijbeginNode(1:iLevel1)=TransitionArray(1:iLevel1)
      ijbeginNode(iLevel)=ijbeginNode(iLevel1)+NumberOfNodesMG(iLevel1)
      deallocate(TransitionArray)
c
      allocate(TransitionArray(iLevel1))
      TransitionArray=ijbeginENode
      deallocate(ijbeginENode)
      allocate(ijbeginENode(iLevel))
      ijbeginENode(1:iLevel1)=TransitionArray(1:iLevel1)
      ijbeginENode(iLevel)=size(ListOfElementNodesMG)
      deallocate(TransitionArray)
c
      ij1=size(ListOfElementNodesMG)
      allocate(TransitionArray(ij1))
      TransitionArray=ListOfElementNodesMG
      deallocate(ListOfElementNodesMG)
      allocate(ListOfElementNodesMG(ij1+NumberOfParents*nmax))
      ListOfElementNodesMG(1:ij1)=TransitionArray(1:ij1)
      ListOfElementNodesMG(ij1+1:ij1+NumberOfParents*nmax)=
     *       ListOfElementNodesMGTemp(1:NumberOfParents*nmax)
      deallocate(ListOfElementNodesMGTemp)
      deallocate(TransitionArray)
c
      ij1=size(NumbOfElementNodesMG)
      allocate(TransitionArray(ij1))
      TransitionArray=NumbOfElementNodesMG
      deallocate(NumbOfElementNodesMG)
      allocate(NumbOfElementNodesMG(ij1+NumberOfParents))
      NumbOfElementNodesMG(1:ij1)=TransitionArray(1:ij1)
      NumbOfElementNodesMG(ij1+1:ij1+NumberOfParents)=
     *       NumbOfElementNodesMGTemp(1:NumberOfParents)
      deallocate(NumbOfElementNodesMGTemp)
      deallocate(TransitionArray)
c
      ij1=size(VolumeMG)
      allocate(TransitionArray(ij1))
      TransitionArray=VolumeMG
      deallocate(VolumeMG)
      allocate(VolumeMG(ij1+NumberOfParents))
      VolumeMG(1:ij1)=TransitionArray(1:ij1)
      VolumeMG(ij1+1:ij1+NumberOfParents)=
     *          VolumeMGTemp(1:NumberOfParents)
      deallocate(VolumeMGTemp)
      deallocate(TransitionArray)
c
c--- Calculate the current minimum and maximum number of children
c
      maxChildren=-1
      minChildren=100
c
      do i=1,NumberOfParents
c
        ij=ijbeginE(iLevel1)+i
        j=NumberOfChildren(ij)
        maxChildren=max(maxChildren,j)
        minChildren=min(minChildren,j)
c
      enddo
c
c--- Calculate the exact maximum number of children and save it
c
      maxChildren=nmax
      NBChildrenMaxMG(iLevel1)=maxChildren
c
c---------------------------------------------------------------------------
c
c
c--- Allocate for connectivity of the coarse grid
c
      allocate(TransitionArray(iLevel1))
c
      do i=1,iLevel1 
c
        TransitionArray(i)=NumberOfElementsMG(i)
c
      enddo
c
      deallocate(NumberOfElementsMG)
c
      allocate(NumberOfElementsMG(iLevel))
c
      do i=1,iLevel1 
c
        NumberOfElementsMG(i)=TransitionArray(i)
c
      enddo
c
      deallocate(TransitionArray)      
c
c-----------------------------------------------------------------------------
      NumberOfElementsMG(iLevel)=NumberOfParents
c      NElementFacesMG(iLevel1)=max(nmax,mmax)
      call SetLevelIndices(iLevel)
c-----------------------------------------------------------------------------
c
      ij=0
c
      do i=1,iLevel1
c
        ij=ij+NumberOfElementsMG(i)*NElementFacesMG(i)
c
      enddo
c
      allocate(TransitionArray(ij))
      TransitionArray=ElementNeighborMG   
      deallocate(ElementNeighborMG)
      ij1=ij+NumberOfElementsMG(iLevel)*3*NElementFacesMG(iLevel1)
      allocate(ElementNeighborMG(ij1))
      ElementNeighborMG=0
c      
      do i=1,iLevel1    
        do j=1,NumberOfElementsMG(i)
c
          ij=ijbeginE(i)+j          
          ij2=ijbeginN(i)+j          
c
          do k=1,NumberofElementNeighborsMG(ij)
c
            ij1=ij2+(k-1)*NumberOfElementsMG(i)
            ElementNeighborMG(ij1)=TransitionArray(ij1)     
c
          enddo
c
        enddo
      enddo
c
      deallocate(TransitionArray)
c
c-----------------------------------------------------------------------------
c
      ij=0
      do i=1,iLevel1    
c
        ij=ij+NumberOfElementsMG(i)           
c
      enddo
c
      allocate(TransitionArray(ij))      
c
      TransitionArray=NumberofElementNeighborsMG
      deallocate(NumberofElementNeighborsMG)
      ij1=ij+NumberOfElementsMG(iLevel)
      allocate(NumberofElementNeighborsMG(ij1))
      NumberofElementNeighborsMG=0
c
      do i=1,iLevel1
        do j=1,NumberOfElementsMG(i)
c
          ij=ijbeginE(i)+j
          NumberofElementNeighborsMG(ij)=TransitionArray(ij)
c
        enddo
      enddo
c
      deallocate(TransitionArray)      
c
c--- Setup connectivity of the coarse grid
c
      do iElement=1,NumberOfFineElements
c
        ij=ijbeginE(iLevel-1)+iElement
c
        do iNBlocal=1,NumberofElementNeighborsMG(ij)
c
          ij1=ijbeginN(iLevel-1)+iElement+
     *                   (iNBlocal-1)*NumberOfFineElements
          iNB = ElementNeighborMG(ij1)
c
          if(iNB.ne.0) then
c
            iNB1=ijbeginE(iLevel-1)+iNB
c
            if (Parents(ij).ne.Parents(iNB1)) then
c
              iParent=Parents(ij)
              ij2=ijbeginE(iLevel)+iParent
              iMax=NumberofElementNeighborsMG(ij2)
              DoIt=.true.
              i2=Parents(iNB1)
c
              do i=1,iMax
c
                ij3=ijbeginN(iLevel)+iParent+
     *                 (i-1)*NumberOfElementsMG(iLevel)
                i1=ElementNeighborMG(ij3)
c
                if(i1.eq.i2) then
c
                  DoIt=.false.
c
                endif
c
              enddo
c
              if(DoIt) then
c
                NumberofElementNeighborsMG(ij2)= 
     *                    NumberofElementNeighborsMG(ij2)+1
                j=NumberofElementNeighborsMG(ij2)
                ij3=ijbeginN(iLevel)+iParent+
     *                    (j-1)*NumberOfElementsMG(iLevel)
                ElementNeighborMG(ij3)=Parents(iNB1)
c
              endif
c
            endif
c
          endif
c
        enddo
c
      enddo
c
c--- Correct memory allocation
c
      maxNeighbors=-1
      minNeighbors=100
c
      do i=1,NumberOfParents
c
        ij=ijbeginE(iLevel)+i
        j=NumberofElementNeighborsMG(ij)
        maxNeighbors=max(maxNeighbors,j)
        minNeighbors=min(minNeighbors,j)
c
      enddo
c
      maxNeighbors=max(maxNeighbors,NBChildrenMaxMG(iLevel1))
c
      allocate(TransitionArray(iLevel))
      TransitionArray(1:iLevel-1)=ijbeginNOld(1:iLevel-1)
      deallocate(ijbeginNOld)
      allocate(ijbeginNOld(iLevel))
      ijbeginNOld(1:iLevel-1)=TransitionArray(1:iLevel-1)
      deallocate(TransitionArray)
c
      do i=1,iLevel
        ijbeginNOld(i)=ijbeginN(i)
      enddo
c
      call SetLevelIndices(iLevel)
c----------------------------------------------------------
      allocate(TransitionArray(iLevel1))
      TransitionArray=NElementFacesMG
      deallocate(NElementFacesMG)
      allocate(NElementFacesMG(iLevel))
      do i=1,iLevel1
        NElementFacesMG(i)=TransitionArray(i)
      enddo
      NElementFacesMG(iLevel)=3*maxNeighbors     
      deallocate(TransitionArray)
c
      deallocate(NBChildrenMaxMG)
      allocate(NBChildrenMaxMG(iLevel))
      NBChildrenMaxMG=NElementFacesMG
c
c--- Reallocate and reindex "Children"
c
      allocate(TransitionArray(size(Children)))
      TransitionArray=Children
      deallocate(Children)
c
      ij=0
      do i=1,iLevel
c
        ij=ij+NumberOfElementsMG(i)*NBChildrenMaxMG(i)
c
      enddo
      allocate(Children(ij))
c      
      do i=1,iLevel-1
        do j=1,NumberOfElementsMG(i)
c
          ij2=ijbeginE(i)+j
c
          do k=1,NumberOfChildren(ij2)
c
            ij1=ijbeginN(i)+j+(k-1)*NumberOfElementsMG(i)
            ij3=ijbeginNOld(i)+j+(k-1)*NumberOfElementsMG(i)
            Children(ij1)=TransitionArray(ij3) 
c
          enddo
c
        enddo
      enddo          
c      
      deallocate(TransitionArray)
c
c--- Reallocate and reindex "ElementNeighborMG"
c
      allocate(TransitionArray(size(ElementNeighborMG)))
      TransitionArray=ElementNeighborMG
      deallocate(ElementNeighborMG)
c
      ij=0
      do i=1,iLevel
c
        ij=ij+NumberOfElementsMG(i)*NBChildrenMaxMG(i)
c
      enddo       
      allocate(ElementNeighborMG(ij))      
c      
      do i=1,iLevel
        do j=1,NumberOfElementsMG(i)
c
          ij2=ijbeginE(i)+j
c
          do k=1,NumberofElementNeighborsMG(ij2)
c
            ij1=ijbeginN(i)+j+(k-1)*NumberOfElementsMG(i)
            ij3=ijbeginNOld(i)+j+(k-1)*NumberOfElementsMG(i)
            ElementNeighborMG(ij1)=TransitionArray(ij3) 
c
          enddo
c
        enddo
      enddo          
c
      deallocate(TransitionArray)
c
c--- Reallocate and reindex "anbMG"
c
      allocate(TransitionArrayR(size(anbMG)))
      TransitionArrayR=anbMG
      deallocate(anbMG)
c
      ij=0
      do i=1,iLevel
c
        ij=ij+NumberOfElementsMG(i)*NBChildrenMaxMG(i)
c
      enddo
      allocate(anbMG(ij))
      anbMG=0.
c      
      do i=1,iLevel-1
        do j=1,NumberOfElementsMG(i)
c
          ij2=ijbeginE(i)+j
c
          do k=1,NumberofElementNeighborsMG(ij2)
c
            ij1=ijbeginN(i)+j+(k-1)*NumberOfElementsMG(i)
            ij3=ijbeginNOld(i)+j+(k-1)*NumberOfElementsMG(i)
            anbMG(ij1)=TransitionArrayR(ij3) 
c
          enddo
c
        enddo
      enddo          
c      
      deallocate(TransitionArrayR)
c
c--- Deallocate and reallocate maxanb
c
      deallocate(maxanb)
      allocate(maxanb(NumberOfElementsMG(iLevel)))
c
c--- Deallocate and reallocate Parents
c
      ij3=0
      do i=1,iLevel
c
      ij3=ij3+NumberOfElementsMG(i)
c
      enddo
      allocate(TransitionArray(ij3))
c      
      do i=1,iLevel1    
        do j=1,NumberOfElementsMG(i)
c
          ij=ijbeginE(i)+j           
          TransitionArray(ij)=parents(ij)    
c
        enddo
      enddo
c
      deallocate(parents)
      allocate(parents(ij3))
      parents=TransitionArray
c
c--- Deallocate and reallocate NumberOfChildren
c
      do i=1,iLevel1    
        do j=1,NumberOfElementsMG(i)
          ij=ijbeginE(i)+j           
          TransitionArray(ij)=NumberOfChildren(ij)    
        enddo
      enddo
c
      deallocate(NumberOfChildren)
      allocate(NumberOfChildren(ij3))
      NumberOfChildren=TransitionArray
c
c--- Deallocate and reallocate ElementParent
c
      do i=1,iLevel1    
        do j=1,NumberOfElementsMG(i)
          ij=ijbeginE(i)+j           
          TransitionArray(ij)=ElementParent(ij)    
        enddo
      enddo
c
      deallocate(ElementParent)
      allocate(ElementParent(ij3))
      ElementParent=TransitionArray
c
      deallocate(TransitionArray)
c
      return
      end