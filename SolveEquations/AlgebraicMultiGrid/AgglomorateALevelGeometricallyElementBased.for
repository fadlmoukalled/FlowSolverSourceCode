c
C#############################################################################################
c
      SUBROUTINE AgglomorateLevelGeometricallyElement(iLevel)
c
C#############################################################################################
c
      use User0
      use Constants1, only:tiny
      use Geometry1, only:NumberOfElements
      use MultiGrid1
      use reallocateData
      use Variables2, only:anb
      use Geometry3, only:NElementFaces,ElementNeighbor,
     *                    NumberofElementNeighbors
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
c*********************************************************************************************
c
      iLevel1=iLevel-1
c
      NumberOfFineElements=NumberOfElementsMG(iLevel1)
c
      do i=1,NumberOfFineElements
c
        ij=ijbeginE(iLevel1)+i
        parents(ij)=0
        NumberOfChildren(ij)=0
        ElementParent(ij)=0
c
        do j=1,NBChildrenMaxMG(iLevel1)
          ij1=ijbeginN(iLevel1)+i+(j-1)*NumberOfFineElements
          children(ij1)=0
        enddo
c
      enddo
c
c--- Set the minimum and maximum agglomerated element per level
c
c
      if(iLevel.eq.2) then
c
        minimumNumberOfChildren=3
        maximumNumberOfChildren=8
c
      endif
c
      if(iLevel.eq.3) then
c
        minimumNumberOfChildren=2
        maximumNumberOfChildren=4
c
      endif
c
      if(iLevel.eq.4) then
c
        minimumNumberOfChildren=1
        maximumNumberOfChildren=3
c
      endif
c
      if(iLevel.ge.5) then
c
        minimumNumberOfChildren=1
        maximumNumberOfChildren=2
c
      endif      
c        
      iParent = 1
      NumberOfChildren0=0
c
c--- Step 1 Agglomeration: immediate neighbors
c
      do iSeed1=1,NumberOfFineElements
c
        iSeed=ijbeginE(iLevel1)+iSeed1
        iNBchildren=0
c
        if(parents(iSeed).eq.0) then
c
          parents(iSeed) = iParent
          iNBlocal1=0
c
          do iNBlocal=1,NumberofElementNeighborsMG(iSeed)
c
            iSeed2=ijbeginN(iLevel1)+iSeed1+
     *             (iNBlocal-1)*NumberOfFineElements
            iNB=ElementNeighborMG(iSeed2)
c
            if(iNB.ne.0) then
c
              if(parents(ijbeginE(iLevel1)+iNB).eq.0) then
c

                if(iNBchildren.ge.maximumNumberOfChildren) exit

                iNBlocal1=iNBlocal1+1
                parents(ijbeginE(iLevel1)+iNB) = iParent
                ichildren=ijbeginN(iLevel1)+iParent+
     *                    (iNBlocal1-1)*NumberOfFineElements
                children(ichildren) = iNB
                iNBchildren=iNBchildren+1
                NumberOfChildren0=NumberOfChildren0+1
c
              endif
c
            endif
c
          enddo
c
c--- Step 2 Agglomeration: neighbors of neighbors
c
          iNBchildren2=iNBchildren
c
          if(iNBchildren.lt.maximumNumberOfChildren) then
c
            do iChildlocal=1,iNBchildren2
c
              if(iNBchildren.ge.maximumNumberOfChildren) exit
c
              iParentT=ijbeginN(iLevel1)+iParent+
     *              (iChildlocal-1)*NumberOfFineElements
              iChild = children(iParentT)
c
              if(iChild.ne.0) then
c
                iChildT=ijbeginE(iLevel1)+iChild
c
                do iChildNBlocal=
     *              1,NumberofElementNeighborsMG(iChildT)
c
                  if(iNBchildren.ge.maximumNumberOfChildren) exit
c
                  iChildNBT=ijbeginN(iLevel1)+iChild+
     *                   (iChildNBlocal-1)*NumberOfFineElements
                  iChildNB=ElementNeighborMG(iChildNBT)
c
                  if(iChildNB.ne.0) then
c
                    if(parents(ijbeginE(iLevel1)+iChildNB).eq.0) then
c
                      parents(ijbeginE(iLevel1)+iChildNB)=iParent
                      iNBchildren=iNBchildren+1
                      ichildren=ijbeginN(iLevel1)+iParent+
     *                    (iNBchildren-1)*NumberOfFineElements
                      children(ichildren)=iChildNB
                      NumberOfChildren0=NumberOfChildren0+1
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
          endif
c
          if(iNBchildren.eq.0) then
c
            parents(iSeed) = 0
c
          else
c
            ij=ijbeginE(iLevel1)+iParent
            ElementParent(ij)=iSeed
            NumberOfChildren(ij)=iNBchildren
            iParent = iParent + 1
c
          endif
c
        endif
c
      enddo
c
c--- Last Step Agglomeration
c
      do iOrphan=1,NumberOfFineElements
        ij=ijbeginE(iLevel1)+iOrphan
        if(parents(ij).eq.0) then
          index=0
          do iNBlocal=1,NumberofElementNeighborsMG(ij)
c
            ij1=ijbeginN(iLevel1)+iOrphan+
     *                (iNBlocal-1)*NumberOfFineElements
            iNB=ElementNeighborMG(ij1)
            if(iNB.ne.0) then
              ij2=ijbeginE(iLevel1)+iNB
              if(parents(ij2).ne.0.and.index.ne.1) then
                index=1
                ij2=ijbeginE(iLevel1)+iNB
                ij3=ijbeginE(iLevel1)+parents(ij2)
c
                if(NumberOfChildren(ij3).ge.
     *                  maximumNumberOfChildren) cycle
c
                parents(ij) = parents(ij2)
                itemp=NumberOfChildren(ij3)
                ij4=ijbeginN(iLevel1)+Parents(ijbeginE(iLevel1)+iNB)+
     *                        itemp*NumberOfFineElements
                children(ij4)=iOrphan
                NumberOfChildren(ij3)=itemp+1
                NumberOfChildren0=NumberOfChildren0+1
              endif
            endif
          enddo
        endif
c
      enddo
c
c---- Repeat without constraints
c
      do iOrphan=1,NumberOfFineElements
        ij=ijbeginE(iLevel1)+iOrphan
        if(parents(ij).eq.0) then
          index=0
          do iNBlocal=1,NumberofElementNeighborsMG(ij)
c
            ij1=ijbeginN(iLevel1)+iOrphan+
     *                (iNBlocal-1)*NumberOfFineElements
            iNB=ElementNeighborMG(ij1)
            if(iNB.ne.0) then
              ij2=ijbeginE(iLevel1)+iNB
              if(parents(ij2).ne.0.and.index.ne.1) then
                index=1
                ij2=ijbeginE(iLevel1)+iNB
                ij3=ijbeginE(iLevel1)+parents(ij2)
                parents(ij) = parents(ij2)
                itemp=NumberOfChildren(ij3)
                ij4=ijbeginN(iLevel1)+Parents(ijbeginE(iLevel1)+iNB)+
     *                        itemp*NumberOfFineElements
                children(ij4)=iOrphan
                NumberOfChildren(ij3)=itemp+1
                NumberOfChildren0=NumberOfChildren0+1
              endif
            endif
          enddo
        endif
c
        if(parents(ijbeginE(iLevel1)+iOrphan).eq.0) then
          parents(ijbeginE(iLevel1)+iOrphan) = iParent
          iParent = iParent + 1
          write(*,*) 'the orphan could not find a parent'
        endif
      enddo
c
c-------------------------------------------------------------------------------------
c
      NumberOfParents = iParent-1
c
c--- Calculate the current minimum and maximum number of children
c
      maxChildren=-1
      minChildren=100
      do i=1,NumberOfParents
        ij=ijbeginE(iLevel1)+i
        j=NumberOfChildren(ij)
        maxChildren=max(maxChildren,j)
        minChildren=min(minChildren,j)
      enddo
c
c--- Calculate the exact maximum number of children and save it
c
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
