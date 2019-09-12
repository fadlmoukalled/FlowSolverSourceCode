c
C#############################################################################################
c
      SUBROUTINE AgglomorateLevelAlgebraicly(iLevel)
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
c
      integer, dimension(:), allocatable :: TransitionArray
      double precision, dimension(:), allocatable :: TransitionArrayR
c
      data ratiomin/0.5/
c*********************************************************************************************
c
      iLevel1=iLevel-1
c
      NumberOfFineElements=NumberOfElementsMG(iLevel1)
c
      do i=1,NumberOfFineElements
        ij=ijbeginE(iLevel1)+i
          maxanb(i)=-1.
        do j=1,NumberofElementNeighborsMG(ij)
          ij1=ijbeginN(iLevel1)+i+(j-1)*NumberOfFineElements
c
          maxanb(i)=dmax1(maxanb(i),-anbMG(ij1))
c
        enddo
      enddo
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
                ratio=-anbMG(iSeed2)/dmax1(maxanb(iSeed1),tiny)
c
                  if(ratio.gt.ratiomin.and.iNBchildren
     *                            .lt.maximumNumberOfChildren) then
c
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
            endif
c
          enddo
c
c--- Step 2 Agglomeration: neighbors of neighbors
c
          iNBchildren2=iNBchildren
c
          do iChildlocal=1,iNBchildren2
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
                iChildNBT=ijbeginN(iLevel1)+iChild+
     *                   (iChildNBlocal-1)*NumberOfFineElements
                iChildNB=ElementNeighborMG(iChildNBT)
c
                if(iChildNB.ne.0) then
c
                  if(parents(ijbeginE(iLevel1)+iChildNB).eq.0) then
c
                    ratio=-anbMG(iChildNBT)/dmax1(maxanb(iChild),tiny)
c
                    if(ratio.gt.ratiomin.and.iNBchildren
     *                               .lt.maximumNumberOfChildren) then
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
                endif
c
              enddo
c
            endif
c
          enddo
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
c--- Assign orphan elements to Parents
c
!      do iOrphan=1,NumberOfFineElements
!c
!        ij=ijbeginE(iLevel1)+iOrphan
!c
!        if(parents(ij).eq.0) then
!c
!          NChildren=1000000000
!c
!          do iNBlocal=1,NumberofElementNeighborsMG(ij)
!c
!            ij1=ijbeginN(iLevel1)+iOrphan+
!     *                (iNBlocal-1)*NumberOfFineElements
!            iNB=ElementNeighborMG(ij1)
!c
!            if(iNB.ne.0) then
!c
!              ij2=ijbeginE(iLevel1)+iNB
!c
!              if(parents(ij2).ne.0) then
!c
!                ij3=ijbeginE(iLevel1)+parents(ij2)
!c
!                if(NChildren.gt.NumberOfChildren(ij3)) then
!c
!                 NChildren=NumberOfChildren(ij3)
!                 iNBuse=iNB
!c
!                endif
!c
!              endif
!c
!            endif
!c
!          enddo
!c
!          iNB=iNBuse
!          iOfiNB=ijbeginE(iLevel1)+iNB
!          iOfiOrphan=ijbeginE(iLevel1)+iOrphan
!c
!          if(iNB.ne.0) then
!          if(parents(iOfiNB).ne.0) then
!c
!            parents(iOfiOrphan)=parents(iOfiNB)
!            ij3=ijbeginE(iLevel1)+Parents(iOfiNB)
!            itemp=NumberOfChildren(ij3)
!            ij4=ijbeginN(iLevel1)+Parents(iOfiNB)+
!     *                        itemp*NumberOfFineElements
!            children(ij4)=iOrphan
!            NumberOfChildren(ij3)=itemp+1
!            NumberOfChildren0=NumberOfChildren0+1
!c
!          endif
!          endif
!c
!        endif
!c
!        if(parents(ijbeginE(iLevel1)+iOrphan).eq.0) then
!c
!          iOfiOrphan=ijbeginE(iLevel1)+iOrphan
!          parents(iOfiOrphan) = iParent
!          ElementParent(ijbeginE(iLevel1)+iParent)=iOfiOrphan
!          iParent = iParent + 1
!          write(*,*) 'the orphan could not find a parent'
!c
!        endif
!c
!      enddo
c
c--- Last Step Agglomeration
c
      do iOrphan=1,NumberOfFineElements
        ij=ijbeginE(iLevel1)+iOrphan
        if(parents(ij).eq.0) then
          strength = 0.
          index=0
          do iNBlocal=1,NumberofElementNeighborsMG(ij)
c
            ij1=ijbeginN(iLevel1)+iOrphan+
     *                (iNBlocal-1)*NumberOfFineElements
            iNB=ElementNeighborMG(ij1)
            if(iNB.ne.0) then
              ij2=ijbeginE(iLevel1)+iNB
              if(parents(ij2).ne.0.and.index.ne.1) then
c     *        .and.
c     *   NumberOfChildren(ijbeginE(iLevel1)+parents(ij2)).lt.
c     *      maximumNumberOfChildren) then
                ratio=-anbMG(ij1)/maxanb(iNB)
                if(strength.lt.ratio) then
                  index=1
                  strength = ratio
                  ij2=ijbeginE(iLevel1)+iNB
                  parents(ij) = parents(ij2)
                  ij3=ijbeginE(iLevel1)+parents(ij2)
                  itemp=NumberOfChildren(ij3)
                  ij4=ijbeginN(iLevel1)+Parents(ijbeginE(iLevel1)+iNB)+
     *                        itemp*NumberOfFineElements
                  children(ij4)=iOrphan
                  NumberOfChildren(ij3)=itemp+1
                  NumberOfChildren0=NumberOfChildren0+1
                endif
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
!c
!c--- Set the minimum and maximum agglomerated element per level
!c
!      OldNumberOfParents=NumberOfParents
!      minChildrenNew=0
!      minChildrenOld=minChildren
!c          
!c----Start diffusing parents with low number of children           
!c
!      do while(minChildrenNew.ne.minChildrenOld)  
!c
!        DroppedParents=1
!c
!        do while(DroppedParents.ne.0) 
!c
!          do iParent=1,NumberOfParents
!c
!            ij1=ijbeginE(iLevel1)+iParent
!            k=NumberOfChildren(ij1)
!c
!            if(k.ne.0) then
!c
!              if(k.le.minChildren.and.k.lt.minimumNumberOfChildren) then
!c
!                j=1
!                do while(j.le.k)
!c
!                  ij2=ijbeginN(iLevel1)+
!     *                   iParent+(j-1)*NumberOfFineElements
!                  ichildren=Children(ij2)
!                  agglomorated=.false.
!                  ij3=ijbeginE(iLevel1)+ichildren
!                  iNB1=NumberofElementNeighborsMG(ij3)
!c
!                  if(iNB1.gt.0) then
!c
!                    do iNBlocal=1,iNB1
!c
!                      if(.not.agglomorated) then
!c
!                        ij4=ijbeginN(iLevel1)+ichildren+
!     *                         (iNBlocal-1)*NumberOfFineElements
!                        iNB=ElementNeighborMG(ij4)
!c
!                        if(iNB.gt.0) then
!c
!                          iNBParent=Parents(ijbeginE(iLevel1)+iNB)
!c
!                          if(iNBParent.gt.0) then
!c
!                            if(iNBParent.ne.iParent) then
!c
!                              ij5=ijbeginE(iLevel1)+iNBParent
!                              iNBChildren=NumberOfChildren(ij5)
!c
!                              if(iNBChildren.gt.0.and.iNBChildren.lt.
!     *                                  maximumNumberOfChildren+2) then
!c
!                                indexn=0
!                                do while(indexn.eq.0)
!                                  agglomorated=.true.
!                                  Parents(ij3)=iNBParent
!                                  k1=NumberOfChildren(ij5)
!                                  NumberOfChildren(ij5)=k1+1
!                                  NumberOfChildren(ij1)=
!     *                                   NumberOfChildren(ij1)-1
!                                  ij6=ijbeginN(iLevel1)+iNBParent+
!     *                                         k1*NumberOfFineElements
!                                  Children(ij6)=ichildren
!c
!c--- Reindex the remaining children---------------------------------------------------------                               
!c                                 
!                                  do jn=j+1,k
!                                    ij2n1=ijbeginN(iLevel1)+
!     *                               iParent+(jn-2)*NumberOfFineElements
!                                    ij2n=ijbeginN(iLevel1)+
!     *                               iParent+(jn-1)*NumberOfFineElements
!                                    ichildren=Children(ij2)
!                                    ij3=ijbeginE(iLevel1)+ichildren
!                                    children(ij2n1)=children(ij2n)
!                                  enddo
!                                  ij2n=ijbeginN(iLevel1)+
!     *                               iParent+(k-1)*NumberOfFineElements
!                                  children(ij2n)=0
!                                  j=j-1
!                                  k=k-1
!c
!c--- Make sure no connection between the parent and any of the children is lost ------------
!c
!                                  n1loop: do n1=1,k
!                                    ijn1=ijbeginN(iLevel1)+
!     *                               iParent+(n1-1)*NumberOfFineElements
!                                    ichildren=Children(ijn1)
!                                    n2=NumberofElementNeighborsMG(ij1)
!                                    
!                                    indexn=0
!                                    n3loop: do n3=1,n2
!                                      ijn2=ijbeginN(iLevel1)+iParent+
!     *                                      (n3-1)*NumberOfFineElements
!                                      n4=ElementNeighborMG(ijn2)
!                                      if(n4.ne.0) then
!                                        ijn3=ijbeginE(iLevel1)+n4
!                                        if(ijn3.eq.ijn1) then
!                                          indexn=1
!                                          exit n3loop
!                                        endif
!                                      endif
!                                    enddo n3loop
!c                                  
!                                    if(indexn.eq.0) then 
!c                                  
!                                      n5loop: do n5=1,k
!                                        if(n5.ne.n1) then
!                                          ijn4=ijbeginN(iLevel1)+iParent
!     *                                      +(n5-1)*NumberOfFineElements
!                                          n6=ElementNeighborMG(ijn4)
!                                          if(n6.ne.0) then
!                                            ijn5=ijbeginE(iLevel1)+n6
!                                            n7=
!     *                                       NumberofElementNeighborsMG
!     *                                                     (ijn5)
!                                            n8loop: do n8=1,n7
!                                              ijn6=ijbeginN(iLevel1)+n6+
!     *                                       (n8-1)*NumberOfFineElements
!                                              n9=ElementNeighborMG(ijn6)
!                                              if(n9.ne.0) then
!                                                ijn7=
!     *                                             ijbeginE(iLevel1)+n9
!                                                if(ijn7.eq.ijn1) then
!                                                  indexn=1
!                                                  exit n8loop 
!                                                endif
!                                              endif
!                                            enddo n8loop
!                                          endif
!                                        endif        
!                                    
!                                      enddo n5loop
!c                                  
!                                    endif
!c                                  
!                                    if(indexn.eq.0) then
!                                      exit n1loop
!                                    endif
!c                                  
!                                  enddo n1loop
!                                  indexn=1
!                                enddo  
!c
!c-------------------------------------------------------------------------------------------                                
!c
!                                if(NumberOfChildren(ij1).eq.0) then
!c
!                                  k3=NumberOfChildren(ij5)
!                                  NumberOfChildren(ij5)=k3+1
!                                  i=ElementParent(ij1)
!                                  ElementParent(ij1)=0
!                                  ij7=ijbeginN(iLevel1)+iNBParent+
!     *                                       k3*NumberOfFineElements
!                                  Children(ij7)=i-ijbeginE(iLevel1)
!                                  Parents(i)=iNBParent
!                                  NumberOfParents=NumberOfParents-1
!c
!                                endif
!c
!                              endif  
!c
!                            endif
!c
!                          endif
!c
!                        endif
!c
!                      endif
!c
!                    enddo
!c
!                  endif
!c
!                  j=j+1
!                enddo  
!c
!              endif    
!c
!            endif    
!c
!          enddo
!c
!          index=1
!c
!          do while (index.eq.1)
!c
!            do i=1,OldNumberOfParents
!c
!              if(ElementParent(ijbeginE(iLevel1)+i).le.0) then
!c
!                do k3=1,NumberOfFineElements
!c
!                  ij=ijbeginE(iLevel1)+k3
!c
!                  if(Parents(ij).gt.i) then
!c
!                    iParent=Parents(ij)
!                    Parents(ij)=iParent-1
!c
!                  endif
!c
!                enddo
!c        
!                do j=i+1,OldNumberOfParents      
!c
!                  ElementParent(ijbeginE(iLevel1)+j-1)=
!     *                            ElementParent(ijbeginE(iLevel1)+j)
!c
!                  k=NumberOfChildren(ijbeginE(iLevel1)+j)
!c
!                  do k1=1,k
!c
!                    ij1=ijbeginN(iLevel1)+j+(k1-1)*NumberOfFineElements
!                    k2=Children(ij1)
!                    ij2=ijbeginN(iLevel1)+
!     *                                 j-1+(k1-1)*NumberOfFineElements
!                    Children(ij2)=k2
!c
!                  enddo
!c
!                  NumberOfChildren(ijbeginE(iLevel1)+j-1)=k
!c
!                enddo
!c
!              endif
!c
!            enddo
!c        
!            index=0
!c
!            do k=1,NumberOfParents
!c
!              if(ElementParent(ijbeginE(iLevel1)+k).le.0) index=1
!c
!            enddo
!c
!          enddo
!c
!c--- Calculate the number of dropped Parents
!c
!          droppedParents=OldNumberOfParents-NumberOfParents
!c
!c----------------------------------------------------------------------------
!          do i=NumberOfParents+1,OldNumberOfParents
!c
!            ij=ijbeginE(iLevel1)+i
!c
!            k=NumberOfChildren(ij)
!c
!            do k1=1,k
!c
!              ij1=ijbeginN(iLevel1)+i+(k1-1)*NumberOfFineElements
!              Children(ij1)=0
!c
!            enddo
!            NumberOfChildren(ij)=0
!c
!          enddo
!c
!c-----------------------------------------------------------------------------       
!c
!          OldNumberOfParents=NumberOfParents
!c
!        enddo
!c
!c--- Check for the new number of minimum and maximum children
!c
!        maxChildren=-1
!        minChildren=100
!        do i=1,NumberOfParents
!          ij=ijbeginE(iLevel1)+i
!          j=NumberOfChildren(ij)
!          maxChildren=max(maxChildren,j)
!          minChildren=min(minChildren,j)
!        enddo
!        minChildrenNew=minChildrenOld
!        minChildrenOld=minChildren
!c
!      enddo
c          
c----End of diffusing parents with low number of children           
c
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
c      pause
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