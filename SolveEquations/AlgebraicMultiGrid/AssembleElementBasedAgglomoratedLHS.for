c
C#############################################################################################
c
      SUBROUTINE AssembleAgglomeratedLHS(iLevel)
c
C#############################################################################################
c
      use Multigrid1, only:ijbeginE,ijbeginN,acMG,anbMG,ElementParent,
     *                    NumberofElementNeighborsMG,NumberOfElementsMG,
     *                    ElementNeighborMG,Parents,Children,
     *                    NumberOfChildren
c
      use Variables2, only:ac,anb
      use Geometry3, only:NumberofElementNeighbors
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer iLevel,iLevel1
      integer iParent,iChildren,ijChildren,ijChildLocal,iNBParent
      integer ij1,ij2,ij3,ij4,j,ijC,ijFine1,ijFine2,ijFine3,ijFine
      integer iNBlocal,iNB1,iNB,iNBP,k1,k2,k3,k3local
      integer i,k,ijC1,ijCN,iNC1
c*********************************************************************************************
c
      iLevel1=iLevel-1
c
c--- Coarse grid coefficient calculation
c
      do iParent=1,NumberOfElementsMG(iLevel)
c
        ij1=ijbeginE(iLevel)+iParent
        ij2=ijbeginN(iLevel)+iParent
c
        acMG(ij1)=0.
c
        do j=1,NumberofElementNeighborsMG(ij1)
c
          ij3=ij2+(j-1)*NumberOfElementsMG(iLevel)
          anbMG(ij3)=0.
c
        enddo
c
      enddo
c      
      do iParent=1,NumberOfElementsMG(iLevel)
c
        ijC=ijbeginE(iLevel)+iParent
c
        ijFine1=ijbeginE(iLevel1)+iParent
        ijFine2=ijbeginN(iLevel1)+iParent
        ijFine=ElementParent(ijFine1)
c
        acMG(ijC)=acMG(ijC)+acMG(ijFine)
c
        do i=1,NumberofElementNeighborsMG(ijFine)
c
          iNB1=ijbeginN(iLevel1)-ijbeginE(iLevel1)+ijFine+
     *                     (i-1)*NumberOfElementsMG(iLevel1)
          iNB=ElementNeighborMG(iNB1)
c
          if(iNB.gt.0) then
c
            iNBP=ijbeginE(iLevel1)+iNB
            iNBParent=Parents(iNBP)
c
            if(iNBParent.gt.0) then
c
              if(iNBParent.eq.iParent) then
c
                acMG(ijC)=acMG(ijC)+anbMG(iNB1)
c                           
              elseif(iNBParent.ne.iParent) then
c                      
                k1=NumberofElementNeighborsMG(ijC)
c
                do k2=1,k1
c
                  k3local=ijbeginN(iLevel)+iParent+
     *                        (k2-1)*NumberOfElementsMG(iLevel)
                  k3=ElementNeighborMG(k3local)
c
                  if(k3.eq.iNBParent) then
c
                    anbMG(k3local)=anbMG(k3local)+anbMG(iNB1) 
c
                  endif
                enddo
              endif
            endif
          endif
        enddo
c
        iNC1=NumberOfChildren(ijFine1)
        do j=1,iNC1
c
          ijFine3=ijFine2+(j-1)*NumberOfElementsMG(iLevel1)
          iChildren=Children(ijFine3)
c
          if(iChildren.gt.0) then
c
            ijChildren=ijbeginE(iLevel1)+iChildren
c
            acMG(ijC)=acMG(ijC)+acMG(ijChildren)
c   
            iNB1=NumberofElementNeighborsMG(ijChildren)
c         
            if(iNB1.gt.0) then
c
              do iNBlocal=1,iNB1
c
                ijChildLocal=ijbeginN(iLevel1)+iChildren+
     *                   (iNBlocal-1)*NumberOfElementsMG(iLevel1)
                iNB=ElementNeighborMG(ijChildLocal)
c
                if(iNB.gt.0) then
c
                  iNBP=ijbeginE(iLevel1)+iNB
                  iNBParent=Parents(iNBP)
c
                  if(iNBParent.gt.0) then
c
                    if(iNBParent.eq.iParent) then
c
                      acMG(ijC)=acMG(ijC)+anbMG(ijChildLocal)
c                           
                    elseif(iNBParent.ne.iParent) then
c                      
                      k1=NumberofElementNeighborsMG(ijC)
c
                      do k2=1,k1
c
                        k3local=ijbeginN(iLevel)+iParent+
     *                        (k2-1)*NumberOfElementsMG(iLevel)
                        k3=ElementNeighborMG(k3local)
c
                        if(k3.eq.iNBParent) then
c
                          anbMG(k3local)=anbMG(k3local)+
     *                                         anbMG(ijChildLocal)
c
                        endif
                      enddo
                    endif
                  endif
                endif
              enddo
            endif
          endif
        enddo
      enddo
c
      return
      end