c
C#############################################################################################
c
      SUBROUTINE SolveEquationsUsingSORMG(iLevel)
c
C#############################################################################################
c
      use MultiGrid1, only:ijbeginE,ijbeginN,acMG,anbMG,bcMG,dphiMG,
     *                     NumberOfElementsMG,
     *                     NumberofElementNeighborsMG,ElementNeighborMG
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,ij,ij1,ij2,iLevel
      double precision localdphi
c********************************************************************************************
c
c-- Start computations from first to the last element
c
      do i=1,NumberOfElementsMG(iLevel)
c
        ij=ijbeginE(iLevel)+i
        localdphi = bcMG(ij)
c        
        do j = 1,NumberofElementNeighborsMG(ij)
c
          ij1=ijbeginN(iLevel)+i+(j-1)*NumberOfElementsMG(iLevel)
          k = ElementNeighborMG(ij1)
c
          if(k.ne.0) then
c
            ij2=ijbeginE(iLevel)+k
c
            localdphi = localdphi - anbMG(ij1)*dphiMG(ij2)
c
          endif
c
        enddo
c
        dphiMG(ij) = localdphi/acMG(ij)
c
      enddo
c
c-- Reverse the direction of computations from the last to the first element
c
      do i=NumberOfElementsMG(iLevel),1,-1
c
        ij=ijbeginE(iLevel)+i
        localdphi = bcMG(ij)
c        
        do j = 1,NumberofElementNeighborsMG(ij)
c
          ij1=ijbeginN(iLevel)+i+(j-1)*NumberOfElementsMG(iLevel)
          k = ElementNeighborMG(ij1)
c
          if(k.ne.0) then
c
            ij2=ijbeginE(iLevel)+k
c
            localdphi = localdphi - anbMG(ij1)*dphiMG(ij2)
c
          endif
c
        enddo
c
        dphiMG(ij) = localdphi/acMG(ij)
c
      enddo
c
      return
      end
