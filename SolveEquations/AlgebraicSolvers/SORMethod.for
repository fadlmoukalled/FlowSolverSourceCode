c
C#############################################################################################
c
      SUBROUTINE SolveEquationsUsingSOR
c
C#############################################################################################
c
      use Variables2, only: ac,anb,bc,dphi
      use Geometry1, only: NumberOfElements
      use Geometry3, only: ElementNeighbor,NumberofElementNeighbors
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k
      double precision localdphi
c********************************************************************************************
c
c-- Start computations from first to the last element
c
      do i=1,NumberOfElements
c
        localdphi = bc(i)
c        
        do j = 1,NumberofElementNeighbors(i)
c
          k = ElementNeighbor(i,j)
c
          if(k.ne.0) then
c
            localdphi = localdphi - anb(i,j)*dphi(k)
c
          endif
c
        enddo
c
        dphi(i) = localdphi/ac(i)
c
      enddo
c
c-- Reverse the direction of computations from the last to the first element
c
      do i=NumberOfElements,1,-1
c
        localdphi = bc(i)
c
        do j = 1,NumberofElementNeighbors(i)
c
          k = ElementNeighbor(i,j)
c
          if(k.ne.0) then
c
            localdphi = localdphi - anb(i,j)*dphi(k)
c
          endif
c
        enddo
c
        dphi(i) = localdphi/ac(i)
c
      enddo
c
      return
      end