C
C############################################################################################
      SUBROUTINE CheckConvergence(NF,AbsoluteResidual)
C############################################################################################
C
      use User0
      use Geometry1, only: NumberOfElements
      use Geometry3, only: ElementNeighbor,NumberofElementNeighbors
      use Variables2, only: ac,anb,bc,dphi
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,NF
      double precision AbsoluteResidual,ElementResidual
c********************************************************************************************
       AbsoluteResidual=0.
       ElementResidual=0.    
C
C---- Calculate all types of residuals	
C
	do i=1,NumberOfElements
C
        ElementResidual=bc(i)-ac(i)*dphi(i)
c
        do j = 1,NumberofElementNeighbors(i)
c
          k = ElementNeighbor(i,j)
          if(k.ne.0) then
c
            ElementResidual=ElementResidual-anb(i,j)*dphi(k)
c
          endif

        enddo
C
        AbsoluteResidual=AbsoluteResidual+dabs(ElementResidual)
C
	enddo
C	
      return
      end