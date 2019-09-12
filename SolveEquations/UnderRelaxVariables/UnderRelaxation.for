c
C#############################################################################################
c
      SUBROUTINE UnderRelaxEquation(urf)
c
C#############################################################################################
c
      use Variables2, only: ac
      use Variables4
      use Geometry1, only: NumberOfElements
      use Geometry3, only: NumberofElementNeighbors
c********************************************************************************************
      implicit none
c********************************************************************************************
      double precision :: urf      
      integer i,j
c********************************************************************************************
c
      do I=1,NumberOfElements
c
        ac(i)=ac(i)/urf
c
      enddo
c
      return
      end
