c
C#############################################################################################
c
      SUBROUTINE storeDvalues(dphi1,dphi2,dphi3,dphi4)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry3, only: NumberofElementNeighbors 
      use Geometry4, only: Volume
      use Variables2, only: ac,anb
      use User0, only: Algorithm
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision, dimension(:) :: dphi1,dphi2,dphi3,dphi4
c********************************************************************************************
c
	if(algorithm.eq.'simple') then
c
	  do i=1,NumberOfElements
c
	    dphi1(i) = Volume(i)/ac(i)
	    dphi2(i)=dphi1(i)
	    dphi4(i)=dphi3(i)
c
	  enddo
c
	elseif(algorithm.eq.'simplec') THEN
c
	  do i=1,NumberOfElements
c
	    dphi1(i) = Volume(i)/ac(i)
          dphi2(i)=ac(i)
c
          do j=1,NumberofElementNeighbors(i)
c
            dphi2(i)=dphi2(i)+anb(i,j)
c
          enddo
c
	    dphi2(i)=Volume(i)/dphi2(i)
	    dphi4(i)=dphi3(i)
c
	  enddo
c
	endif
c
      return
      end
