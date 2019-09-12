c
C#############################################################################################
c
      SUBROUTINE AssembleCoriolisSourceTerm(phi1,phi2,phi3,phi4)
c
C#############################################################################################
c
      use PhysicalProperties1, only:Density
      use Geometry4, only: Volume
      use Geometry1, only: NumberOfElements
      use Variables3, only: FluxTE
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i
      double precision :: term
      double precision :: phi2,phi4
      double precision, dimension(:) :: phi1
      double precision, dimension(:) :: phi3
c********************************************************************************************
c
      do i=1,NumberOfElements     
c 
        term=2.*Density(i)*Volume(i)
        FluxTE(i)=FluxTE(i)+term*(phi1(i)*phi2-phi3(i)*phi4)
c
      enddo
c	
	return
	end
