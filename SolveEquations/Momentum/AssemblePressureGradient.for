c
C#############################################################################################
c
      SUBROUTINE AssemblePressureGradientTerm(dPdxy)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables3, only: FluxTE
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i
      double precision, dimension(:) :: dPdxy
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        FluxTE(i)=FluxTE(i)+dPdxy(i)*Volume(i)
c
      enddo
c	
	return
	end
