c
C#############################################################################################
c
      SUBROUTINE AssembleFalseTransientrField(Variable,dtfalse)
c
C#############################################################################################
c
      use Variables3, only: FluxCE
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i,j
      character*10 Variable
      double precision dtfalse
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        FluxCE(i)=FluxCE(i)+Volume(i)/dtfalse
c
      enddo
c
      return
      end
