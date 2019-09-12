c
C#############################################################################################
c
      SUBROUTINE AssembleSourceTerm(Variable,FiT,Sc,Sb)
c
C#############################################################################################
c
      use Variables3
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i
      character*10 Variable
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: Sc
      double precision, dimension(:) :: Sb
c********************************************************************************************
c
c--- Assemble element fluxes
c      
      do i=1,NumberOfElements
c      
        if(SC(i).lt.0.) then
c
          FluxCE(i)=FluxCE(i)-Sc(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-(Sb(i)+Sc(i)*FiT(i))*Volume(i)
c
        else
c
          FluxCE(i)=FluxCE(i)+0.
          FluxTE(i)=FluxTE(i)-(Sb(i)+Sc(i)*FiT(i))*Volume(i)
c
        endif
c
      enddo
c
      return
      end