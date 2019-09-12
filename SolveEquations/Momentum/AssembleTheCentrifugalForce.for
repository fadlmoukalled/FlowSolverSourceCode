c
C#############################################################################################
c
      SUBROUTINE AssembleCentrifugalForceSourceTerm(Variable)
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Variables3, only: FluxTE
      use Centrifugal1
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer i
c********************************************************************************************
c
      if(Variable.eq.'velx') then
c
        do i=1,NumberOfElements
c
          FluxTE(i)=FluxTE(i)+xCentrifugal(i)
c
        enddo
c
      elseif(Variable.eq.'vely') then
c
        do i=1,NumberOfElements
c
          FluxTE(i)=FluxTE(i)+yCentrifugal(i)
c
        enddo
c
      elseif(Variable.eq.'velz') then
c
        do i=1,NumberOfElements
c
          FluxTE(i)=FluxTE(i)+zCentrifugal(i)
c
        enddo
c
      endif
c
      return
      end
