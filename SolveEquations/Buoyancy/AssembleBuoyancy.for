c
c#############################################################################################
c
      SUBROUTINE AssembleBuoyancyTerm(Variable)
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables3, only: FluxTE
      use Variables1, only: Buoyancyx,Buoyancyy,Buoyancyz,
     *                      uVelocity,vVelocity,wVelocity
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i
      double precision :: FluxTELocal
c********************************************************************************************
c
      if(Variable.eq.'velx') then
c
        do i=1,NumberOfElements
c
          FluxTELocal=Buoyancyx(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
        enddo
c
      elseif(Variable.eq.'vely') then
c
        do i=1,NumberOfElements
c
          FluxTELocal=Buoyancyy(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
        enddo
c
      elseif(Variable.eq.'velz') then
c
        do i=1,NumberOfElements
c
          FluxTELocal=Buoyancyz(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
        enddo
c
      elseif(Variable.eq.'htotal'.or.Variable.eq.'temperature') then
c
        do i=1,NumberOfElements
c
          FluxTELocal=(Buoyancyx(i)*uVelocity(i)+
     *                   Buoyancyy(i)*vVelocity(i)+
     *                     Buoyancyz(i)*wVelocity(i))*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
        enddo
c
      endif
c
      return
      end