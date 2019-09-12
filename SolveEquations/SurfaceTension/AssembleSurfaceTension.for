c
c#############################################################################################
c
      SUBROUTINE AssembleSurfaceTensionSourceTerm(Variable)
c
c#############################################################################################
c
      use User0, only: SurfaceTension
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables3, only: FluxTE
      use Variables1, only: uVelocity,vVelocity,wVelocity
      use VolumeOfFluid1, only: Curvature
      use TransferrField1, only: rFieldGradxT,rFieldGradyT,rFieldGradzT
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
          FluxTELocal=
     *        SurfaceTension*Curvature(i)*rFieldGradxT(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
        enddo
c
      elseif(Variable.eq.'vely') then
c
        do i=1,NumberOfElements
c
          FluxTELocal=
     *        SurfaceTension*Curvature(i)*rFieldGradyT(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
        enddo
c
      elseif(Variable.eq.'velz') then
c
        do i=1,NumberOfElements
c
c
          FluxTELocal=
     *        SurfaceTension*Curvature(i)*rFieldGradzT(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
c
        enddo
c
      elseif(Variable.eq.'htotal'.or.Variable.eq.'temperature') then
c
        do i=1,NumberOfElements
c
          FluxTELocal=SurfaceTension*Curvature(i)*
     *                 (rFieldGradxT(i)*uVelocity(i)+
     *                   rFieldGradyT(i)*vVelocity(i)+
     *                     rFieldGradzT(i)*wVelocity(i))*Volume(i)
          FluxTE(i)=FluxTE(i)-FluxTELocal
c
        enddo
c
      endif
c
      return
      end
