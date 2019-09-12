c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentV2Sources
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TurbulentED,
     *                      TfRelaxation,TurbulentV2
      use Variables3, only: FLuxCE,FLuxTE
      use PhysicalProperties1, only: Density
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        term1=Density(i)*dmax1(TurbulentKE(i)*TfRelaxation(i),0.)
        FLuxTE(i)=FLuxTE(i)-term1*Volume(i)
c
        term1=dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
        term1=6.*term1*Density(i)*Volume(i)
        FLuxCE(i)=FLuxCE(i)+term1
        FLuxTE(i)=FLuxTE(i)+term1*dmax1(TurbulentV2(i),0.)
c
      enddo
c
      return
      end