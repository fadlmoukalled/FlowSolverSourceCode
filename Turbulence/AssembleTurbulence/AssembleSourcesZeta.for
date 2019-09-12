c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentZetaSources
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TfRelaxation,
     *                      TurbulentZeta,TurbulentZetaGradx,
     *                      TurbulentZetaGrady,TurbulentZetaGradz,
     *                      TKEGradx,TKEGrady,TKEGradz,
     *                      TurbulenceProduction
      use Variables3, only: FLuxCE,FLuxTE
      use PhysicalProperties1, only: Density,TurbulentViscosity
      use Turbulence1, only: sigTKE
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
        term1=Density(i)*dmax1(TfRelaxation(i),0.)
        FLuxTE(i)=FLuxTE(i)-term1*Volume(i)
c
        term1=TurbulenceProduction(i)*Volume(i)/
     *                         (dmax1(TurbulentKE(i),tiny))
        FLuxCE(i)=FLuxCE(i)+term1
        FLuxTE(i)=FLuxTE(i)+term1*dmax1(TurbulentZeta(i),0.)
c
        term1=2.*TurbulentViscosity(i)/
     *             (sigTKE*dmax1(TurbulentKE(i),tiny))
        term1=term1*Volume(i)*(
     *           TurbulentZetaGradx(i)*TKEGradx(i)+
     *           TurbulentZetaGrady(i)*TKEGrady(i)+
     *           TurbulentZetaGradz(i)*TKEGradz(i))
        FLuxCE(i)=FLuxCE(i)-
     *              dmin1(term1/dmax1(TurbulentZeta(i),tiny),0.)
        FLuxTE(i)=FLuxTE(i)-term1
c
      enddo
c
      return
      end