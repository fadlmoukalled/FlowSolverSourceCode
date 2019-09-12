c
c#############################################################################################
c
      SUBROUTINE CalculateDiffusionCoefficientScalar(iScalar)
c
C#############################################################################################
c
      Use User0, only: urfConductivity,Linviscid,
     *                 ConstantDiffusionCoefficientScalar
      Use PhysicalProperties1, only: Viscosity,BViscosity,
     *                               PrLaminar,SpecificHeatScalar,
     *                               BSpecificHeatScalar,
     *                               DiffusionCoefficient,
     *                               BDiffusionCoefficient,
     *                               LSoutherlandScalar
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,iScalar
      double precision :: kOld
c********************************************************************************************
c      
      if(Linviscid) return
c
      if(.not.LSoutherlandScalar(iScalar)) then
c      
        if(nIter.gt.1) return
        do i=1,NumberOfElements
c
          DiffusionCoefficient(i,iScalar)=
     *               ConstantDiffusionCoefficientScalar(iScalar)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BDiffusionCoefficient(i,j,iScalar)=
     *                  ConstantDiffusionCoefficientScalar(iScalar)
c
          enddo 
        enddo 
c      
      else
c
        do i=1,NumberOfElements
c
          kOld=DiffusionCoefficient(i,iScalar)
          DiffusionCoefficient(i,iScalar)=(1.-urfConductivity)*kOld+
     *       urfConductivity*Viscosity(i)*
     *          SpecificHeatScalar(i,iScalar)/PrLaminar
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            kOld=BDiffusionCoefficient(i,j,iScalar)
            BDiffusionCoefficient(i,j,iScalar)=
     *         (1.-urfConductivity)*kOld+
     *           urfConductivity*BViscosity(i,j)*
     *             BSpecificHeatScalar(i,j,iScalar)/PrLaminar
c
          enddo 
        enddo 
c
      endif
c
      return
      end
