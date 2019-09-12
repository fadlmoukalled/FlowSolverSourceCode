c
c#############################################################################################
c
      SUBROUTINE CalculateLaminarThermalConductivity
c
C#############################################################################################
c
      Use User0, only: urfConductivity,Linviscid,LSolveEnergy,
     *                 ConstantConductivity
      Use PhysicalProperties1, only: Viscosity,BViscosity,
     *                               PrLaminar,SpecificHeat,
     *                               BSpecificHeat,Conductivity,
     *                               BConductivity,LSoutherland
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: kOld
c********************************************************************************************
c
      if(Linviscid) return
      if(.not.LSolveEnergy) return
c
      if(.not.LSoutherland) then
c      
        if(nIter.gt.1) return
        do i=1,NumberOfElements
c
          Conductivity(i)=ConstantConductivity
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BConductivity(i,j)=ConstantConductivity
c
          enddo 
        enddo 
c      
      else
c
        do i=1,NumberOfElements
c
          kOld=Conductivity(i)
          Conductivity(i)=(1.-urfConductivity)*kOld+
     *       urfConductivity*Viscosity(i)*SpecificHeat(i)/PrLaminar
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            kOld=BConductivity(i,j)
            BConductivity(i,j)=(1.-urfConductivity)*kOld+
     *                     urfConductivity*BViscosity(i,j)*
     *                              BSpecificHeat(i,j)/PrLaminar
c
          enddo 
        enddo 
c
      endif
c
      return
      end