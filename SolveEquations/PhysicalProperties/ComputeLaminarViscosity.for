c
c#############################################################################################
c
      SUBROUTINE CalculateLaminarViscosity
c
C#############################################################################################

      Use User0, only: urfViscosity,LSolveEnergy,Lcompressible,
     *                 Linviscid,ConstantViscosity,LSolveMomentum,
     *                 LConvectScalar
      Use PhysicalProperties1, only: SoutherlandLambda,SoutherlandC,
     *                               Viscosity,BViscosity,LSoutherland
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: Temperature,BTemperature
      use MultiGrid2, only: nIter
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: VisOld
c********************************************************************************************
c
      if(Linviscid) return
      if(.not.LSolveMomentum.and..not.LConvectScalar) return
c
      if(.not.LSolveEnergy.or..not.LSoutherland) then
c
        if(nIter.gt.1) return
        do i=1,NumberOfElements
c
          Viscosity(i)=ConstantViscosity
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BViscosity(i,j)=ConstantViscosity
c
          enddo 
        enddo 
c
      else
c
        do i=1,NumberOfElements
c
          VisOld=Viscosity(i)
          Viscosity(i)=(1.-urfViscosity)*VisOld+
     *       urfViscosity*SoutherlandLambda*
     *                (dmax1(Temperature(i),tiny))**1.5/
     *                         (Temperature(i)+SoutherlandC)
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            VisOld=BViscosity(i,j)
            BViscosity(i,j)=(1.-urfViscosity)*VisOld+
     *       urfViscosity*SoutherlandLambda*
     *                  (dmax1(BTemperature(i,j),tiny))**1.5/
     *                          (BTemperature(i,j)+SoutherlandC)
c
          enddo 
        enddo 
c
      endif
c
      return
      end