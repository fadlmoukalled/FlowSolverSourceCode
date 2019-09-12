c
c#############################################################################################
c
      SUBROUTINE CalculateLaminarViscosityrField(irField)
c
C#############################################################################################

      Use User0, only: urfViscosity,LSolveEnergy,Lcompressible,
     *                 Linviscid,ConstantViscosityrField
      Use PhysicalProperties1, only: SoutherlandLambdarField,
     *                               SoutherlandCrField,
     *                               LSoutherlandrField,
     *                               ViscosityrField,BViscosityrField
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: Temperature,BTemperature
      use MultiGrid2, only: nIter
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,irField
      double precision :: VisOld
c********************************************************************************************
c      
      if(Linviscid) return
c
      if(.not.LSolveEnergy.or..not.LSoutherlandrField(irField)) then
c      
        if(nIter.gt.1) return
        do i=1,NumberOfElements
c
          ViscosityrField(i,irField)=ConstantViscosityrField(irField)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BViscosityrField(i,j,irField)=
     *                  ConstantViscosityrField(irField)
c
          enddo 
        enddo 
c      
      else
c
        do i=1,NumberOfElements
c
          VisOld=ViscosityrField(i,irField)
          ViscosityrField(i,irField)=(1.-urfViscosity)*VisOld+
     *       urfViscosity*SoutherlandLambdarField(irField)*
     *           (dmax1(Temperature(i),tiny))**1.5/
     *            (Temperature(i)+SoutherlandCrField(irField))
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            VisOld=BViscosityrField(i,j,irField)
            BViscosityrField(i,j,irField)=(1.-urfViscosity)*VisOld+
     *       urfViscosity*SoutherlandLambdarField(irField)*
     *         (dmax1(BTemperature(i,j),tiny))**1.5/
     *            (BTemperature(i,j)+SoutherlandCrField(irField))
c
          enddo 
        enddo 
c
      endif
c
      return
      end