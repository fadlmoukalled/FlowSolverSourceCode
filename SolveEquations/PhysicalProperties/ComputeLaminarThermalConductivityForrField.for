c
c#############################################################################################
c
      SUBROUTINE CalculateLaminarThermalConductivityrField(irField)
c
C#############################################################################################
c
      Use User0, only: urfConductivity,Linviscid,LSolveEnergy,
     *                 ConstantConductivityrField
      Use PhysicalProperties1, only: ViscosityrField,BViscosityrField,
     *                               PrLaminarrField,SpecificHeatrField,
     *                               BSpecificHeatrField,
     *                               ConductivityrField,
     *                               BConductivityrField,
     *                               LSoutherlandrField
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,irField
      double precision :: kOld
c********************************************************************************************
c      
      if(Linviscid) return
      if(.not.LSolveEnergy) return
c
      if(.not.LSoutherlandrField(irField)) then
c      
        if(nIter.gt.1) return
        do i=1,NumberOfElements
c
          ConductivityrField(i,irField)=
     *               ConstantConductivityrField(irField)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BConductivityrField(i,j,irField)=
     *                  ConstantConductivityrField(irField)
c
          enddo 
        enddo 
c      
      else
c
        do i=1,NumberOfElements
c
          kOld=ConductivityrField(i,irField)
          ConductivityrField(i,irField)=(1.-urfConductivity)*kOld+
     *       urfConductivity*ViscosityrField(i,irField)*
     *          SpecificHeatrField(i,irField)/PrLaminarrField(irField)
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            kOld=BConductivityrField(i,j,irField)
            BConductivityrField(i,j,irField)=(1.-urfConductivity)*kOld+
     *         urfConductivity*BViscosityrField(i,j,irField)*
     *         BSpecificHeatrField(i,j,irField)/PrLaminarrField(irField)
c
          enddo 
        enddo 
c
      endif
c
      return
      end