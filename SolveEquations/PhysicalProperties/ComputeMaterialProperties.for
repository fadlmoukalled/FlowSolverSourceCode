c
c#############################################################################################
c
      SUBROUTINE CalculateMaterialProperties
c
C#############################################################################################
c
      Use User0, only: NumberOfrFieldsToSolve,LSolveMomentum,
     *                 LSolveEnergy,Linviscid
      Use PhysicalProperties1, only: Density,BDensity,
     *                               Viscosity,BViscosity,
     *                               Conductivity,BConductivity,
     *                               DensityrField,BDensityrField,
     *                               ViscosityrField,BViscosityrField,
     *                               ConductivityrField,SpecificHeat,
     *                               BConductivityrField,BSpecificHeat,
     *                               SpecificHeatrField,
     *                               BSpecificHeatrField
      use VolumeOfFluid1, only: rField,BrField
      use Variables1, only: Pressure,BPressure,Temperature,BTemperature
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
c********************************************************************************************
c
      Density=0.
      Viscosity=0.
      Conductivity=0.
      SpecificHeat=0.
      BDensity=0.
      BViscosity=0.
      BConductivity=0.
      BSpecificHeat=0.
c
      do i=1,NumberOfElements
        do j=1,NumberOfrFieldsToSolve
c
          Density(i)=Density(i)+rField(i,j)*DensityrField(i,j)
c
        enddo
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
          do k=1,NumberOfrFieldsToSolve
c
            BDensity(i,j)=BDensity(i,j)+
     *                        BrField(i,j,k)*BDensityrField(i,j,k)
c
          enddo
        enddo 
      enddo 
c      
      if(LSolveMomentum.and..not.Linviscid) then
c
      do i=1,NumberOfElements
        do j=1,NumberOfrFieldsToSolve
c
          Viscosity(i)=Viscosity(i)+rField(i,j)*ViscosityrField(i,j)
c
        enddo
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
          do k=1,NumberOfrFieldsToSolve
c
            BViscosity(i,j)=BViscosity(i,j)+
     *                        BrField(i,j,k)*BViscosityrField(i,j,k)
c
          enddo
        enddo 
      enddo 
c      
      endif
c
      if(LSolveEnergy) then
c
      do i=1,NumberOfElements
        do j=1,NumberOfrFieldsToSolve
c
          Conductivity(i)=Conductivity(i)+
     *                             rField(i,j)*ConductivityrField(i,j)
          SpecificHeat(i)=SpecificHeat(i)+
     *                             rField(i,j)*SpecificHeatrField(i,j)
c
        enddo
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
          do k=1,NumberOfrFieldsToSolve
c
            BConductivity(i,j)=BConductivity(i,j)+
     *                        BrField(i,j,k)*BConductivityrField(i,j,k)
            BSpecificHeat(i,j)=BSpecificHeat(i,j)+
     *                        BrField(i,j,k)*BSpecificHeatrField(i,j,k)
c
          enddo
        enddo 
      enddo 
c      
      endif
c
      call CalculateFaceDensity       
c
      return
      end
