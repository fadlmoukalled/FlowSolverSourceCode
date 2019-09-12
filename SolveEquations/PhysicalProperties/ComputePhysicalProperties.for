c
c#############################################################################################
c
      SUBROUTINE CalculatePhysicalProperties
c
C#############################################################################################
c
      Use User0, only: LFreeSurfaceFlow,Linviscid,Lcompressible,
     *                 NumberOfrFieldsToSolve,NumberOfScalarsToSolve,
     *                 LSurfaceTension,LsolverField
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
c********************************************************************************************
c
      if(LFreeSurfaceFlow) then
c
        do i=1,NumberOfrFieldsToSolve
c
          call CalculateSpecificHeatrField(i)
          call CalculateDensityrField(i)
c
        enddo    
c
        if(.not.Linviscid) then
c
          do i=1,NumberOfrFieldsToSolve
c
            call CalculateLaminarViscosityrField(i)
            call CalculateLaminarThermalConductivityrField(i)
c
          enddo    
c
        endif
c
        if(LSurfaceTension) then
c
          do i=1,NumberOfrFieldsToSolve
c
            if(LSolverField(i)) call CalculateSurfaceTensionTerm(i)
c
          enddo    
c
        endif
c
        call CalculateMaterialProperties
c
      else
c
        call CalculateSpecificHeat
        if(.not.Lcompressible) call Setdensity
c
        if(.not.Linviscid) then
c
          call CalculateLaminarViscosity      
          call CalculateLaminarThermalConductivity
c
        endif
c
      endif
c
      do i=1,NumberOfScalarsToSolve
c
        call CalculateSpecificHeatScalar(i)
c
        if(.not.Linviscid) then
c
          call CalculateDiffusionCoefficientScalar(i)
c
        endif
c
      enddo   
c      
      return
      end
