c
c#############################################################################################
c
      SUBROUTINE CalculateSpecificHeatScalar(iScalar)
c
C#############################################################################################
c
      Use User0, only: ConstantSpecificHeatScalar
      Use PhysicalProperties1, only: SpecificHeatScalar,
     *                               BSpecificHeatScalar,
     *                               RGasScalar,GammaGasScalar,
     *                               EquationOfStateScalar
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use MultiGrid2, only: nIter
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,iScalar
c********************************************************************************************
c      
      if(nIter.gt.1) return
c
      if(EquationOfStateScalar(iScalar).eq.'constant') then
c      
        do i=1,NumberOfElements
c
          SpecificHeatScalar(i,iScalar)=
     *                   ConstantSpecificHeatScalar(iScalar)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeatScalar(i,j,iScalar)=
     *                   ConstantSpecificHeatScalar(iScalar)
c
          enddo 
        enddo 
c
      elseif(EquationOfStateScalar(iScalar).eq.'tait') then
c      
        do i=1,NumberOfElements
c
          SpecificHeatScalar(i,iScalar)=
     *                   ConstantSpecificHeatScalar(iScalar)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeatScalar(i,j,iScalar)=
     *                   ConstantSpecificHeatScalar(iScalar)
c
          enddo 
        enddo 
c      
      elseif(EquationOfStateScalar(iScalar).eq.'idealgas') then
c
        do i=1,NumberOfElements
c
          SpecificHeatScalar(i,iScalar)=GammaGasScalar(iScalar)*
     *                  RGasScalar(iScalar)/(GammaGasScalar(iScalar)-1.)
c
        enddo 
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BSpecificHeatScalar(i,j,iScalar)=GammaGasScalar(iScalar)*
     *                  RGasScalar(iScalar)/(GammaGasScalar(iScalar)-1.)
c
          enddo 
        enddo 
c
      endif
c
      return
      end