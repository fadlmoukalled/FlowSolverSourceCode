c
C#############################################################################################
c
      SUBROUTINE AssembleLinearPressureTerm(Variable)
c
C#############################################################################################
c
      use Variables3, only: FluxTE
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry4, only: Volume
      use BoundaryConditions1, only: BoundaryType
      use BoundaryConditions2, only: xTranslationUse,yTranslationUse,
     *                               zTranslationUse,periodicBeta
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i
      character*10 Variable
      double precision beta,LTranslation
c********************************************************************************************
c
      LTranslation=dsqrt(xTranslationUse**2+
     *                 yTranslationUse**2+zTranslationUse**2)
c
      if(Variable.eq.'velx') then
c
        beta=periodicBeta*xTranslationUse/LTranslation
c
      elseif(Variable.eq.'vely') then
c
        beta=periodicBeta*yTranslationUse/LTranslation
c
      elseif(Variable.eq.'velz') then
c
        beta=periodicBeta*zTranslationUse/LTranslation
c
      endif
c
c--- Assemble element fluxes
c      
      do i=1,NumberOfElements
c      
        FluxTE(i)=FluxTE(i)-beta*Volume(i)
c
      enddo
c
      return
      end
