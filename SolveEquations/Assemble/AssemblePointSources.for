c
C#############################################################################################
c
      SUBROUTINE AssemblePointSourceTerm(FiT,ScPointSource,
     *                                             SbPointSource)
c
C#############################################################################################
c
      use BoundaryConditions1, only: NumberOfPointSources,
     *                               iElementPointSource
      use Variables3
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i,j
      double precision, dimension(:)  :: FiT
      double precision, dimension(:)  :: ScPointSource
      double precision, dimension (:) :: SbPointSource
c********************************************************************************************
c
c--- Assemble element fluxes
c      
      do j=1,NumberOfPointSources
c      
        i=iElementPointSource(j)
        if(ScPointSource(j).lt.0.) then
c
          FluxCE(i)=FluxCE(i)-ScPointSource(j)
          FluxTE(i)=FluxTE(i)-(SbPointSource(j)+ScPointSource(j)*FiT(i))
c
        else
c
          FluxCE(i)=FluxCE(i)+0.
          FluxTE(i)=FluxTE(i)-(SbPointSource(j)+ScPointSource(j)*FiT(i))
c
        endif
c
      enddo
c
      return
      end
