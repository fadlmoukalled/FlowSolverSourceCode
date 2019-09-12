c
c#############################################################################################
c
      SUBROUTINE saveDensityForFalseTransient
c
C#############################################################################################
c
      use User0, only: LFalseTransientMomentum
      use PhysicalProperties1, only: Density,BDensity,
     *                               DensityStar,BDensityStar
c
c********************************************************************************************
c
      implicit none
c
c********************************************************************************************
c
      if(LFalseTransientMomentum) then
c
        DensityStar=Density
        BDensityStar=BDensity
c
      endif
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE AssembleFalseTransient(Variable,dtfalse)
c
C#############################################################################################
c
      use PhysicalProperties1, only: SpecificHeat,Density,
     *                               SpecificHeatScalar
      use Variables3
      use Variables2, only: acStar
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Scalar2
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i,j
      character*10 Variable
      double precision dtfalse
c********************************************************************************************
c
      if(Variable.eq.'velx'.or.Variable.eq.'vely'
     *                        .or.Variable.eq.'velz') then
c
        do i=1,NumberOfElements
c
          acStar(i)=Density(i)*Volume(i)/dtfalse
          FluxCE(i)=FluxCE(i)+acStar(i)
c
        enddo
c
      elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *         Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *          Variable.eq.'tkl'.or.Variable.eq.'tgamma'.or.
     *           Variable.eq.'tv2'.or.Variable.eq.'tzeta'.or.
     *            Variable.eq.'tretheta'.or.Variable.eq.'htotal'.or.
     *             Variable.eq.'frelax') then
c
        do i=1,NumberOfElements
c
          FluxCE(i)=FluxCE(i)+Density(i)*Volume(i)/dtfalse
c
        enddo
c
      elseif(Variable.eq.'lambda') then
c
        do i=1,NumberOfElements
c
          FluxCE(i)=FluxCE(i)+Volume(i)/dtfalse
c
        enddo
c
      elseif(Variable.eq.'temp') then
c
        do i=1,NumberOfElements
c
          FluxCE(i)=FluxCE(i)+
     *          SpecificHeat(i)*Density(i)*Volume(i)/dtfalse
c
        enddo
c
      else
c
        j=iScalarVariable
        do i=1,NumberOfElements
c
          FluxCE(i)=FluxCE(i)+SpecificHeatScalar(i,j)*
     *                                  Density(i)*Volume(i)/dtfalse
c
        enddo
c
      endif
c
      return
      end
c