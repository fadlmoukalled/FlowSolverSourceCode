c
C#############################################################################################
      SUBROUTINE InterpolateGamaToFace(Variable,Gam,BGam)
C#############################################################################################
      use User0, only: InterpolationSchemeGamaMomentum,
     *                 InterpolationSchemeGamaEnergy,
     *                 InterpolationSchemeGamaTKE,
     *                 InterpolationSchemeGamaTED,
     *                 InterpolationSchemeGamaTOmega,
     *                 InterpolationSchemeGamaTKL,
     *                 InterpolationSchemeGamaMED,
     *                 InterpolationSchemeGamaTGamma,
     *                 InterpolationSchemeGamaTReTheta,
     *                 InterpolationSchemeGamaTfRelaxation,
     *                 InterpolationSchemeGamaTurbulentV2,
     *                 InterpolationSchemeGamaTurbulentZeta
      use PhysicalProperties1, only: InterpolationSchemeGamaScalar, 
     *                               GamaFace,BGamaFace      
      use Scalar2, only: iScalarVariable
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces,NIFaces,NBFaceOwner
      use Geometry4, only: xc,yc,zc,BFaceCentroidx,
     *                     BFaceCentroidy,BFaceCentroidz
      use BoundaryConditions2, only: Iperiodic,IperiodicOwner,
     *                               IperiodicNumberOfBCSets,
     *                               IperiodicNBFaces,PeriodicPair,
     *                               Icorrespondingface,
     *                               a1r,b1r,c1r,a2r,b2r,c2r,
     *                               a3r,b3r,c3r,xTranslation,
     *                               yTranslation,zTranslation,
     *                               LRotationalPeriodicity,
     *                               LTranslationalPeriodicity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j,j1,j2,j3,j4,k1
      character*10 Variable
      character*16 InterpolationScheme
      double precision :: xF1,yF1,zF1,distance1,distance2,GFactCF
      double precision, dimension(:) :: Gam
      double precision, dimension(:,:) :: Bgam
c********************************************************************************************
c--- Interfaces
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE InterpolateElementToFace
     *                 (InterpolationScheme,Gam,GamF)
c--------------------------------------------------------------
          implicit none
c--------------------------------------------------------------
          character*16 InterpolationScheme
          double precision, dimension(:) :: Gam
          double precision, dimension(:) :: GamF
c--------------------------------------------------------------
        end SUBROUTINE InterpolateElementToFace
c********************************************************************************************
      end interface
c********************************************************************************************

      if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.
     *                             Variable.eq.'velz') then
c
        InterpolationScheme=InterpolationSchemeGamaMomentum
c
      elseif(Variable.eq.'temp') then
c
        InterpolationScheme=InterpolationSchemeGamaEnergy
c
      elseif(Variable.eq.'tke') then
c
        InterpolationScheme=InterpolationSchemeGamaTKE
c
      elseif(Variable.eq.'ted') then
c
        InterpolationScheme=InterpolationSchemeGamaTED
c
      elseif(Variable.eq.'tomega') then
c
        InterpolationScheme=InterpolationSchemeGamaTOmega
c
      elseif(Variable.eq.'tkl') then
c
        InterpolationScheme=InterpolationSchemeGamaTKL
c
      elseif(Variable.eq.'med') then
c
        InterpolationScheme=InterpolationSchemeGamaMED
c
      elseif(Variable.eq.'tgamma') then
c
        InterpolationScheme=InterpolationSchemeGamaTGamma
c
      elseif(Variable.eq.'tretheta') then
c
        InterpolationScheme=InterpolationSchemeGamaTReTheta
c
      elseif(Variable.eq.'tv2') then
c
        InterpolationScheme=InterpolationSchemeGamaTurbulentV2
c
      elseif(Variable.eq.'tzeta') then
c
        InterpolationScheme=InterpolationSchemeGamaTurbulentZeta
c
      elseif(Variable.eq.'frelax') then
c
        InterpolationScheme=InterpolationSchemeGamaTfRelaxation
c
      else
c
        InterpolationScheme=
     *       InterpolationSchemeGamaScalar(iScalarVariable)
c
      endif
c
      call InterpolateElementToFace(InterpolationScheme,Gam,GamaFace)
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BGamaFace(i,j)=BGam(i,j)
c
        enddo
      enddo
c
c----------------------------------------------------------------------
c
      do i=1,Iperiodic
c
        i1=IperiodicOwner(i)
        i2=IperiodicNumberOfBCSets(i)
        i3=IperiodicNBFaces(i)
        i4=NIFaces
c
        do j=1,i2-1
c
          i4=i4+NBFaces(j)
c
        enddo
c
        i4=i4+i3
c
        j2=PeriodicPair(i2)         
        j3=Icorrespondingface(i2,i3)
        j1=NBFaceOwner(j2,j3)
c
        do k1=1,j2-1
c
          j4=j4+NBFaces(k1)
c
        enddo
c
        j4=j4+j3
c
        if(LRotationalPeriodicity) then
c
          xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
          yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
          zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
        elseif(LTranslationalPeriodicity) then
c
          xF1=xc(j1)+xTranslation(j2)
          yF1=yc(j1)+yTranslation(j2)
          zF1=zc(j1)+zTranslation(j2)
c
        endif
c
        distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
        distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
        GFactCF=distance2/(distance1+distance2)
c
        BGamaFace(i2,i3)=GFactCF*Gam(i1)+(1.-GFactCF)*Gam(j1)
c
      enddo
c
      return
      end