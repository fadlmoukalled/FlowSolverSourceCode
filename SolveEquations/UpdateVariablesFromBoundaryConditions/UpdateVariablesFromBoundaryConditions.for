c
C#############################################################################################
c
      SUBROUTINE updateValuesFromBoundaryConditions
     *               (Variable,Gam,BGam,FiT,BFiT,BdfidxT,
     *                  BdfidyT,BdfidzT,dfidxfT,dfidyfT,dfidzfT)
c
C#############################################################################################
      use User0
      use Variables1
      use Scalar1
      use Scalar2
      use BoundaryConditions2
      use BoundaryConditionsScalar2
      use BoundaryConditionsTurbulence2
      use BoundaryFluxes
      use Geometry4, only: gDiff,BgDiff,BFaceTx,BFaceTy,BFaceTz,
     *                     BFaceArea,xc,yc,zc,BFaceCentroidx,
     *                     BFaceCentroidy,BFaceCentroidz,
     *                     gDiffp,BgDiffp,BFaceTxp,BFaceTyp,BFaceTzp
      use Geometry3, only: NBFaceOwner
      use Constants1, only: tiny
      use PhysicalProperties1, only: BSpecificHeat,ReferenceTemperature
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i5,j1,j2,j3
      double precision :: term1,term2,gamf,dfidxTf,dfidyTf,dfidzTf,
     *                    xF1,yF1,zF1,distance1,distance2,GFactCF, Tref
      character*10 Variable
c
      double precision, dimension(:) :: FiT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
      double precision, dimension(:) :: dfidxfT
      double precision, dimension(:) :: dfidyfT
      double precision, dimension(:) :: dfidzfT
      double precision, dimension(:) :: Gam
      double precision, dimension(:,:) :: Bgam
c********************************************************************************************
c
      Tref=ReferenceTemperature
c
      if(Variable.eq.'velx') then
c



c
      elseif(Variable.eq.'vely') then
c




c
      elseif(Variable.eq.'velz') then
c




c
      elseif(Variable.eq.'pressc') then
c


c
      elseif(Variable.eq.'press') then
c


c
      elseif(Variable.eq.'tke') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'tv2') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'tzeta') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'frelax') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'med') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'tomega') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'ted') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'tkl') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'tgamma') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'tretheta') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'temp') then
c
        if(LanisotropicDiffusion) then
c
          do i=1,IWallVonNeumann
c
            i1=IWallVonNeumannOwner(i)
            i2=IWallVonNeumannNumberOfBCSets(i)
            i3=IWallVonNeumannNBFaces(i)
c
c--- Update the boundary value           
c
	      BFiT(i2,i3)=FiT(i1)-
     *                 HeatFlux(i2,i3)/(BgDiffp(i2,i3)+tiny)
c
          enddo
c
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
c--- Update the boundary value           
c
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+BgDiffp(i2,i3)
            term2=dfidxTf*BFaceTxp(i2,i3)+
     *              dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3)
c
	      BFiT(i2,i3)=
     *        (HinfinityRobin(i)*BFaceArea(i2,i3)*TinfinityRobin(i)+
     *            BgDiffp(i2,i3)*FiT(i1)-term2)/(term1+tiny)
c
          enddo
c
        else
c
          do i=1,IWallVonNeumann
c
            i1=IWallVonNeumannOwner(i)
            i2=IWallVonNeumannNumberOfBCSets(i)
            i3=IWallVonNeumannNBFaces(i)
c
c--- Update the boundary value           
c
            gamf=BGam(i2,i3)
	      BFiT(i2,i3)=FiT(i1)-
     *                 HeatFlux(i2,i3)/(gamf*BgDiff(i2,i3)+tiny)
c
          enddo
c
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
c--- Update the boundary value           
c
            gamf=BGam(i2,i3)
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+gamf*BgDiff(i2,i3)
            term2=gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
	      BFiT(i2,i3)=
     *        (HinfinityRobin(i)*BFaceArea(i2,i3)*TinfinityRobin(i)+
     *            gamf*BgDiff(i2,i3)*FiT(i1)-term2)/(term1+tiny)
c
          enddo
c
        endif
c
        do i=1,IinletSpecifiedStagnationTemperature
        enddo
c
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,IoutletFullyDevelopedEnergy
c
          i1=IoutletFullyDevelopedEnergyOwner(i)
          i2=IoutletFullyDevelopedEnergyNumberOfBCSets(i)
          i3=IoutletFullyDevelopedEnergyNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'htotal') then
c
        if(LanisotropicDiffusion) then
c
          do i=1,IWallDirichlet
c
            i2=IWallDirichletNumberOfBCSets(i)
            i3=IWallDirichletNBFaces(i)
c
            BFiT(i2,i3)=BSpecificHeat(i2,i3)*(BTemperature(i2,i3)-Tref)+
     *           0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                                        BwVelocity(i2,i3)**2)
c
          enddo
c
          do i=1,IWallVonNeumann
c
            i1=IWallVonNeumannOwner(i)
            i2=IWallVonNeumannNumberOfBCSets(i)
            i3=IWallVonNeumannNBFaces(i)
c
c--- Update the boundary value           
c
	      BTemperature(i2,i3)=Temperature(i1)-
     *                 HeatFlux(i2,i3)/(BgDiffp(i2,i3)+tiny)
            BFiT(i2,i3)=BSpecificHeat(i2,i3)*(BTemperature(i2,i3)-Tref)+
     *           0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                                        BwVelocity(i2,i3)**2)
c
          enddo
c
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
c--- Update the boundary value           
c
            dfidxTf=BTempGradx(i2,i3)
            dfidyTf=BTempGrady(i2,i3)
            dfidzTf=BTempGradz(i2,i3)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+BgDiffp(i2,i3)
            term2=dfidxTf*BFaceTxp(i2,i3)+
     *              dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3)
c
	      BTemperature(i2,i3)=
     *        (HinfinityRobin(i)*BFaceArea(i2,i3)*TinfinityRobin(i)+
     *            BgDiffp(i2,i3)*Temperature(i1)-term2)/(term1+tiny)
c
            BFiT(i2,i3)=BSpecificHeat(i2,i3)*(BTemperature(i2,i3)-Tref)+
     *           0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                                        BwVelocity(i2,i3)**2)
c
          enddo
c
        else
c
          do i=1,IWallDirichlet
c
            i1=IWallDirichletOwner(i)
            i2=IWallDirichletNumberOfBCSets(i)
            i3=IWallDirichletNBFaces(i)
c
            BFiT(i2,i3)=BSpecificHeat(i2,i3)*(BTemperature(i2,i3)-Tref)+
     *           0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                                        BwVelocity(i2,i3)**2)
c
          enddo
c
          do i=1,IWallVonNeumann
c
            i1=IWallVonNeumannOwner(i)
            i2=IWallVonNeumannNumberOfBCSets(i)
            i3=IWallVonNeumannNBFaces(i)
c
c--- Update the boundary value           
c
            gamf=BGam(i2,i3)
	      BTemperature(i2,i3)=Temperature(i1)-
     *                 HeatFlux(i2,i3)/(gamf*BgDiff(i2,i3)+tiny)
            BFiT(i2,i3)=BSpecificHeat(i2,i3)*(BTemperature(i2,i3)-Tref)+
     *           0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                                        BwVelocity(i2,i3)**2)
c
          enddo
c
          do i=1,IWallRobin
c
            i1=IWallRobinOwner(i)
            i2=IWallRobinNumberOfBCSets(i)
            i3=IWallRobinNBFaces(i)
c
c--- Update the boundary value           
c
            gamf=BGam(i2,i3)
            dfidxTf=BTempGradx(i2,i3)
            dfidyTf=BTempGrady(i2,i3)
            dfidzTf=BTempGradz(i2,i3)
c
            term1=HinfinityRobin(i)*BFaceArea(i2,i3)+gamf*BgDiff(i2,i3)
            term2=gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
	      BTemperature(i2,i3)=
     *        (HinfinityRobin(i)*BFaceArea(i2,i3)*TinfinityRobin(i)+
     *            gamf*BgDiff(i2,i3)*Temperature(i1)-term2)/(term1+tiny)
c
            BFiT(i2,i3)=BSpecificHeat(i2,i3)*(BTemperature(i2,i3)-Tref)+
     *           0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                                        BwVelocity(i2,i3)**2)
c
          enddo
c
        endif
c
        do i=1,IinletSpecifiedStagnationTemperature
        enddo
c
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
c
c--- Update the boundary value           
c
	    BTemperature(i2,i3)=Temperature(i1)
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,IoutletFullyDevelopedEnergy
c
          i1=IoutletFullyDevelopedEnergyOwner(i)
          i2=IoutletFullyDevelopedEnergyNumberOfBCSets(i)
          i3=IoutletFullyDevelopedEnergyNBFaces(i)
c
c--- Update the boundary value           
c
	    BTemperature(i2,i3)=Temperature(i1)
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
	    BTemperature(i2,i3)=Temperature(i1)
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BTemperature(i2,i3)=GFactCF*Temperature(i1)+
     *                                 (1.-GFactCF)*Temperature(j1)
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
	    BTemperature(i2,i3)=Temperature(i1)
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          BFiT(i2,i3)=BSpecificHeat(i2,i3)*(BTemperature(i2,i3)-Tref)+
     *           0.5*(BuVelocity(i2,i3)**2+BvVelocity(i2,i3)**2+
     *                                        BwVelocity(i2,i3)**2)
c
        enddo









c
      elseif(Variable.eq.'lambda') then
c
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c












c
      else
c
        i5=iScalarVariable
c
        if(LanisotropicDiffusion) then
c
          do i=1,IWallVonNeumannScalar(i5)
c
            i1=IWallVonNeumannScalarOwner(i,i5)
            i2=IWallVonNeumannScalarNumberOfBCSets(i,i5)
            i3=IWallVonNeumannScalarNBFaces(i,i5)
c
c--- Update the boundary value           
c
	      BFiT(i2,i3)=FiT(i1)-
     *             ScalarFlux(i2,i3,i5)/(BgDiffp(i2,i3)+tiny)
c
          enddo
c
          do i=1,IWallRobinScalar(i5)
c
            i1=IWallRobinScalarOwner(i,i5)
            i2=IWallRobinScalarNumberOfBCSets(i,i5)
            i3=IWallRobinScalarNBFaces(i,i5)
c
c--- Update the boundary value           
c
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            term1=ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)+
     *                                          BgDiffp(i2,i3)
            term2=(dfidxTf*BFaceTxp(i2,i3)+
     *              dfidyTf*BFaceTyp(i2,i3)+dfidzTf*BFaceTzp(i2,i3))
c
	      BFiT(i2,i3)=(ConvectionCoefficientRobin(i,i5)*
     *              BFaceArea(i2,i3)*PhiinfinityRobin(i,i5)+
     *                BgDiffp(i2,i3)*FiT(i1)-term2)/(term1+tiny)
c
          enddo
c
        else
c
          do i=1,IWallVonNeumannScalar(i5)
c
            i1=IWallVonNeumannScalarOwner(i,i5)
            i2=IWallVonNeumannScalarNumberOfBCSets(i,i5)
            i3=IWallVonNeumannScalarNBFaces(i,i5)
c
c--- Update the boundary value           
c
            gamf=BGam(i2,i3)
	      BFiT(i2,i3)=FiT(i1)-
     *             ScalarFlux(i2,i3,i5)/(gamf*BgDiff(i2,i3)+tiny)
c
          enddo
c
          do i=1,IWallRobinScalar(i5)
c
            i1=IWallRobinScalarOwner(i,i5)
            i2=IWallRobinScalarNumberOfBCSets(i,i5)
            i3=IWallRobinScalarNBFaces(i,i5)
c
c--- Update the boundary value           
c
            gamf=BGam(i2,i3)
            dfidxTf=BdfidxT(i2,i3)
            dfidyTf=BdfidyT(i2,i3)
            dfidzTf=BdfidzT(i2,i3)
c
            term1=ConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)+
     *                                          gamf*BgDiff(i2,i3)
            term2=gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
	      BFiT(i2,i3)=(ConvectionCoefficientRobin(i,i5)*
     *              BFaceArea(i2,i3)*PhiinfinityRobin(i,i5)+
     *                gamf*BgDiff(i2,i3)*FiT(i1)-term2)/(term1+tiny)
c
          enddo
c
        endif
c
        do i=1,IoutletsupersonicScalar(i5)
c
          i1=IoutletsupersonicScalarOwner(i,i5)
          i2=IoutletsupersonicScalarNumberOfBCSets(i,i5)
          i3=IoutletsupersonicScalarNBFaces(i,i5)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,IoutletFullyDevelopedScalar(i5)
c
          i1=IoutletFullyDevelopedScalarOwner(i,i5)
          i2=IoutletFullyDevelopedScalarNumberOfBCSets(i,i5)
          i3=IoutletFullyDevelopedScalarNBFaces(i,i5)
c
c--- Update the boundary value           
c
	    BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
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
          BFiT(i2,i3)=GFactCF*FiT(i1)+(1.-GFactCF)*FiT(j1)
c
        enddo
c
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      endif
c
      return
      end