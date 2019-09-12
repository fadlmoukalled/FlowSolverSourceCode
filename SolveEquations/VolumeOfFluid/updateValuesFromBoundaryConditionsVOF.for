c
C#############################################################################################
c
      SUBROUTINE updateValuesFromBoundaryConditionsrField
     *               (Variable,Gam,BGam,FiT,BFiT,BdfidxT,
     *                  BdfidyT,BdfidzT,dfidxfT,dfidyfT,dfidzfT)
c
C#############################################################################################
      use User0
      use Variables1
      use VolumeOfFluid1
      use VolumeOfFluid2, only: irFieldVariable
      use BoundaryConditions2
      use BoundaryConditionsrField2
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
      i5=irFieldVariable
c
      do i=1,IWallVonNeumannrField(i5)
c
        i1=IWallVonNeumannrFieldOwner(i,i5)
        i2=IWallVonNeumannrFieldNumberOfBCSets(i,i5)
        i3=IWallVonNeumannrFieldNBFaces(i,i5)
c
c--- Update the boundary value           
c
c        gamf=BGam(i2,i3)
	  BFiT(i2,i3)=FiT(i1) !-
c     *             rFieldFlux(i2,i3,i5)/(gamf*BgDiff(i2,i3)+tiny)
c
      enddo
c
!      do i=1,IWallRobinrField(i5)
!c
!        i1=IWallRobinrFieldOwner(i,i5)
!        i2=IWallRobinrFieldNumberOfBCSets(i,i5)
!        i3=IWallRobinrFieldNBFaces(i,i5)
!c
!c--- Update the boundary value           
!c
!        gamf=BGam(i2,i3)
!        dfidxTf=BdfidxT(i2,i3)
!        dfidyTf=BdfidyT(i2,i3)
!        dfidzTf=BdfidzT(i2,i3)
!c
!        term1=rFieldConvectionCoefficientRobin(i,i5)*BFaceArea(i2,i3)+
!     *                                          gamf*BgDiff(i2,i3)
!        term2=gamf*(dfidxTf*BFaceTx(i2,i3)+
!     *                   dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
!c
!	  BFiT(i2,i3)=(rFieldConvectionCoefficientRobin(i,i5)*
!     *              BFaceArea(i2,i3)*PhiinfinityRobin(i,i5)+
!     *                gamf*BgDiff(i2,i3)*FiT(i1)-term2)/(term1+tiny)
!c
!      enddo
c
      do i=1,IoutletsupersonicrField(i5)
c
        i1=IoutletsupersonicrFieldOwner(i,i5)
        i2=IoutletsupersonicrFieldNumberOfBCSets(i,i5)
        i3=IoutletsupersonicrFieldNBFaces(i,i5)
c
c--- Update the boundary value           
c
	  BFiT(i2,i3)=FiT(i1)
c
      enddo
c
      do i=1,IoutletFullyDevelopedrField(i5)
c
        i1=IoutletFullyDevelopedrFieldOwner(i,i5)
        i2=IoutletFullyDevelopedrFieldNumberOfBCSets(i,i5)
        i3=IoutletFullyDevelopedrFieldNBFaces(i,i5)
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
      return
      end