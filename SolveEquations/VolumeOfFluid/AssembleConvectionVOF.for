c
c#############################################################################################
c
      SUBROUTINE AssembleConvectionTermrField(Variable,Bleed,
     *        ConvectionScheme,HRFramework,FiT,BFiT,dfidxT,
     *              dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c
c#############################################################################################
      use User0, only: nIterStartApplyingHR
      use MultiGrid2, only: nIter
      use Variables1, only: mdot,Bmdot,effdiv,
     *                      uVelocity,vVelocity,wVelocity
      use Variables2
      use Variables3
      use Geometry1
      use Geometry3
      use Geometry4
      use PhysicalProperties1
      use BoundaryConditions2
      use BoundaryConditionsrField2
      use VolumeOfFluid2, only: irFieldVariable
      use BoundaryConditionsTurbulence2
      use Constants1
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      integer :: i1,i2,i3,i4,i5
      character*10 Variable,VariableT
      double precision :: Bleed
      double precision :: Gamf,dfidxTf,dfidyTf,FluxTfLocal
      double precision :: cpf,FluxCflocal,FluxFflocal,FluxVflocal
      double precision :: FluxCElocal,FluxVElocal,FluxTElocal
      double precision :: xF1,yF1,zF1,distance1,distance2,GFactCF
      integer :: k1,j1,j2,j3,j4
      double precision :: PhiB,uF1,vF1,wF1
c
      character*20 ConvectionScheme
      character*4 HRFramework
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE HRSchemesCorrections(Variable,ConvectionScheme,
     *               Bleed,HRFramework,FiT,dfidxT,dfidyT,dfidzT)
c--------------------------------------------------------------------------
          character*20 ConvectionScheme
          character*10 Variable
          character*4 HRFramework
          double precision :: Bleed
          double precision :: phic,phid,phiu,phiTeldaC,diff,phiHR,
     *                        cpf,gf,rf,psirf,phiTeldaf,dcfx,dcfy
          integer i,j,k,iUpwind
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
        end SUBROUTINE HRSchemesCorrections
c--------------------------------------------------------------------------
      end interface
C********************************************************************************************
c
c--- Modify mdot field at interior faces
c
      do k=1,NIFaces
        mdot(k)=mdot(k)/Densityf(k)
      enddo
c
c--- Modify mdot field at boundary faces
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          Bmdot(i,j)=Bmdot(i,j)/BDensity(i,j)
c
        enddo
      enddo
c
c--- Interior faces (base on upwind. HR via deferred correction)
c
      do k=1,NIFaces
c        
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        FluxCflocal=dmax1(mdot(k),0.)
        FluxFflocal=-dmax1(-mdot(k),0.)
c
        FluxCf(k)=FluxCf(k)+FluxCflocal
        FluxFf(k)=FluxFf(k)+FluxFflocal
c
        FluxTf(k)=FluxTf(k)+FluxCflocal*FiT(i)+FluxFflocal*FiT(j)
c
      enddo
c
      variableT=Variable
      Variable='rfield'
c
      if(HRFramework.eq.'nvf'.or.HRFramework.eq.'tvd') then
c
        if(nIter.ge.nIterStartApplyingHR) then
c
          if(ConvectionScheme.eq.'stacs'.or.
     *        ConvectionScheme.eq.'hric'.or.
     *           ConvectionScheme.eq.'cicsam') 
     *                   call CalculateCosineThetaFace
          call HRSchemesCorrections(Variable,ConvectionScheme,
     *               Bleed,HRFramework,FiT,dfidxT,dfidyT,dfidzT)
c
        endif
c
      endif
c-----------------------------------------------------------------------
c
      call ComputeEffectiveDivergence(Variable)
c
      variable=VariableT
c
c---- Assemble convection-divergence term
c
      do i=1,NumberOfElements
c
        FluxCElocal =   dmax1(effdiv(i),0.0) - effdiv(i)
        FluxVElocal = - dmax1(effdiv(i),0.0) * FiT(i)
        FluxTElocal = FluxCElocal * FiT(i) + FluxVElocal
c
        FluxCE(i)=FluxCE(i)+FluxCElocal
        FluxTE(i)=FluxTE(i)+FluxTElocal
c
      enddo
c------------------------------------------------------------------------------
c--- Boundary faces
c------------------------------------------------------------------------------
      i5=irFieldVariable
c
      do i=1,IinletsupersonicrField(i5)
c
        i1=IinletsupersonicrFieldOwner(i,i5)
        i2=IinletsupersonicrFieldNumberOfBCSets(i,i5)
        i3=IinletsupersonicrFieldNBFaces(i,i5)
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
        FluxCflocal= dmax1(Bmdot(i2,i3),0.)
        FluxFflocal=-dmax1(-Bmdot(i2,i3),0.)
          
        FluxCf(i4)=FluxCf(i4)+FluxCflocal
        FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
        FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
      enddo
c-------------------------------------------------------------------
c
      do i=1,IinletSpecifiedValuerField(i5)
c
        i1=IinletSpecifiedValuerFieldOwner(i,i5)
        i2=IinletSpecifiedValuerFieldNumberOfBCSets(i,i5)
        i3=IinletSpecifiedValuerFieldNBFaces(i,i5)
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
        FluxCflocal= dmax1(Bmdot(i2,i3),0.)
        FluxFflocal=-dmax1(-Bmdot(i2,i3),0.)
          
        FluxCf(i4)=FluxCf(i4)+FluxCflocal
        FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
        FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
      enddo
c-------------------------------------------------------------------
      do i=1,IoutletsupersonicrField(i5)
c
        i1=IoutletsupersonicrFieldOwner(i,i5)
        i2=IoutletsupersonicrFieldNumberOfBCSets(i,i5)
        i3=IoutletsupersonicrFieldNBFaces(i,i5)
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
        FluxCflocal= dmax1(Bmdot(i2,i3),0.)
        FluxFflocal=-dmax1(-Bmdot(i2,i3),0.)
          
        FluxCf(i4)=FluxCf(i4)+FluxCflocal
        FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
        FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
      enddo
c
c-------------------------------------------------------------------
      do i=1,IpressureFarField
c
        i1=IpressureFarFieldOwner(i)
        i2=IpressureFarFieldNumberOfBCSets(i)
        i3=IpressureFarFieldNBFaces(i)
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
        FluxCflocal= dmax1(Bmdot(i2,i3),0.)
        FluxFflocal=-dmax1(-Bmdot(i2,i3),0.)
          
        FluxCf(i4)=FluxCf(i4)+FluxCflocal
        FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
        FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
      enddo
c
c-------------------------------------------------------------------
      do i=1,IoutletFullyDevelopedrField(i5)
c
        i1=IoutletFullyDevelopedrFieldOwner(i,i5)
        i2=IoutletFullyDevelopedrFieldNumberOfBCSets(i,i5)
        i3=IoutletFullyDevelopedrFieldNBFaces(i,i5)
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
        FluxCflocal= dmax1(Bmdot(i2,i3),0.)
        FluxFflocal=-dmax1(-Bmdot(i2,i3),0.)
          
        FluxCf(i4)=FluxCf(i4)+FluxCflocal
        FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
        FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
      enddo
c-------------------------------------------------------------------
      do i=1,Isymmetry
c
        i1=IsymmetryOwner(i)
        i2=IsymmetryNumberOfBCSets(i)
        i3=IsymmetryNBFaces(i)
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
        FluxCf(i4)=FluxCf(i4)+0.D0
        FluxFf(i4)=FluxFf(i4)+0.D0
        FluxVf(i4)=FluxVf(i4)+0.D0
        FluxTf(i4)=FluxTf(i4)+0.D0
c
      enddo
c-------------------------------------------------------------------
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
        FluxCflocal= dmax1(Bmdot(i2,i3),0.)
        FluxFflocal=-dmax1(-Bmdot(i2,i3),0.)
          
        FluxCf(i4)=FluxCf(i4)+FluxCflocal
        FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
        FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*FiT(j1)
c
      enddo
c-------------------------------------------------------------------
      do i=1,Iaxis
c
        i1=IaxisOwner(i)
        i2=IaxisNumberOfBCSets(i)
        i3=IaxisNBFaces(i)
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
        FluxCf(i4)=FluxCf(i4)+0.d0
        FluxFf(i4)=FluxFf(i4)+0.d0
        FluxVf(i4)=FluxVf(i4)+0.d0
        FluxTf(i4)=FluxTf(i4)+0.d0
c
      enddo
c-------------------------------------------------------------------
c
c--- Return mdot field at interior faces to its correct values
c
      do k=1,NIFaces
        mdot(k)=mdot(k)*Densityf(k)
      enddo
c
c--- Return mdot field at boundary faces to its correct values
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          Bmdot(i,j)=Bmdot(i,j)*BDensity(i,j)
c
        enddo
      enddo
c
      return
      end