c
c#############################################################################################
c
      SUBROUTINE AssembleConvectionTerm(Variable,Bleed,
     *        ConvectionScheme,NVF,TVD,FiT,BFiT,dfidxT,
     *              dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c
c#############################################################################################
      use User0, only: nIterStartApplyingHR
      use MultiGrid2, only: nIter
      use Variables1, only: mdot,Bmdot,effdiv,
     *                      uVelocity,vVelocity,wVelocity
      use Variables2
      use Variables3
      use Scalar2
      use Geometry1
      use Geometry3
      use Geometry4
      use PhysicalProperties1
      use BoundaryConditions2
      use BoundaryConditionsScalar2
      use BoundaryConditionsTurbulence2
      use Constants1
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      integer :: i1,i2,i3,i4,i5
      character*10 Variable
      double precision :: Bleed
      double precision :: Gamf,dfidxTf,dfidyTf,FluxTfLocal
      double precision :: cpf,FluxCflocal,FluxFflocal,FluxVflocal
      double precision :: FluxCElocal,FluxVElocal,FluxTElocal
      double precision :: xF1,yF1,zF1,distance1,distance2,GFactCF
      integer :: k1,j1,j2,j3,j4
      double precision :: PhiB,uF1,vF1,wF1
c
      character*20 ConvectionScheme
      logical :: NVF,TVD
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
     *               Bleed,NVF,TVD,FiT,dfidxT,dfidyT,dfidzT)
c--------------------------------------------------------------------------
          character*20 ConvectionScheme
          character*10 Variable
          logical NVF,TVD
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
c--- Interior faces (base on upwind. HR via deferred correction)
c
      if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.
     *                             Variable.eq.'velz') then
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
      elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *        Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *          Variable.eq.'tkl'.or.Variable.eq.'htotal'.or.
     *            Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *               Variable.eq.'tv2'.or.Variable.eq.'tzeta') then
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
      elseif(Variable.eq.'temp') then
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          cpf=GFactorCF(k)*SpecificHeat(i)+
     *               (1.-GFactorCF(k))*SpecificHeat(j)
c
          FluxCflocal=dmax1(mdot(k),0.)*cpf
          FluxFflocal=-dmax1(-mdot(k),0.)*cpf
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
c
          FluxTf(k)=FluxTf(k)+FluxCflocal*FiT(i)+FluxFflocal*FiT(j)
c
        enddo
c
      else
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          cpf=GFactorCF(k)*SpecificHeatScalar(i,iScalarVariable)+
     *          (1.-GFactorCF(k))*SpecificHeatScalar(j,iScalarVariable)
c
          FluxCflocal=dmax1(mdot(k),0.)*cpf
          FluxFflocal=-dmax1(-mdot(k),0.)*cpf
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
c
          FluxTf(k)=FluxTf(k)+FluxCflocal*FiT(i)+FluxFflocal*FiT(j)
c
        enddo
c
      endif
c
      if(NVF.or.TVD) then
c
        if(nIter.ge.nIterStartApplyingHR) then
c
          call HRSchemesCorrections(Variable,ConvectionScheme,
     *               Bleed,NVF,TVD,FiT,dfidxT,dfidyT,dfidzT)
c
        endif
c
      endif
c-----------------------------------------------------------------------
c
      call ComputeEffectiveDivergence(Variable)
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
c
c--- Boundary faces
c
      if(variable.eq.'velx'.or.variable.eq.'vely'
     *                          .or.variable.eq.'velz') then
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
c-------------------------------------------------------------------
        do i=1,IoutletTransmissive
c
          i1=IoutletTransmissiveOwner(i)
          i2=IoutletTransmissiveNumberOfBCSets(i)
          i3=IoutletTransmissiveNBFaces(i)
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
        do i=1,Iinletsupersonic
c
          i1=IinletsupersonicOwner(i)
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
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
        do i=1,IinletSpecifiedVelocity
c
          i1=IinletSpecifiedVelocityOwner(i)
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
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
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
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
        do i=1,IinletSpecifiedStaticPressure
c
          i1=IinletSpecifiedStaticPressureOwner(i)
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
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
        do i=1,IinletSpecifiedStagnationPressure
c
          i1=IinletSpecifiedStagnationPressureOwner(i)
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
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
c=====================================================================
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
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
          FluxCflocal= max(Bmdot(i2,i3),0.)
          FluxFflocal=-max(-Bmdot(i2,i3),0.)
          FluxVflocal=Bmdot(i2,i3)*(BdfidxT(i2,i3)*BDistanceCFx(i2,i3)+
     *                         BdfidyT(i2,i3)*BDistanceCFy(i2,i3)+
     *                               BdfidzT(i2,i3)*BDistanceCFz(i2,i3))
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c-------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i1=IoutletspecifiedVelocityOwner(i)
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
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
        do i=1,IoutletSpecifiedStaticPressure
c
          i1=IoutletSpecifiedStaticPressureOwner(i)
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
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
          FluxCflocal= max(Bmdot(i2,i3),0.)
          FluxFflocal=-max(-Bmdot(i2,i3),0.)
          FluxVflocal=Bmdot(i2,i3)*(BdfidxT(i2,i3)*BDistanceCFx(i2,i3)+
     *                     BdfidyT(i2,i3)*BDistanceCFy(i2,i3)+
     *                               BdfidzT(i2,i3)*BDistanceCFz(i2,i3))
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c-------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i1=IoutletSpecifiedAverageStaticPressureOwner(i)
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
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
          FluxCflocal= max(Bmdot(i2,i3),0.)
          FluxFflocal=-max(-Bmdot(i2,i3),0.)
          FluxVflocal=Bmdot(i2,i3)*(BdfidxT(i2,i3)*BDistanceCFx(i2,i3)+
     *                     BdfidyT(i2,i3)*BDistanceCFy(i2,i3)+
     *                               BdfidzT(i2,i3)*BDistanceCFz(i2,i3))
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c-------------------------------------------------------------------
        do i=1,IoutletSpecifiedResistance
c
          i1=IoutletSpecifiedResistanceOwner(i)
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
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
          FluxCflocal= max(Bmdot(i2,i3),0.)
          FluxFflocal=-max(-Bmdot(i2,i3),0.)
          FluxVflocal=Bmdot(i2,i3)*(BdfidxT(i2,i3)*BDistanceCFx(i2,i3)+
     *                     BdfidyT(i2,i3)*BDistanceCFy(i2,i3)+
     *                               BdfidzT(i2,i3)*BDistanceCFz(i2,i3))
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
c
        enddo
c-------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
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
        do i=1,IoutletFullyDeveloped
c
          i1=IoutletFullyDevelopedOwner(i)
          i2=IoutletFullyDevelopedNumberOfBCSets(i)
          i3=IoutletFullyDevelopedNBFaces(i)
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
          FluxCflocal= max(Bmdot(i2,i3),0.)
          FluxFflocal=-max(-Bmdot(i2,i3),0.)
          FluxVflocal=Bmdot(i2,i3)*(BdfidxT(i2,i3)*BDistanceCFx(i2,i3)+
     *                     BdfidyT(i2,i3)*BDistanceCFy(i2,i3)+
     *                               BdfidzT(i2,i3)*BDistanceCFz(i2,i3))
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)+FluxVflocal
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
c
        enddo
c-------------------------------------------------------------------
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
c          
          if(Variable.eq.'velx') then 
c
            PhiB=a1r(j2)*uVelocity(j1)+b1r(j2)*vVelocity(j1)+
     *                                       c1r(j2)*wVelocity(j1)
c
          elseif(Variable.eq.'vely') then
c
            PhiB=a2r(j2)*uVelocity(j1)+b2r(j2)*vVelocity(j1)+
     *                                       c2r(j2)*wVelocity(j1)
c
          elseif(Variable.eq.'velz') then
c
            PhiB=a3r(j2)*uVelocity(j1)+b3r(j2)*vVelocity(j1)+
     *                                       c3r(j2)*wVelocity(j1)
c
          endif
c
          FluxCflocal= dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          FluxTfLocal=FluxCflocal*FiT(i1)+FluxFflocal*PhiB
          FluxTf(i4)=FluxTf(i4)+FluxTfLocal
c
        enddo
c
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
      elseif(Variable.eq.'tgamma'.or.Variable.eq.'tretheta') then
c-------------------------------------------------------------------
c
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
c-------------------------------------------------------------------
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
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
        do i=1,IinletTurbulence
c
          i1=IinletTurbulenceOwner(i)
          i2=IinletTurbulenceNumberOfBCSets(i)
          i3=IinletTurbulenceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
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
c
c-------------------------------------------------------------------
      elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *          Variable.eq.'tomega'.or. Variable.eq.'med'.or. 
     *          Variable.eq.'tkl'.or.Variable.eq.'tv2'.or.
     *                            Variable.eq.'tzeta') then
c-------------------------------------------------------------------
c
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
c-------------------------------------------------------------------
        do i=1,IoutletTurbulence
c
          i1=IoutletTurbulenceOwner(i)
          i2=IoutletTurbulenceNumberOfBCSets(i)
          i3=IoutletTurbulenceNBFaces(i)
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
        do i=1,IinletTurbulence
c
          i1=IinletTurbulenceOwner(i)
          i2=IinletTurbulenceNumberOfBCSets(i)
          i3=IinletTurbulenceNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)+0.d0
          FluxFf(i4)=FluxFf(i4)+0.d0
          FluxVf(i4)=FluxVf(i4)+0.d0
          FluxTf(i4)=FluxTf(i4)+0.d0
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
c
c-------------------------------------------------------------------
      elseif(variable.eq.'temp') then
c-------------------------------------------------------------------
c
        do i=1,Iinletsupersonic
c
          i1=IinletsupersonicOwner(i)
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
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
          cpf=BSpecificHeat(i2,i3)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c-------------------------------------------------------------------
c
        do i=1,IinletSpecifiedStaticTemperature
c
          i1=IinletSpecifiedStaticTemperatureOwner(i)
          i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticTemperatureNBFaces(i)
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
          cpf=BSpecificHeat(i2,i3)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c-------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationTemperature
c
          i1=IinletSpecifiedStagnationTemperatureOwner(i)
          i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
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
          cpf=BSpecificHeat(i2,i3)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
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
          cpf=BSpecificHeat(i2,i3)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
c          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
c          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c
c-------------------------------------------------------------------
        do i=1,IoutletTransmissive
c
          i1=IoutletTransmissiveOwner(i)
          i2=IoutletTransmissiveNumberOfBCSets(i)
          i3=IoutletTransmissiveNBFaces(i)
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
          cpf=BSpecificHeat(i2,i3)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
c          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
c          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c
c-------------------------------------------------------------------        
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
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
          cpf=BSpecificHeat(i2,i3)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c-------------------------------------------------------------------
        do i=1,IoutletFullyDevelopedEnergy
c
          i1=IoutletFullyDevelopedEnergyOwner(i)
          i2=IoutletFullyDevelopedEnergyNumberOfBCSets(i)
          i3=IoutletFullyDevelopedEnergyNBFaces(i)
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
          cpf=BSpecificHeat(i2,i3)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
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
          cpf=GFactCF*SpecificHeat(i1)+(1.-GFactCF)*SpecificHeat(j1)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
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
c
c-------------------------------------------------------------------
      elseif(variable.eq.'htotal') then
c-------------------------------------------------------------------
c
        do i=1,Iinletsupersonic
c
          i1=IinletsupersonicOwner(i)
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
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
        do i=1,IinletSpecifiedStaticTemperature
c
          i1=IinletSpecifiedStaticTemperatureOwner(i)
          i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticTemperatureNBFaces(i)
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
        do i=1,IinletSpecifiedStagnationTemperature
c
          i1=IinletSpecifiedStagnationTemperatureOwner(i)
          i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
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
c          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
c          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c
c-------------------------------------------------------------------
        do i=1,IoutletTransmissive
c
          i1=IoutletTransmissiveOwner(i)
          i2=IoutletTransmissiveNumberOfBCSets(i)
          i3=IoutletTransmissiveNBFaces(i)
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
c          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
c          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c
c-------------------------------------------------------------------        
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
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
        do i=1,IoutletFullyDevelopedEnergy
c
          i1=IoutletFullyDevelopedEnergyOwner(i)
          i2=IoutletFullyDevelopedEnergyNumberOfBCSets(i)
          i3=IoutletFullyDevelopedEnergyNBFaces(i)
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
c
c-------------------------------------------------------------------
      else
c-------------------------------------------------------------------
c
        i5=iScalarVariable
c
        do i=1,IinletsupersonicScalar(i5)
c
          i1=IinletsupersonicScalarOwner(i,i5)
          i2=IinletsupersonicScalarNumberOfBCSets(i,i5)
          i3=IinletsupersonicScalarNBFaces(i,i5)
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
          cpf=BSpecificHeatScalar(i2,i3,iScalarVariable)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c-------------------------------------------------------------------
c
        do i=1,IinletSpecifiedValueScalar(i5)
c
          i1=IinletSpecifiedValueScalarOwner(i,i5)
          i2=IinletSpecifiedValueScalarNumberOfBCSets(i,i5)
          i3=IinletSpecifiedValueScalarNBFaces(i,i5)
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
          cpf=BSpecificHeatScalar(i2,i3,iScalarVariable)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c-------------------------------------------------------------------
        do i=1,IoutletsupersonicScalar(i5)
c
          i1=IoutletsupersonicScalarOwner(i,i5)
          i2=IoutletsupersonicScalarNumberOfBCSets(i,i5)
          i3=IoutletsupersonicScalarNBFaces(i,i5)
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
          cpf=BSpecificHeatScalar(i2,i3,iScalarVariable)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
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
          cpf=BSpecificHeatScalar(i2,i3,iScalarVariable)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
          FluxCf(i4)=FluxCf(i4)+FluxCflocal
          FluxFf(i4)=FluxFf(i4)+FluxFflocal
          
          FluxTf(i4)=FluxTf(i4)+
     *         FluxCflocal*FiT(i1)+FluxFflocal*BFiT(i2,i3)
c
        enddo
c
c-------------------------------------------------------------------
        do i=1,IoutletFullyDevelopedScalar(i5)
c
          i1=IoutletFullyDevelopedScalarOwner(i,i5)
          i2=IoutletFullyDevelopedScalarNumberOfBCSets(i,i5)
          i3=IoutletFullyDevelopedScalarNBFaces(i,i5)
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
          cpf=BSpecificHeatScalar(i2,i3,iScalarVariable)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
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
          cpf=GFactCF*SpecificHeatScalar(i1,iScalarVariable)+
     *              (1.-GFactCF)*SpecificHeatScalar(j1,iScalarVariable)
c
          FluxCflocal= cpf*dmax1(Bmdot(i2,i3),0.)
          FluxFflocal=-cpf*dmax1(-Bmdot(i2,i3),0.)
          
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
      endif
c
      return
      end