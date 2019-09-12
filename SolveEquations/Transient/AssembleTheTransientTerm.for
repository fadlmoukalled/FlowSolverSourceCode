c
C#############################################################################################
c
      SUBROUTINE AssembleTransientTerm(Variable,FiT,FiTold,FiToldold)
c
C#############################################################################################
c
      use User0, only: dt,TransientScheme
      use Transient1
      use PhysicalProperties1, only: SpecificHeat,Density,
     *                               SpecificHeatOld,DensityOld,
     *                               SpecificHeatOldOld,DensityOldOld,
     *                               SpecificHeatScalar,
     *                               SpecificHeatScalarOld,
     *                               SpecificHeatScalarOldOld
      use Variables3
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use MultiGrid2, only: nIter
      use scalar2
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i,j
      character*10 Variable
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: FiTold
      double precision, dimension(:) :: FiToldold
      double precision :: FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
      double precision :: FluxVElocal,FluxTElocal
      double precision :: dtCN
      double precision :: factor1,factor2,factor3
c
      double precision, save, dimension(:), allocatable :: drhodt
      double precision, save, dimension(:), allocatable :: rhocp
      double precision, save, dimension(:), allocatable :: rhocpOld
      double precision, save, dimension(:), allocatable :: rhocpOldOld
c
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE ComputeTransientGradient(phi,phiOld,phiOldOld,dphidt)
c--------------------------------------------------------------
          double precision, dimension(:) :: phi,phiOld,phiOldOld,dphidt
c--------------------------------------------------------------
        end SUBROUTINE ComputeTransientGradient
c--------------------------------------------------------------
      end interface
c
c********************************************************************************************
c
      if(TransientScheme.eq.'euler') then
c
        if(Variable.eq.'velx'.or.Variable.eq.'vely'
     *                         .or.Variable.eq.'velz') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*Volume(i)/dt
            FluxCEoldlocal=-DensityOld(i)*Volume(i)/dt
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *          Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *           Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *             Variable.eq.'tv2'.or.Variable.eq.'tzeta'.or.
     *                  Variable.eq.'tkl'.or.Variable.eq.'htotal') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*Volume(i)/dt
            FluxCEoldlocal=-DensityOld(i)*Volume(i)/dt
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        elseif(Variable.eq.'temp') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*SpecificHeat(i)*Volume(i)/dt
            FluxCEoldlocal=-DensityOld(i)*
     *                          SpecificHeatOld(i)*Volume(i)/dt
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        else
c
          j=iScalarVariable
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*SpecificHeatScalar(i,j)*Volume(i)/dt
            FluxCEoldlocal=-DensityOld(i)*
     *                          SpecificHeatScalarOld(i,j)*Volume(i)/dt
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        endif
c
      elseif(TransientScheme.eq.'cranknicolson1') then
c
        dtCN=dt/2.
c
        if(Variable.eq.'velx'.or.Variable.eq.'vely'
     *                            .or.Variable.eq.'velz') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*Volume(i)/dtCN
            FluxCEoldlocal=-DensityOld(i)*Volume(i)/dtCN
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *          Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *           Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *             Variable.eq.'tv2'.or.Variable.eq.'tzeta'.or.
     *                  Variable.eq.'tkl'.or.Variable.eq.'htotal') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*Volume(i)/dtCN
            FluxCEoldlocal=-DensityOld(i)*Volume(i)/dtCN
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        elseif(Variable.eq.'temp') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*SpecificHeat(i)*Volume(i)/dtCN
            FluxCEoldlocal=-DensityOld(i)*
     *                          SpecificHeatOld(i)*Volume(i)/dtCN
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        else
c
          j=iScalarVariable
          do i=1,NumberOfElements
c
            FluxCElocal=Density(i)*
     *                     SpecificHeatScalar(i,j)*Volume(i)/dtCN
            FluxCEoldlocal=-DensityOld(i)*
     *                  SpecificHeatScalarOld(i,j)*Volume(i)/dtCN
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
          enddo
c
        endif
c
      elseif(TransientScheme.eq.'cranknicolson2') then
c
        factor1=(dtOld/dt)/(dt+dtOld)
        factor2=1./(dt+dtOld)-(dtOldOld/dt)/(dtOld+dtOldOld)
        factor3=-(dtOld/dt)/(dtOld+dtOldOld)
c
        if(Variable.eq.'velx'.or.Variable.eq.'vely'
     *                           .or.Variable.eq.'velz') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*Volume(i)
            FluxCEoldlocal=factor2*DensityOld(i)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *             FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *          Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *           Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *             Variable.eq.'tv2'.or.Variable.eq.'tzeta'.or.
     *                  Variable.eq.'tkl'.or.Variable.eq.'htotal') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*Volume(i)
            FluxCEoldlocal=factor2*DensityOld(i)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *             FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        elseif(Variable.eq.'temp') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*SpecificHeat(i)*Volume(i)
            FluxCEoldlocal=factor2*DensityOld(i)*
     *                          SpecificHeatOld(i)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*
     *                          SpecificHeatOldOld(i)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *            FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        else
c
          j=iScalarVariable
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*
     *                   SpecificHeatScalar(i,j)*Volume(i)
            FluxCEoldlocal=factor2*DensityOld(i)*
     *                   SpecificHeatScalarOld(i,j)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*
     *                   SpecificHeatScalarOldOld(i,j)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *            FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        endif
c
      elseif(TransientScheme.eq.'adamsmoulton') then
c
        factor1=1./dt+1./(dt+dtOld)
        factor3=(dtOld/dt)/(dtOld+dtOldOld)
        factor2=factor1+factor3
c
        if(Variable.eq.'velx'.or.Variable.eq.'vely'
     *                             .or.Variable.eq.'velz') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*Volume(i)
            FluxCEoldlocal=-factor2*DensityOld(i)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *             FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *          Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *           Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *             Variable.eq.'tv2'.or.Variable.eq.'tzeta'.or.
     *                  Variable.eq.'tkl'.or.Variable.eq.'htotal') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*Volume(i)
            FluxCEoldlocal=-factor2*DensityOld(i)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *             FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        elseif(Variable.eq.'temp') then
c
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*SpecificHeat(i)*Volume(i)
            FluxCEoldlocal=-factor2*DensityOld(i)*
     *                          SpecificHeatOld(i)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*
     *                          SpecificHeatOldOld(i)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *            FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        else
c
          j=iScalarVariable
          do i=1,NumberOfElements
c
            FluxCElocal=factor1*Density(i)*
     *                     SpecificHeatScalar(i,j)*Volume(i)
            FluxCEoldlocal=-factor2*DensityOld(i)*
     *                     SpecificHeatScalarOld(i,j)*Volume(i)
            FluxCEoldoldlocal=factor3*DensityOldOld(i)*
     *                     SpecificHeatScalarOldOld(i,j)*Volume(i)
c
            FluxCE(i)=FluxCE(i)+FluxCElocal
            FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
            FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
            FluxTE(i)=FluxTE(i)+
     *            FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
          enddo
c
        endif
c
      endif
c
c--- Assemble the transient part of the effective diffusion term
c
      allocate(drhodt(NumberOfElements))
c
      if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.
     *    Variable.eq.'velz'.or.Variable.eq.'tke'.or.
     *     Variable.eq.'ted'.or.Variable.eq.'tomega'.or.
     *       Variable.eq.'med'.or.Variable.eq.'htotal'.or.
     *         Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *             Variable.eq.'tv2'.or.Variable.eq.'tzeta'.or.
     *                                     Variable.eq.'tkl') then
c
        call ComputeTransientGradient
     *            (Density,DensityOld,DensityOldOld,drhodt)
c
        do i=1,NumberOfElements
c
          FluxCElocal =   dmax1(drhodt(i),0.0) - drhodt(i)
          FluxVElocal = - dmax1(drhodt(i),0.0) * FiT(i)
          FluxTElocal = FluxCElocal * FiT(i) + FluxVElocal
c
          FluxCE(i)=FluxCE(i)+FluxCElocal*Volume(i)
          FluxTE(i)=FluxTE(i)+FluxTElocal*Volume(i)
c
        enddo
c
      elseif(Variable.eq.'temp') then
c
        allocate(rhocp(NumberOfElements))
        allocate(rhocpOld(NumberOfElements))
        allocate(rhocpOldOld(NumberOfElements))
c
        rhocp=Density*SpecificHeat
        rhocpOld=DensityOld*SpecificHeatOld
        rhocpOldOld=DensityOldOld*SpecificHeatOldOld
c
        call ComputeTransientGradient
     *            (rhocp,rhocpOld,rhocpOldOld,drhodt)
c
        do i=1,NumberOfElements
c
          FluxCElocal =   dmax1(drhodt(i),0.0) - drhodt(i)
          FluxVElocal = - dmax1(drhodt(i),0.0) * FiT(i)
          FluxTElocal = FluxCElocal * FiT(i) + FluxVElocal
c
          FluxCE(i)=FluxCE(i)+FluxCElocal*Volume(i)
          FluxTE(i)=FluxTE(i)+FluxTElocal*Volume(i)
c
        enddo
c
        deallocate(rhocp)
        deallocate(rhocpOld)
        deallocate(rhocpOldOld)
c
      else
c
        j=iScalarVariable
c
        allocate(rhocp(NumberOfElements))
        allocate(rhocpOld(NumberOfElements))
        allocate(rhocpOldOld(NumberOfElements))
c
        do i=1,NumberOfElements
c
          rhocp(i)=Density(i)*SpecificHeatScalar(i,j)
          rhocpOld(i)=DensityOld(i)*SpecificHeatScalarOld(i,j)
          rhocpOldOld(i)=DensityOldOld(i)*SpecificHeatScalarOldOld(i,j)
c
        enddo
c
        call ComputeTransientGradient
     *            (rhocp,rhocpOld,rhocpOldOld,drhodt)
c
        do i=1,NumberOfElements
c
          FluxCElocal =   dmax1(drhodt(i),0.0) - drhodt(i)
          FluxVElocal = - dmax1(drhodt(i),0.0) * FiT(i)
          FluxTElocal = FluxCElocal * FiT(i) + FluxVElocal
c
          FluxCE(i)=FluxCE(i)+FluxCElocal*Volume(i)
          FluxTE(i)=FluxTE(i)+FluxTElocal*Volume(i)
c
        enddo
c
        deallocate(rhocp)
        deallocate(rhocpOld)
        deallocate(rhocpOldOld)
c
      endif
c
      deallocate(drhodt)
c
      return
      end