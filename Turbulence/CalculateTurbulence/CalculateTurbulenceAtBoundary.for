c
c#############################################################################################
c
      SUBROUTINE CalculateInletTurbulence
c
C#############################################################################################
      use User0, only: urfTViscosity,LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate
      use Turbulence1, only: cmu,cmu75
      use BoundaryConditionsTurbulence1
      use BoundaryConditionsTurbulence2
      use Geometry1, only: NumberOfBCSets     
      use Geometry3, only: NBFaces
      use Variables1, only: BTurbulentKE,BTurbulentED,BTurbulentOmega,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: BTurbulentViscosity,BViscosity,
     *                               BDensity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: Bvelocity2,intensity2,length,TVisOld,ratio   !ratio=mut/mu
c********************************************************************************************
c
      do i=1,NumberOfBCSets
c
        if(inletTypeT(i).eq.'specifiedil') then
c
          do j=1,NBFaces(i)
c
            Bvelocity2=BuVelocity(i,j)**2+
     *                    BvVelocity(i,j)**2+BwVelocity(i,j)**2
            intensity2=inletTurbulenceIntensity(i)**2
            BTurbulentKE(i,j)=1.5*Bvelocity2*intensity2
            length=inletTurbulenceLengthScale(i)
            if(LSolveTurbulenceDissipationRate) then
              BTurbulentED(i,j)=
     *           cmu75*(dmax1(BTurbulentKE(i,j),0.)**1.5)/length
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *        cmu*urfTViscosity*BDensity(i,j)*BTurbulentKE(i,j)**2/
     *                                 dmax1(BTurbulentED(i,j),tiny)
            elseif(LSolveTurbulenceSpecificDissipationRate) then
              BTurbulentOmega(i,j)=
     *           cmu75*(dmax1(BTurbulentKE(i,j),0.)**1.5)/length
              BTurbulentOmega(i,j)=BTurbulentOmega(i,j)/
     *                        (cmu*dmax1(BTurbulentKE(i,j),tiny))
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*BDensity(i,j)*BTurbulentKE(i,j)/
     *                                 dmax1(BTurbulentOmega(i,j),tiny)
            endif
c
          enddo
c
        elseif(inletTypeT(i).eq.'specifiedkr') then
c
          do j=1,NBFaces(i)
c
            ratio=inletTurbulenceViscosityRatio(i)
            if(LSolveTurbulenceDissipationRate) then
              BTurbulentED(i,j)=cmu*BDensity(i,j)*
     *           (BTurbulentKE(i,j)**2)/(ratio*BViscosity(i,j))
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *           urfTViscosity*BViscosity(i,j)*ratio
            elseif(LSolveTurbulenceSpecificDissipationRate) then
              BTurbulentED(i,j)=BDensity(i,j)*
     *                  BTurbulentKE(i,j)/(BViscosity(i,j)*ratio)
              TVisOld=BTurbulentViscosity(i,j)
              BTurbulentViscosity(i,j)=(1.-urfTViscosity)*TVisOld+
     *           urfTViscosity*BViscosity(i,j)*ratio
            endif
c
          enddo
c
        endif
c
      enddo
c      
      return
      end
c
C#############################################################################################
c
      SUBROUTINE UpdateTurbulentZetaValueAtInlet(BFiT)
c
C#############################################################################################
      use BoundaryConditionsTurbulence2, only: IinletTurbulence,
     *         IinletTurbulenceNumberOfBCSets,IinletTurbulenceNBFaces
      use Constants1, only: twothird
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i2,i3
c********************************************************************************************
      double precision, dimension(:,:) :: BFiT
c********************************************************************************************
c
      do i=1,IinletTurbulence
c
        i2=IinletTurbulenceNumberOfBCSets(i)
        i3=IinletTurbulenceNBFaces(i)
c
        BFiT(i2,i3)=twothird
c
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE UpdateTurbulentV2ValueAtInlet(BFiT)
c
C#############################################################################################
      use Variables1, only: BTurbulentKE
      use BoundaryConditionsTurbulence2, only: IinletTurbulence,
     *         IinletTurbulenceNumberOfBCSets,IinletTurbulenceNBFaces
      use Constants1, only: twothird
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i2,i3
c********************************************************************************************
      double precision, dimension(:,:) :: BFiT
c********************************************************************************************
c
      do i=1,IinletTurbulence
c
        i2=IinletTurbulenceNumberOfBCSets(i)
        i3=IinletTurbulenceNBFaces(i)
c
        BFiT(i2,i3)=twothird*BTurbulentKE(i2,i3)
c
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE UpdateTfRelaxationValueAtInlet(BFiT)
c
C#############################################################################################
      use Variables1, only: TfRelaxation
      use BoundaryConditionsTurbulence2, only: IinletTurbulence,
     *         IinletTurbulenceOwner,IinletTurbulenceNumberOfBCSets,
     *         IinletTurbulenceNBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
c********************************************************************************************
      double precision, dimension(:,:) :: BFiT
c********************************************************************************************
c
      do i=1,IinletTurbulence
c
        i1=IinletTurbulenceOwner(i)
        i2=IinletTurbulenceNumberOfBCSets(i)
        i3=IinletTurbulenceNBFaces(i)
c
        BFiT(i2,i3)=TfRelaxation(i1)
c
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE UpdateTGammaValueAtInlet(BFiT)
c
C#############################################################################################
      use BoundaryConditionsTurbulence2, only: IinletTurbulence,
     *         IinletTurbulenceNumberOfBCSets,IinletTurbulenceNBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i2,i3
c********************************************************************************************
      double precision, dimension(:,:) :: BFiT
c********************************************************************************************
c
      do i=1,IinletTurbulence
c
        i2=IinletTurbulenceNumberOfBCSets(i)
        i3=IinletTurbulenceNBFaces(i)
c
        BFiT(i2,i3)=1.
c
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE UpdateTReThetaValueAtInlet(BFiT)
c
C#############################################################################################
      use BoundaryConditionsTurbulence2, only: IinletTurbulence,
     *         IinletTurbulenceNumberOfBCSets,IinletTurbulenceNBFaces
      use Variables1, only: BuVelocity,BvVelocity,
     *                      BwVelocity,BTurbulentKE
      use Constants1, only: twothird,tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i2,i3
      double precision :: Uvel,TU
c********************************************************************************************
      double precision, dimension(:,:) :: BFiT
c********************************************************************************************
c
      do i=1,IinletTurbulence
c
        i2=IinletTurbulenceNumberOfBCSets(i)
        i3=IinletTurbulenceNBFaces(i)
c
        Uvel=dsqrt(BuVelocity(i2,i3)**2+
     *                BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2)
        TU=100.*dsqrt(twothird*BTurbulentKE(i2,i3))/dmax1(Uvel,tiny)
        TU=dmax1(TU,0.027)
        if(TU.le.1.3) then
          BFiT(i2,i3)=1173.51-589.428*TU+0.2196/(TU**2)
        elseif(TU.gt.1.3) then
          BFiT(i2,i3)=331.5*(TU-0.5658)**(-0.671)
        endif
c
        enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE updateFarFieldTurbulenceValues(Variable,FiT,BFiT)
c
C#############################################################################################
c
      use Variables1, only: Bmdot,BTurbulentKE
      use BoundaryConditions2, only: TKEFarField,TEDFarField,
     *                               TOmegaFarField,MEDFarField,
     *        IpressureFarField,IpressureFarFieldOwner,
     *        IpressureFarFieldNumberOfBCSets,IpressureFarFieldNBFaces,
     *        TurbulentKLFarField,MachFarField,TemperatureFarField
      use PhysicalProperties1, only: RGas,GammaGas
      use Constants1, only: twothird
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      character*10 Variable
      double precision :: cInfinity,vInfinity,TU
c
      double precision, dimension(:) :: FiT
      double precision, dimension(:,:) :: BFiT
c********************************************************************************************
c
      if(Variable.eq.'tke') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=TKEFarField(i2)
c
          endif
c
        enddo
c
      elseif(Variable.eq.'med') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=MEDFarField(i2)
c
          endif
c
        enddo
c
      elseif(Variable.eq.'tomega') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=TOmegaFarField(i2)
c
          endif
c
        enddo
c
      elseif(Variable.eq.'ted') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=TEDFarField(i2)
c
          endif
c
        enddo
c
      elseif(Variable.eq.'frelax') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          BFiT(i2,i3)=FiT(i1)
c
        enddo
c
      elseif(Variable.eq.'tv2') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=twothird*BTurbulentKE(i2,i3)
c
          endif
c
        enddo
c
      elseif(Variable.eq.'tzeta') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=twothird
c
          endif
c
        enddo
c
      elseif(Variable.eq.'tkl') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=TurbulentKLFarField(i2)
c
          endif
c
        enddo
c
      elseif(Variable.eq.'tgamma') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            BFiT(i2,i3)=1.
c
          endif
c
        enddo
c
      elseif(Variable.eq.'tretheta') then
c
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          if(Bmdot(i2,i3).gt.0.) then
c
            BFiT(i2,i3)=FiT(i1)
c
          else
c
            cInfinity=dsqrt(GammaGas*RGas*TemperatureFarField(i2))
            vInfinity=MachFarField(i2)*cInfinity
            TU=100.*dsqrt(twothird*TKEFarField(i2))/vInfinity
            TU=dmax1(TU,0.027)
            if(TU.le.1.3) then
              BFiT(i2,i3)=1173.51-589.428*TU+0.2196/(TU**2)
            elseif(TU.gt.1.3) then
              BFiT(i2,i3)=331.5*(TU-0.5658)**(-0.671)
            endif
c
          endif
c
        enddo
c
      endif
c
      return
      end
c