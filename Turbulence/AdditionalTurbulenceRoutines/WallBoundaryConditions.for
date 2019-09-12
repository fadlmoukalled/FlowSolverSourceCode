c
c#############################################################################################
c
      SUBROUTINE SetTOmegaWallBC
c
C#############################################################################################
c
      use User0, only: LRough,WallBCModel,GrainSize
      use Turbulence1, only: cmu50,betta,cappa,uTau,KsPlus
      use Variables1, only: BTurbulentOmega
      use PhysicalProperties1, only: BDensity,BViscosity
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use Constants1, only:tiny,big,twothird
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: wWplus,ksmin,dNorm,WallVelocity,
     *                    uTaulocal,ks,SR,d0,KsLocal
c********************************************************************************************
      interface
c********************************************************************************************
        FUNCTION TangentialVelocity(i1)
c********************************************************************************************
          integer :: i1
          double precision :: TangentialVelocity
c********************************************************************************************
        end FUNCTION TangentialVelocity
c********************************************************************************************
      end interface
c********************************************************************************************
c
      if(.not.LRough) then
c
        if(WallBCModel.eq.'wilcox') then
c
          ksmin=big
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            WallVelocity=TangentialVelocity(i1)
            uTaulocal=dsqrt(BViscosity(i2,i3)*WallVelocity/
     *                            (dNorm*BDensity(i2,i3)))
            ks=3.*BViscosity(i2,i3)/(BDensity(i2,i3)*uTaulocal)
            ksmin=dmin1(ksmin,ks)
c
          enddo
c
          ksmin=ksmin*ksmin
          ksmin=dmax1(ksmin,tiny)
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            ksmin=dmin1(ksmin,ks)
            BTurbulentOmega(i2,i3)=40000.*BViscosity(i2,i3)/
     *                                     (BDensity(i2,i3)*ksmin)
c
          enddo
c
        elseif(WallBCModel.eq.'menter') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            dNorm=WallDistance(i1)
            BTurbulentOmega(i2,i3)=10.*6.*BViscosity(i2,i3)/
     *                          (BDensity(i2,i3)*betta*dNorm**2)
c
          enddo
c
        endif
c
      elseif(LRough) then
c
        if(WallBCModel.eq.'wilcox') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            if(KsPlus(i).le.25) then
              SR=(50./dmax1(KsPlus(i),tiny))**2
            elseif(KsPlus(i).gt.25) then
              SR=100./KsPlus(i)
            endif
c
c            if(KsPlus(i).le.5) then
c              SR=(200./dmax1(KsPlus(i),tiny))**2
c            elseif(KsPlus(i).gt.5) then
c              SR=100./KsPlus(i)+
c     *         (200./KsPlus(i)-100./KsPlus(i))*dexp(5.-KsPlus(i))
c            endif
c
            BTurbulentOmega(i2,i3)=uTau(i)*uTau(i)*
     *                   SR*BDensity(i2,i3)/BViscosity(i2,i3)
c
          enddo
c
        elseif(WallBCModel.eq.'knopp') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            d0=0.03*GrainSize(i2)*
     *             dmin1(1.,(KsPlus(i)/30.)**twothird)*
     *                  dmin1(1.,(KsPlus(i)/45.)**0.25)*
     *                      dmin1(1.,(KsPlus(i)/60.)**0.25)
            BTurbulentOmega(i2,i3)=uTau(i)/(cmu50*cappa*d0)
c
          enddo
c
        elseif(WallBCModel.eq.'nikuradse') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            KsLocal=dmax1(KsPlus(i),tiny)
            wWplus=400000./((KsLocal**4)*
     *           dtanh(10000./(3.*KsLocal**3)))+
     *              70.*(1.-dexp(-KsPlus(i)/300.))/KsLocal
     *      
            BTurbulentOmega(i2,i3)=
     *        uTau(i)*uTau(i)*wWplus*BDensity(i2,i3)/BViscosity(i2,i3)
c
          enddo
c
        elseif(WallBCModel.eq.'colebrook') then
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            KsLocal=dmax1(KsPlus(i),tiny)
            wWplus=300./((KsLocal**2)*dtanh(15./(4.*KsLocal)))+
     *                          191.*(1.-dexp(-KsPlus(i)/250.))/KsLocal
     *      
            BTurbulentOmega(i2,i3)=
     *        uTau(i)*uTau(i)*wWplus*BDensity(i2,i3)/BViscosity(i2,i3)
c
          enddo
c
        endif
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE SetTKEWallBC
c
C#############################################################################################
c
      use User0, only: LRough,WallBCModel
      use Turbulence1, only: cmu50,uTau,KsPlus
      use Variables1, only: BTurbulentKE 
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces

c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: k0Plus,kwPlus
c********************************************************************************************
c
      if(.not.LRough) then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          BTurbulentKE(i2,i3)=0.
c
        enddo
c
      elseif(LRough) then
c
        if(WallBCModel.eq.'wilcox') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          BTurbulentKE(i2,i3)=0.
c
        enddo
c
        elseif(WallBCModel.eq.'knopp') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          BTurbulentKE(i2,i3)=uTau(i)*uTau(i)*
     *                          dmin1(1.,KsPlus(i)/90.)/cmu50
c
        enddo
c
        elseif(WallBCModel.eq.'nikuradse') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          k0Plus=dtanh((dlog(KsPlus(i)/30.)/dlog(8.d0)+
     *      0.5*(1.-dtanh(KsPlus(i)/100.)))*dtanh(KsPlus(i)/75.))/cmu50
          kwPlus=dmax1(0.,k0Plus)
c
          BTurbulentKE(i2,i3)=uTau(i)*uTau(i)*kwPlus
c
        enddo
c
        elseif(WallBCModel.eq.'colebrook') then
c
        do i=1,IwallTurbulence
c
          i1=IWallTurbulenceOwner(i)
          i2=IWallTurbulenceNumberOfBCSets(i)
          i3=IWallTurbulenceNBFaces(i)
c
          k0Plus=dtanh((dlog(KsPlus(i)/30.)/dlog(10.d0)+
     *          1.-dtanh(KsPlus(i)/125.))*dtanh(KsPlus(i)/125.))/cmu50
          kwPlus=dmax1(0.,k0Plus)
c
          BTurbulentKE(i2,i3)=uTau(i)*uTau(i)*kwPlus
c
        enddo
c
        endif
c
      endif
c
      return
      end