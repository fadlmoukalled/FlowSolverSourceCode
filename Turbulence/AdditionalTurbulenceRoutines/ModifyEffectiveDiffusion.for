c
c#############################################################################################
c
      SUBROUTINE ModifyEffectiveDiffusionCoefficient(Variable)
c
c#############################################################################################
      use User0, only: TurbulenceModel,urfTViscosity,WallTreatment
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use PhysicalProperties1, only: TurbulentViscosity1,
     *                               BTurbulentViscosity1,
     *                               BTurbulentViscosity,
     *                               Density,BDensity,
     *                               Viscosity,BViscosity,
     *                               eDiffCoefficient,BeDiffCoefficient
      use Variables1, only: TurbulentKE,BTurbulentKE,
     *                      TurbulentED,BTurbulentED,
     *                      TurbulentOmega,BTurbulentOmega,
     *                      uVelGradx,vVelGradx,wVelGradx,
     *                      BuVelGradx,BvVelGradx,BwVelGradx,
     *                      uVelGrady,vVelGrady,wVelGrady,
     *                      BuVelGrady,BvVelGrady,BwVelGrady,
     *                      uVelGradz,vVelGradz,wVelGradz,
     *                      BuVelGradz,BvVelGradz,BwVelGradz
      use Turbulence1, only: sigTKE,sigTED,alfastar,Balfastar
      use BoundaryConditionsTurbulence1, only: inletTypeT
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                   IWallTurbulenceNumberOfBCSets,
     *                                         IWallTurbulenceNBFaces
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,i2,i3,j,k
      double precision :: TVisOld
c********************************************************************************************
c
      if(Variable.eq.'tke') then
c
        if(TurbulenceModel.eq.'komega2006') then
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity1(i)
            TurbulentViscosity1(i)=(1.-urfTViscosity)*TVisOld+
     *         urfTViscosity*Density(i)*dmax1(TurbulentKE(i),0.)/
     *                                 dmax1(TurbulentOmega(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity1(i,j)
              BTurbulentViscosity1(i,j)=(1.-urfTViscosity)*TVisOld+
     *        urfTViscosity*BDensity(i,j)*dmax1(BTurbulentKE(i,j),0.)/
     *                                  dmax1(BTurbulentOmega(i,j),tiny)
c
            enddo 
          enddo 
c
          do i=1,NumberOfBCSets
c
            if(inletTypeT(i).eq.'specifiedkew') then
c
              do j=1,NBFaces(i)
c
                k=NBFaceOwner(i,j)
                BTurbulentViscosity1(i,j)=TurbulentViscosity1(k)
c
              enddo
c
            endif
c
          enddo
c
          do i=1,IwallTurbulence
c
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            BTurbulentViscosity1(i2,i3)=BTurbulentViscosity(i2,i3)
c
          enddo
c
        elseif(TurbulenceModel.eq.'komega2006lrn') then
c
          do i=1,NumberOfElements
c
            TVisOld=TurbulentViscosity1(i)
            TurbulentViscosity1(i)=(1.-urfTViscosity)*TVisOld+
     *         urfTViscosity*alfaStar(i)*Density(i)*
     *           dmax1(TurbulentKE(i),0.)/dmax1(TurbulentOmega(i),tiny)
c
          enddo 
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              TVisOld=BTurbulentViscosity1(i,j)
              BTurbulentViscosity1(i,j)=(1.-urfTViscosity)*TVisOld+
     *           urfTViscosity*BalfaStar(i,j)*BDensity(i,j)*
     *                      dmax1(BTurbulentKE(i,j),0.)/
     *                          dmax1(BTurbulentOmega(i,j),tiny)
c
            enddo 
          enddo 
c
          do i=1,NumberOfBCSets
c
            if(inletTypeT(i).eq.'specifiedkew') then
c
              do j=1,NBFaces(i)
c
                k=NBFaceOwner(i,j)
                BTurbulentViscosity1(i,j)=TurbulentViscosity1(k)
c
              enddo
c
            endif
c
          enddo
c
          do i=1,IwallTurbulence
c
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            BTurbulentViscosity1(i2,i3)=BTurbulentViscosity(i2,i3)
c
          enddo
c
        endif
c
      endif
c
      if(Variable.eq.'tke') then
c
        do i=1,NumberOfElements
c
          eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity1(i)/sigTKE 
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                            BTurbulentViscosity1(i,j)/sigTKE
c
          enddo
        enddo
c
      elseif(Variable.eq.'tomega') then
c
        do i=1,NumberOfElements
c
          eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity1(i)/sigTED
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                             BTurbulentViscosity1(i,j)/sigTED
c
          enddo
        enddo
c
      endif
c
      return
      end