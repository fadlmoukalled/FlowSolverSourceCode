c
c#############################################################################################
c
      SUBROUTINE CalculateOneEquationModelTurbulentKE
c
c#############################################################################################
c
      use User0, only: Lcompressible
      use Variables1, only: TurbulentKE,BTurbulentKE,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz
      use Turbulence1, only: cmu,S11,S12,S13,S22,S23,S33,
     *                       BS11,BS12,BS13,BS22,BS23,BS33,
     *                       StrainRate,BStrainRate
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFacesMax
      use PhysicalProperties1, only: Density,BDensity,
     *                      TurbulentViscosity,BTurbulentViscosity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: term1,Bterm1,factor
      double precision :: S11T,S12T,S13T,S22T,S23T,S33T
      double precision :: BS11T,BS12T,BS13T,BS22T,BS23T,BS33T
      double precision, save, dimension(:), allocatable :: Sstar
      double precision, save, dimension(:,:), allocatable :: BSstar
c********************************************************************************************
c
      if(LCompressible) Then
c
        allocate(Sstar(NumberOfElements))
        allocate(BSstar(NumberOfBCSets,NBFacesMax))
c 
        factor=dsqrt(cmu)
c
        do i=1,NumberOfElements    
c
          term1=(uVelGradx(i)+vVelGrady(i)+wVelGradz(i))/3.
c
          S11(i)=uVelGradx(i)
          S12(i)=0.5*(uVelGrady(i)+vVelGradx(i))
          S13(i)=0.5*(uVelGradz(i)+wVelGradx(i))
          S22(i)=vVelGrady(i)
          S23(i)=0.5*(vVelGradz(i)+wVelGrady(i))
          S33(i)=wVelGradz(i)
c
          S11T=S11(i)-term1
          S12T=S12(i)
          S13T=S13(i)
          S22T=S22(i)-term1
          S23T=S23(i)
          S33T=S33(i)-term1
c          
          Sstar(i)=dsqrt(2.*
     *         (S11T*S11T+S12T*S12T+S13T*S13T+
     *          S12T*S12T+S22T*S22T+S23T*S23T+
     *          S13T*S13T+S23T*S23T+S33T*S33T))
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            Bterm1=(BuVelGradx(i,j)+BvVelGrady(i,j)+BwVelGradz(i,j))/3.
c
            BS11(i,j)=BuVelGradx(i,j)
            BS12(i,j)=0.5*(BuVelGrady(i,j)+BvVelGradx(i,j))
            BS13(i,j)=0.5*(BuVelGradz(i,j)+BwVelGradx(i,j))
            BS22(i,j)=BvVelGrady(i,j)
            BS23(i,j)=0.5*(BvVelGradz(i,j)+BwVelGrady(i,j))
            BS33(i,j)=BwVelGradz(i,j)
c
            BS11T=BS11(i,j)-Bterm1
            BS12T=BS12(i,j)
            BS13T=BS13(i,j)
            BS22T=BS22(i,j)-Bterm1
            BS23T=BS23(i,j)
            BS33T=BS33(i,j)-Bterm1
c          
            BSstar(i,j)=dsqrt(2.*
     *            (BS11T*BS11T+BS12T*BS12T+BS13T*BS13T+
     *             BS12T*BS12T+BS22T*BS22T+BS23T*BS23T+
     *             BS13T*BS13T+BS23T*BS23T+BS33T*BS33T))
c
          enddo
        enddo
c 
        do i=1,NumberOfElements    
c
          TurbulentKE(i)=TurbulentViscosity(i)*
     *                    Sstar(i)/(factor*Density(i))
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BTurbulentKE(i,j)=BTurbulentViscosity(i,j)*
     *                    BSstar(i,j)/(factor*BDensity(i,j))
c
          enddo
        enddo
c
        deallocate(Sstar)
        deallocate(BSstar)
c
      else
c
        factor=dsqrt(cmu)
c
        do i=1,NumberOfElements    
c
          TurbulentKE(i)=TurbulentViscosity(i)*
     *                    StrainRate(i)/(factor*Density(i))
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BTurbulentKE(i,j)=BTurbulentViscosity(i,j)*
     *                    BStrainRate(i,j)/(factor*BDensity(i,j))
c
          enddo
        enddo
c
      endif

      return
      end
c
c#############################################################################################
c
      SUBROUTINE ExtractTurbulentKE
c
c#############################################################################################
c
      use Variables1, only: TurbulentKE,BTurbulentKE
      use Turbulence1, only: cmu50,StrainRate,BStrainRate
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use PhysicalProperties1, only: Density,BDensity,
     *                      TurbulentViscosity,BTurbulentViscosity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements    
c
        TurbulentKE(i)=TurbulentViscosity(i)*
     *                    StrainRate(i)/(cmu50*Density(i))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentKE(i,j)=BTurbulentViscosity(i,j)*
     *                    BStrainRate(i,j)/(cmu50*BDensity(i,j))
c
        enddo
      enddo
c
      return
      end
