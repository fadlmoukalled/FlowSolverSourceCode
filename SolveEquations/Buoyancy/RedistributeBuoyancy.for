c
c#############################################################################################
c
      SUBROUTINE reDistributeTheBuoyancyTerm
c
c#############################################################################################
c
      use User0, only: BuoyancyModel,nIterStartApplyingHR,
     *                 ConvectionSchemeEnergy,BleedEnergy,
     *                 LNVFEnergy,LTVDEnergy
      use PhysicalProperties1, only:CoefficientOfThermalExpansion,
     *                              ReferenceDensity,GravityX,GravityY,
     *                              GravityZ,ReferenceTemperature,
     *                              Densityf,BDensity
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor,
     *                     NBFaces,NBFaceOwner
      use Geometry4, only: FaceAreax,FaceAreay,FaceAreaz,
     *                     BFaceAreax,BFaceAreay,BFaceAreaz,
     *                     xc,yc,zc,Volume,
     *                     FaceCentroidx,FaceCentroidy,FaceCentroidz,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
      use MultiGrid2, only: nIter
      use Variables1, only: Temperature,BTemperature,
     *                      TempGradx,TempGrady,TempGradz,
     *                      Buoyancyx,Buoyancyy,Buoyancyz,
     *                      BBuoyancyx,BBuoyancyy,BBuoyancyz,
     *                      Buoyancyfx,Buoyancyfy,Buoyancyfz
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      character*20, save :: ConvectionScheme1
      double precision, save :: dxi,dyi,dzi,dxj,dyj,dzj,BuoyancyDotdi,
     *                          BuoyancyDotdj
      double precision, save, dimension(:), allocatable :: Temperaturef
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE InterpolateToFaceUsingHRscheme(ConvectionScheme,
     *               Bleed,NVF,TVD,FiT,dfidxT,dfidyT,dfidzT,FiTf)
c--------------------------------------------------------------
          character*20 ConvectionScheme
          logical NVF,TVD
          double precision :: Bleed
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:) :: FiTf
c--------------------------------------------------------------
        end SUBROUTINE InterpolateToFaceUsingHRscheme
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      if(BuoyancyModel.eq.'boussinesq') then
c
        allocate(Temperaturef(NIFaces))
c
        if(nIter.ge.nIterStartApplyingHR) then
c
          call InterpolateToFaceUsingHRscheme(ConvectionSchemeEnergy,
     *     BleedEnergy,LNVFEnergy,LTVDEnergy,Temperature,TempGradx,
     *             TempGrady,TempGradz,Temperaturef)
c
        else
c
          ConvectionScheme1='upwind'
          call InterpolateToFaceUsingHRscheme(ConvectionScheme1,
     *     BleedEnergy,LNVFEnergy,LTVDEnergy,Temperature,TempGradx,
     *             TempGrady,TempGradz,Temperaturef)
c
        endif
c
        do i=1,NIFaces
c
          Buoyancyfx(i)=-ReferenceDensity*CoefficientOfThermalExpansion*
     *                   (Temperaturef(i)-ReferenceTemperature)*GravityX
          Buoyancyfy(i)=-ReferenceDensity*CoefficientOfThermalExpansion*
     *                   (Temperaturef(i)-ReferenceTemperature)*GravityY
          Buoyancyfz(i)=-ReferenceDensity*CoefficientOfThermalExpansion*
     *                   (Temperaturef(i)-ReferenceTemperature)*GravityZ
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BBuoyancyx(i,j)=-ReferenceDensity*
     *                         CoefficientOfThermalExpansion*
     *               (BTemperature(i,j)-ReferenceTemperature)*GravityX
            BBuoyancyy(i,j)=-ReferenceDensity*
     *                         CoefficientOfThermalExpansion*
     *               (BTemperature(i,j)-ReferenceTemperature)*GravityY
            BBuoyancyz(i,j)=-ReferenceDensity*
     *                         CoefficientOfThermalExpansion*
     *               (BTemperature(i,j)-ReferenceTemperature)*GravityZ
c
          enddo
        enddo
c
        Buoyancyx=0.
        Buoyancyy=0.
        Buoyancyz=0.
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          dxi=FaceCentroidx(k)-xc(i)
          dyi=FaceCentroidy(k)-yc(i)
          dzi=FaceCentroidz(k)-zc(i)
c
          dxj=FaceCentroidx(k)-xc(j)
          dyj=FaceCentroidy(k)-yc(j)
          dzj=FaceCentroidz(k)-zc(j)
c
          BuoyancyDotdi=Buoyancyfx(k)*dxi+Buoyancyfy(k)*dyi+
     *                                        Buoyancyfz(k)*dzi         
          BuoyancyDotdj=Buoyancyfx(k)*dxj+Buoyancyfy(k)*dyj+
     *                                        Buoyancyfz(k)*dzj        
c
          Buoyancyx(i)=Buoyancyx(i)+BuoyancyDotdi*FaceAreax(k)
          Buoyancyy(i)=Buoyancyy(i)+BuoyancyDotdi*FaceAreay(k)
          Buoyancyz(i)=Buoyancyz(i)+BuoyancyDotdi*FaceAreaz(k)
c
          Buoyancyx(j)=Buoyancyx(j)-BuoyancyDotdj*FaceAreax(k)
          Buoyancyy(j)=Buoyancyy(j)-BuoyancyDotdj*FaceAreay(k)
          Buoyancyz(j)=Buoyancyz(j)-BuoyancyDotdj*FaceAreaz(k)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            dxi=BFaceCentroidx(i,j)-xc(k)
            dyi=BFaceCentroidy(i,j)-yc(k)
            dzi=BFaceCentroidz(i,j)-zc(k)
c
            BuoyancyDotdi=BBuoyancyx(i,j)*dxi+BBuoyancyy(i,j)*dyi+
     *                                           BBuoyancyz(i,j)*dzi
c
            Buoyancyx(k)=Buoyancyx(k)+BuoyancyDotdi*BFaceAreax(i,j)
            Buoyancyy(k)=Buoyancyy(k)+BuoyancyDotdi*BFaceAreay(i,j)
            Buoyancyz(k)=Buoyancyz(k)+BuoyancyDotdi*BFaceAreaz(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          Buoyancyx(i)=Buoyancyx(i)/Volume(i)
          Buoyancyy(i)=Buoyancyy(i)/Volume(i)
          Buoyancyz(i)=Buoyancyz(i)/Volume(i)
c
        enddo
c
        deallocate(Temperaturef)
c
      elseif(BuoyancyModel.eq.'rhog') then
c
        do i=1,NIFaces
c
          Buoyancyfx(i)=Densityf(i)*GravityX
          Buoyancyfy(i)=Densityf(i)*GravityY
          Buoyancyfz(i)=Densityf(i)*GravityZ
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BBuoyancyx(i,j)=BDensity(i,j)*GravityX
            BBuoyancyy(i,j)=BDensity(i,j)*GravityY
            BBuoyancyz(i,j)=BDensity(i,j)*GravityZ
c
          enddo
        enddo
c
        Buoyancyx=0.
        Buoyancyy=0.
        Buoyancyz=0.
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          dxi=FaceCentroidx(k)-xc(i)
          dyi=FaceCentroidy(k)-yc(i)
          dzi=FaceCentroidz(k)-zc(i)
c
          dxj=FaceCentroidx(k)-xc(j)
          dyj=FaceCentroidy(k)-yc(j)
          dzj=FaceCentroidz(k)-zc(j)
c
          BuoyancyDotdi=Buoyancyfx(k)*dxi+Buoyancyfy(k)*dyi+         
     *                                       Buoyancyfz(k)*dzi         
          BuoyancyDotdj=Buoyancyfx(k)*dxj+Buoyancyfy(k)*dyj+        
     *                                       Buoyancyfz(k)*dzj         
c
          Buoyancyx(i)=Buoyancyx(i)+BuoyancyDotdi*FaceAreax(k)
          Buoyancyy(i)=Buoyancyy(i)+BuoyancyDotdi*FaceAreay(k)
          Buoyancyz(i)=Buoyancyz(i)+BuoyancyDotdi*FaceAreaz(k)
c
          Buoyancyx(j)=Buoyancyx(j)-BuoyancyDotdj*FaceAreax(k)
          Buoyancyy(j)=Buoyancyy(j)-BuoyancyDotdj*FaceAreay(k)
          Buoyancyz(j)=Buoyancyz(j)-BuoyancyDotdj*FaceAreaz(k)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            dxi=BFaceCentroidx(i,j)-xc(k)
            dyi=BFaceCentroidy(i,j)-yc(k)
            dzi=BFaceCentroidz(i,j)-zc(k)
c
            BuoyancyDotdi=BBuoyancyx(i,j)*dxi+
     *                     BBuoyancyy(i,j)*dyi+BBuoyancyz(i,j)*dzi
c
            Buoyancyx(k)=Buoyancyx(k)+BuoyancyDotdi*BFaceAreax(i,j)
            Buoyancyy(k)=Buoyancyy(k)+BuoyancyDotdi*BFaceAreay(i,j)
            Buoyancyz(k)=Buoyancyz(k)+BuoyancyDotdi*BFaceAreaz(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          Buoyancyx(i)=Buoyancyx(i)/Volume(i)
          Buoyancyy(i)=Buoyancyy(i)/Volume(i)
          Buoyancyz(i)=Buoyancyz(i)/Volume(i)
c
        enddo
c
      endif
c
      return
      end