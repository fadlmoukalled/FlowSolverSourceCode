c
c#############################################################################################
c
      SUBROUTINE CalculateTurbulentReynoldsNumber
c
c#############################################################################################
      use User0, only: TurbulenceModel
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Turbulence1, only: cmu,ReT,BReT,coefficientFW,BcoefficientFW
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,TurbulentOmega,BTurbulentOmega
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      if(TurbulenceModel.eq.'kepsilon'.or.
     *    TurbulenceModel.eq.'kepsilonrt'.or.
     *     TurbulenceModel.eq.'kepsilonv2f'.or.
     *      TurbulenceModel.eq.'kepsilonzetaf'.or.
     *       TurbulenceModel.eq.'kepsilonsharma'.or.
     *         TurbulenceModel.eq.'kepsilonchc'.or.
     *          TurbulenceModel.eq.'kepsilonchien'.or.
     *           TurbulenceModel.eq.'kepsilonkasagi'.or.
     *             TurbulenceModel.eq.'kepsilontagawa'.or.
     *               TurbulenceModel.eq.'kepsilonhishida'.or.
     *                 TurbulenceModel.eq.'kepsilonrng') then
c
        do i=1,NumberOfElements
c        
          ReT(i)=Density(i)*TurbulentKE(i)**2/
     *             (Viscosity(i)*dmax1(TurbulentED(i),tiny))
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c        
            BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)**2/
     *             (BViscosity(i,j)*dmax1(BTurbulentED(i,j),tiny))
c
          enddo
        enddo
c
      elseif(TurbulenceModel.eq.'komega'.or.
     *     TurbulenceModel.eq.'komegaepsilon'.or.
     *       TurbulenceModel.eq.'komega2006'.or.
     *        TurbulenceModel.eq.'komega2006lrn'.or.
     *              TurbulenceModel.eq.'komegabsl'.or.
     *                 TurbulenceModel.eq.'komegasst'.or.
     *                   TurbulenceModel.eq.'sstgama'.or.
     *                   TurbulenceModel.eq.'sstgamaretheta') then
c
        do i=1,NumberOfElements
c        
          ReT(i)=Density(i)*TurbulentKE(i)/
     *           (Viscosity(i)*dmax1(TurbulentOmega(i),tiny))
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c        
            BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)/
     *          (BViscosity(i,j)*dmax1(BTurbulentOmega(i,j),tiny))
c
          enddo
        enddo
c
      elseif(TurbulenceModel.eq.'kklomega') then
c
        do i=1,NumberOfElements
c        
          ReT(i)=Density(i)*TurbulentKE(i)*coefficientFW(i)**2/
     *                (Viscosity(i)*dmax1(TurbulentOmega(i),tiny))
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c        
            BReT(i,j)=BDensity(i,j)*BTurbulentKE(i,j)*
     *              BcoefficientFW(i,j)**2/(BViscosity(i,j)*
     *                           dmax1(BTurbulentOmega(i,j),tiny))
c
          enddo
        enddo
c
      endif
c
      return
      end