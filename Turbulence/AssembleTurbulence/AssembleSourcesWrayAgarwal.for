c
c#############################################################################################
c
      SUBROUTINE AssembleWrayAgarwalSources
c
C#############################################################################################
c
      use User0, only: Lcompressible,WallTreatment,
     *                 MethodCalcGradientModifiedED,LWallDistanceFreeWA,
     *                 nIterGradientModifiedED,LimitGradientModifiedED,
     *                 LimitGradientModifiedEDMethod,LRough,GrainSize,
     *                 MethodCalcGradientDensity,nIterGradientDensity,
     *                 LimitGradientDensity,LimitGradientDensityMethod
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BFaceAreanx,BFaceAreany,
     *                     BFaceAreanz,BDistanceCFx,BDistanceCFy,
     *                     BDistanceCFz,BgDiff,BFaceTx,BFaceTy,BFaceTz
      use Variables1, only: ModifiedED,ModifiedEDGradx,
     *                      ModifiedEDGrady,ModifiedEDGradz,BModifiedED,
     *                      BModifiedEDGradx,BModifiedEDGrady,
     *                      BModifiedEDGradz
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density, BDensity,
     *                               eDiffCoefficient,BeDiffCoefficient,
     *                               DensGradx,DensGrady,DensGradz,
     *                               BDensGradx,BDensGrady,BDensGradz
      use Turbulence1, only: StrainRate,BStrainRate,SRateGradx,
     *                       BSRateGradx,SRateGrady,
     *                       BSRateGrady,SRateGradz,
     *                       BSRateGradz,f1WA,c1KW,c1KE,
     *                       c2KW,c2KE,cm,fr1Coefficient
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use WallDistance1, only: WallDistance,iTau
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,i1,i2,i3,i4,j
      double precision :: term1,term2,term3,term4,sum,c1,gamf,
     *                    dfidxTf,dfidyTf,dfidzTf,FluxCflocal,
     *                    FluxFflocal,FluxVflocal,ErTerm,c2KWLocal
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------
          character*10 Variable
          integer :: MethodCalcGradient,nIterGradientPhi
          logical :: LimitGradient
          integer :: LimitGradientMethod
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------
      end interface
c-------------------------------------------------------------------------------------------
c
      do i=1,NumberOfElements     
c
        StrainRate(i)=dmax1(StrainRate(i),1.d-16)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BStrainRate(i,j)=dmax1(BStrainRate(i,j),1.d-16)
c
        enddo
      enddo
c
c---- Calculate Gradients of the phi variable
c
      Variable='strainrate'
      call Gradient(Variable,MethodCalcGradientModifiedED,
     *      StrainRate,SRateGradx,SRateGrady,SRateGradz,
     *        BStrainRate,BSRateGradx,BSRateGrady,
     *          BSRateGradz,nIterGradientModifiedED,
     *            LimitGradientModifiedED,LimitGradientModifiedEDMethod)
c
c--- Calculate turbulence production
c
      if(LWallDistanceFreeWA) then
c
        do i=1,NumberOfElements     
c
          c1=f1WA(i)*(c1KW-c1KE)+c1KE
          term1=c1*fr1Coefficient(i)*Density(i)*
     *             StrainRate(i)*ModifiedED(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-term1
c
          term1=f1WA(i)*c2KW*Density(i)/StrainRate(i)
          term2=ModifiedEDGradx(i)*SRateGradx(i)+
     *               ModifiedEDGrady(i)*SRateGrady(i)+
     *                  ModifiedEDGradz(i)*SRateGradz(i)
          term1=term1*term2*Volume(i)
c        
          FluxCE(i)=FluxCE(i)-dmin1(term1,0.)
          FluxTE(i)=FluxTE(i)-term1*ModifiedED(i)
c
          term1=SRateGradx(i)**2+SRateGrady(i)**2+SRateGradz(i)**2
          term1=term1*c2KE*((ModifiedED(i)/StrainRate(i))**2)
          term2=cm*c2KE*(ModifiedEDGradx(i)**2+
     *              ModifiedEDGrady(i)**2+ModifiedEDGradz(i)**2)
c
          term3=(1.-f1WA(i))*Density(i)*dmin1(term1,term2)*Volume(i)
          FluxCE(i)=FluxCE(i)+term3/ModifiedED(i)
          FluxTE(i)=FluxTE(i)+term3
c
        enddo
c
        if(Lcompressible) then
c
c---- Calculate Density Gradients
c
          Variable='density'
          call Gradient(Variable,MethodCalcGradientDensity,
     *      Density,DensGradx,DensGrady,DensGradz,
     *        BDensity,BDensGradx,BDensGrady,
     *          BDensGradz,nIterGradientDensity,
     *            LimitGradientDensity,LimitGradientDensityMethod)
c
          do i=1,NumberOfElements    
c
            term1=-(DensGradx(i)*ModifiedEDGradx(i)+
     *              DensGrady(i)*ModifiedEDGrady(i)+
     *                 DensGradz(i)*ModifiedEDGradz(i))
            term1=term1*Volume(i)*eDiffCoefficient(i)/Density(i)
c
            FluxTE(i)=FluxTE(i)-term1
c
          enddo
c
        endif
c
        if(WallTreatment.eq.'lowreynoldsnumber') then
c      
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
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
            BModifiedED(i2,i3)=0.
            if(LRough) BModifiedED(i2,i3)=
     *         ((18133.*GrainSize(i2)-58.4)*GrainSize(i2)+
     *                              0.0999)*GrainSize(i2)+0.0000354
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BModifiedEDGradx(i2,i3)
            dfidyTf=BModifiedEDGrady(i2,i3)
            dfidzTf=BModifiedEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                    dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*ModifiedED(i1)+
     *                      FluxFflocal*BModifiedED(i2,i3)+FluxVflocal
c
          enddo
c
        elseif(WallTreatment.eq.'wallfunctions') then
c      
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            BModifiedED(i2,i3)=ModifiedED(i1)
c
          enddo            
c
        endif
c 
      else
c
        do i=1,NumberOfElements     
c
          c1=f1WA(i)*(c1KW-c1KE)+c1KE
          term1=c1*fr1Coefficient(i)*Density(i)*
     *             StrainRate(i)*ModifiedED(i)*Volume(i)
          FluxTE(i)=FluxTE(i)-term1
c
          c2KWLocal=c2KW
          if(LRough) c2KWLocal=c2KW*WallDistance(i)/
     *                     (WallDistance(i)+0.03*GrainSize(iTau(i)))
          term1=f1WA(i)*c2KWLocal*Density(i)/StrainRate(i)
          term2=ModifiedEDGradx(i)*SRateGradx(i)+
     *               ModifiedEDGrady(i)*SRateGrady(i)+
     *                  ModifiedEDGradz(i)*SRateGradz(i)
          term1=term1*term2*Volume(i)
c        
          FluxCE(i)=FluxCE(i)-dmin1(term1,0.)
          FluxTE(i)=FluxTE(i)-term1*ModifiedED(i)
c
          term1=SRateGradx(i)**2+SRateGrady(i)**2+SRateGradz(i)**2
c
c----- A modified implementation of wray agarwal
c
!          term1=term1*((ModifiedED(i)/StrainRate(i))**2)
!          term2=ModifiedEDGradx(i)**2+
!     *              ModifiedEDGrady(i)**2+ModifiedEDGradz(i)**2
!c
!c          ErTerm=cm*term2*dtanh(term1/(cm*term2))  !either can be used
!          ErTerm=dmin1(term1,cm*term2)
!          term3=(1.-f1WA(i))*c2KE*Density(i)*ErTerm*Volume(i)
!          FluxCE(i)=FluxCE(i)+term3/dmax1(ModifiedED(i),tiny)
!          FluxTE(i)=FluxTE(i)+term3
c
c---- Original implementation (NASA website)
c
          term1=term1*ModifiedED(i)/(StrainRate(i)**2)
          term3=(1.-f1WA(i))*c2KE*Density(i)*term1*Volume(i)
          FluxCE(i)=FluxCE(i)+term3
          FluxTE(i)=FluxTE(i)+term3*dmax1(ModifiedED(i),0.)
c
        enddo
c
        if(Lcompressible) then
c
c---- Calculate Density Gradients
c
          Variable='density'
          call Gradient(Variable,MethodCalcGradientDensity,
     *      Density,DensGradx,DensGrady,DensGradz,
     *        BDensity,BDensGradx,BDensGrady,
     *          BDensGradz,nIterGradientDensity,
     *            LimitGradientDensity,LimitGradientDensityMethod)
c
          do i=1,NumberOfElements    
c
            term1=-(DensGradx(i)*ModifiedEDGradx(i)+
     *              DensGrady(i)*ModifiedEDGrady(i)+
     *                 DensGradz(i)*ModifiedEDGradz(i))
            term1=term1*Volume(i)*eDiffCoefficient(i)/Density(i)
c
            FluxTE(i)=FluxTE(i)-term1
c
          enddo
c
        endif
c
        if(WallTreatment.eq.'lowreynoldsnumber') then
c      
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
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
            BModifiedED(i2,i3)=0.
c
            gamf=BeDiffCoefficient(i2,i3)
            dfidxTf=BModifiedEDGradx(i2,i3)
            dfidyTf=BModifiedEDGrady(i2,i3)
            dfidzTf=BModifiedEDGradz(i2,i3)
c
            FluxCflocal= gamf*BgDiff(i2,i3)
            FluxFflocal=-gamf*BgDiff(i2,i3)
            FluxVflocal=-gamf*(dfidxTf*BFaceTx(i2,i3)+
     *                    dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
            FluxCf(i4)=FluxCf(i4)+FluxCflocal
            FluxFf(i4)=FluxFf(i4)+FluxFflocal
            FluxVf(i4)=FluxVf(i4)+FluxVflocal
            FluxTf(i4)=FluxTf(i4)+FluxCflocal*ModifiedED(i1)+
     *                      FluxFflocal*BModifiedED(i2,i3)+FluxVflocal
c
          enddo
c
        elseif(WallTreatment.eq.'wallfunctions') then
c      
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            BModifiedED(i2,i3)=ModifiedED(i1)
c
          enddo            
c
        endif
c
      endif
c
      return
      end