c
C#############################################################################################
      SUBROUTINE InterpolateFromElementsToNodes(NF,index,LGradient)
C#############################################################################################
c
      use User0, only:LsolveMomentum,NumberOfrFieldsToSolve,
     *                TurbulenceModel
      use ReferenceValues1, only: UInfinity
      use Geometry1
      use Geometry3, only:NBFaces,NBFacesMax,NBFaceNodes,NIFaces,
     *                    NumberOfElementFaceNodes,
     *                    GlobalFaceNumberOfNodes,NBFaceOwner
      use Geometry4, only:BFaceCentroidx,BFaceCentroidy,
     *                    BFaceCentroidz,xc,yc,zc
      use Variables1
      use Scalar1, only: Scalar,BScalar,
     *                   ScalarGradx,ScalarGrady,ScalarGradz,
     *                   BScalarGradx,BScalarGrady,BScalarGradz
      use VolumeOfFluid1, only: rField,BrField,
     *                          rFieldGradx,rFieldGradY,rFieldGradZ,
     *                          BrFieldGradx,BrFieldGradY,BrFieldGradZ
      use PhysicalProperties1, only: Density,BDensity,
     *                     TurbulentViscosity,BTurbulentViscosity,
     *                     Viscosity,BViscosity,SpecificHeat,
     *                     BSpecificHeat,ReferenceTemperature,
     *                     eDiffCoefficient,BeDiffCoefficient
      use Turbulence1, only: ReT,BReT,
     *                       TurbulentViscosityTs,BTurbulentViscosityTs,
     *                       TurbulentViscosityTl,BTurbulentViscosityTl,
     *                       ProductionKT,ProductionKL,StrainRate,
     *                       BStrainRate,Vorticity,BVorticity
      use tecplot1
      use WallDistance1, only: WallDistance,BWallDistance
      use Constants1, only: tiny,twothird
c*********************************************************************************************
      implicit none
C********************************************************************************************
      character*10 Variable
      character*15 zonetype
      logical :: LGradient
      integer i,j,k,k1,NF,j1,j2,index,NFTemp
      double precision d,Tref
C********************************************************************************************
c
      if(.not.LGradient.and.index.eq.1) then     
c
        Tref=ReferenceTemperature
c 
        if(NF.eq.1) then
c
          dummyE=uVelocity
          BdummyE=BuVelocity
c
        elseif(NF.eq.2) then
c
          dummyE=vVelocity
          BdummyE=BvVelocity
c
        elseif(NF.eq.3) then
c
          dummyE=wVelocity
          BdummyE=BwVelocity
c
        elseif(NF.eq.14) then
c
          dummyE=dsqrt(uVelocity**2+vVelocity**2+wVelocity**2)
          BdummyE=dsqrt(BuVelocity**2+BvVelocity**2+BwVelocity**2)
c
        elseif(NF.eq.4) then
c
          dummyE=Pressure
          BdummyE=BPressure
c
        elseif(NF.eq.5) then
c
          dummyE=Temperature
          BdummyE=BTemperature
c
        elseif(NF.eq.6) then
c
          dummyE=Density
          BdummyE=BDensity
c
        elseif(NF.eq.7) then
c
          dummyE=MachNumber
          BdummyE=BMachNumber
c
        elseif(NF.eq.8) then
c
          dummyE=TurbulentKE
          BdummyE=BTurbulentKE
c
        elseif(NF.eq.9) then
c
          dummyE=TurbulentED
          BdummyE=BTurbulentED
c
        elseif(NF.eq.37) then
c
          dummyE=TurbulentV2
          BdummyE=BTurbulentV2
c
        elseif(NF.eq.38) then
c
          dummyE=TurbulentZeta
          BdummyE=BTurbulentZeta
c
        elseif(NF.eq.39) then
c
          dummyE=TfRelaxation
          BdummyE=BTfRelaxation
c
        elseif(NF.eq.10) then
c
          dummyE=TurbulentOmega
          BdummyE=BTurbulentOmega
c
        elseif(NF.eq.18) then
c
          dummyE=ModifiedED
          BdummyE=BModifiedED
c
        elseif(NF.eq.19) then
c
          dummyE=Viscosity
          BdummyE=BViscosity
c
        elseif(NF.eq.11) then
c
          dummyE=TurbulentViscosity
          BdummyE=BTurbulentViscosity
c
        elseif(NF.eq.12) then
c
          do i=1,NumberOfElements
            dummyE(i)=TurbulenceProduction(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              k=NBFaceOwner(i,j)
              BdummyE(i,j)=TurbulenceProduction(k)
            enddo
          enddo
c
        elseif(NF.eq.13) then
c
          dummyE=TurbulenceProductionB
          BdummyE=BTurbulenceProductionB
c
        elseif(NF.eq.15) then
c
          dummyE=ReT
          BdummyE=BReT
c
        elseif(NF.eq.16) then
c
          dummyE=yplusPlot
          BdummyE=ByplusPlot
          deallocate(yplusPlot)
          deallocate(ByplusPlot)
c
        elseif(NF.eq.17) then
c
          dummyE=uplusPlot
          BdummyE=BuplusPlot
          deallocate(uplusPlot)
          deallocate(BuplusPlot)
c
        elseif(NF.eq.20) then
c
          dummyE=TurbulentKL
          BdummyE=BTurbulentKL
c
        elseif(NF.eq.21) then
c
          dummyE=WallDistance
          BdummyE=BWallDistance
c
        elseif(NF.eq.22) then
c
          do i=1,NumberOfElements
            dummyE(i)=SpecificHeat(i)*(Temperature(i)-Tref)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BSpecificHeat(i,j)*
     *                 (BTemperature(i,j)-Tref)
            enddo
          enddo
c
        elseif(NF.eq.23) then
c
          do i=1,NumberOfElements
            dummyE(i)=Htotal(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BHtotal(i,j)
            enddo
          enddo
c
        elseif(NF.eq.24) then
c
          do i=1,NumberOfElements
            dummyE(i)=Temperature(i)
            if(LsolveMomentum) dummyE(i)=dummyE(i)+0.5*(uVelocity(i)**2+
     *                  vVelocity(i)**2+wVelocity(i)**2)/SpecificHeat(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BTemperature(i,j)
              if(LsolveMomentum) BdummyE(i,j)=BdummyE(i,j)+
     *             0.5*(BuVelocity(i,j)**2+BvVelocity(i,j)**2+
     *                            BwVelocity(i,j)**2)/BSpecificHeat(i,j)
            enddo
          enddo
c
        elseif(NF.eq.25) then
c
          do i=1,NumberOfElements
            dummyE(i)=TGamma(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BTGamma(i,j)
            enddo
          enddo
c
        elseif(NF.eq.26) then
c
          do i=1,NumberOfElements
            dummyE(i)=TReTheta(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BTReTheta(i,j)
            enddo
          enddo
c
        elseif(NF.eq.27) then
c
          do i=1,NumberOfElements
            dummyE(i)=TurbulentKE(i)+TurbulentKL(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BTurbulentKE(i,j)+BTurbulentKL(i,j)
            enddo
          enddo
c
        elseif(NF.eq.28) then
c
          do i=1,NumberOfElements
            dummyE(i)=Density(i)*TurbulentViscosityTl(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BDensity(i,j)*BTurbulentViscosityTl(i,j)
            enddo
          enddo
c
        elseif(NF.eq.29) then
c
          do i=1,NumberOfElements
            dummyE(i)=Density(i)*TurbulentViscosityTs(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BDensity(i,j)*BTurbulentViscosityTs(i,j)
            enddo
          enddo
c
        elseif(NF.eq.30) then
c
          do i=1,NumberOfElements
            dummyE(i)=Viscosity(i)+TurbulentViscosity(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BViscosity(i,j)+BTurbulentViscosity(i,j)
            enddo
          enddo
c
        elseif(NF.eq.31) then
c
          do i=1,NumberOfElements
            dummyE(i)=TurbulentViscosity(i)/Viscosity(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BTurbulentViscosity(i,j)/BViscosity(i,j)
            enddo
          enddo
c
        elseif(NF.eq.32) then
c
          if(TurbulenceModel.eq.'kklomega') then         
c
            do i=1,NumberOfElements
              dummyE(i)=100.*dsqrt(twothird*(TurbulentKE(i)+
     *                    TurbulentKL(i)))/dmax1(Uinfinity,tiny)
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
                BdummyE(i,j)=100.*dsqrt(twothird*(BTurbulentKE(i,j)+
     *                     BTurbulentKL(i,j)))/dmax1(Uinfinity,tiny)
              enddo
            enddo
c
          else
c
            do i=1,NumberOfElements
              dummyE(i)=100.*dsqrt(twothird*TurbulentKE(i))
     *                                   /dmax1(Uinfinity,tiny)
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
                BdummyE(i,j)=100.*dsqrt(twothird*BTurbulentKE(i,j))/
     *                                        dmax1(Uinfinity,tiny)
              enddo
            enddo
c
          endif
c
        elseif(NF.eq.33) then
c
          do i=1,NumberOfElements
            dummyE(i)=Density(i)*ProductionKT(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              k=NBFaceOwner(i,j)
              BdummyE(i,j)=Density(k)*ProductionKT(k)
            enddo
          enddo
c
        elseif(NF.eq.34) then
c
          do i=1,NumberOfElements
            dummyE(i)=Density(i)*ProductionKL(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              k=NBFaceOwner(i,j)
              BdummyE(i,j)=Density(k)*ProductionKL(k)
            enddo
          enddo
c
        elseif(NF.eq.35) then
c
          Variable='temp'
          call CalculateEffectiveDiffusionCoefficient(Variable)
          do i=1,NumberOfElements
            dummyE(i)=eDiffCoefficient(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BeDiffCoefficient(i,j)
            enddo
          enddo
c
        elseif(NF.eq.36) then
c
          do i=1,NumberOfElements
            dummyE(i)=TGammaEff(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              k=NBFaceOwner(i,j)
              BdummyE(i,j)=TGammaEff(k)
            enddo
          enddo
c
        elseif(NF.eq.40) then
c
          do i=1,NumberOfElements
            dummyE(i)=LambdaELE(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              BdummyE(i,j)=BLambdaELE(i,j)
            enddo
          enddo
c
        elseif(NF.eq.41) then
c
          do i=1,NumberOfElements
            dummyE(i)=InitialVelDivergence(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              k=NBFaceOwner(i,j)
              BdummyE(i,j)=InitialVelDivergence(k)
            enddo
          enddo
c
        elseif(NF.eq.42) then
c
          do i=1,NumberOfElements
            dummyE(i)=FinalVelDivergence(i)
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
              k=NBFaceOwner(i,j)
              BdummyE(i,j)=FinalVelDivergence(k)
            enddo
          enddo
c
        elseif(NF.gt.42) then
c
          if(NF.gt.(42+NumberOfrFieldsToSolve))then
c
            do i=1,NumberOfElements
              dummyE(i)=Scalar(i,NF-NumberOfrFieldsToSolve-42)
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
                BdummyE(i,j)=BScalar(i,j,NF-NumberOfrFieldsToSolve-42)
              enddo
            enddo
c
          else
c
            do i=1,NumberOfElements
              dummyE(i)=rField(i,NF-42)
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
                BdummyE(i,j)=BrField(i,j,NF-42)
              enddo
            enddo
c
          endif
c
        endif
c
      elseif(LGradient.and.index.eq.2) then
c 
        if(NF.eq.1) then
c
          dummyE=uVelGradx
          BdummyE=BuVelGradx
c
        elseif(NF.eq.2) then
c
          dummyE=uVelGrady
          BdummyE=BuVelGrady
c
        elseif(NF.eq.3) then
c
          dummyE=uVelGradz
          BdummyE=BuVelGradz
c 
        elseif(NF.eq.4) then
c
          dummyE=vVelGradx
          BdummyE=BvVelGradx
c
        elseif(NF.eq.5) then
c
          dummyE=vVelGrady
          BdummyE=BvVelGrady
c
        elseif(NF.eq.6) then
c
          dummyE=vVelGradz
          BdummyE=BvVelGradz
c 
        elseif(NF.eq.7) then
c
          dummyE=wVelGradx
          BdummyE=BwVelGradx
c
        elseif(NF.eq.8) then
c
          dummyE=wVelGrady
          BdummyE=BwVelGrady
c
        elseif(NF.eq.9) then
c
          dummyE=wVelGradz
          BdummyE=BwVelGradz
c 
        elseif(NF.eq.10) then
c
          dummyE=PressGradx
          BdummyE=BPressGradx
c
        elseif(NF.eq.11) then
c
          dummyE=PressGrady
          BdummyE=BPressGrady
c
        elseif(NF.eq.12) then
c
          dummyE=PressGradz
          BdummyE=BPressGradz
c 
        elseif(NF.eq.13) then
c
          dummyE=TKEGradx
          BdummyE=BTKEGradx
c
        elseif(NF.eq.14) then
c
          dummyE=TKEGrady
          BdummyE=BTKEGrady
c
        elseif(NF.eq.15) then
c
          dummyE=TKEGradz
          BdummyE=BTKEGradz
c 
        elseif(NF.eq.16) then
c
          dummyE=TEDGradx
          BdummyE=BTEDGradx
c
        elseif(NF.eq.17) then
c
          dummyE=TEDGrady
          BdummyE=BTEDGrady
c
        elseif(NF.eq.18) then
c
          dummyE=TEDGradz
          BdummyE=BTEDGradz
c 
        elseif(NF.eq.19) then
c
          dummyE=TOmegaGradx
          BdummyE=BTOmegaGradx
c
        elseif(NF.eq.20) then
c
          dummyE=TOmegaGrady
          BdummyE=BTOmegaGrady
c
        elseif(NF.eq.21) then
c
          dummyE=TOmegaGradz
          BdummyE=BTOmegaGradz
c 
        elseif(NF.eq.22) then
c
          dummyE=TGammaGradx
          BdummyE=BTGammaGradx
c
        elseif(NF.eq.23) then
c
          dummyE=TGammaGrady
          BdummyE=BTGammaGrady
c
        elseif(NF.eq.24) then
c
          dummyE=TGammaGradz
          BdummyE=BTGammaGradz
c 
        elseif(NF.eq.25) then
c
          dummyE=TReThetaGradx
          BdummyE=BTReThetaGradx
c
        elseif(NF.eq.26) then
c
          dummyE=TReThetaGrady
          BdummyE=BTReThetaGrady
c
        elseif(NF.eq.27) then
c
          dummyE=TReThetaGradz
          BdummyE=BTReThetaGradz
c 
        elseif(NF.eq.28) then
c
          dummyE=TurbulentKLGradx
          BdummyE=BTurbulentKLGradx
c
        elseif(NF.eq.29) then
c
          dummyE=TurbulentKLGrady
          BdummyE=BTurbulentKLGrady
c
        elseif(NF.eq.30) then
c
          dummyE=TurbulentKLGradz
          BdummyE=BTurbulentKLGradz
c 
        elseif(NF.eq.31) then
c
          dummyE=ModifiedEDGradx
          BdummyE=BModifiedEDGradx
c
        elseif(NF.eq.32) then
c
          dummyE=ModifiedEDGrady
          BdummyE=BModifiedEDGrady
c
        elseif(NF.eq.33) then
c
          dummyE=ModifiedEDGradz
          BdummyE=BModifiedEDGradz
c 
        elseif(NF.eq.34) then
c
          dummyE=TempGradx
          BdummyE=BTempGradx
c
        elseif(NF.eq.35) then
c
          dummyE=TempGrady
          BdummyE=BTempGrady
c
        elseif(NF.eq.36) then
c
          dummyE=TempGradz
          BdummyE=BTempGradz
c 
        elseif(NF.eq.37) then
c
          dummyE=HtotalGradx
          BdummyE=BHtotalGradx
c
        elseif(NF.eq.38) then
c
          dummyE=HtotalGrady
          BdummyE=BHtotalGrady
c
        elseif(NF.eq.39) then
c
          dummyE=HtotalGradz
          BdummyE=BHtotalGradz
c
        elseif(NF.eq.40) then
c
          dummyE=StrainRate
          BdummyE=BStrainRate
c
        elseif(NF.eq.41) then
c
          dummyE=Vorticity
          BdummyE=BVorticity
c
        elseif(NF.eq.42) then
c
          dummyE=TurbulentV2Gradx
          BdummyE=BTurbulentV2Gradx
c
        elseif(NF.eq.43) then
c
          dummyE=TurbulentV2Grady
          BdummyE=BTurbulentV2Grady
c
        elseif(NF.eq.44) then
c
          dummyE=TurbulentV2Gradz
          BdummyE=BTurbulentV2Gradz
c
        elseif(NF.eq.45) then
c
          dummyE=TurbulentZetaGradx
          BdummyE=BTurbulentZetaGradx
c
        elseif(NF.eq.46) then
c
          dummyE=TurbulentZetaGrady
          BdummyE=BTurbulentZetaGrady
c
        elseif(NF.eq.47) then
c
          dummyE=TurbulentZetaGradz
          BdummyE=BTurbulentZetaGradz
c
        elseif(NF.eq.48) then
c
          dummyE=TfRelaxationGradx
          BdummyE=BTfRelaxationGradx
c
        elseif(NF.eq.49) then
c
          dummyE=TfRelaxationGrady
          BdummyE=BTfRelaxationGrady
c
        elseif(NF.eq.50) then
c
          dummyE=TfRelaxationGradz
          BdummyE=BTfRelaxationGradz
c
        elseif(NF.eq.51) then
c
          dummyE=LambdaELEGradx
          BdummyE=BLambdaELEGradx
c
        elseif(NF.eq.52) then
c
          dummyE=LambdaELEGrady
          BdummyE=BLambdaELEGrady
c
        elseif(NF.eq.53) then
c
          dummyE=LambdaELEGradz
          BdummyE=BLambdaELEGradz
c
        elseif(NF.gt.53) then
c
          if(NF.gt.(53+3*NumberOfrFieldsToSolve))then
c
            NFTemp=mod(NF,3)
c
            if(NFTemp.eq.0) then
              do i=1,NumberOfElements
                dummyE(i)=ScalarGradx
     *                   (i,(NF-3*NumberOfrFieldsToSolve-51)/3)
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
                  BdummyE(i,j)=BScalarGradx
     *                   (i,j,(NF-3*NumberOfrFieldsToSolve-51)/3)
                enddo
              enddo
c
            elseif(NFTemp.eq.1) then
              do i=1,NumberOfElements
                dummyE(i)=ScalarGrady
     *                   (i,(NF-3*NumberOfrFieldsToSolve-52)/3)
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
                  BdummyE(i,j)=BScalarGrady
     *                   (i,j,(NF-3*NumberOfrFieldsToSolve-52)/3)
                enddo
              enddo
c
            elseif(NFTemp.eq.2) then
              do i=1,NumberOfElements
                dummyE(i)=ScalarGradz
     *                   (i,(NF-3*NumberOfrFieldsToSolve-53)/3)
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
                  BdummyE(i,j)=BScalarGradz
     *                   (i,j,(NF-3*NumberOfrFieldsToSolve-53)/3)
                enddo
              enddo
c
            endif
c
          else
c
            NFTemp=mod(NF,3)
c
            if(NFTemp.eq.0) then
              do i=1,NumberOfElements
                dummyE(i)=rFieldGradx(i,(NF-51)/3)
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
                  BdummyE(i,j)=BrFieldGradx(i,j,(NF-51)/3)
                enddo
              enddo
c
            elseif(NFTemp.eq.1) then
              do i=1,NumberOfElements
                dummyE(i)=rFieldGrady(i,(NF-52)/3)
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
                  BdummyE(i,j)=BrFieldGrady(i,j,(NF-52)/3)
                enddo
              enddo
c
            elseif(NFTemp.eq.2) then
              do i=1,NumberOfElements
                dummyE(i)=rFieldGradz(i,(NF-53)/3)
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
                  BdummyE(i,j)=BrFieldGradz(i,j,(NF-53)/3)
                enddo
              enddo
c
            endif
c
          endif
c
        endif
c
      endif
c
c--- Start
c
      allocate(SumInvd(NumberOfNodes))
      allocate(SumPhiInvd(NumberOfNodes))
      allocate(dummyN(NumberOfNodes))
c
c--- First find for interior nodes
c
      SumInvd=0.
      SumPhiInvd=0.
c
      do i=1,NumberOfElements
        do j=1,NumbOfElementNodes(i)
c
          k=ListOfElementNodes(i,j)

          if(NodeFlag(k).eq.1) then
c
            d=dsqrt((xc(i)-x(k))**2+(yc(i)-y(k))**2+(zc(i)-z(k))**2)
            SumInvd(k)=SumInvd(k)+1./d
            SumPhiInvd(k)=SumPhiInvd(k)+dummyE(i)/d
c
          endif
c
        enddo
      enddo
c
      do i=1,NumberOfElements
        do j=1,NumbOfElementNodes(i)
c
          k=ListOfElementNodes(i,j)

          if(NodeFlag(k).eq.1) then
c
            dummyN(k)=SumPhiInvd(k)/SumInvd(k)
c
          endif
c
        enddo
      enddo
c
c--- Second find for boundary nodes
c
      j1=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
c            j1=NElementBC(i,j)
c            j2=NElementBCFace(i,j)
c
          j1=j1+1
          do k=1,GlobalFaceNumberOfNodes(j1)
c            do k=1,NumberOfElementFaceNodes(j1,j2)
c
            k1=NBFaceNodes(i,j,k)
            d=dsqrt((BFaceCentroidx(i,j)-x(k1))**2+
     *                   (BFaceCentroidy(i,j)-y(k1))**2+
     *                         (BFaceCentroidz(i,j)-z(k1))**2)
            SumInvd(k1)=SumInvd(k1)+1./d
            SumPhiInvd(k1)=SumPhiInvd(k1)+BdummyE(i,j)/d        
c               
          enddo
c
        enddo
      enddo
c
      j1=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
c          j1=NElementBC(i,j)
c          j2=NElementBCFace(i,j)
c
          j1=j1+1
          do k=1,GlobalFaceNumberOfNodes(j1)
c          do k=1,NumberOfElementFaceNodes(j1,j2)
c
            k1=NBFaceNodes(i,j,k)
            dummyN(k1)=SumPhiInvd(k1)/SumInvd(k1)
c               
          enddo
        enddo
      enddo
c
      deallocate(SumInvd)
      deallocate(SumPhiInvd)
c
      return
      end