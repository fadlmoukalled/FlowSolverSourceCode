c
c#############################################################################################
c
      SUBROUTINE AssembleMdot
c
c#############################################################################################
c
      use User0, only: MethodDecomposeS,Lcompressible,LUnsteady,
     *                 LRelaxMomentum,MethodCalcGradientDensity,
     *                 nIterGradientDensity,LimitGradientDensity,
     *                 LimitGradientDensityMethod,LTVDDensity,
     *                 ConvectionSchemeDensity,LNVFDensity,
     *                 urfMomentum,LFalseTransientMomentum,
     *                 LBuoyancy,dt,LFreeSurfaceFlow,LSurfaceTension,
     *                 SurfaceTension,LHarmonic,nIterStartApplyingHR,
     *                 LSolveTurbulenceKineticEnergy,BleedDensity
      use VolumeOfFluid1, only: Curvature,BCurvature
      use TransferrField1, only: rFieldGradxT,rFieldGradyT,rFieldGradzT,
     *                           BrFieldGradxT,BrFieldGradyT,
     *                           BrFieldGradzT,rFieldGradfxT,
     *                           rFieldGradfyT,rFieldGradfzT,
     *                           rFieldT,BrFieldT
      use Multigrid2, only: nIter
      use BoundaryConditions2
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces,NIFaceOwner,NIFaceNeighbor,
     *                     NBFacesMax,NBFaceOwner
      use Geometry4, only: xc,yc,zc,Volume,
     *                     FaceAreax,FaceAreay,FaceAreaz,
     *                  DistanceCF,DistanceCFx,DistanceCFy,DistanceCFz,
     *                  DistanceCFux,DistanceCFuy,DistanceCFuz,
     *                     BFaceAreax,BFaceAreay,BFaceAreaz,
     *              BDistanceCF,BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     BDistanceCFux,BDistanceCFuy,BDistanceCFuz,
     *         GFactorCF,BFaceAreanx,BFaceAreany,BFaceAreanz,BFaceArea,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
      use Variables1, only: Du1Velocity,Dv1Velocity,Dw1Velocity,
     *                      uVelocityStar,vVelocityStar,wVelocityStar,
     *                      Pressure,PressGradx,PressGrady,PressGradz,
     *                      BPressure,BPressGradx,BPressGrady,
     *                      BPressGradz,uVelocity,vVelocity,wVelocity,
     *                      uVelocityOld,vVelocityOld,wVelocityOld,
     *                      uVelocityOldOld,vVelocityOldOld,
     *                      wVelocityOldOld,Temperature,mdot,mdotOld,
     *                      mdotOldOld,Bmdot,BmdotOld,BmdotOldOld,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      PressGradfx,PressGradfy,PressGradfz,
     *                      BTemperature,xVeldirection,yVeldirection,
     *                      zVeldirection,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
     *                      Buoyancyx,Buoyancyy,Buoyancyz,
     *                      Buoyancyfx,Buoyancyfy,Buoyancyfz,
     *                      BBuoyancyx,BBuoyancyy,BBuoyancyz
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf,acold,
     *                      acoldold,acStar
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only:Density,BDensity,
     *                              DensityOld,BDensityOld,
     *                              DensityOldOld,BDensityOldOld,
     *                              drhodp,Bdrhodp,Densityf,
     *                              BSpecificHeat,RGas,
     *                              DensityStar,BDensityStar
      use constants1, only: tiny
      use turbulence1, only: rhok,Brhok,drhokdx,drhokdy,drhokdz,
     *                       Bdrhokdx,Bdrhokdy,Bdrhokdz
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j,j1,j2,j3,k
      Character*16 Interpolation
      Character*20 ConvectionScheme1
      character*10 Variable
      double precision :: rhof,sfx,sfy,sfz,ex,ey,ez,cf,mdotBarf,
     *                    DuSf,DvSf,DwSf,DuEf,DvEf,DwEf,
     *                    DuTf,DvTf,DwTf,geoDiff,FluxCfLocal,
     *                    FluxFfLocal,FluxVfLocal,Magnitude,eDuSf,
     *                    eDvSf,eDwSf,mdotOldBar,mdotOldOldBar,mdotstar,
     *                    velMagnitude,xdirection,ydirection,zdirection,
     *                    TotalTemperature,dpdub,c1,diff,
     *                    velocity,velocity2,BuoyancyxBar,BuoyancyyBar,
     *                    BuoyancyzBar,BuoyancyBar,BuoyancyFace,
     *                    SurfaceTensionxBar,SurfaceTensionyBar,
     *                    SurfaceTensionzBar,CurvactureF
c
      double precision :: nx,ny,nz,vdotInfinity,vdotin,cinfinity,cin,
     *                    vdotb,cb,Tb,rhoInfinity,sbEntropy,rhob,
     *                    vel2,xdir,ydir,zdir,uTangential,vTangential,
     *                    wTangential,ratio,xF1,yF1,zF1
      double precision :: distance1,distance2,GFactCF,
     *                    DistCFx,DistCFy,DistCFz,DistCF,
     *                    DistCFux,DistCFuy,DistCFuz,uBar,vBar,wBar,
     *                    dfdxBar,dfdyBar,dfdzBar,correction,
     *                    dfdxf,dfdyf,dfdzf,rhoOf,uOf,vOf,wOf,rhoOOf,
     *                    uOOf,vOOf,wOOf,DVfO,DVfOO,uStarf,vStarf,
     *                    wStarf,uF1,vF1,wF1,uF1star,vF1star,wF1star,
     *                    pF1x,pF1y,pF1z,DuF1,DvF1,DwF1,value,dotproduct
c
      double precision, save, dimension(:), allocatable :: Du1f
      double precision, save, dimension(:), allocatable :: Dv1f
      double precision, save, dimension(:), allocatable :: Dw1f
      double precision, save, dimension(:), allocatable :: uVelBar
      double precision, save, dimension(:), allocatable :: vVelBar
      double precision, save, dimension(:), allocatable :: wVelBar
      double precision, save, dimension(:), allocatable :: PressGradxBar
      double precision, save, dimension(:), allocatable :: PressGradyBar
      double precision, save, dimension(:), allocatable :: PressGradzBar
      double precision, save, dimension(:), allocatable :: uVelOldf
      double precision, save, dimension(:), allocatable :: vVelOldf
      double precision, save, dimension(:), allocatable :: wVelOldf
      double precision, save, dimension(:), allocatable :: rhoOldf
      double precision, save, dimension(:), allocatable :: rhoOldOldf
      double precision,save, dimension(:), allocatable :: uVelocityStarf
      double precision,save, dimension(:), allocatable :: vVelocityStarf
      double precision,save, dimension(:), allocatable :: wVelocityStarf
c
      double precision, save, dimension(:), allocatable :: drhodpf
      double precision, save, dimension(:), allocatable :: drhodt
c
      double precision, save, dimension(:), allocatable :: drhoxOld
      double precision, save, dimension(:), allocatable :: drhoyOld
      double precision, save, dimension(:), allocatable :: drhozOld
      double precision, save, dimension(:,:), allocatable :: BdrhoxOld
      double precision, save, dimension(:,:), allocatable :: BdrhoyOld
      double precision, save, dimension(:,:), allocatable :: BdrhozOld
      double precision, save, dimension(:), allocatable :: drhoxOldOld
      double precision, save, dimension(:), allocatable :: drhoyOldOld
      double precision, save, dimension(:), allocatable :: drhozOldOld
      double precision, save, dimension(:,:),allocatable :: BdrhoxOldOld
      double precision, save, dimension(:,:),allocatable :: BdrhoyOldOld
      double precision, save, dimension(:,:),allocatable :: BdrhozOldOld
      double precision, save, dimension(:), allocatable :: DutOldT
      double precision, save, dimension(:), allocatable :: DutOldOldT
c
      double precision, save, dimension(:), allocatable :: uVelOldOldf
      double precision, save, dimension(:), allocatable :: vVelOldOldf
      double precision, save, dimension(:), allocatable :: wVelOldOldf
      double precision, save, dimension(:), allocatable :: DutOld
      double precision, save, dimension(:), allocatable :: DvtOld
      double precision, save, dimension(:), allocatable :: DwtOld
      double precision, save, dimension(:), allocatable :: DutOldOld
      double precision, save, dimension(:), allocatable :: DvtOldOld
      double precision, save, dimension(:), allocatable :: DwtOldOld
      double precision, save, dimension(:), allocatable :: DVtfO
      double precision, save, dimension(:), allocatable :: DVtfOO
c
      double precision, save, dimension(:), allocatable :: rhoStarf
      double precision, save, dimension(:), allocatable :: DutStar
      double precision, save, dimension(:), allocatable :: DvtStar
      double precision, save, dimension(:), allocatable :: DwtStar
      double precision, save, dimension(:), allocatable :: DVtfStar
      double precision, save, dimension(:), allocatable :: DutStarT
      double precision, save, dimension(:), allocatable :: drhoxStar
      double precision, save, dimension(:), allocatable :: drhoyStar
      double precision, save, dimension(:), allocatable :: drhozStar
      double precision, save, dimension(:,:), allocatable :: BdrhoxStar
      double precision, save, dimension(:,:), allocatable :: BdrhoyStar
      double precision, save, dimension(:,:), allocatable :: BdrhozStar
c
      double precision, save, dimension(:), allocatable :: drhokdxBar
      double precision, save, dimension(:), allocatable :: drhokdyBar
      double precision, save, dimension(:), allocatable :: drhokdzBar
      double precision, save, dimension(:), allocatable :: drhokdfx
      double precision, save, dimension(:), allocatable :: drhokdfy
      double precision, save, dimension(:), allocatable :: drhokdfz
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
       SUBROUTINE InterpolateElementToFace(InterpolationScheme,FiT,FiTf)
c--------------------------------------------------------------
          character*16 InterpolationScheme
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: FiTf
        end SUBROUTINE InterpolateElementToFace
c--------------------------------------------------------------
        SUBROUTINE InterpolateGradientToFace
     *   (GradientInterpolationScheme,FiT,dfidxT,dfidyT,dfidzT,
     *                                  dfidxfT,dfidyfT,dfidzfT)
c--------------------------------------------------------------
          character*16 GradientInterpolationScheme
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:) :: dfidxfT
          double precision, dimension(:) :: dfidyfT
          double precision, dimension(:) :: dfidzfT
c--------------------------------------------------------------
        end SUBROUTINE InterpolateGradientToFace
c--------------------------------------------------------------------------------
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
      allocate(Du1f(NIFaces))
      allocate(Dv1f(NIFaces))
      allocate(Dw1f(NIFaces))
      allocate(uVelBar(NIFaces))
      allocate(vVelBar(NIFaces))
      allocate(wVelBar(NIFaces))
      allocate(PressGradxBar(NIFaces))
      allocate(PressGradyBar(NIFaces))
      allocate(PressGradzBar(NIFaces))
      allocate(uVelocityStarf(NIFaces))
      allocate(vVelocityStarf(NIFaces))
      allocate(wVelocityStarf(NIFaces))
c
      if(LSolveTurbulenceKineticEnergy) then
c
        allocate(drhokdxBar(NIFaces))
        allocate(drhokdyBar(NIFaces))
        allocate(drhokdzBar(NIFaces))
        allocate(drhokdfx(NIFaces))
        allocate(drhokdfy(NIFaces))
        allocate(drhokdfz(NIFaces))
c
      endif
c
      interpolation='average'
c
      call InterpolateElementToFace(interpolation,Du1Velocity,Du1f)
      call InterpolateElementToFace(interpolation,Dv1Velocity,Dv1f)
      call InterpolateElementToFace(interpolation,Dw1Velocity,Dw1f)
      call InterpolateElementToFace(interpolation,uVelocity,uVelBar)
      call InterpolateElementToFace(interpolation,vVelocity,vVelBar)
      call InterpolateElementToFace(interpolation,wVelocity,wVelBar)
c
      call InterpolateElementToFace
     *                  (interpolation,uVelocityStar,uVelocityStarf)
      call InterpolateElementToFace
     *                  (interpolation,vVelocityStar,vVelocityStarf)
      call InterpolateElementToFace
     *                  (interpolation,wVelocityStar,wVelocityStarf)
c
      if(LFreeSurfaceFlow.and.LHarmonic) then

        interpolation='harmonic'
c
        call InterpolateGradientToFace(interpolation,Pressure,
     *           PressGradx,PressGrady,PressGradz,
     *               PressGradxBar,PressGradyBar,PressGradzBar)
c
      else
c
        call InterpolateGradientToFace(interpolation,Pressure,
     *           PressGradx,PressGrady,PressGradz,
     *               PressGradxBar,PressGradyBar,PressGradzBar)
c
      endif
c
      if(LSolveTurbulenceKineticEnergy) then
c
        interpolation='average'
c
        call InterpolateGradientToFace(interpolation,rhok,
     *        drhokdx,drhokdy,drhokdz,drhokdxBar,drhokdyBar,drhokdzBar)
c
        interpolation='averagecorrected'
c
        call InterpolateGradientToFace(interpolation,rhok,
     *           drhokdx,drhokdy,drhokdz,drhokdfx,drhokdfy,drhokdfz)
c
      endif
c
      if(LUnsteady) then
c
        allocate(uVelOldf(NIFaces))
        allocate(vVelOldf(NIFaces))
        allocate(wVelOldf(NIFaces))
        allocate(uVelOldOldf(NIFaces))
        allocate(vVelOldOldf(NIFaces))
        allocate(wVelOldOldf(NIFaces))
        allocate(rhoOldf(NIFaces))
        allocate(rhoOldOldf(NIFaces))
        allocate(DutOld(NumberOfElements))
        allocate(DvtOld(NumberOfElements))
        allocate(DwtOld(NumberOfElements))
        allocate(DutOldOld(NumberOfElements))
        allocate(DvtOldOld(NumberOfElements))
        allocate(DwtOldOld(NumberOfElements))
        allocate(DVtfO(NIFaces))
        allocate(DVtfOO(NIFaces))
        allocate(DutOldT(NumberOfElements))
        allocate(DutOldOldT(NumberOfElements))
c
        allocate(drhoxOld(NumberOfElements))
        allocate(drhoyOld(NumberOfElements))
        allocate(drhozOld(NumberOfElements))
        allocate(BdrhoxOld(NumberOfBCSets,NBFacesMax))
        allocate(BdrhoyOld(NumberOfBCSets,NBFacesMax))
        allocate(BdrhozOld(NumberOfBCSets,NBFacesMax))
        allocate(drhoxOldOld(NumberOfElements))
        allocate(drhoyOldOld(NumberOfElements))
        allocate(drhozOldOld(NumberOfElements))
        allocate(BdrhoxOldOld(NumberOfBCSets,NBFacesMax))
        allocate(BdrhoyOldOld(NumberOfBCSets,NBFacesMax))
        allocate(BdrhozOldOld(NumberOfBCSets,NBFacesMax))
c
        interpolation='average'
c
        call InterpolateElementToFace
     *                  (interpolation,uVelocityOld,uVelOldf)
        call InterpolateElementToFace
     *                  (interpolation,vVelocityOld,vVelOldf)
        call InterpolateElementToFace
     *                  (interpolation,wVelocityOld,wVelOldf)
        call InterpolateElementToFace
     *                  (interpolation,uVelocityOldOld,uVelOldOldf)
        call InterpolateElementToFace
     *                  (interpolation,vVelocityOldOld,vVelOldOldf)
        call InterpolateElementToFace
     *                  (interpolation,wVelocityOldOld,wVelOldOldf)
c
        if(LFreeSurfaceFlow.and.LHarmonic) then 
c
          Interpolation='harmonic'
          call InterpolateElementToFace
     *            (Interpolation,DensityOld,rhoOldf)
          call InterpolateElementToFace
     *            (Interpolation,DensityOldOld,rhoOldOldf)
c
        elseif(LNVFDensity.or.LTVDDensity) then
c
          Variable='DensO'
          call Gradient(Variable,MethodCalcGradientDensity,
     *      DensityOld,drhoxOld,drhoyOld,drhozOld,BDensityOld,
     *            BdrhoxOld,BdrhoyOld,BdrhozOld,nIterGradientDensity,
     *             LimitGradientDensity,LimitGradientDensityMethod)
c
          Variable='DensOO'
          call Gradient(Variable,MethodCalcGradientDensity,
     *      DensityOldOld,drhoxOldOld,drhoyOldOld,drhozOldOld,
     *            BDensityOldOld,BdrhoxOldOld,BdrhoyOldOld,BdrhozOldOld,
     *             nIterGradientDensity,LimitGradientDensity,
     *                              LimitGradientDensityMethod)
c
          if(nIter.ge.nIterStartApplyingHR) then
c
            call InterpolateToFaceUsingHRscheme(ConvectionSchemeDensity,
     *        BleedDensity,LNVFDensity,LTVDDensity,DensityOld,
     *                            drhoxOld,drhoyOld,drhozOld,rhoOldf)
c
            call InterpolateToFaceUsingHRscheme(ConvectionSchemeDensity,
     *       BleedDensity,LNVFDensity,LTVDDensity,DensityOldOld,
     *                   drhoxOldOld,drhoyOldOld,drhozOldOld,rhoOldOldf)
c
          else
c
            ConvectionScheme1='upwind'
c
            call InterpolateToFaceUsingHRscheme(ConvectionScheme1,
     *        BleedDensity,LNVFDensity,LTVDDensity,DensityOld,
     *                            drhoxOld,drhoyOld,drhozOld,rhoOldf)
c
            call InterpolateToFaceUsingHRscheme(ConvectionScheme1,
     *       BleedDensity,LNVFDensity,LTVDDensity,DensityOldOld,
     *                   drhoxOldOld,drhoyOldOld,drhozOldOld,rhoOldOldf)
c
          endif
c
        else
c
          Interpolation='average'
          call InterpolateElementToFace
     *            (Interpolation,DensityOld,rhoOldf)
          call InterpolateElementToFace
     *            (Interpolation,DensityOldOld,rhoOldOldf)
c
        endif
c
        deallocate(drhoxOld)
        deallocate(drhoyOld)
        deallocate(drhozOld)
        deallocate(BdrhoxOld)
        deallocate(BdrhoyOld)
        deallocate(BdrhozOld)
        deallocate(drhoxOldOld)
        deallocate(drhoyOldOld)
        deallocate(drhozOldOld)
        deallocate(BdrhoxOldOld)
        deallocate(BdrhoyOldOld)
        deallocate(BdrhozOldOld)
c
        do i=1,NumberOfElements
c         
          value=acold(i)/Volume(i)
          DutOld(i)=Du1Velocity(i)*value      
          DvtOld(i)=Dv1Velocity(i)*value        
          DwtOld(i)=Dw1Velocity(i)*value        
c
          value=acoldold(i)/Volume(i)
          DutOldOld(i)=Du1Velocity(i)*value      
          DvtOldOld(i)=Dv1Velocity(i)*value        
          DwtOldOld(i)=Dw1Velocity(i)*value        
c
        enddo
c
        do i=1,NumberOfElements
c
          DutOldT(i)=(DutOld(i)+DvtOld(i)+DwtOld(i))/3.
          DutOldOldT(i)=(DutOldOld(i)+DvtOldOld(i)+DwtOldOld(i))/3.
c
        enddo
c
        Interpolation='average'
c
        call InterpolateElementToFace
     *                  (interpolation,DutOld,DVtfO)
        call InterpolateElementToFace
     *                  (interpolation,DutOldOld,DVtfOO)
c
        deallocate(DutOld)
        deallocate(DvtOld)
        deallocate(DwtOld)
        deallocate(DutOldOld)
        deallocate(DvtOldOld)
        deallocate(DwtOldOld)
c
      endif
c
      if(LFalseTransientMomentum) then
c
        allocate(rhoStarf(NIFaces))
        allocate(DutStar(NumberOfElements))
        allocate(DvtStar(NumberOfElements))
        allocate(DwtStar(NumberOfElements))
        allocate(DVtfStar(NIFaces))
        allocate(DutStarT(NumberOfElements))
c
        allocate(drhoxStar(NumberOfElements))
        allocate(drhoyStar(NumberOfElements))
        allocate(drhozStar(NumberOfElements))
        allocate(BdrhoxStar(NumberOfBCSets,NBFacesMax))
        allocate(BdrhoyStar(NumberOfBCSets,NBFacesMax))
        allocate(BdrhozStar(NumberOfBCSets,NBFacesMax))
c
        if(LFreeSurfaceFlow.and.LHarmonic) then 
c
          Interpolation='harmonic'
          call InterpolateElementToFace
     *            (Interpolation,DensityStar,rhoStarf)
c
        elseif(LNVFDensity.or.LTVDDensity) then
c
          Variable='DensO'
          call Gradient(Variable,MethodCalcGradientDensity,
     *      DensityStar,drhoxStar,drhoyStar,drhozStar,BDensityStar,
     *            BdrhoxStar,BdrhoyStar,BdrhozStar,nIterGradientDensity,
     *             LimitGradientDensity,LimitGradientDensityMethod)
c
          if(nIter.ge.nIterStartApplyingHR) then
c
            call InterpolateToFaceUsingHRscheme(ConvectionSchemeDensity,
     *        BleedDensity,LNVFDensity,LTVDDensity,DensityStar,
     *                          drhoxStar,drhoyStar,drhozStar,rhoStarf)
c
          else
c
            ConvectionScheme1='upwind'
c
            call InterpolateToFaceUsingHRscheme(ConvectionScheme1,
     *        BleedDensity,LNVFDensity,LTVDDensity,DensityStar,
     *                         drhoxStar,drhoyStar,drhozStar,rhoStarf)
c
          endif
c
        else
c
          Interpolation='average'
          call InterpolateElementToFace
     *            (Interpolation,DensityStar,rhoStarf)
c
        endif
c
        deallocate(drhoxStar)
        deallocate(drhoyStar)
        deallocate(drhozStar)
        deallocate(BdrhoxStar)
        deallocate(BdrhoyStar)
        deallocate(BdrhozStar)
c
        do i=1,NumberOfElements
c         
          value=acStar(i)/Volume(i)
          DutStar(i)=Du1Velocity(i)*value      
          DvtStar(i)=Dv1Velocity(i)*value        
          DwtStar(i)=Dw1Velocity(i)*value        
c
        enddo
c
        do i=1,NumberOfElements
c
          DutStarT(i)=(DutStar(i)+DvtStar(i)+DwtStar(i))/3.
c
        enddo
c
        interpolation='average'
        call InterpolateElementToFace
     *                  (interpolation,DutStarT,DVtfStar)
c
        deallocate(DutStar)
        deallocate(DvtStar)
        deallocate(DwtStar)
c
      endif
c
c--- Compute the mass flow rate at interior faces using Rhie-Chow intepolation
c
      do k=1,NIFaces
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        rhof=Densityf(k)
        sfx=FaceAreax(k)
        sfy=FaceAreay(k)
        sfz=FaceAreaz(k)
        ex=DistanceCFux(k)
        ey=DistanceCFuy(k)
        ez=DistanceCFuz(k)
        cf=DistanceCF(k)
c
c---  mdot bar term
c
        mdotBarf=uVelBar(k)*sfx+vVelBar(k)*sfy+wVelBar(k)*sfz
        FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
        DuSf=Du1f(k)*sfx
        DvSf=Dv1f(k)*sfy
        DwSf=Dw1f(k)*sfz
c
        if(MethodDecomposeS.eq.1) then
c
          dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
          DuEf=dotproduct*ex
          DvEf=dotproduct*ey
          DwEf=dotproduct*ez
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
        elseif(MethodDecomposeS.eq.2) then
c          
          Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
          DuEf=ex*Magnitude
          DvEf=ey*Magnitude
          DwEf=ez*Magnitude
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
        elseif(MethodDecomposeS.eq.3) then
c
          Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
          eDuSf=DuSf/Magnitude
          eDvSf=DvSf/Magnitude
          eDwSf=DwSf/Magnitude
c
          dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
          DuEf=Magnitude*ex/dotproduct
          DvEf=Magnitude*ey/dotproduct
          DwEf=Magnitude*ez/dotproduct
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
        endif
c          
        DuTf=DuSf-DuEf         
        DvTf=DvSf-DvEf
        DwTf=DwSf-DwEf
c
        FluxCfLocal=rhof*geoDiff
        FluxFfLocal=-rhof*geoDiff
        FluxVfLocal=FluxVfLocal-rhof*(PressGradfx(k)*DuTf+
     *                   PressGradfy(k)*DvTf+PressGradfz(k)*DwTf)
c
        FluxVfLocal=FluxVfLocal+rhof*(PressGradxBar(k)*DuSf+
     *               PressGradyBar(k)*DvSf+PressGradzBar(k)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
        if(LSolveTurbulenceKineticEnergy) then
c
          FluxVfLocal=FluxVfLocal-rhof*(drhokdfx(k)*DuTf+
     *                   drhokdfy(k)*DvTf+drhokdfz(k)*DwTf)
          FluxVfLocal=FluxVfLocal+rhof*(drhokdxBar(k)*DuSf+
     *               drhokdyBar(k)*DvSf+drhokdzBar(k)*DwSf)
          FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i)+FluxFfLocal*rhok(j)
c
        endif
c
c--- Transient term correction
c
        if(LUnsteady) then
c
          mdotOldBar=rhoOldf(k)*(uVelOldf(k)*sfx+
     *                           vVelOldf(k)*sfy+wVelOldf(k)*sfz)   
          mdotOldOldBar=rhoOldOldf(k)*(uVelOldOldf(k)*sfx+
     *                      vVelOldOldf(k)*sfy+wVelOldOldf(k)*sfz)
          FluxVfLocal=FluxVfLocal+
     *                      DVtfO(k)*(mdotOld(k)-mdotOldBar)+
     *                          DVtfOO(k)*(mdotOldOld(k)-mdotOldOldBar)
c
        endif
c
c--- False transient term correction
c
        if(LFalseTransientMomentum) then
c
          mdotOldBar=rhoStarf(k)*(uVelocityStarf(k)*sfx+
     *               vVelocityStarf(k)*sfy+wVelocityStarf(k)*sfz)   
          FluxVfLocal=FluxVfLocal+DVtfStar(k)*(mdot(k)-mdotOldBar)
c
        endif
c
c--- underrelaxation correction
c
        if(LRelaxMomentum) then
c
          mdotstar=rhof*(uVelocityStarf(k)*sfx+
     *               vVelocityStarf(k)*sfy+wVelocityStarf(k)*sfz) 
          FluxVfLocal=FluxVfLocal+(1.-urfMomentum)*(mdot(k)-mdotstar)
c
        endif
c
c--- Surface tension correction
c
        if(LFreeSurfaceFlow.and.LSurfaceTension) then
c
          if(LHarmonic) then
c
            SurfaceTensionxBar=rhof*(
     *        GFactorCF(k)*Curvature(i)*rFieldGradxT(i)/Density(i)+
     *        (1.-GFactorCF(k))*Curvature(j)*rFieldGradxT(j)/Density(j))
            SurfaceTensionyBar=rhof*(
     *        GFactorCF(k)*Curvature(i)*rFieldGradyT(i)/Density(i)+
     *        (1.-GFactorCF(k))*Curvature(j)*rFieldGradyT(j)/Density(j))
            SurfaceTensionzBar=rhof*(
     *        GFactorCF(k)*Curvature(i)*rFieldGradzT(i)/Density(i)+
     *        (1.-GFactorCF(k))*Curvature(j)*rFieldGradzT(j)/Density(j))
c
          else
c
            SurfaceTensionxBar=
     *          GFactorCF(k)*Curvature(i)*rFieldGradxT(i)+
     *               (1.-GFactorCF(k))*Curvature(j)*rFieldGradxT(j)
            SurfaceTensionyBar=
     *          GFactorCF(k)*Curvature(i)*rFieldGradyT(i)+
     *               (1.-GFactorCF(k))*Curvature(j)*rFieldGradyT(j)
            SurfaceTensionzBar=
     *          GFactorCF(k)*Curvature(i)*rFieldGradzT(i)+
     *               (1.-GFactorCF(k))*Curvature(j)*rFieldGradzT(j)
c
          endif
c
          CurvactureF=GFactorCF(k)*Curvature(i)+
     *                           (1.-GFactorCF(k))*Curvature(j)
          FluxVfLocal=FluxVfLocal+rhof*geoDiff*
     *            SurfaceTension*CurvactureF*(rFieldT(j)-rFieldT(i))
          FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
     *                (rFieldGradfxT(k)*DuTf+rFieldGradfyT(k)*DvTf+
     *                                         rFieldGradfzT(k)*DwTf)
          FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
     *        (SurfaceTensionxBar*DuSf+
     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
c
        endif
c
c--- Buoyancy correction
c
        if(LBuoyancy) then
c
          if(LFreeSurfaceFlow.and.LHarmonic) then
c
            BuoyancyxBar=rhof*(
     *         GFactorCF(k)*Buoyancyx(i)/Density(i)+
     *         (1.-GFactorCF(k))*Buoyancyx(j)/Density(j))
            BuoyancyyBar=rhof*(
     *         GFactorCF(k)*Buoyancyy(i)/Density(i)+
     *         (1.-GFactorCF(k))*Buoyancyy(j)/Density(j))
            BuoyancyzBar=rhof*(
     *         GFactorCF(k)*Buoyancyz(i)/Density(i)+
     *        (1.-GFactorCF(k))*Buoyancyz(j)/Density(j))
c
          else
c
            BuoyancyxBar=GFactorCF(k)*Buoyancyx(i)+
     *                   (1.-GFactorCF(k))*Buoyancyx(j)
            BuoyancyyBar=GFactorCF(k)*Buoyancyy(i)+
     *                   (1.-GFactorCF(k))*Buoyancyy(j)
            BuoyancyzBar=GFactorCF(k)*Buoyancyz(i)+
     *                    (1.-GFactorCF(k))*Buoyancyz(j)
c
          endif
c
!          BuoyancyBar=DistanceCFx(k)*BuoyancyxBar+
!     *         DistanceCFy(k)*BuoyancyyBar+DistanceCFz(k)*BuoyancyzBar
!          BuoyancyFace=Buoyancyfx(k)*DistanceCFx(k)+
!     *         Buoyancyfy(k)*DistanceCFy(k)+Buoyancyfz(k)*DistanceCFz(k)
!          FluxVfLocal=FluxVfLocal+
!     *                     rhof*geoDiff*(BuoyancyFace-BuoyancyBar)
c
          FluxVfLocal=FluxVfLocal+rhof*(
     *                 DuSf*(Buoyancyfx(k)-BuoyancyxBar)+
     *                    DvSf*(Buoyancyfy(k)-BuoyancyyBar)+
     *                       DwSf*(Buoyancyfz(k)-BuoyancyzBar))
c
        endif
c
        mdot(k)=FluxCfLocal*Pressure(i)+
     *                   FluxFfLocal*Pressure(j)+FluxVfLocal
c
      enddo
c
c----------------------------------------------------------------------
c--- Calculate mdot along periodic boundaries for in/compressible flows
c----------------------------------------------------------------------
c
      if(LRotationalPeriodicity) then

        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
            yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
            zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            rhof=BDensity(i2,i3)
            sfx=BFaceAreax(i2,i3)
            sfy=BFaceAreay(i2,i3)
            sfz=BFaceAreaz(i2,i3)
c
            DistCFx=xF1-xc(i1)
            DistCFy=yF1-yc(i1)
            DistCFz=zF1-zc(i1)
c
            DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
            DistCFux=DistCFx/DistCF
            DistCFuy=DistCFy/DistCF
            DistCFuz=DistCFz/DistCF
c
            ex=DistCFux
            ey=DistCFuy
            ez=DistCFuz
            cf=DistCF
c
c---  mdot bar term
c
            uF1=a1r(j2)*uVelocity(j1)+b1r(j2)*vVelocity(j1)+
     *                                         c1r(j2)*wVelocity(j1)
            vF1=a2r(j2)*uVelocity(j1)+b2r(j2)*vVelocity(j1)+
     *                                         c2r(j2)*wVelocity(j1)
            wF1=a3r(j2)*uVelocity(j1)+b3r(j2)*vVelocity(j1)+
     *                                         c3r(j2)*wVelocity(j1)
            uBar=GFactCF*uVelocity(i1)+(1.-GFactCF)*uF1
            vBar=GFactCF*vVelocity(i1)+(1.-GFactCF)*vF1
            wBar=GFactCF*wVelocity(i1)+(1.-GFactCF)*wF1
            mdotBarf=uBar*sfx+vBar*sfy+wBar*sfz
            FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
            DuSf=GFactCF*Du1Velocity(i1)+(1.-GFactCF)*Du1Velocity(j1)
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv1Velocity(i1)+(1.-GFactCF)*Dv1Velocity(j1)
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw1Velocity(i1)+(1.-GFactCF)*Dw1Velocity(j1)
            DwSf=DwSf*sfz
c
            if(MethodDecomposeS.eq.1) then
c
              DuEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ex
              DvEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ey
              DwEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ez
c          
              geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
            elseif(MethodDecomposeS.eq.2) then
c          
              DuEf=ex*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
              DvEf=ey*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
              DwEf=ez*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
c          
              geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
            elseif(MethodDecomposeS.eq.3) then
c
              Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
              eDuSf=DuSf/Magnitude
              eDvSf=DvSf/Magnitude
              eDwSf=DwSf/Magnitude
c
              DuEf=Magnitude*ex/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
              DvEf=Magnitude*ey/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
              DwEf=Magnitude*ez/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
c          
              geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
            endif
c          
            DuTf=DuSf-DuEf         
            DvTf=DvSf-DvEf
            DwTf=DwSf-DwEf
c
            dfdxBar=BPressGradx(i2,i3)
            dfdyBar=BPressGrady(i2,i3)
            dfdzBar=BPressGradz(i2,i3)
c
            correction=(Pressure(j1)-Pressure(i1))/DistCF
            correction=correction-(dfdxBar*DistCFux+
     *                         dfdyBar*DistCFuy+dfdzBar*DistCFuz)
c
            dfdxf=dfdxBar+correction*DistCFux
            dfdyf=dfdyBar+correction*DistCFuy
            dfdzf=dfdzBar+correction*DistCFuz
c
            FluxCfLocal=rhof*geoDiff
            FluxFfLocal=-rhof*geoDiff
            FluxVfLocal=FluxVfLocal-
     *                   rhof*(dfdxf*DuTf+dfdyf*DvTf+dfdzf*DwTf)
            FluxVfLocal=FluxVfLocal+
     *                  rhof*(dfdxBar*DuSf+dfdyBar*DvSf+dfdzBar*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
            if(LSolveTurbulenceKineticEnergy) then
c
              dfdxBar=Bdrhokdx(i2,i3)
              dfdyBar=Bdrhokdy(i2,i3)
              dfdzBar=Bdrhokdz(i2,i3)
c
              correction=(rhok(j1)-rhok(i1))/DistCF
              correction=correction-(dfdxBar*DistCFux+
     *                         dfdyBar*DistCFuy+dfdzBar*DistCFuz)
c
              dfdxf=dfdxBar+correction*DistCFux
              dfdyf=dfdyBar+correction*DistCFuy
              dfdzf=dfdzBar+correction*DistCFuz
c
              FluxVfLocal=FluxVfLocal-
     *                 rhof*(dfdxf*DuTf+dfdyf*DvTf+dfdzf*DwTf)
              FluxVfLocal=FluxVfLocal+
     *                  rhof*(dfdxBar*DuSf+dfdyBar*DvSf+dfdzBar*DwSf)
              FluxVfLocal=FluxVfLocal+
     *                FluxCfLocal*rhok(i1)+FluxFfLocal*rhok(j1)
c
            endif
c
c--- Transient term correction
c
            if(LUnsteady) then
c
              rhoOf=GFactCF*DensityOld(i1)+(1.-GFactCF)*DensityOld(j1)
              uOf=GFactCF*uVelocityOld(i1)+(1.-GFactCF)*uVelocityOld(j1)
              vOf=GFactCF*vVelocityOld(i1)+(1.-GFactCF)*vVelocityOld(j1)
              wOf=GFactCF*wVelocityOld(i1)+(1.-GFactCF)*wVelocityOld(j1)
c
              rhoOOf=GFactCF*DensityOldOld(i1)+
     *                                 (1.-GFactCF)*DensityOldOld(j1)
              uOOf=GFactCF*uVelocityOldOld(i1)+
     *                                 (1.-GFactCF)*uVelocityOldOld(j1)
              vOOf=GFactCF*vVelocityOldOld(i1)+
     *                                 (1.-GFactCF)*vVelocityOldOld(j1)
              wOOf=GFactCF*wVelocityOldOld(i1)+
     *                                 (1.-GFactCF)*wVelocityOldOld(j1)
              mdotOldBar=rhoOf*(uOf*sfx+vOf*sfy+wOf*sfz)   
              mdotOldOldBar=rhoOOf*(uOOf*sfx+vOOf*sfy+wOOf*sfz)
              DVfO=GFactCF*DutOldT(i1)+(1.-GFactCF)*DutOldT(j1)
              DVfOO=GFactCF*DutOldOldT(i1)+(1.-GFactCF)*DutOldOldT(j1)
              FluxVfLocal=FluxVfLocal+
     *                DVfO*(BmdotOld(i2,i3)-mdotOldBar)+
     *                    DVfOO*(BmdotOldOld(i2,i3)-mdotOldOldBar)
c
            endif
c
c--- False transient term correction
c
            if(LFalseTransientMomentum) then
c
              rhoOf=GFactCF*DensityStar(i1)+(1.-GFactCF)*DensityStar(j1)
              uStarf=GFactCF*uVelocityStarf(i1)+
     *                           (1.-GFactCF)*uVelocityStarf(j1)
              vStarf=GFactCF*vVelocityStarf(i1)+
     *                           (1.-GFactCF)*vVelocityStarf(j1)
              wStarf=GFactCF*wVelocityStarf(i1)+
     *                           (1.-GFactCF)*wVelocityStarf(j1)
              mdotOldBar=rhoOf*(uStarf*sfx+vStarf*sfy+wStarf*sfz)   
              DVfO=GFactCF*DutStarT(i1)+(1.-GFactCF)*DutStarT(j1)
              FluxVfLocal=FluxVfLocal+DVfO*(Bmdot(i2,i3)-mdotOldBar)
c
            endif
c
c--- underrelaxation correction
c
            if(LRelaxMomentum) then
c
              uF1star=a1r(j2)*uVelocitystar(j1)+
     *                        b1r(j2)*vVelocitystar(j1)+
     *                            c1r(j2)*wVelocitystar(j1)
              vF1star=a2r(j2)*uVelocitystar(j1)+
     *                        b2r(j2)*vVelocitystar(j1)+
     *                            c2r(j2)*wVelocitystar(j1)
              wF1star=a3r(j2)*uVelocitystar(j1)+
     *                        b3r(j2)*vVelocitystar(j1)+
     *                            c3r(j2)*wVelocitystar(j1)
              uStarf=GFactCF*uVelocityStar(i1)+
     *                           (1.-GFactCF)*uF1star
              vStarf=GFactCF*vVelocityStar(i1)+
     *                           (1.-GFactCF)*vF1star
              wStarf=GFactCF*wVelocityStar(i1)+
     *                           (1.-GFactCF)*wF1star
              mdotstar=rhof*(uStarf*sfx+vStarf*sfy+wStarf*sfz) 
              FluxVfLocal=FluxVfLocal+
     *                     (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
            endif
c
c--- Surface tension correction
c
            if(LFreeSurfaceFlow.and.LSurfaceTension) then
c
              if(LHarmonic) then
c
                SurfaceTensionxBar=rhof*(
     *            GFactCF*Curvature(i1)*rFieldGradxT(i1)/Density(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*
     *                                    rFieldGradxT(j1)/Density(j1))
                SurfaceTensionyBar=rhof*(
     *           GFactCF*Curvature(i1)*rFieldGradyT(i1)/Density(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*
     *                                    rFieldGradyT(j1)/Density(j1))
                SurfaceTensionzBar=rhof*(
     *            GFactCF*Curvature(i1)*rFieldGradzT(i1)/Density(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*
     *                                    rFieldGradzT(j1)/Density(j1))
c
              else
c
                SurfaceTensionxBar=
     *            GFactCF*Curvature(i1)*rFieldGradxT(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*rFieldGradxT(j1)
                SurfaceTensionyBar=
     *           GFactCF*Curvature(i1)*rFieldGradyT(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*rFieldGradyT(j1)
                SurfaceTensionzBar=
     *            GFactCF*Curvature(i1)*rFieldGradzT(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*rFieldGradzT(j1)
c
              endif
c
              CurvactureF=GFactCF*Curvature(i1)+
     *                                     (1.-GFactCF)*Curvature(j1)
              FluxVfLocal=FluxVfLocal+rhof*geoDiff*
     *            SurfaceTension*CurvactureF*(rFieldT(j1)-rFieldT(i1))
              FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
     *                                      BrFieldGradzT(i2,i3)*DwTf)
              FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
     *        (SurfaceTensionxBar*DuSf+
     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
c
            endif
c
c--- Buoyancy correction
c
            if(LBuoyancy) then
c
              if(LFreeSurfaceFlow.and.LHarmonic) then
c
                BuoyancyxBar=rhof*(GFactCF*Buoyancyx(i1)/Density(i1)+
     *                          (1.-GFactCF)*Buoyancyx(j1)/Density(j1))
                BuoyancyyBar=rhof*(GFactCF*Buoyancyy(i1)/Density(i1)+
     *                          (1.-GFactCF)*Buoyancyy(j1)/Density(j1))
                BuoyancyzBar=rhof*(GFactCF*Buoyancyz(i1)/Density(i1)+
     *                          (1.-GFactCF)*Buoyancyz(j1)/Density(j1))
c
              else
c
                BuoyancyxBar=GFactCF*Buoyancyx(i1)+
     *                             (1.-GFactCF)*Buoyancyx(j1)      
                BuoyancyyBar=GFactCF*Buoyancyy(i1)+
     *                             (1.-GFactCF)*Buoyancyy(j1)      
                BuoyancyzBar=GFactCF*Buoyancyz(i1)+
     *                             (1.-GFactCF)*Buoyancyz(j1)
c
              endif      
c
!              BuoyancyBar=DistCFx*BuoyancyxBar+
!     *                   DistCFy*BuoyancyyBar+DistCFz*BuoyancyzBar
!              BuoyancyFace=BBuoyancyx(i2,i3)*DistCFx+
!     *              BBuoyancyy(i2,i3)*DistCFy+BBuoyancyz(i2,i3)*DistCFz
!              FluxVfLocal=FluxVfLocal+
!     *                     rhof*geoDiff*(BuoyancyFace-BuoyancyBar)
c
               FluxVfLocal=FluxVfLocal+rhof*(
     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
c
            endif
c
            Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*Pressure(j1)+FluxVfLocal
c
            Bmdot(j2,j3)=-Bmdot(i2,i3)
c
          endif
c
        enddo
c
      elseif(LTranslationalPeriodicity) then
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=xc(j1)+xTranslation(j2)
            yF1=yc(j1)+yTranslation(j2)
            zF1=zc(j1)+zTranslation(j2)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            rhof=BDensity(i2,i3)
            sfx=BFaceAreax(i2,i3)
            sfy=BFaceAreay(i2,i3)
            sfz=BFaceAreaz(i2,i3)
c
            DistCFx=xF1-xc(i1)
            DistCFy=yF1-yc(i1)
            DistCFz=zF1-zc(i1)
c
            DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
            DistCFux=DistCFx/DistCF
            DistCFuy=DistCFy/DistCF
            DistCFuz=DistCFz/DistCF
c
            ex=DistCFux
            ey=DistCFuy
            ez=DistCFuz
            cf=DistCF
c
c---  mdot bar term
c
            uF1=a1r(j2)*uVelocity(j1)+b1r(j2)*vVelocity(j1)+
     *                                          c1r(j2)*wVelocity(j1)
            vF1=a2r(j2)*uVelocity(j1)+b2r(j2)*vVelocity(j1)+
     *                                          c2r(j2)*wVelocity(j1)
            wF1=a3r(j2)*uVelocity(j1)+b3r(j2)*vVelocity(j1)+
     *                                          c3r(j2)*wVelocity(j1)
            uBar=GFactCF*uVelocity(i1)+(1.-GFactCF)*uF1
            vBar=GFactCF*vVelocity(i1)+(1.-GFactCF)*vF1
            wBar=GFactCF*wVelocity(i1)+(1.-GFactCF)*wF1
            mdotBarf=uBar*sfx+vBar*sfy+wBar*sfz
            FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
            DuSf=GFactCF*Du1Velocity(i1)+(1.-GFactCF)*Du1Velocity(j1)
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv1Velocity(i1)+(1.-GFactCF)*Dv1Velocity(j1)
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw1Velocity(i1)+(1.-GFactCF)*Dw1Velocity(j1)
            DwSf=DwSf*sfz
c
            if(MethodDecomposeS.eq.1) then
c
              DuEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ex
              DvEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ey
              DwEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ez
c          
              geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
            elseif(MethodDecomposeS.eq.2) then
c          
              DuEf=ex*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
              DvEf=ey*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
              DwEf=ez*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
c          
              geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
            elseif(MethodDecomposeS.eq.3) then
c
              Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
              eDuSf=DuSf/Magnitude
              eDvSf=DvSf/Magnitude
              eDwSf=DwSf/Magnitude
c
              DuEf=Magnitude*ex/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
              DvEf=Magnitude*ey/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
              DwEf=Magnitude*ez/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
c          
              geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
            endif
c          
            DuTf=DuSf-DuEf         
            DvTf=DvSf-DvEf
            DwTf=DwSf-DwEf
c
            dfdxBar=BPressGradx(i2,i3)
            dfdyBar=BPressGrady(i2,i3)
            dfdzBar=BPressGradz(i2,i3)
c
            correction=(Pressure(j1)-Pressure(i1))/DistCF
            correction=correction-
     *              (dfdxBar*DistCFux+dfdyBar*DistCFuy+dfdzBar*DistCFuz)
c
            dfdxf=dfdxBar+correction*DistCFux
            dfdyf=dfdyBar+correction*DistCFuy
            dfdzf=dfdzBar+correction*DistCFuz
c
            FluxCfLocal=rhof*geoDiff
            FluxFfLocal=-rhof*geoDiff
            FluxVfLocal=FluxVfLocal-
     *                        rhof*(dfdxf*DuTf+dfdyf*DvTf+dfdzf*DwTf)
            FluxVfLocal=FluxVfLocal+
     *                  rhof*(dfdxBar*DuSf+dfdyBar*DvSf+dfdzBar*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
            if(LSolveTurbulenceKineticEnergy) then
c
              dfdxBar=Bdrhokdx(i2,i3)
              dfdyBar=Bdrhokdy(i2,i3)
              dfdzBar=Bdrhokdz(i2,i3)
c
              correction=(rhok(j1)-rhok(i1))/DistCF
              correction=correction-(dfdxBar*DistCFux+
     *                         dfdyBar*DistCFuy+dfdzBar*DistCFuz)
c
              dfdxf=dfdxBar+correction*DistCFux
              dfdyf=dfdyBar+correction*DistCFuy
              dfdzf=dfdzBar+correction*DistCFuz
c
              FluxVfLocal=FluxVfLocal-
     *                 rhof*(dfdxf*DuTf+dfdyf*DvTf+dfdzf*DwTf)
              FluxVfLocal=FluxVfLocal+
     *                  rhof*(dfdxBar*DuSf+dfdyBar*DvSf+dfdzBar*DwSf)
              FluxVfLocal=FluxVfLocal+
     *                FluxCfLocal*rhok(i1)+FluxFfLocal*rhok(j1)
c
            endif
c
c--- Transient term correction
c
            if(LUnsteady) then
c
              rhoOf=GFactCF*DensityOld(i1)+(1.-GFactCF)*DensityOld(j1)
              uOf=GFactCF*uVelocityOld(i1)+(1.-GFactCF)*uVelocityOld(j1)
              vOf=GFactCF*vVelocityOld(i1)+(1.-GFactCF)*vVelocityOld(j1)
              wOf=GFactCF*wVelocityOld(i1)+(1.-GFactCF)*wVelocityOld(j1)
c
              rhoOOf=GFactCF*DensityOldOld(i1)+
     *                                 (1.-GFactCF)*DensityOldOld(j1)
              uOOf=GFactCF*uVelocityOldOld(i1)+
     *                                 (1.-GFactCF)*uVelocityOldOld(j1)
              vOOf=GFactCF*vVelocityOldOld(i1)+
     *                                 (1.-GFactCF)*vVelocityOldOld(j1)
              wOOf=GFactCF*wVelocityOldOld(i1)+
     *                                 (1.-GFactCF)*wVelocityOldOld(j1)
          
              mdotOldBar=rhoOf*(uOf*sfx+vOf*sfy+wOf*sfz)   
              mdotOldOldBar=rhoOOf*(uOOf*sfx+vOOf*sfy+wOOf*sfz)
              DVfO=GFactCF*DutOldT(i1)+(1.-GFactCF)*DutOldT(j1)
              DVfOO=GFactCF*DutOldOldT(i1)+(1.-GFactCF)*DutOldOldT(j1)
              FluxVfLocal=FluxVfLocal+
     *                DVfO*(BmdotOld(i2,i3)-mdotOldBar)+
     *                    DVfOO*(BmdotOldOld(i2,i3)-mdotOldOldBar)
c
            endif
c
c--- False transient term correction
c
            if(LFalseTransientMomentum) then
c
              rhoOf=GFactCF*DensityStar(i1)+(1.-GFactCF)*DensityStar(j1)
              uStarf=GFactCF*uVelocityStarf(i1)+
     *                           (1.-GFactCF)*uVelocityStarf(j1)
              vStarf=GFactCF*vVelocityStarf(i1)+
     *                           (1.-GFactCF)*vVelocityStarf(j1)
              wStarf=GFactCF*wVelocityStarf(i1)+
     *                           (1.-GFactCF)*wVelocityStarf(j1)
              mdotOldBar=rhoOf*(uStarf*sfx+vStarf*sfy+wStarf*sfz)   
              DVfO=GFactCF*DutStarT(i1)+(1.-GFactCF)*DutStarT(j1)
              FluxVfLocal=FluxVfLocal+DVfO*(Bmdot(i2,i3)-mdotOldBar)
c
            endif
c
c--- underrelaxation correction
c
            if(LRelaxMomentum) then
c
              uF1star=a1r(j2)*uVelocitystar(j1)+
     *                           b1r(j2)*vVelocitystar(j1)+
     *                              c1r(j2)*wVelocitystar(j1)
              vF1star=a2r(j2)*uVelocitystar(j1)+
     *                           b2r(j2)*vVelocitystar(j1)+
     *                              c2r(j2)*wVelocitystar(j1)
              wF1star=a3r(j2)*uVelocitystar(j1)+
     *                           b3r(j2)*vVelocitystar(j1)+
     *                              c3r(j2)*wVelocitystar(j1)
              uStarf=GFactCF*uVelocityStar(i1)+
     *                           (1.-GFactCF)*uF1star
              vStarf=GFactCF*vVelocityStar(i1)+
     *                           (1.-GFactCF)*vF1star
              wStarf=GFactCF*wVelocityStar(i1)+
     *                           (1.-GFactCF)*wF1star
              mdotstar=rhof*(uStarf*sfx+vStarf*sfy+wStarf*sfz) 
              FluxVfLocal=FluxVfLocal+
     *                     (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
            endif
c
c--- Surface tension correction
c
            if(LFreeSurfaceFlow.and.LSurfaceTension) then
c
              if(LHarmonic) then
c
                SurfaceTensionxBar=rhof*(
     *            GFactCF*Curvature(i1)*rFieldGradxT(i1)/Density(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*
     *                                    rFieldGradxT(j1)/Density(j1))
                SurfaceTensionyBar=rhof*(
     *           GFactCF*Curvature(i1)*rFieldGradyT(i1)/Density(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*
     *                                    rFieldGradyT(j1)/Density(j1))
                SurfaceTensionzBar=rhof*(
     *            GFactCF*Curvature(i1)*rFieldGradzT(i1)/Density(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*
     *                                    rFieldGradzT(j1)/Density(j1))
c
              else
c
                SurfaceTensionxBar=
     *            GFactCF*Curvature(i1)*rFieldGradxT(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*rFieldGradxT(j1)
                SurfaceTensionyBar=
     *           GFactCF*Curvature(i1)*rFieldGradyT(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*rFieldGradyT(j1)
                SurfaceTensionzBar=
     *            GFactCF*Curvature(i1)*rFieldGradzT(i1)+
     *                      (1.-GFactCF)*Curvature(j1)*rFieldGradzT(j1)
c
              endif
c
              CurvactureF=GFactCF*Curvature(i1)+
     *                                     (1.-GFactCF)*Curvature(j1)
              FluxVfLocal=FluxVfLocal+rhof*geoDiff*
     *            SurfaceTension*CurvactureF*(rFieldT(j1)-rFieldT(i1))
              FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
     *                                      BrFieldGradzT(i2,i3)*DwTf)
              FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
     *        (SurfaceTensionxBar*DuSf+
     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
c
            endif
c
c--- Buoyancy correction
c
            if(LBuoyancy.and.LHarmonic) then
c
              if(LFreeSurfaceFlow) then
c
                BuoyancyxBar=rhof*(GFactCF*Buoyancyx(i1)/Density(i1)+
     *                          (1.-GFactCF)*Buoyancyx(j1)/Density(j1))
                BuoyancyyBar=rhof*(GFactCF*Buoyancyy(i1)/Density(i1)+
     *                          (1.-GFactCF)*Buoyancyy(j1)/Density(j1))
                BuoyancyzBar=rhof*(GFactCF*Buoyancyz(i1)/Density(i1)+
     *                          (1.-GFactCF)*Buoyancyz(j1)/Density(j1))
c
              else
c
                BuoyancyxBar=GFactCF*Buoyancyx(i1)+
     *                             (1.-GFactCF)*Buoyancyx(j1)      
                BuoyancyyBar=GFactCF*Buoyancyy(i1)+
     *                             (1.-GFactCF)*Buoyancyy(j1)      
                BuoyancyzBar=GFactCF*Buoyancyz(i1)+
     *                             (1.-GFactCF)*Buoyancyz(j1)
c
              endif      
c
!              BuoyancyBar=DistCFx*BuoyancyxBar+
!     *                     DistCFy*BuoyancyyBar+DistCFz*BuoyancyzBar
!              BuoyancyFace=BBuoyancyx(i2,i3)*DistCFx+
!     *                           BBuoyancyy(i2,i3)*DistCFy+
!     *                                   BBuoyancyz(i2,i3)*DistCFz
!              FluxVfLocal=FluxVfLocal+
!     *                     rhof*geoDiff*(BuoyancyFace-BuoyancyBar)
c
              FluxVfLocal=FluxVfLocal+rhof*(
     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
c
            endif
c
            Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*Pressure(j1)+FluxVfLocal
c
            Bmdot(j2,j3)=-Bmdot(i2,i3)
c
          endif
c
        enddo
c
      endif
c
c----------------------------------------------------------------------
c--- Modify along inlet and oulet boundaries for incompressible and compressible flows
c----------------------------------------------------------------------
c
      if(.not.Lcompressible) then
c
c----------------------------------------------------------------------
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          Bmdot(i2,i3)=BDensity(i2,i3)*(BuVelocity(i2,i3)*sfx+
     *                  BvVelocity(i2,i3)*sfy+BwVelocity(i2,i3)*sfz)
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
          BuVelocity(i2,i3)=uVelocity(i1)+
     *           (BuVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BuVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BuVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BvVelocity(i2,i3)=vVelocity(i1)+
     *           (BvVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BvVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BvVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BwVelocity(i2,i3)=wVelocity(i1)+
     *           (BwVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BwVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BwVelGradz(i2,i3)*BDistanceCFz(i2,i3))
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *             BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
          BuVelocity(i2,i3)=uVelocity(i1)+
     *           (BuVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BuVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BuVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BvVelocity(i2,i3)=vVelocity(i1)+
     *           (BvVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BvVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BvVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BwVelocity(i2,i3)=wVelocity(i1)+
     *           (BwVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BwVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BwVelGradz(i2,i3)*BDistanceCFz(i2,i3))
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *            BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
          BuVelocity(i2,i3)=uVelocity(i1)+
     *           (BuVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BuVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BuVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BvVelocity(i2,i3)=vVelocity(i1)+
     *           (BvVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BvVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BvVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BwVelocity(i2,i3)=wVelocity(i1)+
     *           (BwVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BwVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BwVelGradz(i2,i3)*BDistanceCFz(i2,i3))
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *            BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          velocity=dsqrt(uVelocity(i1)**2+
     *                    vVelocity(i1)**2+wVelocity(i1)**2)
          xdirection=uVelocity(i1)/(velocity+tiny)
          ydirection=vVelocity(i1)/(velocity+tiny)
          zdirection=wVelocity(i1)/(velocity+tiny)
c
          velMagnitude=Bmdot(i2,i3)/
     *         (rhof*(xdirection*sfx+ydirection*sfy+zdirection*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xdirection
	    BvVelocity(i2,i3)=velMagnitude*ydirection
	    BwVelocity(i2,i3)=velMagnitude*zdirection
c
        enddo
c----------------------------------------------------------------------
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
	    BuVelocity(i2,i3)=uVelocity(i1)
	    BvVelocity(i2,i3)=vVelocity(i1)
	    BwVelocity(i2,i3)=wVelocity(i1)
c
        enddo
c----------------------------------------------------------------------
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          Bmdot(i2,i3)=BDensity(i2,i3)*(BuVelocity(i2,i3)*sfx+
     *                    BvVelocity(i2,i3)*sfy+BwVelocity(i2,i3)*sfz)
c
        enddo
c----------------------------------------------------------------------
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          rhof=BDensity(i2,i3)
c
          velMagnitude=Bmdot(i2,i3)/(rhof*(xVeldirection(i2,i3)*sfx+
     *               yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*magnitude
            DvEf=ey*magnitude
            DwEf=ez*magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *             BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
          velMagnitude=Bmdot(i2,i3)/(rhof*(xVeldirection(i2,i3)*sfx+
     *               yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *            BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
          velMagnitude=Bmdot(i2,i3)/(rhof*(xVeldirection(i2,i3)*sfx+
     *             yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
        enddo
c----------------------------------------------------------------------
c
      elseif(Lcompressible) then
c
c----------------------------------------------------------------------
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
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
          ratio=BSpecificHeat(i2,i3)/(BSpecificHeat(i2,i3)-RGas)
c
          vdotInfinity=uVelocityFarField(i2)*nx+
     *                vVelocityFarField(i2)*ny+wVelocityFarField(i2)*nz
          vdotin=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
          cinfinity=dsqrt(ratio*RGas*TemperatureFarField(i2))
          cin=dsqrt(ratio*RGas*Temperature(i1))
c
          vdotb=0.5*(vdotInfinity+vdotin)-(cinfinity-cin)/(ratio-1.)
          cb=-(vdotInfinity-vdotin)*(ratio-1.)/4.+0.5*(cinfinity+cin)
          Tb=(cb*cb)/(ratio*RGas)
          rhoInfinity=pressureFarField(i2)/
     *                           (RGas*TemperatureFarField(i2))
c
          if(Bmdot(i2,i3).lt.0.) then
c
            sbEntropy=cinfinity**2/(ratio*(rhoInfinity**(ratio-1)))
c
            uTangential=uVelocityFarField(i2)-
     *          (uVelocityFarField(i2)*nx+vVelocityFarField(i2)*ny+
     *                                   wVelocityFarField(i2)*nz)*nx
            vTangential=vVelocityFarField(i2)-
     *          (uVelocityFarField(i2)*nx+vVelocityFarField(i2)*ny+
     *                                   wVelocityFarField(i2)*nz)*ny
            wTangential=wVelocityFarField(i2)-
     *          (uVelocityFarField(i2)*nx+vVelocityFarField(i2)*ny+
     *                                   wVelocityFarField(i2)*nz)*nz
c
          else
c
            sbEntropy=cin**2/(ratio*(Density(i1)**(ratio-1)))
c                  
            uTangential=uVelocity(i1)-
     *          (uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz)*nx
            vTangential=vVelocity(i1)-
     *          (uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz)*ny
            wTangential=wVelocity(i1)-
     *          (uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz)*nz
c
          endif           
c
          rhob=(cb*cb/(ratio*sbEntropy))**(1./(ratio-1.))
c          
          BPressure(i2,i3)=RGas*rhob*Tb
          BDensity(i2,i3)=rhob
          BTemperature(i2,i3)=Tb
c
          BuVelocity(i2,i3)=vdotb*nx+uTangential
          BvVelocity(i2,i3)=vdotb*ny+vTangential
          BwVelocity(i2,i3)=vdotb*nz+wTangential
          Bmdot(i2,i3)=rhob*vdotb*BFaceArea(i2,i3)         
c
        enddo
c
c----------------------------------------------------------------------
c
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
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
          ratio=BSpecificHeat(i2,i3)/(BSpecificHeat(i2,i3)-RGas)
c
          vdotin=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
          cin=dsqrt(ratio*RGas*Temperature(i1))
c
          vdotb=vdotin
          cb=cin
          Tb=(cb*cb)/(ratio*RGas)
c
          sbEntropy=cin**2/(ratio*(Density(i1)**(ratio-1)))
c                  
          uTangential=uVelocity(i1)-
     *          (uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz)*nx
          vTangential=vVelocity(i1)-
     *          (uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz)*ny
          wTangential=wVelocity(i1)-
     *          (uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz)*nz
c
          rhob=(cb*cb/(ratio*sbEntropy))**(1./(ratio-1.))
c          
          BPressure(i2,i3)=RGas*rhob*Tb
          BDensity(i2,i3)=rhob
          BTemperature(i2,i3)=Tb
c
          BuVelocity(i2,i3)=vdotb*nx+uTangential
          BvVelocity(i2,i3)=vdotb*ny+vTangential
          BwVelocity(i2,i3)=vdotb*nz+wTangential
          Bmdot(i2,i3)=rhob*vdotb*BFaceArea(i2,i3)         
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
          BuVelocity(i2,i3)=uVelocity(i1)
          BvVelocity(i2,i3)=vVelocity(i1)
          BwVelocity(i2,i3)=wVelocity(i1)
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *           BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!          FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
        enddo
c----------------------------------------------------------------------
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          Bmdot(i2,i3)=BDensity(i2,i3)*(BuVelocity(i2,i3)*sfx+
     *                    BvVelocity(i2,i3)*sfy+BwVelocity(i2,i3)*sfz)
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
          BuVelocity(i2,i3)=uVelocity(i1)+
     *           (BuVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BuVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BuVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BvVelocity(i2,i3)=vVelocity(i1)+
     *           (BvVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BvVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BvVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BwVelocity(i2,i3)=wVelocity(i1)+
     *           (BwVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BwVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BwVelGradz(i2,i3)*BDistanceCFz(i2,i3))
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *           BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
          BuVelocity(i2,i3)=uVelocity(i1)+
     *           (BuVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BuVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BuVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BvVelocity(i2,i3)=vVelocity(i1)+
     *           (BvVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BvVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BvVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BwVelocity(i2,i3)=wVelocity(i1)+
     *           (BwVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BwVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BwVelGradz(i2,i3)*BDistanceCFz(i2,i3))
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *            BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
          BuVelocity(i2,i3)=uVelocity(i1)+
     *           (BuVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BuVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BuVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BvVelocity(i2,i3)=vVelocity(i1)+
     *           (BvVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BvVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BvVelGradz(i2,i3)*BDistanceCFz(i2,i3))
          BwVelocity(i2,i3)=wVelocity(i1)+
     *           (BwVelGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BwVelGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                    BwVelGradz(i2,i3)*BDistanceCFz(i2,i3))
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *             BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          velocity=dsqrt(uVelocity(i1)**2+
     *                    vVelocity(i1)**2+wVelocity(i1)**2)
          xdirection=uVelocity(i1)/(velocity+tiny)
          ydirection=vVelocity(i1)/(velocity+tiny)
          zdirection=wVelocity(i1)/(velocity+tiny)
c
          velMagnitude=Bmdot(i2,i3)/
     *         (rhof*(xdirection*sfx+ydirection*sfy+zdirection*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xdirection
	    BvVelocity(i2,i3)=velMagnitude*ydirection
	    BwVelocity(i2,i3)=velMagnitude*zdirection
c
        enddo
c----------------------------------------------------------------------
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          Bmdot(i2,i3)=BDensity(i2,i3)*(BuVelocity(i2,i3)*sfx+
     *                    BvVelocity(i2,i3)*sfy+BwVelocity(i2,i3)*sfz)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
c
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          Bmdot(i2,i3)=BDensity(i2,i3)*(BuVelocity(i2,i3)*sfx+
     *                    BvVelocity(i2,i3)*sfy+BwVelocity(i2,i3)*sfz)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          velMagnitude=Bmdot(i2,i3)/(rhof*(xVeldirection(i2,i3)*sfx+
     *               yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *           BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
          velMagnitude=Bmdot(i2,i3)/(rhof*(xVeldirection(i2,i3)*sfx+
     *               yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
        enddo
c----------------------------------------------------------------------
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c---  mdot bar term
c
          mdotBarf=uVelocity(i1)*sfx+vVelocity(i1)*sfy+wVelocity(i1)*sfz
          FluxVfLocal=rhof*mdotBarf
c
c--- difference in pressure gradients term
c
          DuSf=Du1Velocity(i1)*sfx
          DvSf=Dv1Velocity(i1)*sfy
          DwSf=Dw1Velocity(i1)*sfz
c
          if(MethodDecomposeS.eq.1) then
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c          
          DuTf=DuSf-DuEf         
          DvTf=DvSf-DvEf
          DwTf=DwSf-DwEf
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
          FluxVfLocal=FluxVfLocal-rhof*(BPressGradx(i2,i3)*DuTf+
     *              BPressGrady(i2,i3)*DvTf+BPressGradz(i2,i3)*DwTf)
c
          FluxVfLocal=FluxVfLocal+rhof*(PressGradx(i1)*DuSf+
     *                     PressGrady(i1)*DvSf+PressGradz(i1)*DwSf)
c
c--- difference in -2/3*rho*k gradient term
c
          if(LSolveTurbulenceKineticEnergy) then
c
            FluxVfLocal=FluxVfLocal-rhof*(Bdrhokdx(i2,i3)*DuTf+
     *                   Bdrhokdy(i2,i3)*DvTf+Bdrhokdz(i2,i3)*DwTf)
c
            FluxVfLocal=FluxVfLocal+rhof*(drhokdx(i1)*DuSf+
     *                         drhokdy(i1)*DvSf+drhokdz(i1)*DwSf)

            FluxVfLocal=FluxVfLocal+
     *                 FluxCfLocal*rhok(i1)+FluxFfLocal*Brhok(i2,i3)
c
          endif
c
c--- Transient term correction
c
          if(LUnsteady) then
c
            mdotOldBar=BDensityOld(i2,i3)*(uVelocityOld(i1)*sfx+
     *                    vVelocityOld(i1)*sfy+wVelocityOld(i1)*sfz)   
            mdotOldOldBar=
     *           BDensityOldOld(i2,i3)*(uVelocityOldOld(i1)*sfx+
     *                  vVelocityOldOld(i1)*sfy+wVelocityOldOld(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *           (DutOldT(i1)*(BmdotOld(i2,i3)-mdotOldBar)+
     *                DutOldOldT(i1)*(BmdotOldOld(i2,i3)-mdotOldOldBar))
c
          endif
c
c--- False transient term correction
c
          if(LFalseTransientMomentum) then
c
            mdotOldBar=BDensityStar(i2,i3)*(uVelocityStar(i1)*sfx+
     *                     vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz)
            FluxVfLocal=FluxVfLocal+
     *                     DutStarT(i1)*(Bmdot(i2,i3)-mdotOldBar)
c
          endif
c
c--- underrelaxation correction
c
          if(LRelaxMomentum) then
c
            mdotstar=rhof*(uVelocityStar(i1)*sfx+
     *                    vVelocityStar(i1)*sfy+wVelocityStar(i1)*sfz) 
            FluxVfLocal=FluxVfLocal+
     *                 (1.-urfMomentum)*(Bmdot(i2,i3)-mdotstar)
c
          endif
c
c--- Surface tension correction
c
!          if(LFreeSurfaceFlow.and.LSurfaceTension) then
!c
!            SurfaceTensionxBar=BCurvature(i2,i3)*BrFieldGradxT(i2,i3)
!            SurfaceTensionyBar=BCurvature(i2,i3)*BrFieldGradyT(i2,i3)
!            SurfaceTensionzBar=BCurvature(i2,i3)*BrFieldGradzT(i2,i3)
!c
!            CurvactureF=BCurvature(i2,i3)
!            FluxVfLocal=FluxVfLocal+rhof*geoDiff*
!     *         SurfaceTension*CurvactureF*(BrFieldT(i2,i3)-rFieldT(i1))
!            FluxVfLocal=FluxVfLocal+rhof*CurvactureF*SurfaceTension*
!     *         (BrFieldGradxT(i2,i3)*DuTf+BrFieldGradyT(i2,i3)*DvTf+
!     *                                      BrFieldGradzT(i2,i3)*DwTf)
!            FluxVfLocal=FluxVfLocal-rhof*SurfaceTension*
!     *        (SurfaceTensionxBar*DuSf+
!     *               SurfaceTensionyBar*DvSf+SurfaceTensionzBar*DwSf)
!c
!          endif
c
c--- Buoyancy correction
c
!          if(LBuoyancy) then
!c
!            BuoyancyxBar=Buoyancyx(i1)      
!            BuoyancyyBar=Buoyancyy(i1)     
!            BuoyancyzBar=Buoyancyz(i1)     
!c
!            FluxVfLocal=FluxVfLocal+rhof*(
!     *                 DuSf*(BBuoyancyx(i2,i3)-BuoyancyxBar)+
!     *                    DvSf*(BBuoyancyy(i2,i3)-BuoyancyyBar)+
!     *                       DwSf*(BBuoyancyz(i2,i3)-BuoyancyzBar))
!c
!          endif
c
          Bmdot(i2,i3)=FluxCfLocal*Pressure(i1)+
     *                   FluxFfLocal*BPressure(i2,i3)+FluxVfLocal
c
          velMagnitude=Bmdot(i2,i3)/(rhof*(xVeldirection(i2,i3)*sfx+
     *             yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
        enddo
c
      endif
c
c--- Deallocate storage
c
      deallocate(Du1f)
      deallocate(Dv1f)
      deallocate(Dw1f)
      deallocate(uVelBar)
      deallocate(vVelBar)
      deallocate(wVelBar)
      deallocate(PressGradxBar)
      deallocate(PressGradyBar)
      deallocate(PressGradzBar)
      deallocate(uVelocityStarf)
      deallocate(vVelocityStarf)
      deallocate(wVelocityStarf)
c
      if(LUnsteady) then
c
        deallocate(uVelOldf)
        deallocate(vVelOldf)
        deallocate(wVelOldf)
        deallocate(uVelOldOldf)
        deallocate(vVelOldOldf)
        deallocate(wVelOldOldf)
        deallocate(rhoOldf)
        deallocate(rhoOldOldf)
        deallocate(DVtfO)
        deallocate(DVtfOO)
        deallocate(DutOldT)
        deallocate(DutOldOldT)
c
      endif
c
      if(LFalseTransientMomentum) then
c
        deallocate(rhoStarf)
        deallocate(DVtfStar)
        deallocate(DutStarT)
c
      endif
c
      if(LSolveTurbulenceKineticEnergy) then
c
        deallocate(drhokdxBar)
        deallocate(drhokdyBar)
        deallocate(drhokdzBar)
        deallocate(drhokdfx)
        deallocate(drhokdfy)
        deallocate(drhokdfz)
c
      endif
c
	return
	end