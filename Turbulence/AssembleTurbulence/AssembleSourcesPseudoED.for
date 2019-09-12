c
c#############################################################################################
c
      SUBROUTINE AssemblePseudoEddyViscositySources
c
C#############################################################################################
c
      use User0, only: MethodCalcGradientMomentum,nIterGradientMomentum,
     *                 LimitGradientMomentum,LimitGradientMomentumMethod
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: Volume,BFaceAreanx,BFaceAreany,
     *                     BDistanceCFx,BDistanceCFy,BgDiff,
     *                     BFaceTx,BFaceTy,BFaceTz
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      TurbulenceProduction,
     *                      ModifiedED,ModifiedEDGradx,
     *                      ModifiedEDGrady,ModifiedEDGradz,BModifiedED,
     *                      BModifiedEDGradx,BModifiedEDGrady,
     *                      BModifiedEDGradz,TGamma
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density,BeDiffCoefficient
      use Turbulence1, only: c1,c2,c3,fr1Coefficient,f2RT,velRT,BvelRT,
     *                       velRTGradx,velRTGrady,velRTGradz,
     *                       BvelRTGradx,BvelRTGrady,BvelRTGradz
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use WallDistance1, only: WallDistance
      use ReferenceValues1, only: uVelFrameOfReference,
     *                            vVelFrameOfReference,
     *                            wVelFrameOfReference
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,i1,i2,i3,i4,j
      double precision :: uo,vo,wo,u1,v1,w1,term1,term2,LambdaRT,DRT,
     *                    gamf,dfidxTf,dfidyTf,dfidzTf,FluxCflocal,
     *                    FluxFflocal,FluxVflocal
c********************************************************************************************
c--------------------------------------------------------------
c--- Interfaces
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------------------
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
c--------------------------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------------------------
c
c--- Start by calculating the gradient of the velocity magnitude
c
      uo=uVelFrameOfReference      
      vo=vVelFrameOfReference      
      wo=wVelFrameOfReference      
c
      do i=1,NumberOfElements     
c
        u1=uVelocity(i)-uo
        v1=vVelocity(i)-vo
        w1=wVelocity(i)-wo
        velRT(i)=dsqrt(u1*u1+v1*v1+w1*w1)
c
      enddo        
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
        u1=BuVelocity(i,j)-uo
        v1=BvVelocity(i,j)-vo
        w1=BwVelocity(i,j)-wo
        BvelRT(i,j)=dsqrt(u1*u1+v1*v1+w1*w1)
c
        enddo
      enddo
c
c---- Calculate Gradient of VelRT
c
      Variable='velrt'
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      velRT,velRTGradx,velRTGrady,velRTGradz,
     *        BvelRT,BvelRTGradx,BvelRTGrady,
     *          BvelRTGradz,nIterGradientMomentum,
     *            LimitGradientMomentum,LimitGradientMomentumMethod)
c
      do i=1,NumberOfElements     
c
        term1=Density(i)*dmax1(ModifiedED(i),0.)*TurbulenceProduction(i)
        term2=fr1Coefficient(i)*(c1-c2*f2RT(i))*dsqrt(term1)
c
        FluxTE(i)=FluxTE(i)-term2*Volume(i)
c
        LambdaRT=velRTGradx(i)*ModifiedEDGradx(i)+ 
     *             velRTGrady(i)*ModifiedEDGrady(i)+            
     *               velRTGradz(i)*ModifiedEDGradz(i)           
c
        DRT=0.
        if(LambdaRT.gt.0.) then
          DRT=ModifiedEDGradx(i)*ModifiedEDGradx(i)+ 
     *             ModifiedEDGrady(i)*ModifiedEDGrady(i)+            
     *               ModifiedEDGradz(i)*ModifiedEDGradz(i)
        endif
c
        term1=Density(i)*c3*DRT
        FluxCE(i)=FluxCE(i)+(term1/dmax1(ModifiedED(i),tiny))*Volume(i)
        FluxTE(i)=FluxTE(i)+term1*Volume(i)
c
      enddo
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
     *               dfidyTf*BFaceTy(i2,i3)+dfidzTf*BFaceTz(i2,i3))
c
        FluxCf(i4)=FluxCf(i4)+FluxCflocal
        FluxFf(i4)=FluxFf(i4)+FluxFflocal
        FluxVf(i4)=FluxVf(i4)+FluxVflocal
        FluxTf(i4)=FluxTf(i4)+FluxCflocal*dmax1(ModifiedED(i1),0.)+
     *                      FluxFflocal*BModifiedED(i2,i3)+FluxVflocal
c
      enddo
c
      return
      end