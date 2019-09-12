c
c#############################################################################################
c
      SUBROUTINE CalculateNormalDistanceToWall
c
C#############################################################################################
c
      use User0, only: TurbulenceModel,WallTreatment,
     *                 LReadSavedWallDistance
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFacesMax,NBFaces,NIFaceOwner,
     *                     NIFaceNeighbor,NBFaceOwner
      use Geometry4, only: gDiff,BgDiff,FaceTx,BFaceTx,FaceTy,BFaceTy,
     *                     FaceTz,BFaceTz,BFaceAreanx,BFaceAreany,
     *                     BFaceAreanz,BDistanceCFx,BDistanceCFy,
     *                     BDistanceCFz,xc,yc,zc,BFaceCentroidx,
     *                     BFaceCentroidy,BFaceCentroidz
      use Variables2
      use BoundaryConditionsTurbulence2, only: IwallTurbulence
      use BoundaryConditions1, only: BoundaryType,wallTypeM
      use WallDistance1     
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Constants1, only: big
      use DirectSolver1, only: b
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      character*16, save :: GradientInterpolationSchemeWall='average'
      logical, save :: lstop=.false.
c
      double precision, save, dimension(:), allocatable :: PhiWall
      double precision, save, dimension(:,:), allocatable :: BPhiWall
      double precision, save, dimension(:), allocatable :: PhiWGradx
      double precision, save, dimension(:), allocatable :: PhiWGrady
      double precision, save, dimension(:), allocatable :: PhiWGradz
      double precision, save, dimension(:,:), allocatable :: BPhiWGradx
      double precision, save, dimension(:,:), allocatable :: BPhiWGrady
      double precision, save, dimension(:,:), allocatable :: BPhiWGradz
      double precision, save, dimension(:), allocatable :: PhiWGradfx
      double precision, save, dimension(:), allocatable :: PhiWGradfy
      double precision, save, dimension(:), allocatable :: PhiWGradfz
      double precision, save, dimension(:), allocatable :: ScWall
      double precision, save, dimension(:), allocatable :: SbWall
c
      integer :: i,j,k,i1,i2,i3,j1,NF,OuterIter,iter,itmax
      double precision :: Residuals,IResiduals,FResiduals,term1,term2,
     *                    FluxCflocal,FluxFflocal,FluxVflocal,nx,ny,nz,
     *                    distance,WallD,BWallD

c********************************************************************************************
c--- Interfaces
c--------------------------------------------------------------------------------
      interface
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
c--------------------------------------------------------------------------------
        SUBROUTINE InterpolateGradientToFace
     *   (GradientInterpolationScheme,FiT,dfidxT,dfidyT,dfidzT,
     *                                    dfidxfT,dfidyfT,dfidzfT)
c--------------------------------------------------------------------------------
          character*16 GradientInterpolationScheme
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:) :: dfidxfT
          double precision, dimension(:) :: dfidyfT
          double precision, dimension(:) :: dfidzfT
c--------------------------------------------------------------------------------
        end SUBROUTINE InterpolateGradientToFace
c--------------------------------------------------------------------------------
        SUBROUTINE AssembleSourceTerm(Variable,FiT,Sc,Sb)
          character*10 Variable
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: Sc
          double precision, dimension(:) :: Sb
c--------------------------------------------------------------------------------
        end SUBROUTINE AssembleSourceTerm
c--------------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------------
c
      if(IwallTurbulence.eq.0) return

      allocate(WallDistance(NumberOfElements))
      allocate(BWallDistance(NumberOfBCSets,NBFacesMax))
      allocate(iTau(NumberOfElements))
      allocate(jTau(NumberOfElements))
      allocate(BiTau(NumberOfBCSets,NBFacesMax))
      allocate(BjTau(NumberOfBCSets,NBFacesMax))
c
      if(LReadSavedWallDistance) then
c
        print*,'Reading saved normal distance to wall...Please wait'
c
          read(15,*) WallDistance,iTau,jTau
          read(15,*)  BWallDistance,BiTau,BjTau
c           
          return
c
      endif
c
      print*,'Calculating normal distance to wall...Please wait'
c
      if(MethodCalculateNormalDistance.eq.1) then
c
        WallDistance=big
c
        do i=1,NumberOfElements
c
          do j=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(j)
            i2=IWallTurbulenceNumberOfBCSets(j)
            i3=IWallTurbulenceNBFaces(j)
c
            distance=dsqrt((xc(i)-BFaceCentroidx(i2,i3))**2+
     *                       (yc(i)-BFaceCentroidy(i2,i3))**2+
     *                          (zc(i)-BFaceCentroidz(i2,i3))**2)
c
            if(distance.lt.WallDistance(i)) then
c
              WallDistance(i)=distance
              iTau(i)=i2
              jTau(i)=i3
c
            endif
c
          enddo
c
        enddo            
c
        BWallDistance=big
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            do k=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(k)
              i2=IWallTurbulenceNumberOfBCSets(k)
              i3=IWallTurbulenceNBFaces(k)
c
              distance=dsqrt(
     *             (BFaceCentroidx(i,j)-BFaceCentroidx(i2,i3))**2+
     *             (BFaceCentroidy(i,j)-BFaceCentroidy(i2,i3))**2+
     *             (BFaceCentroidz(i,j)-BFaceCentroidz(i2,i3))**2)
c
              if(distance.lt.BWallDistance(i,j)) then
c
                BWallDistance(i,j)=distance
                BiTau(i,j)=i2
                BjTau(i,j)=i3
c
              endif
c
            enddo
c           
          enddo            
c
        enddo            
c
      endif 
c
      if(MethodCalculateNormalDistance.eq.2) then
c
        allocate(PhiWall(NumberOfElements))
        allocate(BPhiWall(NumberOfBCSets,NBFacesMax))
        allocate(PhiWGradx(NumberOfElements))
        allocate(PhiWGrady(NumberOfElements))
        allocate(PhiWGradz(NumberOfElements))
        allocate(BPhiWGradx(NumberOfBCSets,NBFacesMax))
        allocate(BPhiWGrady(NumberOfBCSets,NBFacesMax))
        allocate(BPhiWGradz(NumberOfBCSets,NBFacesMax))
        allocate(PhiWGradfx(NIFaces))
        allocate(PhiWGradfy(NIFaces))
        allocate(PhiWGradfz(NIFaces))
        allocate(ScWall(NumberOfElements))
        allocate(SbWall(NumberOfElements))
c        
        print*,'Calculating normal distance to wall'
c
        ScWall=0.
        SbWall=1. 
c
        PhiWall=1.
        itmax=ASIterWallDistance
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'wall'.and.
     *                    wallTypeM(i).eq.'noslip') then

            do j=1,NBFaces(i)
c
              BPhiWall(i,j)=0.
c
            enddo
c
          else

            do j=1,NBFaces(i)
c
              BPhiWall(i,j)=1.
c
            enddo
c      
          endif
c            
        enddo
c
        Variable='PhiW'
c
        Residuals=1.e10
        outerIter=1
        do while(OuterIter.LT.IterWDMax.and.Residuals.gt.WallResiduals)
c
          dphi=0.
c
          call InitializeFluxes
          call InitializeCoefficients
c
c---- Calculate Gradients of the phi variable
c
          call Gradient(Variable,2,PhiWall,PhiWGradx,PhiWGrady,
     *                  PhiWGradz,BPhiWall,BPhiWGradx,BPhiWGrady,
     *                  BPhiWGradz,2,.false.,1)
c
          call InterpolateGradientToFace
     *        (GradientInterpolationSchemeWall,PhiWall,
     *               PhiWGradx,PhiWGrady,PhiWGradz,PhiWGradfx,
     *                                   PhiWGradfy,PhiWGradfz)
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            FluxCflocal= gDiff(k)
            FluxFflocal=-gDiff(k)
            FluxVflocal=-PhiWGradfx(k)*FaceTx(k)-
     *                PhiWGradfy(k)*FaceTy(k)-PhiWGradfz(k)*FaceTz(k)
c
            FluxCf(k)=FluxCf(k)+FluxCflocal
            FluxFf(k)=FluxFf(k)+FluxFflocal
            FluxVf(k)=FluxVf(k)+FluxVflocal
            FluxTf(k)= FluxTf(k)+FluxCflocal*PhiWall(i)+
     *                        FluxFflocal*PhiWall(j)+FluxVflocal
c
          enddo
c
          do i=1,NumberOfBCSets
c
            if(BoundaryType(i).eq.'wall'.and.
     *                    wallTypeM(i).eq.'noslip') then

              do j=1,NBFaces(i)
c
                BPhiWall(i,j)=0.
c
                i1=NBFaceOwner(i,j)
                k=NIFaces
c
                do j1=1,i-1
c
                  k=k+NBFaces(j1)
c
                enddo
c
                k=k+j
c
                FluxCflocal= BgDiff(i,j)
                FluxFflocal=-BgDiff(i,j)
                FluxVflocal=-(BPhiWGradx(i,j)*BFaceTx(i,j)+
     *                            BPhiWGrady(i,j)*BFaceTy(i,j)+
     *                               BPhiWGradz(i,j)*BFaceTz(i,j))
c
                FluxCf(k)=FluxCf(k)+FluxCflocal
                FluxFf(k)=FluxFf(k)+FluxFflocal
                FluxVf(k)=FluxVf(k)+FluxVflocal
                FluxTf(k)=FluxTf(k)+FluxCflocal*PhiWall(i1)
     *                     +FluxFflocal*BPhiWall(i,j)+FluxVflocal
c
              enddo
c
            else
c
              do j=1,NBFaces(i)
c
                i1=NBFaceOwner(i,j)
                BPhiWall(i,j)=PhiWall(i1)
                k=NIFaces
c
                do j1=1,i-1
c
                  k=k+NBFaces(j1)
c
                enddo
c
                k=k+j
c
                FluxCflocal=0.
                FluxFflocal=0.
                FluxVflocal=0.
c
                FluxCf(k)=FluxCf(k)+FluxCflocal
                FluxFf(k)=FluxFf(k)+FluxFflocal
                FluxVf(k)=FluxVf(k)+FluxVflocal
                FluxTf(k)=FluxTf(k)+FluxVflocal
c
              enddo
c
            endif
c
          enddo
c
          call AssembleSourceTerm(Variable,PhiWall,ScWall,SbWall)
c
          call AssembleGlobalMatrixFaceFluxes
          call AssembleGlobalMatrixElementFluxes
c
          call UnderRelaxEquation(urfWall)
c
          NF=-1
          call CheckConvergence(NF,IResiduals)
          FResiduals=IResiduals
c
          iter=0
c
          if(LMultigridWD) then
c
            call AlgebraicMultigridWallDistance
     *          (NF,rrfWD,itmax,ASSolverWD,IResiduals,FResiduals)
c
          else
c
            if(ASSolverWD.eq.'sor') then      
c
              do while(FResiduals.GT.rrfWD*IResiduals.and.iter.lt.itmax)
c
                call SolveEquationsUsingSOR
                call CheckConvergence(NF,FResiduals)
c          
                iter=iter+1
              enddo
c
            elseif(ASSolverWD.eq.'ilu') then      
c
              call ILUFactorization
              do while(FResiduals.GT.rrfWD*IResiduals.and.iter.lt.itmax)
c
                call SolveEquationsUsingILU
                call CheckConvergence(NF,FResiduals)
c
                iter=iter+1
              enddo
c
            elseif(ASSolverWD.eq.'pbcg') then      
c
              call allocateForPBCG
              call sprsin
              call ILUFactorization
c
              do while(FResiduals.GT.rrfWD*IResiduals.and.iter.lt.itmax)
c
                call SolveEquationsUsingPBCG(iter)
                call CheckConvergence(NF,FResiduals)
c
                iter=iter+1
              enddo
c
            elseif(ASSolverWD.eq.'direct') then 
c
              call SetupMatrixforDirectSolver(NumberOfElements)
              call ludcmpF77
              call lubksbF77
              dphi=b
              call deallocateMatrixStorage
c
            endif
c
          endif
c
          PhiWall=PhiWall+dphi
          Residuals=FResiduals
          outerIter=outerIter+1
c
          print*,outerIter,'   Residuals= ',FResiduals
c        
        enddo
c
        if(Residuals.GT.WallResiduals) then
c        
          Print*,'Solutions of normal distance to wall diverged'
          Print*,'Program has stopped'
c
          pause
          lstop=.true.
c
        else          
c        
          Print*,'Normal distance to wall successfully generated'
c        
        endif
c
        Do i=1,NumberOfElements
c
          term1=PhiWGradx(i)**2+PhiWGrady(i)**2+PhiWGradz(i)**2
          term2=term1+2.*PhiWall(i)
c
          WallDistance(i)=-dsqrt(term1)+dsqrt(term2)        
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
          term1=BPhiWGradx(i,j)**2+BPhiWGrady(i,j)**2+BPhiWGradz(i,j)**2
          term2=term1+2.*BPhiWall(i,j)
c
          BWallDistance(i,j)=-dsqrt(term1)+dsqrt(term2)        
c
          enddo
        enddo
c
        deallocate(PhiWall)
        deallocate(BPhiWall)
        deallocate(PhiWGradx)
        deallocate(PhiWGrady)
        deallocate(PhiWGradz)
        deallocate(BPhiWGradx)
        deallocate(BPhiWGrady)
        deallocate(BPhiWGradz)
        deallocate(PhiWGradfx)
        deallocate(PhiWGradfy)
        deallocate(PhiWGradfz)
        deallocate(ScWall)
        deallocate(SbWall)
c
c---- Calculate iTau and jTau for printing
c
        do i=1,NumberOfElements
c
          WallD=big
c
          do j=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(j)
            i2=IWallTurbulenceNumberOfBCSets(j)
            i3=IWallTurbulenceNBFaces(j)
c
            distance=dsqrt((xc(i)-BFaceCentroidx(i2,i3))**2+
     *                       (yc(i)-BFaceCentroidy(i2,i3))**2+
     *                          (zc(i)-BFaceCentroidz(i2,i3))**2)
c
            if(distance.lt.WallD) then
c
              WallD=distance
              iTau(i)=i2
              jTau(i)=i3
c
            endif
c
          enddo
c
        enddo            
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BWallD=big
c
            do k=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(k)
              i2=IWallTurbulenceNumberOfBCSets(k)
              i3=IWallTurbulenceNBFaces(k)
c
              distance=dsqrt(
     *             (BFaceCentroidx(i,j)-BFaceCentroidx(i2,i3))**2+
     *             (BFaceCentroidy(i,j)-BFaceCentroidy(i2,i3))**2+
     *             (BFaceCentroidz(i,j)-BFaceCentroidz(i2,i3))**2)
c
              if(distance.lt.BWallD) then
c
                BWallD=distance
                BiTau(i,j)=i2
                BjTau(i,j)=i3
c
              endif
c
            enddo
c           
          enddo            
c
        enddo            
c
      endif
c
c--- Over write distance to wall for the first adjacent cell to wall 
c
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        nx=BFaceAreanx(i2,i3)
        ny=BFaceAreany(i2,i3)
        nz=BFaceAreanz(i2,i3)
        WallDistance(i1)=BDistanceCFx(i2,i3)*nx+
     *              BDistanceCFy(i2,i3)*ny+BDistanceCFz(i2,i3)*nz
c
        iTau(i1)=i2
        jTau(i1)=i3
c
      enddo
c
c--- Save WallDistance for later runs
c
      write(15,*) WallDistance,iTau,jTau
      write(15,*)  BWallDistance,BiTau,BjTau
c
      if(lstop) stop
c      
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CheckPositivityOfNormalDistanceToWall
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use WallDistance1, only: WallDistance,BWallDistance   
      use Geometry3, only: NBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
        Do i=1,NumberOfElements
c
          if(WallDistance(i).lt.0.) then
c
            print*,'Normal distance to wall has negative values'
            print*,'Program stopped'
            print*,WallDistance(i)
            pause
            stop
c
          endif        
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            if(BWallDistance(i,j).lt.0.) then
c
              print*,'Normal distance to wall has negative values'
              print*,'Program stopped'
              print*,WallDistance(i)
              pause
              stop
c
            endif        
c
          enddo
        enddo
c
      return
      end
 