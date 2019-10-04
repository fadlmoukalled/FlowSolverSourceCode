c
c#############################################################################################
c
      SUBROUTINE SolveLambdaEulerLagrangeEquation
c
C#############################################################################################
c
      use User0
      use PhysicalProperties1
      use Variables1
      use Variables4
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      character*10 Variable
c********************************************************************************************
c--- Interfaces
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
c--------------------------------------------------------------
        SUBROUTINE AssembleAnisotropicDiffusionTerm
     *       (Variable,FiT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *                            dfidxfT,dfidyfT,dfidzfT)
c--------------------------------------------------------------
          character*10 Variable
          double precision, dimension(:) :: FiT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
          double precision, dimension(:) :: dfidxfT
          double precision, dimension(:) :: dfidyfT
          double precision, dimension(:) :: dfidzfT
c--------------------------------------------------------------
        end SUBROUTINE AssembleAnisotropicDiffusionTerm
c--------------------------------------------------------------
        SUBROUTINE AssembleSourceTerm(Variable,FiT,Sc,Sb)
          character*10 Variable
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: Sc
          double precision, dimension(:) :: Sb
c--------------------------------------------------------------
        end SUBROUTINE AssembleSourceTerm
c--------------------------------------------------------------
        SUBROUTINE SolveEquation(Variable,FiT,rrF,itmax,
     *                                       solver,LMultigrid)
c--------------------------------------------------------------
          character*10 Variable
          character*6 solver
          logical LMultigrid
          integer itmax  
          double precision rrF
          double precision, dimension(:)  :: FiT
c--------------------------------------------------------------
        end SUBROUTINE SolveEquation
c********************************************************************************************
      end interface
C*********************************************************************************************
c
      Variable='lambda'
c
      call InitializeFluxes
      call InitializeCoefficients
c
c---- Calculate Gradients of the phi variable
c
      call Gradient(Variable,MethodCalcGradientLambdaELE,
     *    LambdaELE,LambdaELEGradx,LambdaELEGrady,
     *     LambdaELEGradz,BLambdaELE,BLambdaELEGradx,
     *      BLambdaELEGrady,BLambdaELEGradz,
     *        nIterGradientLambdaELE,LimitGradientLambdaELE,
     *                              LimitGradientLambdaELEMethod)
c
      call InterpolateGradientToFace
     *      (GradientInterpolationSchemeLambdaELE,
     *        LambdaELE,LambdaELEGradx,LambdaELEGrady,
     *         LambdaELEGradz,LambdaELEGradfx,
     *                 LambdaELEGradfy,LambdaELEGradfz)
c
c---- Assemble terms
c
      call AssembleAnisotropicDiffusionTerm(Variable,LambdaELE,
     *   BLambdaELE,BLambdaELEGradx,BLambdaELEGrady,BLambdaELEGradz,
     *                LambdaELEGradfx,LambdaELEGradfy,LambdaELEGradfz)
c
      call AssembleSourceTerm
     *           (Variable,LambdaELE,ScLambdaELE,SbLambdaELE)
c
c---- Underrelax using false transient
c
      if(LFalseTransientLambdaELE) 
     *      call AssembleFalseTransient(Variable,FalsedtLambdaELE)
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxLambdaELE) call UnderRelaxEquation(urfLambdaELE)
c
c---- Solve equations
c
	call SolveEquation(Variable,LambdaELE,rrfLambdaELE,
     *   ASIterLambdaELE,ASSolverLambdaELE,LMultigridLambdaELE)
c
	return
	end
c
c#############################################################################################
c
      SUBROUTINE CalculateLambdaELESources
c
C#############################################################################################
c
      use User0, only: MethodCalcGradientMomentum,nIterGradientMomentum,
     *                LimitGradientMomentum,LimitGradientMomentumMethod,
     *                alpha1Lambda,IAbsDivergence,IMaxDivergence,
     *                IrmsDivergence
      use Variables1, only: uVelocity,uVelGradx,uVelGrady,uVelGradz,
     *                    BuVelocity,BuVelGradx,BuVelGrady,BuVelGradz,
     *                    vVelocity,vVelGradx,vVelGrady,vVelGradz,
     *                    BvVelocity,BvVelGradx,BvVelGrady,BvVelGradz,
     *                    wVelocity,wVelGradx,wVelGrady,wVelGradz,
     *                    BwVelocity,BwVelGradx,BwVelGrady,BwVelGradz,
     *                    LambdaELE,initialVelDivergence,
     *                    BinitialVelDivergence
      use Variables4, only: ScLambdaELE,SbLambdaELE
      use Geometry1, only: NumberOfElements
      use constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i
      double precision :: Term,Term1
c********************************************************************************************
c--- Interfaces
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
C*********************************************************************************************
c
c---- Calculate Gradients of the u,v,and w variable
c
      Variable='velx'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelocity,uVelGradx,uVelGrady,uVelGradz,BuVelocity,
     *       BuVelGradx,BuVelGrady,BuVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      Variable='vely'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelocity,vVelGradx,vVelGrady,vVelGradz,BvVelocity,
     *        BvVelGradx,BvVelGrady,BvVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      Variable='velz'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelocity,wVelGradx,wVelGrady,wVelGradz,BwVelocity,
     *        BwVelGradx,BwVelGrady,BwVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      do i=1,NumberOfElements
c
        Term=2.*alpha1Lambda*alpha1Lambda*
     *            (uVelGradx(i)+vVelGrady(i)+wVelGradz(i))
        ScLambdaELE(i)=0.
        SbLambdaELE(i)=Term
c
      enddo
c
c--- Save Initial velocity field divergence
c
      InitialVelDivergence=uVelGradx+vVelGrady+wVelGradz
      BInitialVelDivergence=BuVelGradx+BvVelGrady+BwVelGradz
c
c--- Calculate the sum of the absolute divergence values
c
      IAbsDivergence=0.    
      do i=1,NumberOfElements
        IAbsDivergence=IAbsDivergence+dabs(InitialVelDivergence(i))
      enddo
c
c--- Calculate the maximum divergence value
c
      IMaxDivergence=-1    
      do i=1,NumberOfElements
        IMaxDivergence=
     *         dmax1(IMaxDivergence,dabs(InitialVelDivergence(i)))
      enddo
c
c--- Calculate the rms divergence value
c
      IrmsDivergence=0.    
      do i=1,NumberOfElements
        IrmsDivergence=IrmsDivergence+InitialVelDivergence(i)**2
      enddo
      IrmsDivergence=dsqrt(IrmsDivergence/NumberOfElements)
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE UpdateEulerLagrangeVelocity
c
C#############################################################################################
c
      use User0, only: MethodCalcGradientMomentum,nIterGradientMomentum,
     *                LimitGradientMomentum,LimitGradientMomentumMethod,
     *                FAbsDivergence,FMaxDivergence,FrmsDivergence,
     *                alpha1Lambda,alpha2Lambda
      use Variables1, only: uVelocity,uVelGradx,uVelGrady,uVelGradz,
     *                    BuVelocity,BuVelGradx,BuVelGrady,BuVelGradz,
     *                    vVelocity,vVelGradx,vVelGrady,vVelGradz,
     *                    BvVelocity,BvVelGradx,BvVelGrady,BvVelGradz,
     *                    wVelocity,wVelGradx,wVelGrady,wVelGradz,
     *                    BwVelocity,BwVelGradx,BwVelGrady,BwVelGradz,
     *                    LambdaELE,LambdaELEGradx,LambdaELEGrady,
     *                    LambdaELEGradz,BLambdaELEGradx,
     *                    BLambdaELEGrady,BLambdaELEGradz,
     *                    FinalVelDivergence,BFinalVelDivergence
      use BoundaryConditions1, only: wallTypeL
      use Geometry1, only: NumberOfBCSets,NumberOfElements
      use Geometry3, only: NBFaces
C*********************************************************************************************
      implicit none
C*********************************************************************************************
      integer :: i,j
      double precision :: factor1,factor2
      character*10 :: Variable
c********************************************************************************************
c--- Interfaces
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
C*********************************************************************************************
c
      factor1=1./(2.*alpha1Lambda*alpha1Lambda)
      factor2=1./(2.*alpha2Lambda*alpha2Lambda)
      uVelocity=uVelocity+factor1*LambdaELEGradx
      vVelocity=vVelocity+factor1*LambdaELEGrady
      wVelocity=wVelocity+factor2*LambdaELEGradz
c
      do i=1,NumberOfBCSets
        if(wallTypeL(i).eq.'dirichlet') then
          do j=1,NBFaces(i)
c
            BuVelocity(i,j)=BuVelocity(i,j)+factor1*BLambdaELEGradx(i,j)
            BvVelocity(i,j)=BvVelocity(i,j)+factor1*BLambdaELEGrady(i,j)
            BwVelocity(i,j)=BwVelocity(i,j)+factor2*BLambdaELEGradz(i,j)
c
          enddo  
        endif 
      enddo  
c
c---- Calculate Gradients of the u,v,and w variable (for printing)
c
      Variable='velx'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      uVelocity,uVelGradx,uVelGrady,uVelGradz,BuVelocity,
     *       BuVelGradx,BuVelGrady,BuVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      Variable='vely'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      vVelocity,vVelGradx,vVelGrady,vVelGradz,BvVelocity,
     *        BvVelGradx,BvVelGrady,BvVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
      Variable='velz'
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *      wVelocity,wVelGradx,wVelGrady,wVelGradz,BwVelocity,
     *        BwVelGradx,BwVelGrady,BwVelGradz,nIterGradientMomentum,
     *             LimitGradientMomentum,LimitGradientMomentumMethod)
c
c--- Save Final velocity field divergence
c
      FinalVelDivergence=uVelGradx+vVelGrady+wVelGradz
      BFinalVelDivergence=BuVelGradx+BvVelGrady+BwVelGradz
c
c--- Calculate the sum of the absolute divergence values
c
      FAbsDivergence=0.    
      do i=1,NumberOfElements
        FAbsDivergence=FAbsDivergence+dabs(FinalVelDivergence(i))
      enddo
c
c--- Calculate the maximum divergence value
c
      FMaxDivergence=-1    
      do i=1,NumberOfElements
        FMaxDivergence=dmax1(FMaxDivergence,dabs(FinalVelDivergence(i)))
      enddo
c
c--- Calculate the rms divergence value
c
      FrmsDivergence=0.    
      do i=1,NumberOfElements
        FrmsDivergence=FrmsDivergence+FinalVelDivergence(i)**2
      enddo
      FrmsDivergence=DSQRT(FrmsDivergence/NumberOfElements)
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateLambdaDiffusionTensor
c
C#############################################################################################
c
      use User0, only: alpha1Lambda,alpha2Lambda
      use PhysicalProperties1, only: Conductivity11,Conductivity12,
     *                               Conductivity13,Conductivity22,
     *                               Conductivity23,Conductivity33,
     *                               BConductivity11,BConductivity12,
     *                               BConductivity13,BConductivity22,
     *                               BConductivity23,BConductivity33
C*********************************************************************************************
      implicit none
C*********************************************************************************************
c 
      Conductivity11=1.
      Conductivity12=0.
      Conductivity13=0.
      Conductivity22=1.
      Conductivity23=0.
      Conductivity33=(alpha1Lambda/alpha2Lambda)**2
      BConductivity11=1.
      BConductivity12=0.
      BConductivity13=0.
      BConductivity22=1.
      BConductivity23=0.
      BConductivity33=(alpha1Lambda/alpha2Lambda)**2
c      
      return
      end
c
c#############################################################################################
c
      SUBROUTINE IntializeLambdaVelocityField
c
C#############################################################################################
c
      use ReadpolyMesh
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use BoundaryConditions1, only: wallTypeL
      use Variables1, only: uVelocity,vVelocity,wVelocity,        
     *                      BuVelocity,BvVelocity,BwVelocity
C*********************************************************************************************
      implicit none
C*********************************************************************************************
      integer :: i,j,k
C*********************************************************************************************
c
      call InitializeV
c
        do i=1,NumberOfBCSets
c
          if(wallTypeL(i).eq.'dirichlet') then
c
            do j=1,NBFaces(i)      
c
              k=NBFaceOwner(i,j)
c              
              BuVelocity(i,j)=uVelocity(k)
              BvVelocity(i,j)=vVelocity(k)
              BwVelocity(i,j)=wVelocity(k)
c             
            enddo
c
          elseif(wallTypeL(i).eq.'vonneumann') then
c
            do j=1,NBFaces(i)      
c
              BuVelocity(i,j)=0.
              BvVelocity(i,j)=0.
              BwVelocity(i,j)=0.
c             
            enddo
c
          endif
        enddo
c
      
      
      
      
      
      return
      end