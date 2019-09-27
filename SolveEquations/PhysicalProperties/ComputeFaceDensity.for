c
c#############################################################################################
c
      SUBROUTINE CalculateFaceDensity
c
c#############################################################################################
c
      use User0, only: MethodCalcGradientDensity,BleedDensity,
     *                 nIterGradientDensity,LimitGradientDensity,
     *                 LimitGradientDensityMethod,HRFrameworkDensity,
     *                 ConvectionSchemeDensity,
     *                 nIterStartApplyingHR,LConvectScalar,
     *                 LSolveMomentum,LFreeSurfaceFlow,LHarmonic
      use MultiGrid2, only: nIter
      use PhysicalProperties1, only: Density,DensGradx,DensGrady,
     *                               DensGradz,BDensity,BDensGradx,
     *                               BDensGrady,BDensGradz,Densityf
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      character*10 Variable 
      character*20 ConvectionScheme1 
      character*20 InterpolationScheme 
c
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
        SUBROUTINE InterpolateToFaceUsingHRscheme(ConvectionScheme,
     *               Bleed,HRFramework,FiT,dfidxT,dfidyT,dfidzT,FiTf)
c--------------------------------------------------------------
          character*20 ConvectionScheme
          character*4 HRFramework
          double precision :: Bleed
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:) :: FiTf
c--------------------------------------------------------------
        end SUBROUTINE InterpolateToFaceUsingHRscheme
c--------------------------------------------------------------
       SUBROUTINE InterpolateElementToFace(InterpolationScheme,FiT,FiTf)
c--------------------------------------------------------------
          character*16 InterpolationScheme
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: FiTf
        end SUBROUTINE InterpolateElementToFace
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      If(.not.LSolveMomentum.and..not.LConvectScalar) then
c
        InterpolationScheme='average'
        call InterpolateElementToFace
     *            (InterpolationScheme,Density,Densityf)
c
      elseif(LFreeSurfaceFlow.and.LHarmonic) then 
c
        InterpolationScheme='harmonic'
        call InterpolateElementToFace
     *            (InterpolationScheme,Density,Densityf)
c
      elseif(HRFrameworkDensity.eq.'nvf'.or.
     *                  HRFrameworkDensity.eq.'tvd') then
c
        Variable='Dens'
c
        call Gradient(Variable,MethodCalcGradientDensity,Density,
     *          DensGradx,DensGrady,DensGradz,BDensity,BDensGradx,
     *            BDensGrady,BDensGradz,nIterGradientDensity,
     *              LimitGradientDensity,LimitGradientDensityMethod)
c
        if(nIter.ge.nIterStartApplyingHR) then
c
          call InterpolateToFaceUsingHRscheme(ConvectionSchemeDensity,
     *     BleedDensity,HRFrameworkDensity,Density,DensGradx,
     *             DensGrady,DensGradz,Densityf)
c
        else
c
          ConvectionScheme1='upwind'
          call InterpolateToFaceUsingHRscheme(ConvectionScheme1,
     *     BleedDensity,HRFrameworkDensity,Density,DensGradx,
     *             DensGrady,DensGradz,Densityf)
c
        endif
c
      else
c
        InterpolationScheme='average'
        call InterpolateElementToFace
     *                  (InterpolationScheme,Density,Densityf)
c
      endif
c
      return
      end
