c
C#############################################################################################
      SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
C#############################################################################################
        implicit none
c*********************************************************************************************
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
c*********************************************************************************************
      interface
c*********************************************************************************************
        SUBROUTINE StandardGreenGauss
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c*********************************************************************************************
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE StandardGreenGauss
c--------------------------------------------------------------------------
        SUBROUTINE StandardGreenGaussFaceCorrection1
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,
     *                               BdfidzT,nIterGradientPhi)
c--------------------------------------------------------------------------
          integer :: nIterGradientPhi
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE StandardGreenGaussFaceCorrection1
c--------------------------------------------------------------------------
        SUBROUTINE StandardGreenGaussFaceCorrection2
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,
     *                               BdfidzT,nIterGradientPhi)
c--------------------------------------------------------------------------
          integer :: nIterGradientPhi
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE StandardGreenGaussFaceCorrection2
c--------------------------------------------------------------------------
        SUBROUTINE StandardGreenGaussFaceCorrection3
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,
     *                               BdfidzT,nIterGradientPhi)
c--------------------------------------------------------------------------
          integer :: nIterGradientPhi
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE StandardGreenGaussFaceCorrection3
c--------------------------------------------------------------------------
        SUBROUTINE NodeBasedGreenGauss
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE NodeBasedGreenGauss
c--------------------------------------------------------------------------
        SUBROUTINE LeastSquares
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE LeastSquares
c--------------------------------------------------------------------------
        SUBROUTINE LimitGradientVariable(FiT,dfidxT,dfidyT,dfidzT,
     *                 BFiT,BdfidxT,BdfidyT,BdfidzT,LimitGradientMethod)
c--------------------------------------------------------------------------
          integer :: LimitGradientMethod
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE LimitGradientVariable
c--------------------------------------------------------------------------
        SUBROUTINE CorrectGradientAtOulet(dfidxT,dfidyT,dfidzT,
     *                                    BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE CorrectGradientAtOulet
c--------------------------------------------------------------------------
        SUBROUTINE CorrectGradientAtPeriodicBoundary
     *              (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE CorrectGradientAtPeriodicBoundary
c--------------------------------------------------------------------------
        SUBROUTINE CorrectPresureGradientAtWalls
     *              (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE CorrectPresureGradientAtWalls
c--------------------------------------------------------------------------
        SUBROUTINE CorrectPressureGradientAtSymmetryBoundary
     *               (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE CorrectPressureGradientAtSymmetryBoundary
c--------------------------------------------------------------------------
      end interface
c--------------------------------------------------------------------------
c
C---   METHOD_CALC_GRADIENT=1 ----->   Standard Green-Gauss Gradient
C---   METHOD_CALC_GRADIENT=2 ----->   Green-Gauss Gradient with face correction (f' at intersection beween CF and face f) option 1 in book
C---   METHOD_CALC_GRADIENT=3 ----->   Green-Gauss Gradient with face correction (f' at center of CF)                      option 2 in book
C---   METHOD_CALC_GRADIENT=4 ----->   Green-Gauss Gradient with face correction (f' is chosen such that ff' is shortest)  option 3 in book
C---   METHOD_CALC_GRADIENT=5 ----->   Node-based Green-Gauss Gradient 
C---   METHOD_CALC_GRADIENT=6 ----->   Least Square Gradient 
c
c--- Start by initializing
c
      if(MethodCalcGradient.eq.1) then
c
        call StandardGreenGauss
     *    (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c
      elseif(MethodCalcGradient.eq.2) then
c
        call StandardGreenGaussFaceCorrection1
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,
     *                               BdfidzT,nIterGradientPhi)
c
      elseif(MethodCalcGradient.eq.3) then
c
        call StandardGreenGaussFaceCorrection2
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,
     *                               BdfidzT,nIterGradientPhi)
c
      elseif(MethodCalcGradient.eq.4) then
c
        call StandardGreenGaussFaceCorrection3
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,
     *                               BdfidzT,nIterGradientPhi)
c
      elseif(MethodCalcGradient.eq.5) then
c
        call NodeBasedGreenGauss
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c
      elseif(MethodCalcGradient.eq.6) then
c
        call LeastSquares
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c
      endif      
c
c---- Limit gradient
c
      if(limitgradient) call LimitGradientVariable(FiT,dfidxT,dfidyT,
     *         dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,LimitGradientMethod)
c
c---- Correct gradient along boundaries
c
      if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.Variable.eq.'velz')
     *      call CorrectGradientAtOulet(dfidxT,dfidyT,dfidzT,
     *                                    BdfidxT,BdfidyT,BdfidzT)
c
      call CorrectGradientAtPeriodicBoundary
     *              (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c
      if(Variable.eq.'press'.or.Variable.eq.'rhok')
     *  call CorrectPressureGradientAtSymmetryBoundary
     *              (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c
      if(Variable.eq.'press'.or.Variable.eq.'rhok')
     *   call CorrectPresureGradientAtWalls
     *              (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c
      return
      end