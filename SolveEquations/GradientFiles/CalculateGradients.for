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
      SUBROUTINE SaveGradient(dfidxT,dfidyT,dfidzT)
c*********************************************************************************************
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
c--------------------------------------------------------------------------
          end SUBROUTINE SaveGradient
c--------------------------------------------------------------------------
        SUBROUTINE RelaxGradient
     *      (Variable,dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c*********************************************************************************************
          character*10 :: Variable
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
          end SUBROUTINE RelaxGradient
c--------------------------------------------------------------------------
        SUBROUTINE StandardGreenGauss
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
c--- Start by saving for under relaxation
c
      call SaveGradient(dfidxT,dfidyT,dfidzT)
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
c---- Under relax gradient
c
      call RelaxGradient
     *      (Variable,dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
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
c
C#############################################################################################
      SUBROUTINE SaveGradient(dfidxT,dfidyT,dfidzT)
C#############################################################################################
      use Variables1, only: dfidxTstar,dfidyTstar,dfidzTstar
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
c*********************************************************************************************
c
      dfidxTstar=dfidxT
      dfidyTstar=dfidyT
      dfidzTstar=dfidzT
c
      return
      end
c
c*********************************************************************************************
      SUBROUTINE RelaxGradient
     *      (Variable,dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c*********************************************************************************************
      use User0, only: ScalarName,rFieldName,LrelaxGradientMomentum,
     *    LrelaxGradientContinuity,LrelaxGradientEnergy,
     *    LrelaxGradientTKE,LrelaxGradientTED,LrelaxGradientTOmega,
     *    LrelaxGradientMED,LrelaxGradientTKL,
     *    LrelaxGradientTGamma,LrelaxGradientTReTheta,
     *    LrelaxGradientTv2,LrelaxGradientTZeta,
     *    LrelaxGradientTfRelaxation,LrelaxGradientLambdaELE,
     *    urfGradientMomentum,urfGradientContinuity,urfGradientEnergy,
     *    urfGradientTKE,urfGradientTED,urfGradientTOmega,
     *    urfGradientMED,urfGradientTKL,
     *    urfGradientTGamma,urfGradientTReTheta,
     *    urfGradientTv2,urfGradientTZeta,
     *    urfGradientTfRelaxation,urfGradientLambdaELE,
     *    LrelaxGradientScalar,urfGradientScalar,
     *    LrelaxGradientrField,urfGradientrField,
     *    LrelaxGradientDensity,urfGradientDensity,
     *    LrelaxGradientOthers,urfGradientOthers
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Scalar2, only: iScalarVariable
      use VolumeOfFluid2, only: irFieldVariable
      use Variables1, only: dfidxTstar,dfidyTstar,dfidzTstar
c*********************************************************************************************
        implicit none
c*********************************************************************************************
        character*10 :: Variable
        integer :: i,j,k
        double precision :: urf,urf1
        double precision, dimension(:) :: dfidxT
        double precision, dimension(:) :: dfidyT
        double precision, dimension(:) :: dfidzT
        double precision, dimension(:,:) :: BdfidxT
        double precision, dimension(:,:) :: BdfidyT
        double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
c
      if(Variable.eq.'velx'.or.Variable.eq.
     *        'vely'.or.Variable.eq.'velz'.or.Variable.eq.'rhok'.or.
     *         Variable.eq.'nvel'.or.Variable.eq.'velrt'.or.
     *                                  Variable.eq.'velgrad') then
c
        if(LrelaxGradientMomentum) then
          urf=urfGradientMomentum
        else
          return
        endif
c
      elseif(Variable.eq.'pressc'.or.Variable.eq.'press') then
c
        if(LrelaxGradientContinuity) then
          urf=urfGradientContinuity
        else
          return
        endif
c
      elseif(Variable.eq.'Dens'.or.Variable.eq.'DensO'.or.
     *       Variable.eq.'DensOO'.or.Variable.eq.'density') then
c
        if(LrelaxGradientDensity) then
          urf=urfGradientDensity
        else
          return
        endif
c
      elseif(Variable.eq.'temp'.or.Variable.eq.'htotal') then
c
        if(LrelaxGradientEnergy) then
          urf=urfGradientEnergy
        else
          return
        endif
c
      elseif(Variable.eq.'tke') then
c
        if(LrelaxGradientTKE) then
          urf=urfGradientTKE
        else
          return
        endif
c
      elseif(Variable.eq.'ted') then
c
        if(LrelaxGradientTED) then
          urf=urfGradientTED
        else
          return
        endif
c
      elseif(Variable.eq.'tomega') then
c
        if(LrelaxGradientTOmega) then
          urf=urfGradientTOmega
        else
          return
        endif
c
      elseif(Variable.eq.'med') then
c
        if(LrelaxGradientMED) then
          urf=urfGradientMED
        else
          return
        endif
c
      elseif(Variable.eq.'tkl') then
c
        if(LrelaxGradientTKL) then
          urf=urfGradientTKL
        else
          return
        endif
c
      elseif(Variable.eq.'tgamma') then
c
        if(LrelaxGradientTGamma) then
          urf=urfGradientTGamma
        else
          return
        endif
c
      elseif(Variable.eq.'tretheta') then
c
        if(LrelaxGradientTReTheta) then
          urf=urfGradientTReTheta
        else
          return
        endif
c
      elseif(Variable.eq.'tv2') then
c
        if(LrelaxGradientTv2) then
          urf=urfGradientTv2
        else
          return
        endif
c
      elseif(Variable.eq.'tzeta') then
c
        if(LrelaxGradientTZeta) then
          urf=urfGradientTZeta
        else
          return
        endif
c
      elseif(Variable.eq.'frelax') then
c
        if(LrelaxGradientTfRelaxation) then
          urf=urfGradientTfRelaxation
        else
          return
        endif
c
      elseif(Variable.eq.'lambda') then
c
        if(LrelaxGradientLambdaELE) then
          urf=urfGradientLambdaELE
        else
          return
        endif
c
      else
c
        if(LrelaxGradientOthers) then
          urf=urfGradientOthers
        else
          return
        endif
c
      endif
c      
      if(irFieldVariable.ne.0) then
        if(Variable.eq.rFieldName(irFieldVariable)) then
c
          if(LrelaxGradientrField(irFieldVariable)) then
            urf=urfGradientrField(irFieldVariable)
          else
            return
          endif
        endif
      endif
c      
      if(iScalarVariable.ne.0) then
        if(Variable.eq.ScalarName(iScalarVariable)) then
c
          if(LrelaxGradientScalar(iScalarVariable)) then
            urf=urfGradientScalar(iScalarVariable)
          else
            return
          endif
        endif
      endif
c
      urf1=1.-urf
      dfidxT=urf*dfidxT+urf1*dfidxTstar
      dfidyT=urf*dfidyT+urf1*dfidyTstar
      dfidzT=urf*dfidzT+urf1*dfidzTstar
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          BdfidxT(i,j)=dfidxT(k)
          BdfidyT(i,j)=dfidyT(k)
          BdfidzT(i,j)=dfidzT(k)
c
        enddo
      enddo
c
      return
      end