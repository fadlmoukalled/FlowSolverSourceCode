c
C#############################################################################################
c
      SUBROUTINE CalculateMaterialDerivative(Variable,FiTemp,
     *                            BFiTemp,FiTempOld,FiTempOldOld)
c
C#############################################################################################
c
      use User0, only:LUnsteady,MethodCalcGradientMomentum,
     *                nIterGradientMomentum,LimitGradientMomentum,
     *                LimitGradientMomentumMethod
      use variables1, only: uVelocity,vVelocity,wVelocity,
     *                      MaterialDerivative
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFacesMax
      use Geometry4, only: Volume
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i
      character*10 Variable
      double precision FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
      double precision, dimension(:) :: FiTemp
      double precision, dimension(:) :: FiTempOld
      double precision, dimension(:) :: FiTempOldOld
      double precision, dimension(:,:) :: BFiTemp
      double precision, dimension(:), allocatable :: FiTempGradx
      double precision, dimension(:), allocatable :: FiTempGrady
      double precision, dimension(:), allocatable :: FiTempGradz
      double precision, dimension(:), allocatable :: dFiTempdt
      double precision, dimension(:,:), allocatable :: BFiTempGradx
      double precision, dimension(:,:), allocatable :: BFiTempGrady
      double precision, dimension(:,:), allocatable :: BFiTempGradz
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
        SUBROUTINE ComputeTransientGradient(phi,phiOld,phiOldOld,dphidt)
c--------------------------------------------------------------
          double precision, dimension(:) :: phi,phiOld,phiOldOld,dphidt
c--------------------------------------------------------------
        end SUBROUTINE ComputeTransientGradient
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      allocate(FiTempGradx(NumberOfElements))
      allocate(FiTempGrady(NumberOfElements))
      allocate(FiTempGradz(NumberOfElements))
      allocate(dFiTempdt(NumberOfElements))
      allocate(BFiTempGradx(NumberOfBCSets,NBFacesMax))
      allocate(BFiTempGrady(NumberOfBCSets,NBFacesMax))
      allocate(BFiTempGradz(NumberOfBCSets,NBFacesMax))
c
      call Gradient(Variable,MethodCalcGradientMomentum,
     *           FiTemp,FiTempGradx,FiTempGrady,FiTempGradz,
     *              BFiTemp,BFiTempGradx,BFiTempGrady,BFiTempGradz,
     *              nIterGradientMomentum,LimitGradientMomentum,
     *                            LimitGradientMomentumMethod)
c
      MaterialDerivative=uVelocity*FiTempGradx+
     *                  vVelocity*FiTempGrady+wVelocity*FiTempGradz
c
      if(LUnsteady) then
c
        dFiTempdt=0.d0
        call ComputeTransientGradient
     *           (FiTemp,FiTempOld,FiTempOldOld,dFiTempdt)
        MaterialDerivative=MaterialDerivative+dFiTempdt
c
      endif
c
      deallocate(FiTempGradx)
      deallocate(FiTempGrady)
      deallocate(FiTempGradz)
      deallocate(dFiTempdt)
      deallocate(BFiTempGradx)
      deallocate(BFiTempGrady)
      deallocate(BFiTempGradz)
c
      return
      end