c
C#############################################################################################
c
      SUBROUTINE AssemblePressureWorkTerm
c
C#############################################################################################
c
      use User0, only:dt,LUnsteady,MethodCalcGradientContinuity,
     *                LSolveMomentum,nIterGradientContinuity,
     *            LimitGradientContinuity,LimitGradientContinuityMethod
      use Variables3, only: FluxTE
      use variables1, only: uVelocity,vVelocity,wVelocity,
     *                      PressureWork,Pressure,PressureOld,
     *                      PressureOldOld,PressGradx,PressGrady,
     *                      PressGradz,BPressure,BPressGradx,
     *                      BPressGrady,BPressGradz
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Transient1, only:dpdt
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i
      character*10 Variable
      double precision FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
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
      if(LSolveMomentum) then
c
        pressureWork=0.d0
c      
        if(LUnsteady) then
c
          dpdt=0.d0
c
          call ComputeTransientGradient
     *            (Pressure,PressureOld,PressureOldOld,dpdt)
c
        endif
c
        call Gradient(Variable,MethodCalcGradientContinuity,
     *           Pressure,PressGradx,PressGrady,PressGradz,BPressure,
     *            BPressGradx,BPressGrady,BPressGradz,
     *            nIterGradientContinuity,LimitGradientContinuity,
     *               LimitGradientContinuityMethod)
c
        do i=1,NumberOfElements
c
          PressureWork(i)=uVelocity(i)*PressGradx(i)+
     *            vVelocity(i)*PressGrady(i)+wVelocity(i)*PressGradz(i)
          if(LUnsteady) PressureWork(i)=PressureWork(i)+dpdt(i)
          PressureWork(i)=Volume(i)*PressureWork(i) 
c
          FluxTE(i)=FluxTE(i)-PressureWork(i)
c
        enddo
c
      endif
c
      return
      end
