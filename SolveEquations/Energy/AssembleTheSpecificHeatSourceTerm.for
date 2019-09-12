c
C#############################################################################################
c
      SUBROUTINE AssembleSpecificHeatSourceTerm
c
C#############################################################################################
c
      use User0, only:dt,LUnsteady,MethodCalcGradientEnergy,
     *                LSolveMomentum,LConvectScalar,nIterGradientEnergy,
     *                LimitGradientEnergy,LimitGradientEnergyMethod
      use Variables3, only: FluxCE,FluxTE
      use variables1, only: uVelocity,vVelocity,wVelocity,Temperature
      use PhysicalProperties1, only: SpecificHeat,SHeatGradx,
     *                                SHeatGrady,SHeatGradz,
     *                                SpecificHeatOld,
     *                                SpecificHeatOldOld,BSpecificHeat,
     *                                BSHeatGradx,BSHeatGrady,
     *                                BSHeatGradz,Density,cpterm
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Transient1, only:dCpdt
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
      cpterm=0.d0
c      
      if(LUnsteady) then
c
        dCpdt=0.d0
        call ComputeTransientGradient
     *    (SpecificHeat,SpecificHeatOld,SpecificHeatOldOld,dCpdt)
c
      endif
c      
      if(LSolveMomentum.or.LConvectScalar) then
c
        call Gradient(Variable,MethodCalcGradientEnergy,
     *           SpecificHeat,SHeatGradx,SHeatGrady,SHeatGradz,
     *              BSpecificHeat,BSHeatGradx,BSHeatGrady,BSHeatGradz,
     *              nIterGradientEnergy,LimitGradientEnergy,
     *                            LimitGradientEnergyMethod)
c
        do i=1,NumberOfElements
c
          cpterm(i)=uVelocity(i)*SHeatGradx(i)+
     *        vVelocity(i)*SHeatGrady(i)+wVelocity(i)*SHeatGradz(i)
c
          if(LUnsteady) cpterm(i)=cpterm(i)+dCpdt(i)
c
          cpterm(i)=Density(i)*Volume(i)*cpterm(i)
c
          if(cpterm(i).lt.0.) then
c          
            FluxCE(i)=FluxCE(i)-cpterm(i)
            FluxTE(i)=FluxTE(i)-cpterm(i)*Temperature(i)
c
          else
c
            FluxTE(i)=FluxTE(i)-cpterm(i)*Temperature(i)
c
          endif        
c
        enddo
c
      endif
c
      return
      end
