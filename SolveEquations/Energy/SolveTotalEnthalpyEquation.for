c
c#############################################################################################
c
      SUBROUTINE SolveTotalEnthalpy
c
C#############################################################################################
c
      use User0
      use Variables1
      use Variables4
      use PhysicalProperties1
      use variables2, only: dphi
      use BoundaryConditions1
c********************************************************************************************
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
        SUBROUTINE AssembleDiffusionTerm(Variable,
     *       Gam,BGam,FiT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *                           dfidxfT,dfidyfT,dfidzfT)
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
          double precision, dimension(:) :: Gam
          double precision, dimension(:,:) :: Bgam
c--------------------------------------------------------------
        end SUBROUTINE AssembleDiffusionTerm
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
        SUBROUTINE AssembleConvectionTerm(Variable,Bleed,
     *        ConvectionScheme,HRFramework,FiT,BFiT,dfidxT,
     *              dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------
          character*10 Variable
          character*20 ConvectionScheme
          character*4 HRFramework
          double precision :: Bleed
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------
        end SUBROUTINE AssembleConvectionTerm
c--------------------------------------------------------------
        SUBROUTINE AssembleSourceTerm(Variable,FiT,Sc,Sb)
          character*10 Variable
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: Sc
          double precision, dimension(:) :: Sb
c--------------------------------------------------------------
        end SUBROUTINE AssembleSourceTerm
c--------------------------------------------------------------
        SUBROUTINE AssemblePointSourceTerm(FiT,ScPointSource,
     *                                                  SbPointSource)
c--------------------------------------------------------------
          double precision, dimension(:)  :: FiT
          double precision, dimension(:)  :: ScPointSource
          double precision, dimension (:) :: SbPointSource
c--------------------------------------------------------------
        end SUBROUTINE AssemblePointSourceTerm
c--------------------------------------------------------------
      SUBROUTINE AssembleTransientTerm(Variable,FiT,FiTold,FiToldold)
c--------------------------------------------------------------
        character*10 Variable
        double precision, dimension(:) :: FiT
        double precision, dimension(:) :: FiTold
        double precision, dimension(:) :: FiToldold
c--------------------------------------------------------------
      end SUBROUTINE AssembleTransientTerm
c--------------------------------------------------------------
        SUBROUTINE updateValuesFromBoundaryConditions
     *               (Variable,Gam,BGam,FiT,BFiT,BdfidxT,
     *                  BdfidyT,BdfidzT,dfidxfT,dfidyfT,dfidzfT)
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
          double precision, dimension(:) :: Gam
          double precision, dimension(:,:) :: Bgam
c--------------------------------------------------------------
        end SUBROUTINE updateValuesFromBoundaryConditions
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
c--------------------------------------------------------------
      end interface
C*********************************************************************************************
c
      Variable='htotal'
c
      call InitializeFluxes
      call InitializeCoefficients
c
c--- Calculate Htotal at boundary and Update static temperature from total temperature
c
      call CalculateBoundaryHtotal
c
c---- Calculate Gradients of the phi variable
c
      call Gradient(Variable,MethodCalcGradientEnergy,
     *    Htotal,HtotalGradx,HtotalGrady,HtotalGradz,BHtotal,
     *      BHtotalGradx,BHtotalGrady,BHtotalGradz,nIterGradientEnergy,
     *             LimitGradientEnergy,LimitGradientEnergyMethod)
c
      call InterpolateGradientToFace(GradientInterpolationSchemeEnergy,
     *            Htotal,HtotalGradx,HtotalGrady,HtotalGradz,
     *                      HtotalGradfx,HtotalGradfy,HtotalGradfz)
c
c---- Assemble terms
c
      if(.not.Linviscid) then
c
        Variable='temp'
        call Gradient(Variable,MethodCalcGradientEnergy,
     *      Temperature,TempGradx,TempGrady,TempGradz,BTemperature,
     *            BTempGradx,BTempGrady,BTempGradz,nIterGradientEnergy,
     *             LimitGradientEnergy,LimitGradientEnergyMethod)
c
        call InterpolateGradientToFace(
     *          GradientInterpolationSchemeEnergy,
     *            Temperature,TempGradx,TempGrady,TempGradz,
     *                           TempGradfx,TempGradfy,TempGradfz)
c
        if(LanisotropicDiffusion) then
c      
          call AssembleAnisotropicDiffusionTerm(Variable,Temperature,
     *           BTemperature,BTempGradx,BTempGrady,BTempGradz,
     *                         TempGradfx,TempGradfy,TempGradfz)
c
        else
c      
          call CalculateEffectiveDiffusionCoefficient(Variable)
          call AssembleDiffusionTerm(Variable,eDiffCoefficient,      
     *           BeDiffCoefficient,Temperature,BTemperature,
     *               BTempGradx,BTempGrady,BTempGradz,
     *                      TempGradfx,TempGradfy,TempGradfz)
c
        endif
c
        Variable='htotal'
c
      endif
c
      if(LConvectScalar.or.LSolveMomentum)
     *   call AssembleConvectionTerm(Variable,BleedEnergy,
     *        ConvectionSchemeEnergy,HRFrameworkEnergy,Htotal,
     *             BHtotal,HtotalGradx,HtotalGrady,HtotalGradz,
     *                       BHtotalGradx,BHtotalGrady,BHtotalGradz)
c
      if(Lcompressible.and.LSolveMomentum) then
c
        if(LUnsteady) call AssembleTransientPressureGradientTerm
c
        if(.not.Linviscid) then
          call AssembleViscousWorkTerm
        endif
c
      endif
c
      if(LUnsteady) then
c
        call AssembleTransientTerm
     *          (Variable,Htotal,HtotalOld,HtotalOldOld)
c
      endif
c
      if(LBuoyancy) call AssembleBuoyancyTerm(Variable) 
c
      call AssembleSourceTerm(Variable,Htotal,ScEnergy,SbEnergy)
      call AssemblePointSourceTerm(Htotal,ScPointSourceEnergy,
     *                                           SbPointSourceEnergy)
c
      if(LFreeSurfaceFlow.and.LSurfaceTension) 
     *           call AssembleSurfaceTensionSourceTerm(Variable)
c
c---- Underrelax using false transient
c
      if(LFalseTransientEnergy) 
     *        call AssembleFalseTransient(Variable,FalsedtEnergy)
c
c---- Assemble global matrix
c
      call AssembleGlobalMatrixFaceFluxes
      call AssembleGlobalMatrixElementFluxes
c
c---- Underrelax equation using a factor
c
      if(LRelaxEnergy) call UnderRelaxEquation(urfEnergy)
c
c---- Solve equations
c
	call SolveEquation(Variable,Htotal,RRFEnergy,
     *                  ASIterEnergy,ASSolverEnergy,LMultigridEnergy)
c
      call updateTemperatureFromHtotal
      if(LimitTemperature) call BoundTemperature
c
      call updateValuesFromBoundaryConditions(Variable,      
     *      eDiffCoefficient,BeDiffCoefficient,Htotal,
     *         BHtotal,BHtotalGradx,BHtotalGrady,BHtotalGradz,
     *                     HtotalGradfx,HtotalGradfy,HtotalGradfz)
c
C----   If the flow is compressible perform the following: 
c
      if(Lcompressible) then
c
        call Calculatedrhodp  
        call CalculateDensity 
        call UpdateBoundaryTemperature
        call CalculateFaceDensity  
        call CalculateMachNumber  
c
      endif
c
C----   Update the buoyancy term 
c
      if(LsolveMomentum.and.LBuoyancy) call reDistributeTheBuoyancyTerm
c
	return
	end
