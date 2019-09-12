c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentfRelaxationSources
c
C#############################################################################################
c
      use User0, only: TurbulenceModel,MethodCalcGradientTurbulentZeta,
     *                 nIterGradientTurbulentZeta,
     *                 LimitGradientTurbulentZeta,
     *                 LimitGradientTurbulentZetaMethod
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentV2,TurbulentKE,TfRelaxation,
     *                      TurbulentZeta,TurbulentZetaGradx,
     *                      TurbulentZetaGrady,TurbulentZetaGradz,
     *                      TKEGradx,TKEGrady,TKEGradz,
     *                      TurbulenceProduction,BTurbulentZetaGradx,
     *                      BTurbulentZetaGrady,BTurbulentZetaGradz
      use Variables3, only: FLuxCE,FLuxTE
      use Turbulence1, only: C1V2f,C2V2f,C1Zeta,C2Zeta,
     *                       TScale,LScale,
     *                       TurbulentZetaGrad2x,TurbulentZetaGrad2y,
     *                       TurbulentZetaGrad2z,TurbulentZetaGradxy,
     *                       TurbulentZetaGradxz,
     *                       BTurbulentZetaGrad2x,BTurbulentZetaGrad2y,
     *                       BTurbulentZetaGrad2z,BTurbulentZetaGradxy,
     *                       BTurbulentZetaGradxz
      use PhysicalProperties1, only: Density,Viscosity
      use Constants1, only: twothird,tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i
      double precision :: term1
c
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c********************************************************************************************
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
c********************************************************************************************
        end SUBROUTINE Gradient
c********************************************************************************************
      end interface
c********************************************************************************************
c
      if(TurbulenceModel.eq.'kepsilonv2f') then
c
        do i=1,NumberOfElements
c
                   
          term1=(C1V2f-1.)*twothird*Volume(i)/dmax1(TScale(i),tiny)
          FluxTE(i)=FluxTE(i)-term1
c
          term1=(C1V2f-1.)*dmax1(TurbulentV2(i),0.)*Volume(i)/
     *                         dmax1(TurbulentKE(i)*TScale(i),tiny)
          FluxCE(i)=FluxCE(i)+term1/dmax1(TfRelaxation(i),tiny)
          FluxTE(i)=FluxTE(i)+term1
c
          term1=C2V2f*TurbulenceProduction(i)*Volume(i)/
     *                    (Density(i)*dmax1(TurbulentKE(i),tiny)) 
          term1=term1+5.*dmax1(TurbulentV2(i),0.)*Volume(i)/
     *                         dmax1(TurbulentKE(i)*TScale(i),tiny)
          FluxTE(i)=FluxTE(i)-term1
c          
          FluxCE(i)=FluxCE(i)+Volume(i)
          FluxTE(i)=FluxTE(i)+dmax1(TfRelaxation(i),0.)*Volume(i)
c
        enddo
c
      elseif(TurbulenceModel.eq.'kepsilonzetaf') then
c
        Variable='tzeta'
        call Gradient(Variable,MethodCalcGradientTurbulentZeta,
     *    TurbulentZetaGradx,TurbulentZetaGrad2x,TurbulentZetaGradxy,
     *     TurbulentZetaGradxz,BTurbulentZetaGradx,BTurbulentZetaGrad2x,
     *        BTurbulentZetaGradxy,BTurbulentZetaGradxz,
     *        nIterGradientTurbulentZeta,LimitGradientTurbulentZeta,
     *                              LimitGradientTurbulentZetaMethod)
        call Gradient(Variable,MethodCalcGradientTurbulentZeta,
     *    TurbulentZetaGrady,TurbulentZetaGradxy,TurbulentZetaGrad2y,
     *     TurbulentZetaGradxz,BTurbulentZetaGrady,BTurbulentZetaGradxy,
     *        BTurbulentZetaGrad2y,BTurbulentZetaGradxz,
     *        nIterGradientTurbulentZeta,LimitGradientTurbulentZeta,
     *                              LimitGradientTurbulentZetaMethod)
        call Gradient(Variable,MethodCalcGradientTurbulentZeta,
     *    TurbulentZetaGradz,TurbulentZetaGradxy,TurbulentZetaGradxz,
     *     TurbulentZetaGrad2z,BTurbulentZetaGradz,BTurbulentZetaGradxy,
     *        BTurbulentZetaGradxz,BTurbulentZetaGrad2z,
     *        nIterGradientTurbulentZeta,LimitGradientTurbulentZeta,
     *                              LimitGradientTurbulentZetaMethod)
c
        do i=1,NumberOfElements
c
          term1=(C1Zeta-1.)*twothird*Volume(i)/dmax1(TScale(i),tiny)
          FluxTE(i)=FluxTE(i)-term1
c
          term1=(C1Zeta-1.)*dmax1(TurbulentZeta(i),0.)*
     *                         Volume(i)/dmax1(TScale(i),tiny)
          FluxCE(i)=FluxCE(i)+term1/dmax1(TfRelaxation(i),tiny)
          FluxTE(i)=FluxTE(i)+term1
c
          term1=C2Zeta*TurbulenceProduction(i)*Volume(i)/
     *                    (Density(i)*dmax1(TurbulentKE(i),tiny)) 
          FluxTE(i)=FluxTE(i)-term1
c          
          FluxCE(i)=FluxCE(i)+Volume(i)
          FluxTE(i)=FluxTE(i)+dmax1(TfRelaxation(i),0.)*Volume(i)
c
          term1=2.*Viscosity(i)/(Density(i)*dmax1(TurbulentKE(i),tiny))
          term1=term1*Volume(i)*(
     *           TurbulentZetaGradx(i)*TKEGradx(i)+
     *           TurbulentZetaGrady(i)*TKEGrady(i)+
     *           TurbulentZetaGradz(i)*TKEGradz(i))
          FLuxCE(i)=FLuxCE(i)-
     *              dmin1(term1/dmax1(TfRelaxation(i),tiny),0.)
          FLuxTE(i)=FLuxTE(i)-term1
c
          term1=(Viscosity(i)*Volume(i)/Density(i))*
     *        (TurbulentZetaGrad2x(i)+TurbulentZetaGrad2y(i)+
     *                                     TurbulentZetaGrad2z(i))  
          FLuxCE(i)=FLuxCE(i)-
     *              dmin1(term1/dmax1(TfRelaxation(i),tiny),0.)
          FLuxTE(i)=FLuxTE(i)-term1
c
        enddo
c
      endif
c
      return
      end