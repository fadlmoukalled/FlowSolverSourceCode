c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentGammaSources
c
C#############################################################################################
c
      use User0, only: MethodCalcGradientMomentum,nIterGradientMomentum,
     *                 LimitGradientMomentum,LimitGradientMomentumMethod
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume,BFaceAreanx,BFaceAreany,BFaceAreanz
      use Variables1, only: TurbulentKE,TurbulentOmega,TGamma,TReTheta,
     *                      NormalVelocity,BNormalVelocity,
     *                      NormalVelocityGradx,BNormalVelocityGradx,
     *                      NormalVelocityGrady,BNormalVelocityGrady,
     *                      NormalVelocityGradz,BNormalVelocityGradz
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density,Viscosity
      use Turbulence1, only: ca1,ca2,ce1,ce2,StrainRate,Vorticity,
     *                       ModelNumber,Cpg1,Cpg2,Cpg3,Cpg1lim,Cpg2lim,
     *                       Ctu1,Ctu2,Ctu3,FlengthGama
      use Constants1, only: tiny,twothird
      use WallDistance1, only: WallDistance,iTau,jTau
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      character*10 :: Variable
      double precision :: ReV,Fonset,Fonset1,Fonset2,Fonset3,
     *                    Fturb,RT,ReW,Fsublayer,FLength,
     *                    PGamma,EGamma,nx,ny,nz,dvdy,Lambdatheta,
     *                    Fpg,tke1,ted1,TuL,ReThetaC1
c********************************************************************************************
      Interface
c********************************************************************************************
        FUNCTION ReThetaC(x)
c*********************************************************************************************
          double precision :: x
          double precision :: ReThetaC
c*********************************************************************************************
        end FUNCTION ReThetaC
c*********************************************************************************************
        FUNCTION Flength1(x)
c*********************************************************************************************
          double precision :: x
          double precision :: FLength1
c*********************************************************************************************
        end FUNCTION FLength1
c*********************************************************************************************
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
      end interface
c--------------------------------------------------------------------------------
c
c-----------------------------------------------------------------------------
      CalculateTGammaSources: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(15) CalculateTGammaSources           !sstgamaretheta
c-----------------------------------------------------------------------------
c
c--- Calculate turbulence production of Gamma
c
          do i=1,NumberOfElements    
c
            ReV=Density(i)*StrainRate(i)*WallDistance(i)*
     *                WallDistance(i)/Viscosity(i)
            Fonset1=ReV/(2.193*dmax1(ReThetaC(TReTheta(i)),tiny))
            Fonset2=
     *        dmin1(dmax1(Fonset1,Fonset1*Fonset1*Fonset1*Fonset1),2.)
            RT=Density(i)*dmax1(TurbulentKE(i),0.)/
     *            (Viscosity(i)*dmax1(TurbulentOmega(i),tiny))
            Fonset3=dmax1(1.-RT*RT*RT/15.625,0.)
            Fonset=dmax1(Fonset2-Fonset3,0.)
            Fturb=dexp(-RT*RT*RT*RT/256.)
            ReW=Density(i)*dmax1(TurbulentOmega(i),0.)*
     *         WallDistance(i)*WallDistance(i)/Viscosity(i)
            Fsublayer=dexp(-ReW*ReW/40000.)
            FLength=FLength1(TReTheta(i))*(1.-Fsublayer)+40.*Fsublayer
c
            PGamma=FLength*ca1*Density(i)*StrainRate(i)*
     *                 dsqrt(dmax1(TGamma(i)*Fonset,0.))*Volume(i)
c
            FluxCE(i)=FluxCE(i)+PGamma*ce1
            FluxTE(i)=FluxTE(i)+PGamma*(ce1*dmax1(TGamma(i),0.)-1.)
c
c--- Calculate turbulence destruction of Gamma
c
            EGamma=ca2*Density(i)*Vorticity(i)*
     *              Fturb*dmax1(TGamma(i),0.)*Volume(i)
            FluxCE(i)=FluxCE(i)+EGamma*ce2
            FluxTE(i)=FluxTE(i)+
     *             EGamma*(ce2*dmax1(TGamma(i),0.)-1.)
c
          enddo
c
c-----------------------------------------------------------------------------
        case(23) CalculateTGammaSources           !sstgamaretheta
c-----------------------------------------------------------------------------
c
          call CalculateNormalVelocity
          variable='nvel'
          call Gradient(Variable,MethodCalcGradientMomentum,
     *      NormalVelocity,NormalVelocityGradx,NormalVelocityGrady,
     *      NormalVelocityGradz,BNormalVelocity,BNormalVelocityGradx,
     *       BNormalVelocityGrady,BNormalVelocityGradz,
     *        nIterGradientMomentum,LimitGradientMomentum,
     *                              LimitGradientMomentumMethod)
c
c--- Calculate turbulence production of Gamma
c
          do i=1,NumberOfElements    
c
            j=itau(i)
            k=jtau(i) 
c
            nx=BFaceAreanx(j,k)
            ny=BFaceAreany(j,k)
            nz=BFaceAreanz(j,k)
c
            dvdy=NormalVelocityGradx(i)*nx+
     *            NormalVelocityGrady(i)*ny+NormalVelocityGradz(i)*nz
            Lambdatheta=0.0128+(-7.57d-3)*Density(i)*dvdy*
     *                   WallDistance(i)*WallDistance(i)/Viscosity(i)  
            Lambdatheta=dmin1(dmax1(Lambdatheta,-1.),1.)
c
            if(Lambdatheta.ge.0.) then
              Fpg=dmin1(1.+Cpg1*Lambdatheta,Cpg1lim)
            else
              Fpg=dmin1(1.+Cpg2*Lambdatheta+
     *                 Cpg3*dmin1(Lambdatheta+0.0681,0.),Cpg2lim)
            endif
            Fpg=dmax1(Fpg,0.)
c            
            tke1=dmax1(TurbulentKE(i),0.)
            ted1=dmax1(TurbulentOmega(i),tiny)
            TuL=dmin1(100.*dsqrt(twothird*tke1)/
     *                        (WallDistance(i)*ted1),100.)
            ReThetaC1=Ctu1+Ctu2*dexp(-Ctu3*TuL*Fpg)        
c                     
            ReV=Density(i)*StrainRate(i)*WallDistance(i)*
     *                WallDistance(i)/Viscosity(i)
            Fonset1=ReV/(2.2*dmax1(ReThetaC1,tiny))
            Fonset2=dmin1(Fonset1,2.)
            RT=Density(i)*tke1/(Viscosity(i)*ted1)
            Fonset3=dmax1(1.-RT*RT*RT/42.875,0.)
            Fonset=dmax1(Fonset2-Fonset3,0.)
            Fturb=dexp(-RT*RT*RT*RT/16.)
c
            PGamma=FlengthGama*Density(i)*StrainRate(i)*
     *                   dmax1(TGamma(i),0.)*Fonset*Volume(i)
c
            FluxCE(i)=FluxCE(i)+PGamma
            FluxTE(i)=FluxTE(i)+PGamma*(dmax1(TGamma(i),0.)-1.)
c
c--- Calculate turbulence destruction of Gamma
c
            EGamma=ca2*Density(i)*Vorticity(i)*
     *              Fturb*dmax1(TGamma(i),0.)*Volume(i)
            FluxCE(i)=FluxCE(i)+EGamma*ce2
            FluxTE(i)=FluxTE(i)+
     *             EGamma*(ce2*dmax1(TGamma(i),0.)-1.)
c
          enddo
c
c-----------------------------------------------------------------------------
      end select CalculateTGammaSources 
c-----------------------------------------------------------------------------
c
      return
      end