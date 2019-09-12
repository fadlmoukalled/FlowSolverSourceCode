c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentReThetaSources
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TurbulentOmega,TGamma,TReTheta,
     *                      uVelocity,vVelocity,wVelocity,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density,Viscosity
      use Turbulence1, only: cot,ce2,Vorticity
      use Constants1, only: tiny,TwoThird
      use WallDistance1, only: WallDistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,iter
      double precision :: Uvel,Uvel2,Tfactor,ReW,Fwake,delta,Fthetat,
     *                    dUveldx,dUveldy,dUveldz,dUds,TU,Thetat,error,
     *                    ThetatOld,LambdaTheta,ReThetaTeqLocal,PThetat,
     *                    FThetatLocal,FPrimeThetatLocal,dLambdadTheta
c********************************************************************************************
      Interface
c********************************************************************************************
        FUNCTION ReThetaTeq(x,y)
!       x=Tu    y=LambdaTheta
c********************************************************************************************
          double precision :: x,y
          double precision :: ReThetaTeq
c*********************************************************************************************
        end FUNCTION ReThetaTeq
c*********************************************************************************************
        FUNCTION ReThetaTeqPrime(x,y,z)
!       x=Tu    y=LambdaTheta
c*********************************************************************************************
          double precision :: x,y,z
          double precision :: ReThetaTeqPrime
c*********************************************************************************************
        end FUNCTION ReThetaTeqPrime
c*********************************************************************************************
      end interface
c*********************************************************************************************
c
c--- Calculate turbulence production of Gamma
c
      do i=1,NumberOfElements    
c
        Uvel2=uVelocity(i)*uVelocity(i)+vVelocity(i)*
     *                    vVelocity(i)+wVelocity(i)*wVelocity(i)
        Uvel2=dmax1(Uvel2,tiny)
        UVel=dsqrt(Uvel2)
        Tfactor=500.*Viscosity(i)/(Density(i)*Uvel2)
        ReW=Density(i)*dmax1(TurbulentOmega(i),0.)*
     *         WallDistance(i)*WallDistance(i)/Viscosity(i)
        Fwake=dexp(-ReW*ReW/1.d10)
        delta=375.*Vorticity(i)*Viscosity(i)*TReTheta(i)*
     *        WallDistance(i)/(Density(i)*Uvel2)
        if(delta.eq.0.) delta=tiny       
        Fthetat=dmin1(dmax1(Fwake*dexp(-(wallDistance(i)/delta)**4),
     *                         1.-((ce2*TGamma(i)-1.)/(ce2-1.))**2),1.)
        dUveldx= uVelocity(i)*uVelGradx(i)+
     *               vVelocity(i)*vVelGradx(i)+
     *                    wVelocity(i)*wVelGradx(i)
        dUveldy= uVelocity(i)*uVelGrady(i)+
     *               vVelocity(i)*vVelGrady(i)+
     *                    wVelocity(i)*wVelGrady(i)
        dUveldz= uVelocity(i)*uVelGradz(i)+
     *               vVelocity(i)*vVelGradz(i)+
     *                    wVelocity(i)*wVelGradz(i)
        dUds=(uVelocity(i)*dUveldx+vVelocity(i)*dUveldy+
     *             wVelocity(i)*dUveldz)/Uvel2
        if(dUds.eq.0.) dUds=tiny       
        TU=100.*dsqrt(TwoThird*dmax1(TurbulentKE(i),0.))/Uvel
        TU=dmax1(TU,0.027)
c
c---- Calculate Reynolds theta equivalent
c
        Thetat=0.
        LambdaTheta=0.
        dLambdadTheta=0.
        error=1.
        iter=1
c
        do while(error.gt.1.e-6.and.iter.le.100)
c
          ReThetaTeqLocal=ReThetaTeq(TU,LambdaTheta)
          FThetatLocal=Density(i)*Uvel*Thetat/
     *                         Viscosity(i)-ReThetaTeqLocal
          FPrimeThetatLocal=Density(i)*Uvel/Viscosity(i)-
     *                   ReThetaTeqPrime(TU,LambdaTheta,dLambdadTheta)
          ThetatOld=Thetat
          Thetat=ThetatOld-FThetatLocal/FPrimeThetatLocal
          LambdaTheta=Density(i)*Thetat*Thetat*dUds/Viscosity(i)
          LambdaTheta=dmax1(dmin1(LambdaTheta,0.1),-0.1)
          ThetaT=dsqrt(LambdaTheta*Viscosity(i)/(dUds*Density(i)))
          dLambdadTheta=2.*Density(i)*Thetat*dUds/Viscosity(i)
c
          error=dabs(Thetat-ThetatOld)
          iter=iter+1
          if(iter.gt.100) then
            print*, 'Thetateq iterations diverged. Error= ',error
          endif
c
        enddo        
c
        PThetat=cot*Density(i)*Volume(i)/Tfactor       
c
        FluxCE(i)=FluxCE(i)+PThetat
        FluxTE(i)=FluxTE(i)-PThetat*
     *               (1.-Fthetat)*(ReThetaTeqLocal-TReTheta(i))
c
      enddo
c
      return
      end