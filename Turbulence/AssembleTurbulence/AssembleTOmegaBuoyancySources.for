c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentOmegaBuoyancySources
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TurbulentOmega,
     *                      uVelocity,vVelocity,wVelocity,
     *                      TurbulenceProductionB
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: GravityX,GravityY,GravityZ
      use Turbulence1, only: ModelNumber,alpha,ce3
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1,vpg,vng,gMagnitude
c********************************************************************************************
c
      CalculateTurbulentOmegaBuoyancySources: select case(ModelNumber)
c-----------------------------------------------------------------------------
        case(11) CalculateTurbulentOmegaBuoyancySources        !komega
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                                 dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+
     *                      term1*dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *           dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(12) CalculateTurbulentOmegaBuoyancySources      !komegaepsilon
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                                dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+
     *                      term1*dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *         dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(13) CalculateTurbulentOmegaBuoyancySources      !komegabsl
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                                dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+
     *                       term1*dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *         dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(14) CalculateTurbulentOmegaBuoyancySources      !komegasst
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                          dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+term1*
     *                             dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *         dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(15) CalculateTurbulentOmegaBuoyancySources     !sstgamaretheta
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                          dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+term1*
     *                             dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *         dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(16) CalculateTurbulentOmegaBuoyancySources      !komega2006
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                               dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+
     *                   term1*dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *          dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(17) CalculateTurbulentOmegaBuoyancySources     !komega2006lrn
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                               dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+
     *                   term1*dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *          dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c
c-----------------------------------------------------------------------------
        case(21) CalculateTurbulentOmegaBuoyancySources       !transitional k-kl-w model
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(23) CalculateTurbulentOmegaBuoyancySources       !sstgama
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
!            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
!            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
!     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
!            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
!     *                              wVelocity(i)**2-vpg**2)
!            ce3=dtanh(vpg/dmax1(vng,tiny))
            ce3=1.
c            
            term1=TurbulentOmega(i)*Volume(i)/
     *                          dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-(1.+alpha)*
     *                term1*ce3*dmax1(TurbulenceProductionB(i),0.)-
     *                       term1*dmax1(-TurbulenceProductionB(i),0.)
            term1=Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxCE(i)=FluxCE(i)+term1*
     *                             dmax1(TurbulenceProductionB(i),0.)
            FluxTE(i)=FluxTE(i)+term1*
     *         dmax1(TurbulenceProductionB(i)*TurbulentOmega(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTurbulentOmegaBuoyancySources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end