c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentEDBuoyancySources
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TurbulentED,
     *                      uVelocity,vVelocity,wVelocity,
     *                      TurbulenceProductionB
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: GravityX,GravityY,GravityZ
      use Turbulence1, only: ce1,ce3,ModelNumber
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1,vpg,vng,gMagnitude
c********************************************************************************************
c
      CalculateTurbulentEDBuoyancySources: select case(ModelNumber)
c
c-------------------------------------------------------------------------------------------
        case(1) CalculateTurbulentEDBuoyancySources      !kepsilon
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(2) CalculateTurbulentEDBuoyancySources      !kepsilonchien
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(3) CalculateTurbulentEDBuoyancySources      !kepsilonsharma
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(4) CalculateTurbulentEDBuoyancySources      !kepsilonchc
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(5) CalculateTurbulentEDBuoyancySources      !kepsilonkasagi
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(6) CalculateTurbulentEDBuoyancySources      !kepsilontagawa
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(7) CalculateTurbulentEDBuoyancySources      !kepsilonhishida
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(8) CalculateTurbulentEDBuoyancySources      !kelambremhorst
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(9) CalculateTurbulentEDBuoyancySources      !kelambremhorstm
c-------------------------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(10) CalculateTurbulentEDBuoyancySources      !realizable
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(22) CalculateTurbulentEDBuoyancySources      !K-Epsilon-Rt
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(25) CalculateTurbulentEDBuoyancySources      !kepsilonrng
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(26) CalculateTurbulentEDBuoyancySources      !kepsilonv2f
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-------------------------------------------------------------------------------------------
        case(27) CalculateTurbulentEDBuoyancySources      !kepsilonzetaf
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements    
c
            gMagnitude=dsqrt(gravityX**2+gravityY**2+gravityZ**2)
            vpg=dabs(gravityX*uVelocity(i)+gravityY*vVelocity(i)+
     *                    gravityZ*wVelocity(i))/dmax1(gMagnitude,tiny)
            vng=dsqrt(uVelocity(i)**2+vVelocity(i)**2+
     *                              wVelocity(i)**2-vpg**2)
            ce3=dtanh(vpg/dmax1(vng,tiny))
c            
            term1=TurbulentED(i)*Volume(i)/dmax1(TurbulentKE(i),tiny)
            FluxTE(i)=FluxTE(i)-
     *                  ce1*term1*ce3*dmax1(TurbulenceProductionB(i),0.)
c        
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTurbulentEDBuoyancySources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end