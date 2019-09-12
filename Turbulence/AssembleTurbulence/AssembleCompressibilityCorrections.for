c
c#############################################################################################
c
      SUBROUTINE AssembleTOmegaCompressibilitySources
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentOmega
      use Variables3, only: FluxTE
      use PhysicalProperties1, only: Density
      use Turbulence1, only: cmu,ModelNumber,XiStarFmt,F1factor
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1,term2
c********************************************************************************************
c
      CalculateTOmegaCompressibilitySources: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(11) CalculateTOmegaCompressibilitySources        !komega
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *               TurbulentOmega(i)*TurbulentOmega(i)*XiStarFmt(i)
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(12) CalculateTOmegaCompressibilitySources        !komegaepsilon
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*TurbulentOmega(i)*
     *                 TurbulentOmega(i)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(13) CalculateTOmegaCompressibilitySources        !komegabsl
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*TurbulentOmega(i)*
     *                 TurbulentOmega(i)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(14) CalculateTOmegaCompressibilitySources        !komegasst
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*TurbulentOmega(i)*
     *                 TurbulentOmega(i)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(15) CalculateTOmegaCompressibilitySources        !sstgamaretheta
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*TurbulentOmega(i)*
     *                 TurbulentOmega(i)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(16) CalculateTOmegaCompressibilitySources        !komega2006
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *               TurbulentOmega(i)*TurbulentOmega(i)*XiStarFmt(i)
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(17) CalculateTOmegaCompressibilitySources        !komega2006lrn
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *               TurbulentOmega(i)*TurbulentOmega(i)*XiStarFmt(i)
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(23) CalculateTOmegaCompressibilitySources        !sstgama
c-----------------------------------------------------------------------------
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*TurbulentOmega(i)*
     *                 TurbulentOmega(i)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxTE(i)=FluxTE(i)-term1
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTOmegaCompressibilitySources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentKECompressibilitySources
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TurbulentED,TurbulentOmega
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: Density
      use Turbulence1, only: cmu,ModelNumber,XiStarFmt,F1factor
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1
c********************************************************************************************
c
      CalculateTKECompressibilitySources: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(1) CalculateTKECompressibilitySources      !kepsilon
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(2) CalculateTKECompressibilitySources      !kepsilonchien
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(3) CalculateTKECompressibilitySources      !kepsilonsharma
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(4) CalculateTKECompressibilitySources      !kepsilonchc
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(5) CalculateTKECompressibilitySources      !kepsilonkasagi
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(6) CalculateTKECompressibilitySources      !kepsilontagawa
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(7) CalculateTKECompressibilitySources      !kepsilonhishida
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(8) CalculateTKECompressibilitySources      !kelambremhorst
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(9) CalculateTKECompressibilitySources      !kelambremhorstm
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(10) CalculateTKECompressibilitySources      !realizable
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(11) CalculateTKECompressibilitySources      !komega
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *          dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(12) CalculateTKECompressibilitySources      !komegaepsilon
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *        dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(13) CalculateTKECompressibilitySources      !komegabsl
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *        dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(14) CalculateTKECompressibilitySources      !komegasst
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *        dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(15) CalculateTKECompressibilitySources      !sstgamaretheta
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *        dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(16) CalculateTKECompressibilitySources      !komega2006
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *          dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(17) CalculateTKECompressibilitySources      !komega2006lrn
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *          dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
          return
c
c-----------------------------------------------------------------------------
        case(18) CalculateTKECompressibilitySources      !kklmodel
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(21) CalculateTKECompressibilitySources      !kklomega
c-----------------------------------------------------------------------------
c
          return
c 
c-----------------------------------------------------------------------------
        case(22) CalculateTKECompressibilitySources      !kepsilonrt
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(23) CalculateTKECompressibilitySources      !sstgama
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements
c
            term1=cmu*Density(i)*Volume(i)*
     *        dmax1(TurbulentOmega(i),0.)*XiStarFmt(i)*(1.-F1factor(i))
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c-----------------------------------------------------------------------------
        case(25) CalculateTKECompressibilitySources      !kepsilonrng
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(26) CalculateTKECompressibilitySources      !kepsilonv2f
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------
        case(27) CalculateTKECompressibilitySources      !kepsilonzetaf
c-----------------------------------------------------------------------------
c
          call CalculateCompressibilityCorrection
c
          do i=1,NumberOfElements    
c
            term1=Density(i)*XiStarFmt(i)*Volume(i)*
     *          dmax1(TurbulentED(i),0.)/dmax1(TurbulentKE(i),tiny)
c
            FluxCE(i)=FluxCE(i)+term1
            FluxTE(i)=FluxTE(i)+term1*dmax1(TurbulentKE(i),tiny)
c        
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTKECompressibilitySources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end
c