c
c#############################################################################################
c
      SUBROUTINE AssembleTurbulentKEBuoyancySources
c
C#############################################################################################
c
      use User0, only: BuoyancyModel
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Variables1, only: TurbulentKE,TempGradx,TempGrady,TempGradz,
     *                      TurbulenceProductionB
      use Variables3, only: FluxCE,FluxTE
      use PhysicalProperties1, only: GravityX,GravityY,GravityZ,
     *                               TurbulentViscosity,Density,
     *                               CoefficientOfThermalExpansion,
     *                               DensGradx,DensGrady,DensGradz
      use Turbulence1, only: SigT,ModelNumber
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1
c********************************************************************************************
c
      CalculateTKEBuoyancySources: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(1) CalculateTKEBuoyancySources           !kepsilon
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *            term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                       GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                          Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                       GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(2) CalculateTKEBuoyancySources           !kepsilonchien
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *                term1*(GravityX*TempGradx(i)+
     *                 GravityY*TempGrady(i)+GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                        GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                             Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(3) CalculateTKEBuoyancySources           !kepsilonsharma
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                             GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                        GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(4) CalculateTKEBuoyancySources           !kepsilonchc
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                            GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                              GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(5) CalculateTKEBuoyancySources           !kepsilonkasagi
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                              GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                       GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(6) CalculateTKEBuoyancySources           !kepsilontagawa
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                      GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                     GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(7) CalculateTKEBuoyancySources           !kepsilonhishida
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                       GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                         GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(8) CalculateTKEBuoyancySources           !kelambremhorst
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                    GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                        GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(9) CalculateTKEBuoyancySources           !kelambremhorstm
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                     GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                             GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(10) CalculateTKEBuoyancySources           !realizable
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *            term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                       GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                          Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                       GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(11) CalculateTKEBuoyancySources           !komega
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                         GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                           GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(12) CalculateTKEBuoyancySources           !komegaepsilon
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                          GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                          GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(13) CalculateTKEBuoyancySources           !komegabsl
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                         GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                             GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(14) CalculateTKEBuoyancySources           !komegasst
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                           GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                               GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(15) CalculateTKEBuoyancySources           !sstgamaretheta
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                           GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                               GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(16) CalculateTKEBuoyancySources           !komega2006
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                           GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *             term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                           GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                             Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(17) CalculateTKEBuoyancySources           !komega2006lrn
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                             GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *             term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                           GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                             Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(18) CalculateTKEBuoyancySources           !kklmodel
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(21) CalculateTKEBuoyancySources           !kklomega
c-----------------------------------------------------------------------------
c
          return
c 
c-----------------------------------------------------------------------------
        case(22) CalculateTKEBuoyancySources           !kepsilonrt
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *            term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                       GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                          Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                       GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(23) CalculateTKEBuoyancySources           !sstgama
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *           term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                           GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                               GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                              Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(25) CalculateTKEBuoyancySources           !kepsilonrng
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *            term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                       GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                          Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                       GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(26) CalculateTKEBuoyancySources           !kepsilonv2f
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *            term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                       GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                          Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                       GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------
        case(27) CalculateTKEBuoyancySources           !kepsilonzetaf
c-----------------------------------------------------------------------------
c
          if(BuoyancyModel.eq.'rhog') then
c
            do i=1,NumberOfElements    
c
              term1=CoefficientOfThermalExpansion*
     *                               TurbulentViscosity(i)/SigT
              TurbulenceProductionB(i)=
     *            term1*(GravityX*TempGradx(i)+GravityY*TempGrady(i)+
     *                       GravityZ*TempGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                          Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          elseif(BuoyancyModel.eq.'boussinesq') then
c
            do i=1,NumberOfElements    
c
              term1=-TurbulentViscosity(i)/(SigT*Density(i))
              TurbulenceProductionB(i)=
     *           term1*(GravityX*DensGradx(i)+GravityY*DensGrady(i)+
     *                       GravityZ*DensGradz(i))
c
              if(TurbulenceProductionB(i).gt.0.) then
c
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c
              else
c
                FluxCE(i)=FluxCE(i)-TurbulenceProductionB(i)*
     *                            Volume(i)/dmax1(TurbulentKE(i),tiny)
                FluxTE(i)=FluxTE(i)-TurbulenceProductionB(i)*Volume(i)
c        
              endif
c        
            enddo
c
          endif
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateTKEBuoyancySources 
c-----------------------------------------------------------------------------------------------      
c
      return
      end