c
c#############################################################################################
c
      SUBROUTINE UpdateBoundaryPressureCorrection
c
c#############################################################################################
c
      use User0, only: Lcompressible,MethodDecomposeS,LFreeSurfaceFlow
      use Variables1, only: PressureC,BPressureC,mdot,Bmdot,
     *                      Du2Velocity,Dv2Velocity,Dw2Velocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      uVelocity,vVelocity,wVelocity,
     *                      Temperature,BTemperature,Pressure,BPressure
      use Variables2, only: FluxCf,FluxFf
      use Geometry1, only: NumberOfBCSets
      use BoundaryConditions2
      use BoundaryConditions1, only: BoundaryType,outletTypeC
      use Geometry3, only: NIFaces,NBFaces,NBFaceOwner
      use PhysicalProperties1, only: BDensity,Density,
     *                               BSpecificHeat,Bdrhodp,RGas
      use Geometry4, only: BFaceAreax,BFaceAreay,BFaceAreaz,
     *                     BDistanceCFux,BDistanceCFuy,BDistanceCFuz,
     *                     BDistanceCF,xc,yc,zc,BFaceArea,
     *                     BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
      use Constants1, only: tiny
      use ArteryResistance1, only: ArteryResistance,geoDiffSum,
     *                             geoDiffB,PressureCResistance
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j,j1,j2,j3
      double precision :: rhof,sfx,sfy,sfz,ex,ey,ez,cf,DuSf,DvSf,DwSf,
     *                    eDuSf,eDvSf,eDwSf,DuEf,DvEf,DwEf,geoDiff,
     *                    Magnitude,velocity2,TotalTemp,dpdub,c1
c
      double precision :: nx,ny,nz,vdotInfinity,vdotin,cinfinity,cin,
     *                    vdotb,cb,Tb,rhoInfinity,sbEntropy,rhob,
     *                    vel2,xdir,ydir,zdir,FluxCfLocal1,ratio,
     *                    GFactCF,FluxCfLocal,c2,xF1,yF1,zF1,distance1,
     *                    distance2,uF1,vF1,Ub1,Ub2,term,term4,
     *                    dotproduct
c********************************************************************************************
c
      if(.not.Lcompressible) then
c
c----------------------------------------------------------------------
c
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i1=IinletSpecifiedVelocityOwner(i)
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
c
          BPressureC(i2,i3)=0.
c
        enddo
c----------------------------------------------------------------------
c
        do i=1,IinletSpecifiedStagnationPressure
c
          i1=IinletSpecifiedStagnationPressureOwner(i)
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
c
          geoDiff=FluxFf(i4)
          if(LFreeSurfaceFlow) geoDiff=rhof*FluxFf(i4)
c
          velocity2=BuVelocity(i2,i3)**2+
     *                 BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2
          BPressureC(i2,i3)=-rhof*velocity2*geoDiff*PressureC(i1)/
     *                             (Bmdot(i2,i3)-geoDiff*rhof*velocity2)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedStaticPressure
c
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
c
          BPressureC(i2,i3)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
c
          BPressureC(i2,i3)=0.
c
        enddo
c----------------------------------------------------------------------
        PressureCResistance=0.
        do i=1,IoutletSpecifiedResistance
c
          i1=IoutletSpecifiedResistanceOwner(i)
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
c
          PressureCResistance(i2)=PressureCResistance(i2)+
     *                             geoDiffB(i2,i3)*PressureC(i1)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'outlet'.and.
     *            outletTypeC(i).eq.'specifiedresistance') then
            PressureCResistance(i)=
     *          ArteryResistance(i)*PressureCResistance(i)/
     *                     (1.+ArteryResistance(i)*geoDiffSum(i))
          endif
c
        enddo        
c 
        do i=1,IoutletSpecifiedResistance
c
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
c
          BPressureC(i2,i3)=PressureCResistance(i2)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedVelocity
c
          i1=IoutletSpecifiedVelocityOwner(i)
          i2=IoutletSpecifiedVelocityNumberOfBCSets(i)
          i3=IoutletSpecifiedVelocityNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletFullyDeveloped
c
          i1=IoutletFullyDevelopedOwner(i)
          i2=IoutletFullyDevelopedNumberOfBCSets(i)
          i3=IoutletFullyDevelopedNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
c
      if(LRotationalPeriodicity) then

        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
            yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
            zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            BPressureC(i2,i3)=GFactCF*PressureC(i1)+
     *                         (1.-GFactCF)*PressureC(j1)
c
            BPressureC(j2,j3)=BPressureC(i2,i3)
c
          endif
c
        enddo
c
      elseif(LTranslationalPeriodicity) then

        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=xc(j1)+xTranslation(j2)
            yF1=yc(j1)+yTranslation(j2)
            zF1=zc(j1)+zTranslation(j2)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            BPressureC(i2,i3)=GFactCF*PressureC(i1)+
     *                         (1.-GFactCF)*PressureC(j1)
c
            BPressureC(j2,j3)=BPressureC(i2,i3)
c
          endif
c
        enddo
c
      endif
c
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
c
      elseif(Lcompressible) then
c
c----------------------------------------------------------------------
        do i=1,IWallSlip
c
          i1=IWallSlipOwner(i)
          i2=IWallSlipNumberOfBCSets(i)
          i3=IWallSlipNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iinletsupersonic
c
          i1=IinletsupersonicOwner(i)
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
c
          BPressureC(i2,i3)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i1=IinletSpecifiedVelocityOwner(i)
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=rhof*dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=rhof*dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=rhof*dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          BPressureC(i2,i3)=geoDiff*PressureC(i1)/
     *            (geoDiff-Bmdot(i2,i3)*Bdrhodp(i2,i3)/rhof)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
c
          BPressureC(i2,i3)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationPressure
c
          i1=IinletSpecifiedStagnationPressureOwner(i)
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
c
          geoDiff=FluxFf(i4)/rhof
          velocity2=BuVelocity(i2,i3)**2+
     *           BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2
          TotalTemp=BTemperature(i2,i3)+
     *                velocity2/(2.*BSpecificHeat(i2,i3))
          dpdub=-BPressure(i2,i3)*Bmdot(i2,i3)/(rhof*TotalTemp*RGas)
          c1=dpdub*geoDiff/(1.+dpdub*geoDiff)
          BPressureC(i2,i3)=c1*PressureC(i1)
c
        enddo
c
c----------------------------------------------------------------------
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
          FluxCfLocal1=FluxFf(i4)/BFaceArea(i2,i3)
c
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
          ratio=BSpecificHeat(i2,i3)/(BSpecificHeat(i2,i3)-RGas)

          vdotInfinity=uVelocityFarField(i2)*nx+
     *           vVelocityFarField(i2)*ny+wVelocityFarField(i2)*nz
          vdotin=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
          cinfinity=dsqrt(ratio*RGas*TemperatureFarField(i2))
          cin=dsqrt(ratio*RGas*Temperature(i1))
c
          vdotb=0.5*(vdotInfinity+vdotin)-(cinfinity-cin)/(ratio-1.)
          cb=-(vdotInfinity-vdotin)*(ratio-1.)/4.+0.5*(cinfinity+cin)
          Tb=(cb*cb)/(ratio*RGas)
          rhoInfinity=pressureFarField(i2)/
     *                  (RGas*TemperatureFarField(i2))
c
          if(Bmdot(i2,i3).lt.0.) then
c
            sbEntropy=cinfinity**2/(ratio*(rhoInfinity**(ratio-1)))
c
          else
c
            sbEntropy=cin**2/(ratio*(Density(i1)**(ratio-1)))
c
          endif           
c
          rhob=(cb*cb/(ratio*sbEntropy))**(1./(ratio-1.))
c
          FluxCfLocal=rhob*(1.+vdotb/cb)*BFaceArea(i2,i3)*
     *     FluxCfLocal1/(1.+FluxCfLocal1*BPressure(i2,i3)/cb)
c          FluxCfLocal=rhob*(1.+vdotb/cb)*BFaceArea(i2,i3)*
c     *         (FluxCfLocal1/2.+cin/(2.*ratio*Pressure(i1)))/
c     *         (1+0.5*FluxCfLocal1*ratio*BPressure(i2,i3)/cb)
c
c          c1=FluxFf(i4)-FluxCfLocal
c          c2=FluxFf(i4)-Bmdot(i2,i3)*Bdrhodp(i2,i3)/rhof
c          
          c1=FluxCfLocal1*BPressure(i2,i3)/cb
          c2=1.+c1
c
          BPressureC(i2,i3)=(c1/c2)*PressureC(i1)
c
        enddo
c
c----------------------------------------------------------------------
        do i=1,IoutletTransmissive
c
          i1=IoutletTransmissiveOwner(i)
          i2=IoutletTransmissiveNumberOfBCSets(i)
          i3=IoutletTransmissiveNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
          FluxCfLocal1=FluxFf(i4)/BFaceArea(i2,i3)
c
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
          ratio=BSpecificHeat(i2,i3)/(BSpecificHeat(i2,i3)-RGas)

          vdotin=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
          cin=dsqrt(ratio*RGas*Temperature(i1))
c
          vdotb=vdotin
          cb=cin
          Tb=(cb*cb)/(ratio*RGas)
c
          sbEntropy=cin**2/(ratio*(Density(i1)**(ratio-1)))
c
          rhob=(cb*cb/(ratio*sbEntropy))**(1./(ratio-1.))
c
          FluxCfLocal=rhob*(1.+vdotb/cb)*BFaceArea(i2,i3)*
     *         FluxCfLocal1/(1.+FluxCfLocal1*BPressure(i2,i3)/cb)
c          FluxCfLocal=rhob*FluxCfLocal1*(1.+vdotb/cb)*BFaceArea(i2,i3)/
c     *                      (1.+FluxCfLocal1*cb/ratio)
c
          c1=FluxFf(i4)-FluxCfLocal
          c2=FluxFf(i4)-Bmdot(i2,i3)*Bdrhodp(i2,i3)/rhof         
          BPressureC(i2,i3)=(c1/c2)*PressureC(i1)
c
        enddo
c
c----------------------------------------------------------------------
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i1=IoutletspecifiedVelocityOwner(i)
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedStaticPressure
c
          i1=IoutletSpecifiedStaticPressureOwner(i)
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
c
          BPressureC(i2,i3)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i1=IoutletSpecifiedAverageStaticPressureOwner(i)
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
c
          BPressureC(i2,i3)=0.
c
        enddo
c----------------------------------------------------------------------
        PressureCResistance=0.
        do i=1,IoutletSpecifiedResistance
c
          i1=IoutletSpecifiedResistanceOwner(i)
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
c
          PressureCResistance(i2)=PressureCResistance(i2)+
     *                             geoDiffB(i2,i3)*PressureC(i1)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'outlet'.and.
     *            outletTypeC(i).eq.'specifiedresistance') then
            PressureCResistance(i)=
     *          ArteryResistance(i)*PressureCResistance(i)/
     *                     (1.+ArteryResistance(i)*geoDiffSum(i))
          endif
c
        enddo        
c 
        do i=1,IoutletSpecifiedResistance
c
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
c
          BPressureC(i2,i3)=PressureCResistance(i2)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
c
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=rhof*dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=rhof*dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c

            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          =
            geoDiff=rhof*dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          BPressureC(i2,i3)=geoDiff*PressureC(i1)/
     *            (geoDiff-Bmdot(i2,i3)*Bdrhodp(i2,i3)/rhof)
c
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
c
      if(LRotationalPeriodicity) then

        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
            yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
            zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            BPressureC(i2,i3)=GFactCF*PressureC(i1)+
     *                         (1.-GFactCF)*PressureC(j1)
c
            BPressureC(j2,j3)=BPressureC(i2,i3)
c
          endif
c
        enddo
c
      elseif(LTranslationalPeriodicity) then

        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=xc(j1)+xTranslation(j2)
            yF1=yc(j1)+yTranslation(j2)
            zF1=zc(j1)+zTranslation(j2)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            BPressureC(i2,i3)=GFactCF*PressureC(i1)+
     *                         (1.-GFactCF)*PressureC(j1)
c
            BPressureC(j2,j3)=BPressureC(i2,i3)
c
          endif
c
        enddo
c
      endif
c
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i1=IaxisOwner(i)
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BPressureC(i2,i3)=PressureC(i1)
c
        enddo
c----------------------------------------------------------------------
c
      endif
c
      return
      end