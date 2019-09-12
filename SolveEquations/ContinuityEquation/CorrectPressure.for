c
c#############################################################################################
c
      SUBROUTINE CorrectPressure
c
c#############################################################################################
c
      use User0, only: ReferencePressureLocation,Lcompressible,
     *                 LfixPressure,FixedPressureValue,FixAtLocation,
     *                 urfPressure
      use ReferencePressure1, only: LSetReferencePressure
      use Variables1, only: Pressure,PressureC
      use Geometry1, only: NumberOfElements
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: RPressureC
c********************************************************************************************
c
      i=ReferencePressureLocation
      if(i.ne.0) then
        RPressureC=PressureC(i)
      else
        RPressureC=0.
      endif
c
      if(LSetReferencePressure.or.Lcompressible) RPressureC=0. 
c
      if(LfixPressure) Then
c
        RPressureC=0.
        i=FixAtLocation
        Pressure(i)=FixedPressureValue
c
      endif
c      
      do i=1,NumberOfElements           
c
        Pressure(i)=Pressure(i)+urfPressure*(PressureC(i)-RPressureC)   
c      
      enddo      
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE UpdateBoundaryPressure
c
c#############################################################################################
c
      use User0, only: ReferencePressureLocation,Lcompressible,
     *                 LfixPressure,MethodDecomposeS,urfPressure
      use ReferencePressure1, only: LSetReferencePressure
      use Variables1, only: PressureC,BPressureC,BPressure,Pressure,
     *                      PressGradx,PressGrady,PressGradz,
     *                      BPressGradx,BPressGrady,BPressGradz
      use BoundaryConditions1, only: BoundaryType,OutletTypeC
      use BoundaryConditions2
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces,NBFaceOwner
      use Geometry4, only: xc,yc,zc,
     *                     BDistanceCF,BDistanceCFx,BDistanceCFy,
     *                     BDistanceCFz,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz,BFaceArea
      use AveragePressure1, only: AverageOutletPressure
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,i1,i2,i3,j1,j2,j3
      double precision :: RPressureC,rhof,xF1,yF1,zF1,
     *                    distance1,distance2,GFactCF,
     *                    AveragePressure,TotalArea
c********************************************************************************************
c
      i=ReferencePressureLocation
      if(i.ne.0) then
        RPressureC=PressureC(i)
      else
        RPressureC=0.
      endif
      if(LSetReferencePressure.or.Lcompressible) RPressureC=0. 
      if(LfixPressure) RPressureC=0. 
c
c----------------------------------------------------------------------
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
c          BPressure(i2,i3)=BPressure(i2,i3)+
c     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
          BPressure(i2,i3)=Pressure(i1)+
     *          BPressGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BPressGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                  BPressGradz(i2,i3)*BDistanceCFz(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
c
c          BPressure(i2,i3)=BPressure(i2,i3)+
c     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
          BPressure(i2,i3)=Pressure(i1)+
     *          BPressGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BPressGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                   BPressGradz(i2,i3)*BDistanceCFz(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationPressure
c
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedStaticPressure
c
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i1=IoutletSpecifiedAverageStaticPressureOwner(i)
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
c
          BPressure(i2,i3)=Pressure(i1)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'outlet'.and.OutletTypeC(i).eq.
     *                   'specifiedaveragestaticpressure') then
c
            AveragePressure=0.
            TotalArea=0.
c
            do j=1,NBFaces(i)
c
              AveragePressure=AveragePressure+
     *                  BFaceArea(i,j)*BPressure(i,j)
              TotalArea=TotalArea+BFaceArea(i,j)
c
            enddo
c
            AveragePressure=AveragePressure/(TotalArea+tiny)
            AveragePressure=AverageOutletPressure(i)-AveragePressure
c
            do j=1,NBFaces(i)
c
              BPressure(i,j)=BPressure(i,j)+AveragePressure
c
            enddo
c
          endif
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedResistance
c
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletFullyDeveloped
c
          i2=IoutletFullyDevelopedNumberOfBCSets(i)
          i3=IoutletFullyDevelopedNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BPressure(i2,i3)=Pressure(i1)+
     *          BPressGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BPressGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                  BPressGradz(i2,i3)*BDistanceCFz(i2,i3)
c          BPressure(i2,i3)=BPressure(i2,i3)+
c     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
c
      if(LRotationalPeriodicity) then
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            BPressure(i2,i3)=GFactCF*Pressure(i1)+
     *                         (1.-GFactCF)*Pressure(j1)
c
            BPressure(j2,j3)=BPressure(i2,i3)
c
          endif
c
        enddo
c
      elseif(LTranslationalPeriodicity) then
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
            BPressure(i2,i3)=GFactCF*Pressure(i1)+
     *                         (1.-GFactCF)*Pressure(j1)
c
            BPressure(j2,j3)=BPressure(i2,i3)
c
          endif
c
        enddo
c
      endif
c----------------------------------------------------------------------
        do i=1,Iaxis
c
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
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
c          BPressure(i2,i3)=BPressure(i2,i3)+
c     *                       urfPressure*BPressureC(i2,i3)
          BPressure(i2,i3)=Pressure(i1)+
     *          BPressGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BPressGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                  BPressGradz(i2,i3)*BDistanceCFz(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
c
c          BPressure(i2,i3)=BPressure(i2,i3)+
c     *                       urfPressure*BPressureC(i2,i3)
          BPressure(i2,i3)=Pressure(i1)+
     *          BPressGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BPressGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                  BPressGradz(i2,i3)*BDistanceCFz(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iinletsupersonic
c
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c

        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationPressure
c
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IpressureFarField
c
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletTransmissive
c
          i2=IoutletTransmissiveNumberOfBCSets(i)
          i3=IoutletTransmissiveNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Ioutletsupersonic
c
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedStaticPressure
c
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i1=IoutletSpecifiedAverageStaticPressureOwner(i)
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
c
          BPressure(i2,i3)=Pressure(i1)
c
        enddo
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'outlet'.and.OutletTypeC(i).eq.
     *                   'specifiedaveragestaticpressure') then
c
            AveragePressure=0.
            TotalArea=0.
c
            do j=1,NBFaces(i)
c
              AveragePressure=AveragePressure+
     *                  BFaceArea(i,j)*BPressure(i,j)
              TotalArea=TotalArea+BFaceArea(i,j)
c
            enddo
c
            AveragePressure=AveragePressure/TotalArea
            AveragePressure=AverageOutletPressure(i)-AveragePressure
c
            do j=1,NBFaces(i)
c
              BPressure(i,j)=BPressure(i,j)+AveragePressure
c
            enddo
c
          endif
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedResistance
c
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                    urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Isymmetry
c
          i1=IsymmetryOwner(i)
          i2=IsymmetryNumberOfBCSets(i)
          i3=IsymmetryNBFaces(i)
c
          BPressure(i2,i3)=Pressure(i1)+
     *          BPressGradx(i2,i3)*BDistanceCFx(i2,i3)+
     *               BPressGrady(i2,i3)*BDistanceCFy(i2,i3)+
     *                  BPressGradz(i2,i3)*BDistanceCFz(i2,i3)
c
c          BPressure(i2,i3)=BPressure(i2,i3)+
c     *                       urfPressure*BPressureC(i2,i3)
c
        enddo
c----------------------------------------------------------------------
c
        if(LRotationalPeriodicity) then
c
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
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
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
              GFactCF=distance2/(distance1+distance2)
c
              BPressure(i2,i3)=GFactCF*Pressure(i1)+
     *                         (1.-GFactCF)*Pressure(j1)
c
              BPressure(j2,j3)=BPressure(i2,i3)
c
            endif
c
          enddo
c
        elseif(LTranslationalPeriodicity) then
c
          do i=1,Iperiodic
c
            i1=IperiodicOwner(i)
            i2=IperiodicNumberOfBCSets(i)
            i3=IperiodicNBFaces(i)
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
              BPressure(i2,i3)=GFactCF*Pressure(i1)+
     *                         (1.-GFactCF)*Pressure(j1)
c
              BPressure(j2,j3)=BPressure(i2,i3)
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
          i2=IaxisNumberOfBCSets(i)
          i3=IaxisNBFaces(i)
c
          BPressure(i2,i3)=BPressure(i2,i3)+
     *                urfPressure*(BPressureC(i2,i3)-RPressureC)
c
        enddo
c----------------------------------------------------------------------
c
      endif
c
      return
      end