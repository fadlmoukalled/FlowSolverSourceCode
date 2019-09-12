c
c#############################################################################################
c
c     IMPLEMENTATION OF THE WINDKESSEL MODELS
c
c#############################################################################################
c
      SUBROUTINE WindKesselModel
c
c#############################################################################################
      use User0, only: WindKesselType
c*********************************************************************************************
      implicit none
c*********************************************************************************************
c
      if(WindKesselType.eq.2) then
        call WindKessel2
      elseif(WindKesselType.eq.3) then
        call WindKessel3
      elseif(WindKesselType.eq.4) then
        call WindKessel4
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE WindKessel2
c
c#############################################################################################
      use User0, only: ConstantDensity,LWindKessel,dt,TransientScheme,
     *                 urfOutletPressure
      use WindKessel1, only: OutletPressure,OutletPressureOld,
     *                       OutletPressureOldOld,ComplianceC,
     *                       TotalPeripheralResistance
      use Transient1, only: dtOld,dtOldOld
      use FlowInOut1, only: MassFlowRate
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces
      use BoundaryConditions1, only: BoundaryType,outletTypeM
      use Variables1, only: BPressure
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j
      double precision :: VolumeFlowRate,dtCN,factor1,factor2,factor3
c*********************************************************************************************
c
      if(TransientScheme.eq.'euler') then
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            OutletPressure(i)=
     *         (OutletPressureOld(i)+VolumeFlowRate*dt/ComplianceC(i))/
     *             (1.+dt/(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      elseif(TransientScheme.eq.'cranknicolson1') then
c
        dtCN=dt/2.
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            OutletPressure(i)=
     *       (OutletPressureOld(i)+VolumeFlowRate*dtCN/ComplianceC(i))/
     *           (1.+dtCN/(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      elseif(TransientScheme.eq.'cranknicolson2') then
c
        factor1=(dtOld/dt)/(dt+dtOld)
        factor2=1./(dt+dtOld)-(dtOldOld/dt)/(dtOld+dtOldOld)
        factor3=-(dtOld/dt)/(dtOld+dtOldOld)
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            OutletPressure(i)=(-factor2*OutletPressureOld(i)-
     *           factor3*OutletPressureOldOld(i)+VolumeFlowRate/
     *             ComplianceC(i))/(factor1+
     *               1./(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      elseif(TransientScheme.eq.'adamsmoulton') then
c
        factor1=1./dt+1./(dt+dtOld)
        factor3=(dtOld/dt)/(dtOld+dtOldOld)
        factor2=factor1+factor3
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            OutletPressure(i)=(factor2*OutletPressureOld(i)-
     *           factor3*OutletPressureOldOld(i)+VolumeFlowRate/
     *             ComplianceC(i))/(factor1+
     *               1./(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE WindKessel3
c
c#############################################################################################
      use User0, only: ConstantDensity,LWindKessel,dt,TransientScheme,
     *                 urfOutletPressure
      use Transient1, only: dtOld,dtOldOld
      use WindKessel1, only: OutletPressure,OutletPressureOld,
     *                       OutletPressureOldOld,ComplianceC,
     *                       ResistanceToBloodFlow,MassFlowrateOld,
     *                       MassFlowrateOldOld,
     *                       TotalPeripheralResistance
      use FlowInOut1, only: MassFlowRate
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces
      use BoundaryConditions1, only: BoundaryType,outletTypeM
      use Variables1, only: BPressure
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j
      double precision :: VolumeFlowRate,VolumeFlowRateOld,
     *                    VolumeFlowRateOldOld,VolumeDiff,
     *                    RightHandSide,dtCN,factor1,factor2,factor3
c*********************************************************************************************
c
      if(TransientScheme.eq.'euler') then
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            VolumeFlowRateOld=MassFlowrateOld(i)/ConstantDensity
            VolumeDiff=VolumeFlowRate-VolumeFlowRateOld
            RightHandSide=VolumeFlowRate*
     *      (1.+ResistanceToBloodFlow(i)/TotalPeripheralResistance(i))+
     *            ComplianceC(i)*ResistanceToBloodFlow(i)*VolumeDiff/dt
            OutletPressure(i)=
     *      (OutletPressureOld(i)+RightHandSide*dt/ComplianceC(i))/
     *            (1.+dt/(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      elseif(TransientScheme.eq.'cranknicolson1') then
c
        dtCN=dt/2.
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            VolumeFlowRateOld=MassFlowrateOld(i)/ConstantDensity
            VolumeDiff=VolumeFlowRate-VolumeFlowRateOld
            RightHandSide=VolumeFlowRate*
     *      (1.+ResistanceToBloodFlow(i)/TotalPeripheralResistance(i))+
     *          ComplianceC(i)*ResistanceToBloodFlow(i)*VolumeDiff/dtCN
            OutletPressure(i)=
     *      (OutletPressureOld(i)+RightHandSide*dtCN/ComplianceC(i))/
     *          (1.+dtCN/(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      elseif(TransientScheme.eq.'cranknicolson2') then
c
        factor1=(dtOld/dt)/(dt+dtOld)
        factor2=1./(dt+dtOld)-(dtOldOld/dt)/(dtOld+dtOldOld)
        factor3=-(dtOld/dt)/(dtOld+dtOldOld)
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            VolumeFlowRateOld=MassFlowrateOld(i)/ConstantDensity
            VolumeFlowRateOldOld=MassFlowrateOldOld(i)/ConstantDensity
            VolumeDiff=factor1*VolumeFlowRate+factor2*
     *            VolumeFlowRateOld+factor3*VolumeFlowRateOldOld
            RightHandSide=VolumeFlowRate*
     *      (1.+ResistanceToBloodFlow(i)/TotalPeripheralResistance(i))+
     *               ComplianceC(i)*ResistanceToBloodFlow(i)*VolumeDiff
            OutletPressure(i)=(-factor2*OutletPressureOld(i)-
     *                 factor3*OutletPressureOldOld(i)+RightHandSide/
     *              ComplianceC(i))/(factor1+
     *               1./(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      elseif(TransientScheme.eq.'adamsmoulton') then
c
        factor1=1./dt+1./(dt+dtOld)
        factor3=(dtOld/dt)/(dtOld+dtOldOld)
        factor2=factor1+factor3
c
        do i=1,NumberOfBCSets
c      
          if(BoundaryType(i).eq.'outlet'.and.outletTypeM(i).eq.
     *          'specifiedstaticpressure'.and.LWindKessel(i)) then
c          
            VolumeFlowRate=MassFlowrate(i)/ConstantDensity
            VolumeFlowRateOld=MassFlowrateOld(i)/ConstantDensity
            VolumeFlowRateOldOld=MassFlowrateOldOld(i)/ConstantDensity
            VolumeDiff=factor1*VolumeFlowRate-factor2*
     *            VolumeFlowRateOld+factor3*VolumeFlowRateOldOld
            RightHandSide=VolumeFlowRate*
     *      (1.+ResistanceToBloodFlow(i)/TotalPeripheralResistance(i))+
     *               ComplianceC(i)*ResistanceToBloodFlow(i)*VolumeDiff
            OutletPressure(i)=(factor2*OutletPressureOld(i)-
     *                 factor3*OutletPressureOldOld(i)+RightHandSide/
     *              ComplianceC(i))/(factor1+
     *               1./(TotalPeripheralResistance(i)*ComplianceC(i)))
c
            do j=1,NBFaces(i)
              BPressure(i,j)=urfOutletPressure*OutletPressure(i)+
     *                            (1.-urfOutletPressure)*BPressure(i,j)
            enddo
c
          endif      
c      
        enddo
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE WindKessel4
c
c#############################################################################################
      use WindKessel1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
c
c
      return
      end
c
c#############################################################################################
c
c     END OF IMPLEMENTATION OF THE WINDKESSEL MODELS
c
c#############################################################################################