c
C#############################################################################################
c
      SUBROUTINE UpdateArteryResistance
c
C#############################################################################################
c
      use Geometry1, only: NumberOfBCSets
      use ArteryResistance1, only: ArteryResistance
      use WindKessel1, only: ResistanceToBloodFlow,
     *                       TotalPeripheralResistance,ComplianceC
      use Transient1, only: dtOld
      use BoundaryConditions1, only: BoundaryType,outletTypeC
      use User0, only: dt,TransientScheme,WindKesselType,LWindKessel
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      Double Precision :: Numerator,Denominator,dtCN,factor1
c********************************************************************************************
c
      CalculateArteryResistance:select case (WindKesselType)
c-----------------------------------------------------------------------------------------------      
        case(1) CalculateArteryResistance   !constant R
c-----------------------------------------------------------------------------------------------      
c
          return
c
c-----------------------------------------------------------------------------------------------      
        case(2) CalculateArteryResistance   !2-element WindKessel model
c-----------------------------------------------------------------------------------------------      
c
          if(TransientScheme.eq.'euler') then
c
            do i=1,NumberOfBCSets
c      
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                         ComplianceC(i)/dt
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
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
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                      ComplianceC(i)/dtCN
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
c
              endif
c
            enddo
c
          elseif(TransientScheme.eq.'cranknicolson2') then
c
            factor1=(dtOld/dt)/(dt+dtOld)
c
            do i=1,NumberOfBCSets
c      
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                      ComplianceC(i)*factor1
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
c
              endif
c
            enddo
c
          elseif(TransientScheme.eq.'adamsmoulton') then
c
            factor1=1./dt+1./(dt+dtOld)
c
            do i=1,NumberOfBCSets
c      
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                      ComplianceC(i)*factor1
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
c
              endif
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(3) CalculateArteryResistance   !3-element WindKessel model
c-----------------------------------------------------------------------------------------------      
c
          if(TransientScheme.eq.'euler') then
c
            do i=1,NumberOfBCSets
c      
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.+ResistanceToBloodFlow(i)/
     *             TotalPeripheralResistance(i)+
     *               ComplianceC(i)*ResistanceToBloodFlow(i)/dt
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                      ComplianceC(i)/dt
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
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
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.+ResistanceToBloodFlow(i)/
     *             TotalPeripheralResistance(i)+
     *             ComplianceC(i)*ResistanceToBloodFlow(i)/dtCN
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                      ComplianceC(i)/dtCN
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
c
              endif
c
            enddo
c
          elseif(TransientScheme.eq.'cranknicolson2') then
c
            factor1=(dtOld/dt)/(dt+dtOld)
c
            do i=1,NumberOfBCSets
c      
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.+ResistanceToBloodFlow(i)/
     *             TotalPeripheralResistance(i)+
     *             ComplianceC(i)*ResistanceToBloodFlow(i)*factor1
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                      ComplianceC(i)*factor1
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
c
              endif
c
            enddo
c
          elseif(TransientScheme.eq.'adamsmoulton') then
c
            factor1=1./dt+1./(dt+dtOld)
c
            do i=1,NumberOfBCSets
c      
              if(LWindKessel(i)) then
c
                if(BoundaryType(i).eq.'outlet'.and.
     *           outletTypeC(i).eq.'specifiedresistance') then
c
                  Numerator=1.+ResistanceToBloodFlow(i)/
     *             TotalPeripheralResistance(i)+
     *             ComplianceC(i)*ResistanceToBloodFlow(i)*factor1
                  Denominator=1./TotalPeripheralResistance(i)+
     *                                      ComplianceC(i)*factor1
                  ArteryResistance(i)=Numerator/Denominator
c
                endif
c
              endif
c
            enddo
c
          endif
c
c-----------------------------------------------------------------------------------------------      
        case(4) CalculateArteryResistance   !4-element WindKessel model
c-----------------------------------------------------------------------------------------------      
c
c-----------------------------------------------------------------------------------------------      
      end select CalculateArteryResistance 
c-----------------------------------------------------------------------------------------------      
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE CalculatePressureAtOutlets
c
C#############################################################################################
c
      use User0, only: urfPressure,LWindKessel
      use BoundaryConditions1, only: BoundaryType,outletTypeC
      use ArteryResistance1, only: LArteryExplicit,ArteryResistance,
     *                             urfPressureResistance
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces
      use Geometry4, only: BFaceArea
      use Variables1, only: BPressure
      use FlowInOut1, only: MassFlowRate
      use PhysicalProperties1, only: BDensity
      use WindKessel1, only: OutletPressure
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
      double precision :: Area
c********************************************************************************************
c
      do i=1,NumberOfBCSets
c      
        if(BoundaryType(i).eq.'outlet'.and.outletTypeC(i).eq.
     *         'specifiedstaticpressure'.and.LArteryExplicit(i)
     *                                .and..not.LWindKessel(i)) then
c
          do j=1,NBFaces(i)
        
            BPressure(i,j)=urfPressureResistance(i)*ArteryResistance(i)*
     *                      MassFlowRate(i)/BDensity(i,j)+
     *                   (1.-urfPressureResistance(i))*BPressure(i,j)
          enddo
c
        elseif(BoundaryType(i).eq.'outlet'.and.outletTypeC(i).eq.
     *         'specifiedresistance'.and.LWindKessel(i)) then
c
          Area=0.
          OutletPressure(i)=0.
c
          do j=1,NBFaces(i)
c
            OutletPressure(i)=OutletPressure(i)+
     *                 BPressure(i,j)*BFaceArea(i,j) 
            Area=Area+BFaceArea(i,j)
c
          enddo
c
          OutletPressure(i)=OutletPressure(i)/Area
c
        elseif(BoundaryType(i).eq.'inlet') then
c
          Area=0.
          OutletPressure(i)=0.
c
          do j=1,NBFaces(i)
c
            OutletPressure(i)=OutletPressure(i)+
     *                 BPressure(i,j)*BFaceArea(i,j) 
            Area=Area+BFaceArea(i,j)
c
          enddo
c
          OutletPressure(i)=OutletPressure(i)/Area
c
        endif
c
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE UpdateSourceResistance(iter,FiT)
c
C#############################################################################################
c
      use PhysicalProperties1, only: BDensity
      use Variables2, only: bc,bcOriginal
      use BoundaryConditions2, only: IoutletSpecifiedResistance,
     *                               IoutletSpecifiedResistanceOwner,
     *                        IoutletSpecifiedResistanceNumberOfBCSets,
     *                          IoutletSpecifiedResistanceNBFaces
      use ArteryResistance1, only: ArteryResistance,PressureCResistance,
     *                             geoDiffB,geoDiffSum
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,iter
      double precision, dimension(:)  :: FiT
c********************************************************************************************
c
      if(iter.eq.1) then
        bcOriginal=bc
      endif
c
      PressureCResistance=0.
      do i=1,IoutletSpecifiedResistance
c
        i1=IoutletSpecifiedResistanceOwner(i)
        i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
        i3=IoutletSpecifiedResistanceNBFaces(i)
c
        PressureCResistance(i2)=PressureCResistance(i2)+
     *                             geoDiffB(i2,i3)*FiT(i1)
c
      enddo
c
      do i=1,IoutletSpecifiedResistance
c
        i1=IoutletSpecifiedResistanceOwner(i)
        i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
        i3=IoutletSpecifiedResistanceNBFaces(i)
c
        bc(i1)=bcOriginal(i1)+BDensity(i2,i3)*geoDiffB(i2,i3)*
     *                  ArteryResistance(i2)*
     *             (PressureCResistance(i2)-geoDiffB(i2,i3)*FiT(i1))/
     *                     (1.+ArteryResistance(i2)*geoDiffSum(i2))
c
      enddo        
c 
      return
      end
c
C#############################################################################################
c
      SUBROUTINE UpdateOutletPressureForResistanceConditions
c
C#############################################################################################
c
      use User0, only: urfPressure,LWindKessel
      use BoundaryConditions1, only: BoundaryType,outletTypeC
      use ArteryResistance1, only: LArteryExplicit,ArteryResistance,
     *                             urfPressureResistance
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: BPressure
      use FlowInOut1, only: MassFlowRate
      use PhysicalProperties1, only: BDensity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfBCSets
c      
        if(BoundaryType(i).eq.'outlet'.and.outletTypeC(i).eq.
     *         'specifiedstaticpressure'.and.LArteryExplicit(i)
     *                                 .and..not.LWindKessel(i)) then
c
          do j=1,NBFaces(i)
        
            BPressure(i,j)=urfPressureResistance(i)*ArteryResistance(i)*
     *                      MassFlowRate(i)/BDensity(i,j)+
     *                   (1.-urfPressureResistance(i))*BPressure(i,j)

          enddo
c
        endif
c
      enddo
c
      return
      end