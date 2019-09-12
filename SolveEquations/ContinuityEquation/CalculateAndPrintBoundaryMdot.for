c
C#############################################################################################
c
      SUBROUTINE CalculateBoundaryMassFlowRates
c
c#############################################################################################
c
      use User0, only: time
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: Bmdot
      use FlowInOut1, only: MassFlowRate,LPrintMassFlowRate
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      MassFlowRate=0.
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c              
          MassFlowRate(i)=MassFlowRate(i)+Bmdot(i,j)
c
        enddo
      enddo
c
      return
c
c-------------------------------------------------------------------------------------------
      entry PrintBoundaryMassFlowRates
c-------------------------------------------------------------------------------------------

      do i=1,NumberOfBCSets
        If(LPrintMassFlowRate(i)) then
c              
          Print*,'Mass flow rate through boundary ',i,' = ',
     *      MassFlowRate(i),'Kg/s'
c
        endif
c
      enddo
c
      return
c
c-------------------------------------------------------------------------------------------
      entry PrintBoundaryMassFlowRatesInTime
c-------------------------------------------------------------------------------------------
c
        write(16,*) time,(MassFlowRate(i),i=1,NumberOfBCSets)
c
      return
      end
c