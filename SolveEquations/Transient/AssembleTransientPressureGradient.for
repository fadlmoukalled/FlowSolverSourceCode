c
C#############################################################################################
c
      SUBROUTINE AssembleTransientPressureGradientTerm
c
C#############################################################################################
c
      use Variables3, only: FluxTE
      use variables1, only: Pressure,PressureOld,PressureOldOld
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use Transient1, only:dpdt
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE ComputeTransientGradient(phi,phiOld,phiOldOld,dphidt)
c********************************************************************************************
          double precision, dimension(:) :: phi,phiOld,phiOldOld,dphidt
c********************************************************************************************
        end SUBROUTINE ComputeTransientGradient
c********************************************************************************************
      end interface
c********************************************************************************************
c
      dpdt=0.d0
c
      call ComputeTransientGradient
     *            (Pressure,PressureOld,PressureOldOld,dpdt)
c
      do i=1,NumberOfElements
c
        FluxTE(i)=FluxTE(i)-dpdt(i)*Volume(i)
c
      enddo
c
      return
      end