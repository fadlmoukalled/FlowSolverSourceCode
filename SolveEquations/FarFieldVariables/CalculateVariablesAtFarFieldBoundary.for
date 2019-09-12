c
c#############################################################################################
c
      SUBROUTINE calculateFarFieldVariables
c
C#############################################################################################
c
      use BoundaryConditions2, only:uVelocityFarField,vVelocityFarField,
     *                             wVelocityFarField,
     *                             pressureFarField,TemperatureFarField,
     *                             MachFarField,xFlowDirectionFarField,
     *                             yFlowDirectionFarField,
     *                             zFlowDirectionFarField,
     *                             SpecificHeatFarField
      use PhysicalProperties1, only: RGas,GammaGas
      use Geometry1, only: NumberOfBCSets
c
c--------------------------------------------------------------
c
      implicit none
      double precision :: ratio,cInfinity,vInfinity
      integer :: i
c
c--------------------------------------------------------------
c
c
      do i=1,NumberOfBCSets
c
        cInfinity=dsqrt(GammaGas*RGas*TemperatureFarField(i))
        vInfinity=MachFarField(i)*cInfinity
        uVelocityFarField(i)=xFlowDirectionFarField(i)*vInfinity
        vVelocityFarField(i)=yFlowDirectionFarField(i)*vInfinity
        wVelocityFarField(i)=zFlowDirectionFarField(i)*vInfinity
c
      enddo
c
      return
      end
