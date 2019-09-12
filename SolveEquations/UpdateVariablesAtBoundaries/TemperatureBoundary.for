c
c#############################################################################################
c
      SUBROUTINE UpdateBoundaryTemperature
c
c#############################################################################################
c
      use User0, only: ConstantDensity
      use PhysicalProperties1, only: RGas,BdrhodP,
     *                               BSpecificHeat,BDensity
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use BoundaryConditions1, only: BCType,inletTypeEnergy,
     *                               outletTypeEnergy
      use Variables1, only: Temperature,BTemperature,
     *                      BStagnationTemperature,BPressure,
     *                      BuVelocity,BvVelocity,BwVelocity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      double precision :: velocity2
c********************************************************************************************
c
c--- Extrapolate temperature along outlet boundaries
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
          if(BCType(i,j).eq.'outlet') then
c
            if(outletTypeEnergy(i,j).ne.'transmissive')
     *                     BTemperature(i,j)=Temperature(k)
            BdrhodP(i,j)=1./(RGas*BTemperature(i,j))
c
          endif
c
        enddo
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
          if(BCType(i,j).eq.'pressurefarfield') then
c
            BdrhodP(i,j)=1./(RGas*BTemperature(i,j))
c
          endif
c
        enddo
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
          if(BCType(i,j).eq.'periodic') then
c
            BdrhodP(i,j)=1./(RGas*BTemperature(i,j))
c
          endif
c
        enddo
      enddo
c
c--- Calculate temperature along inlet boundaries
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          if(BCType(i,j).eq.'inlet') then
c
            if(inletTypeEnergy(i,j).eq.
     *                 'specifiedstagnationtemperature') then
c
              velocity2=BuVelocity(i,j)**2+
     *                  BvVelocity(i,j)**2+BwVelocity(i,j)**2
              BTemperature(i,j)=BStagnationTemperature(i,j)-
     *                            velocity2/(2.*BSpecificHeat(i,j))
              BdrhodP(i,j)=1./(RGas*BTemperature(i,j))
c
            endif
c
          endif
c
        enddo
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BDensity(i,j)=BPressure(i,j)*BdrhodP(i,j)+ConstantDensity
c
        enddo
      enddo
c
      return
      end
