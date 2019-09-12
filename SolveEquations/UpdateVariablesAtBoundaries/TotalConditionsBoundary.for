
c
c#############################################################################################
c
      SUBROUTINE UpdateBoundaryForTotalConditions
c
c#############################################################################################
c
      use User0, only: Lcompressible
      use BoundaryConditions2, only: IinletSpecifiedStagnationPressure,
     *                          IinletSpecifiedStagnationPressureOwner,
     *                 IinletSpecifiedStagnationPressureNumberOfBCSets,
     *                        IinletSpecifiedStagnationPressureNBFaces
c
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                     BPressure,BTemperature,BStagnationPressure,
     *                     BStagnationTemperature,
     *                     uVelGradx,uvelGrady,uvelGradz,
     *                     vVelGradx,vVelGrady,vVelGradz,
     *                     wVelGradx,wVelGrady,wVelGradz
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BDistanceCFx,BDistanceCFy,BDistanceCFz
      use PhysicalProperties1, only: BSpecificHeat,BDensity,RGas
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j
      double precision :: velocitySquared,ratio,ratio1,ratio2,ratio3
c********************************************************************************************
c
      if(.not.Lcompressible) then
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
c          BuVelocity(i2,i3)=uVelocity(i1)
c          BvVelocity(i2,i3)=vVelocity(i1)
c          BwVelocity(i2,i3)=wVelocity(i1)
          velocitySquared=BuVelocity(i2,i3)**2+
     *                BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2
          BPressure(i2,i3)=dmax1(BStagnationPressure(i2,i3)
     *                       -0.5*BDensity(i2,i3)*velocitySquared,0.)
c          BPressure(i2,i3)=BStagnationPressure(i2,i3)
c     *                       -0.5*BDensity(i2,i3)*velocitySquared
c
        enddo
c
      elseif(Lcompressible) then
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
c          BuVelocity(i2,i3)=uVelocity(i1) !+
c     *          uVelGradx(i1)*BDistanceCFx(i2,i3)+
c     *               uVelGrady(i1)*BDistanceCFy(i2,i3)+
c     *                   uVelGradz(i1)*BDistanceCFz(i2,i3)
c          BvVelocity(i2,i3)=vVelocity(i1) !+
c     *          vVelGradx(i1)*BDistanceCFx(i2,i3)+
c     *               vVelGrady(i1)*BDistanceCFy(i2,i3)
c     *                   vVelGradz(i1)*BDistanceCFz(i2,i3)
c          BwVelocity(i2,i3)=wVelocity(i1) !+
c     *          wVelGradx(i1)*BDistanceCFx(i2,i3)+
c     *               wVelGrady(i1)*BDistanceCFy(i2,i3)
c     *                   wVelGradz(i1)*BDistanceCFz(i2,i3)
          velocitySquared=BuVelocity(i2,i3)**2+
     *                 BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2
c           
          ratio=BSpecificHeat(i2,i3)/(BSpecificHeat(i2,i3)-RGas)
          ratio1=ratio-1
          ratio2=ratio/ratio1
          ratio3=velocitySquared/(ratio*RGas*BTemperature(i2,i3))
c
          BPressure(i2,i3)=
     *      BStagnationPressure(i2,i3)/(1+0.5*ratio1*ratio3)**ratio2
          BDensity(i2,i3)=BPressure(i2,i3)/(RGas*BTemperature(i2,i3))
c
        enddo
c
      endif
c
      return
      end
c