c
C#############################################################################################
c
      SUBROUTINE CalculateCentrifugalForce
c
C#############################################################################################
c
      use Geometry4, only: xc,yc,zc,Volume
      use Geometry1, only: NumberOfElements
      use Variables3, only: FluxTE
      use PhysicalProperties1, only: Density
      use Coriolis1
      use Centrifugal1
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i
      double precision :: wrx,wry,wrz,x0,y0,z0,wx,wy,wz
      double precision :: mass
c********************************************************************************************
c
      x0=AxisOfRotationOriginX
      y0=AxisOfRotationOriginY
      z0=AxisOfRotationOriginZ
c
      wx=AngularVelocityX
      wy=AngularVelocityY
      wz=AngularVelocityZ
c
      do i=1,NumberOfElements
c
        wrx=(zc(i)-z0)*wy-(yc(i)-y0)*wz
        wry=(xc(i)-x0)*wz-(zc(i)-z0)*wx
        wrz=(yc(i)-y0)*wx-(xc(i)-x0)*wy
c
        mass=Density(i)*Volume(i)
        xCentrifugal(i)=(wrz*wy-wry*wz)*mass
        yCentrifugal(i)=(wrx*wz-wrz*wx)*mass
        zCentrifugal(i)=(wry*wx-wrx*wy)*mass
c
      enddo
c
      return
      end
