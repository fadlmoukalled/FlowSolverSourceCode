c
c#############################################################################################
c
      FUNCTION BTangentialVelocity(i1,i2)
c
C#############################################################################################
c
      use Variables1, only:BuVelocity,BvVelocity,BwVelocity
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
      use WallDistance1, only: BiTau,BjTau
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i1,i2
      integer :: i3,i4
      double precision :: BTangentialVelocity
      double precision :: nx,ny,nz,uWall,vWall,wWall,dotproduct
c********************************************************************************************
c
      i3=BiTau(i1,i2)
      i4=BjTau(i1,i2)
c
      nx=BFaceAreanx(i3,i4)
      ny=BFaceAreany(i3,i4)
      nz=BFaceAreanz(i3,i4)
c
      uWall=BuVelocity(i3,i4)-BuVelocity(i1,i2)
      vWall=BvVelocity(i3,i4)-BvVelocity(i1,i2)
      wWall=BwVelocity(i3,i4)-BwVelocity(i1,i2)
c
      dotproduct=uWall*nx+vWall*ny+wWall*nz
c
      uWall=uWall-dotproduct*nx
      vWall=vWall-dotproduct*ny
      wWall=wWall-dotproduct*nz
c
      BTangentialVelocity=dsqrt(uWall*uWall+
     *                      vWall*vWall+wWall*wWall)
c
      end FUNCTION BTangentialVelocity