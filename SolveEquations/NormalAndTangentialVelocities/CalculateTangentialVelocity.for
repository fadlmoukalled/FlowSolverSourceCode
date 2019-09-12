c
c#############################################################################################
c
      FUNCTION TangentialVelocity(i1)
c
C#############################################################################################
c
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
      use WallDistance1, only: iTau,jTau
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i1
      integer :: i2,i3
      double precision :: TangentialVelocity
      double precision :: nx,ny,nz,uWall,vWall,wWall,dotproduct
c********************************************************************************************
c
      i2=iTau(i1)
      i3=jTau(i1)
c
      nx=BFaceAreanx(i2,i3)
      ny=BFaceAreany(i2,i3)
      nz=BFaceAreanz(i2,i3)
c
      uWall=BuVelocity(i2,i3)-uVelocity(i1)
      vWall=BvVelocity(i2,i3)-vVelocity(i1)
      wWall=BwVelocity(i2,i3)-wVelocity(i1)
c
      dotproduct=uWall*nx+vWall*ny+wWall*nz
c
      uWall=uWall-dotproduct*nx
      vWall=vWall-dotproduct*ny
      wWall=wWall-dotproduct*nz
c
      TangentialVelocity=dsqrt(uWall*uWall+
     *                      vWall*vWall+wWall*wWall)
c
      end FUNCTION TangentialVelocity
