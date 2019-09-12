c
c#############################################################################################
c
      SUBROUTINE UpdateSlipWallBoundary
c
c#############################################################################################
c
      use BoundaryConditions2, only: IWallSlip,IWallSlipOwner,
     *                               IWallSlipNumberOfBCSets,
     *                               IWallSlipNBFaces
      use VAriables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: nx,ny,nz,Vdotn
c********************************************************************************************
c
      do i=1,IWallSlip
c
        i1=IWallSlipOwner(i)
        i2=IWallSlipNumberOfBCSets(i)
        i3=IWallSlipNBFaces(i)
c
        nx=BFaceAreanx(i2,i3)
        ny=BFaceAreany(i2,i3)
        nz=BFaceAreanz(i2,i3)
c
        Vdotn=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
c
        BuVelocity(i2,i3)=uVelocity(i1)-Vdotn*nx
        BvVelocity(i2,i3)=vVelocity(i1)-Vdotn*ny
        BwVelocity(i2,i3)=wVelocity(i1)-Vdotn*nz
c
      enddo
c
      return
      end
