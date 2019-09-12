c
c#############################################################################################
c
      SUBROUTINE CalculateNormalVelocity
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
      use WallDistance1, only: WallDistance,BWallDistance,
     *                         iTau,jTau,BiTau,BjTau
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      NormalVelocity,BNormalVelocity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,j,k
      double precision :: nx,ny,nz,Nvelx,Nvely,Nvelz
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        j=itau(i)
        k=jtau(i) 
c
        nx=BFaceAreanx(j,k)
        ny=BFaceAreany(j,k)
        nz=BFaceAreanz(j,k)
c
        Nvelx=(uVelocity(i)-BuVelocity(j,k))*nx
        Nvely=(vVelocity(i)-BvVelocity(j,k))*ny
        Nvelz=(wVelocity(i)-BwVelocity(j,k))*nz
        NormalVelocity(i)=Nvelx+Nvely+Nvelz
c
      enddo
c
      do i1=1,NumberOfBCSets
        do i2=1,NBFaces(i1)
c
          j=Bitau(i1,i2)
          k=Bjtau(i1,i2) 
c
          nx=BFaceAreanx(j,k)
          ny=BFaceAreany(j,k)
          nz=BFaceAreanz(j,k)
c
          Nvelx=(BuVelocity(i1,i2)-BuVelocity(j,k))*nx
          Nvely=(BvVelocity(i1,i2)-BvVelocity(j,k))*ny
          Nvelz=(BwVelocity(i1,i2)-BwVelocity(j,k))*nz
          BNormalVelocity(i1,i2)=Nvelx+Nvely+Nvelz
c
        enddo
      enddo
c
      return
      end
