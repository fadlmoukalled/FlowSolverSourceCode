c
c#############################################################################################
c
      SUBROUTINE UpdateSymmetry
c
c#############################################################################################
c
      use BoundaryConditions2, only: Isymmetry,IsymmetryOwner,
     *                               IsymmetryNumberOfBCSets,
     *                               IsymmetryNBFaces,
     *                               Iaxis,IaxisOwner,
     *                               IaxisNumberOfBCSets,
     *                               IaxisNBFaces
      use VAriables1, only: uVelocity,vVelocity,wVelocity,BuVelocity,
     *                      BvVelocity,BwVelocity
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: areanx,areany,areanz,Amdot
c********************************************************************************************
c
      do i=1,Isymmetry
c
        i1=IsymmetryOwner(i)
        i2=IsymmetryNumberOfBCSets(i)
        i3=IsymmetryNBFaces(i)
c
        areanx=BFaceAreanx(i2,i3)
        areany=BFaceAreany(i2,i3)
        areanz=BFaceAreanz(i2,i3)
c
        Amdot=uVelocity(i1)*areanx+vVelocity(i1)*areany+
     *                                wVelocity(i1)*areanz
c
        BuVelocity(i2,i3)=uVelocity(i1)-Amdot*areanx
        BvVelocity(i2,i3)=vVelocity(i1)-Amdot*areany
        BwVelocity(i2,i3)=wVelocity(i1)-Amdot*areanz
c
      enddo
c
      do i=1,Iaxis
c
        i1=IaxisOwner(i)
        i2=IaxisNumberOfBCSets(i)
        i3=IaxisNBFaces(i)
c
        BuVelocity(i2,i3)=uVelocity(i1)
        BvVelocity(i2,i3)=vVelocity(i1)
        BwVelocity(i2,i3)=wVelocity(i1)
c
      enddo
c
      return
      end
