c
c#############################################################################################
c
      SUBROUTINE CalculateCosineThetaFace
c
c#############################################################################################
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor,NBFaces,
     *                     NBFaceOwner
      use Geometry4, only: xc,yc,zc,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz
      use VolumeOfFluid1, only: cosThetaF,BcosThetaF
      use TransferrField1, only: rFieldGradfxT,rFieldGradfyT,
     *                           rFieldGradfzT,BrFieldGradxT,
     *                           BrFieldGradyT,BrFieldGradzT
      use Constants1, only: tiny 
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,k
      double precision :: rpfx,rpfy,rpfz,rpf,rGradf
c********************************************************************************************
c
      do k=1,NIFaces
c        
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        rpfx=xc(j)-xc(i)
        rpfy=yc(j)-yc(i)
        rpfz=zc(j)-zc(i)
        rpf=dsqrt(rpfx*rpfx+rpfy*rpfy+rpfz*rpfz)
c
        rGradf=dsqrt(rFieldGradfxT(k)*rFieldGradfxT(k)+
     *                rFieldGradfyT(k)*rFieldGradfyT(k)+
     *                 rFieldGradfzT(k)*rFieldGradfzT(k))
        cosThetaF(k)=dabs(rFieldGradfxT(k)*rpfx+
     *                    rFieldGradfyT(k)*rpfy+
     *                    rFieldGradfzT(k)*rpfz)/(rpf*rGradf+tiny)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j) 
c
          rpfx=BFaceCentroidx(i,j)-xc(k)
          rpfy=BFaceCentroidy(i,j)-yc(k)
          rpfz=BFaceCentroidz(i,j)-zc(k)
          rpf=dsqrt(rpfx*rpfx+rpfy*rpfy+rpfz*rpfz)
c
          rGradf=dsqrt(BrFieldGradxT(i,j)*BrFieldGradxT(i,j)+
     *                 BrFieldGradyT(i,j)*BrFieldGradyT(i,j)+
     *                 BrFieldGradzT(i,j)*BrFieldGradzT(i,j))
          BcosThetaF(i,j)=dabs(BrFieldGradxT(i,j)*rpfx+
     *                    BrFieldGradyT(i,j)*rpfy+
     *                    BrFieldGradzT(i,j)*rpfz)/(rpf*rGradf+tiny)
c
        enddo
      enddo
c
      return
      end