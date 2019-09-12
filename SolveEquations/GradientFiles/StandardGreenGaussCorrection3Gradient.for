c
C#############################################################################################
      SUBROUTINE StandardGreenGaussFaceCorrection3
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,
     *                               BdfidzT,nIterGradientPhi)
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets,
     *                     ListOfElementNodes,x,y,z,NumbOfElementFaces
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor,NBFaces,
     *                     NBFaceOwner,NElementFaces
      use Geometry4
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k,nIterGradientPhi,Iter
      double precision :: phif,FifSfx,Fifsfy,Fifsfz,gc,gc1,gc2
      double precision :: rCfDOTrCF,rCFSquare,xfPrime,yfPrime,zfPrime
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c***************************************************************************
      double precision, dimension(:), allocatable :: FifPrime
      double precision, dimension(:), allocatable :: dfidxF
      double precision, dimension(:), allocatable :: dfidyF
      double precision, dimension(:), allocatable :: dfidzF
c*********************************************************************************************
c
      dfidxT=0.
      dfidyT=0.
      dfidzT=0.
      BdfidxT=0.
      BdfidyT=0.
      BdfidzT=0.
c
      allocate(FifPrime(NIFaces))
      allocate(dfidxF(NumberOfElements))
      allocate(dfidyF(NumberOfElements))
      allocate(dfidzF(NumberOfElements))
c
      FifPrime=0.
c
c--- Internal faces
c
      do i=1,NIFaces
c
        j=NIFaceOwner(i)
        k=NIFaceNeighbor(i) 
c
        rCfDOTrCF=(FaceCentroidx(i)-xc(j))*(xc(k)-xc(j))+
     *               (FaceCentroidy(i)-yc(j))*(yc(k)-yc(j))+
     *                  (FaceCentroidz(i)-zc(j))*(zc(k)-zc(j))
        rCFSquare=(xc(k)-xc(j))*(xc(k)-xc(j))+
     *                (yc(k)-yc(j))*(yc(k)-yc(j))+
     *                     (zc(k)-zc(j))*(zc(k)-zc(j))
        xfPrime=xc(j)-(rCfDOTrCF/rCFSquare)*(xc(j)-xc(k))
        yfPrime=yc(j)-(rCfDOTrCF/rCFSquare)*(yc(j)-yc(k))
        zfPrime=zc(j)-(rCfDOTrCF/rCFSquare)*(zc(j)-zc(k))
        gc1=(xc(k)-xfPrime)**2+(yc(k)-yfPrime)**2+(zc(k)-zfPrime)**2
        gc2=(xc(k)-xc(j))**2+(yc(k)-yc(j))**2+(zc(k)-zc(j))**2
        gc=dsqrt(gc1/gc2)
c
        FifPrime(i)=FiT(j)*gc+FiT(k)*(1.-gc)
c
        FifSfx=FifPrime(i)*FaceAreax(i)
        Fifsfy=FifPrime(i)*FaceAreay(i)
        Fifsfz=FifPrime(i)*FaceAreaz(i)
c
        dfidxT(j)=dfidxT(j)+FifSfx
        dfidyT(j)=dfidyT(j)+FifSfy
        dfidzT(j)=dfidzT(j)+FifSfz
c
        dfidxT(k)=dfidxT(k)-FifSfx
        dfidyT(k)=dfidyT(k)-FifSfy
        dfidzT(k)=dfidzT(k)-FifSfz
c
      enddo
c
c--- Boundary faces
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          FifSfx=BFiT(i,j)*BFaceAreax(i,j)
          FifSfy=BFiT(i,j)*BFaceAreay(i,j)
          FifSfz=BFiT(i,j)*BFaceAreaz(i,j)
c
          dfidxT(k)=dfidxT(k)+FifSfx
          dfidyT(k)=dfidyT(k)+FifSfy
          dfidzT(k)=dfidzT(k)+FifSfz
c
        enddo
      enddo
c
      do i=1,NumberOfElements
c
        dfidxT(i)=dfidxT(i)/Volume(i)
        dfidyT(i)=dfidyT(i)/Volume(i)
        dfidzT(i)=dfidzT(i)/Volume(i)
c
      enddo
c
      do Iter=1,nIterGradientPhi
c
c--- Initialize storage for updated gradient
c
        dfidxF=0.
        dfidyF=0.
        dfidzF=0.
c
c--- Internal faces
c
        do i=1,NIFaces
c
          j=NIFaceOwner(i)
          k=NIFaceNeighbor(i) 
c
          rCfDOTrCF=(FaceCentroidx(i)-xc(j))*(xc(k)-xc(j))+
     *                (FaceCentroidy(i)-yc(j))*(yc(k)-yc(j))+
     *                    (FaceCentroidz(i)-zc(j))*(zc(k)-zc(j))
          rCFSquare=(xc(k)-xc(j))*(xc(k)-xc(j))+
     *                 (yc(k)-yc(j))*(yc(k)-yc(j))+
     *                      (zc(k)-zc(j))*(zc(k)-zc(j))
          xfPrime=xc(j)-(rCfDOTrCF/rCFSquare)*(xc(j)-xc(k))
          yfPrime=yc(j)-(rCfDOTrCF/rCFSquare)*(yc(j)-yc(k))
          zfPrime=zc(j)-(rCfDOTrCF/rCFSquare)*(zc(j)-zc(k))
          gc1=(xc(k)-xfPrime)**2+(yc(k)-yfPrime)**2+(zc(k)-zfPrime)**2
          gc2=(xc(k)-xc(j))**2+(yc(k)-yc(j))**2+(zc(k)-zc(j))**2
          gc=dsqrt(gc1/gc2)
c
          FifSfx=gc*dfidxT(j)+(1.-gc)*dfidxT(k)
          FifSfy=gc*dfidyT(j)+(1.-gc)*dfidyT(k)
          FifSfz=gc*dfidzT(j)+(1.-gc)*dfidzT(k)
c
          phif=FifPrime(i)+FifSfx*(FaceCentroidx(i)-xfPrime)+
     *                      FifSfy*(FaceCentroidy(i)-yfPrime)+
     *                          FifSfz*(FaceCentroidz(i)-zfPrime)
c
          FifSfx=phif*FaceAreax(i)
          Fifsfy=phif*FaceAreay(i)
          Fifsfz=phif*FaceAreaz(i)
c
          dfidxF(j)=dfidxF(j)+FifSfx
          dfidyF(j)=dfidyF(j)+FifSfy
          dfidzF(j)=dfidzF(j)+FifSfz
c
          dfidxF(k)=dfidxF(k)-FifSfx
          dfidyF(k)=dfidyF(k)-FifSfy
          dfidzF(k)=dfidzF(k)-FifSfz
c
        enddo
c
c--- Boundary faces
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            FifSfx=BFiT(i,j)*BFaceAreax(i,j)
            FifSfy=BFiT(i,j)*BFaceAreay(i,j)
            FifSfz=BFiT(i,j)*BFaceAreaz(i,j)
c
            dfidxF(k)=dfidxF(k)+FifSfx
            dfidyF(k)=dfidyF(k)+FifSfy
            dfidzF(k)=dfidzF(k)+FifSfz
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          dfidxT(i)=dfidxF(i)/Volume(i)
          dfidyT(i)=dfidyF(i)/Volume(i)
          dfidzT(i)=dfidzF(i)/Volume(i)
c
        enddo
c
      enddo
c
c--- Update Gradients along boundaries
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          BdfidxT(i,j)=dfidxT(k)
          BdfidyT(i,j)=dfidyT(k)
          BdfidzT(i,j)=dfidzT(k)
c
        enddo
      enddo
c
      deallocate(FifPrime)
      deallocate(dfidxF)
      deallocate(dfidyF)
      deallocate(dfidzF)
c
      return
      end