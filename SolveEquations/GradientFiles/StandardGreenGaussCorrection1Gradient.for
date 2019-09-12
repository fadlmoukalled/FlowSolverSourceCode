c
C#############################################################################################
      SUBROUTINE StandardGreenGaussFaceCorrection1
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
      double precision :: area,x1,x2,x3,x4,y1,y2,y3,y4
      integer :: i,j,j1,j2,j3,j4,k,nIterGradientPhi,Iter
      double precision :: Fif,FifSfx,Fifsfy,Fifsfz
      double precision :: phif,correctionx,correctiony,correctionz
      double precision :: eDOTn,rfMrCDOTn,xfPrime,yfPrime,zfPrime
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
        FifPrime(i)=FiT(j)*GFactorCF(i)+FiT(k)*(1.-GFactorCF(i))
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
          eDOTn=FaceEx(i)*FaceAreanx(i)+FaceEy(i)*FaceAreany(i)
     *                                    +FaceEz(i)*FaceAreanz(i)
          rfMrCDOTn=(FaceCentroidx(i)-xc(j))*FaceAreanx(i)+
     *                 (FaceCentroidy(i)-yc(j))*FaceAreany(i)+
     *                     (FaceCentroidz(i)-zc(j))*FaceAreanz(i)
c
          xfPrime=xc(j)+(rfMrCDOTn/eDOTn)*FaceEx(i)
          yfPrime=yc(j)+(rfMrCDOTn/eDOTn)*FaceEy(i)
          zfPrime=zc(j)+(rfMrCDOTn/eDOTn)*FaceEz(i)
c
          correctionx=dfidxT(j)*GFactorCF(i)+dfidxT(k)*(1.-GFactorCF(i))
          correctionx=correctionx*(FaceCentroidx(i)-xfPrime)
          correctiony=dfidyT(j)*GFactorCF(i)+dfidyT(k)*(1.-GFactorCF(i))
          correctiony=correctiony*(FaceCentroidy(i)-yfPrime)
          correctionz=dfidzT(j)*GFactorCF(i)+dfidzT(k)*(1.-GFactorCF(i))
          correctionz=correctionz*(FaceCentroidz(i)-zfPrime)
c
          phif=FifPrime(i)+correctionx+correctiony+correctionz
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