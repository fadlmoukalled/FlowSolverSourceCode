c
c#############################################################################################
      SUBROUTINE StandardGreenGauss
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets,
     *                     ListOfElementNodes,x,y,z,NumbOfElementFaces
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor,NBFaces,
     *                     NBFaceOwner,NElementFaces
      use Geometry4
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k
      double precision Fif,FifSfx,Fifsfy,Fifsfz
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c*********************************************************************************************
c
      dfidxT=0.
      dfidyT=0.
      dfidzT=0.
      BdfidxT=0.
      BdfidyT=0.
      BdfidzT=0.
c
c--- Internal faces
c
      do i=1,NIFaces
c
        j=NIFaceOwner(i)
        k=NIFaceNeighbor(i) 
c
        Fif=FiT(j)*GFactorCF(i)+FiT(k)*(1.-GFactorCF(i))
        FifSfx=Fif*FaceAreax(i)
        Fifsfy=Fif*FaceAreay(i)
        Fifsfz=Fif*FaceAreaz(i)
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
      return
      end