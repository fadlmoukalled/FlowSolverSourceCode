c
c#############################################################################################
      SUBROUTINE NodeBasedGreenGauss
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
      use Geometry1, only:NumberOfNodes,ListOfElementNodes,NodeFlag,
     *                    NumberOfElements,NumbOfElementNodes,
     *                    NumbOfElementFaces,x,y,z,NumberOfBCSets,
     *                    NBCDataRecords,NElementBC,NElementBCFace
      use Geometry3, only:NIFaceOwner,NIFaceNeighbor,NIFaceNodes,
     *                    NIFaces,NBFaces,NBFaceOwner,NBFaceNodes,
     *                    NumberOfElementFaceNodes,
     *                    GlobalFaceNumberOfNodes
      use Geometry4, only:xc,yc,zc,FaceAreax,FaceAreay,FaceAreaz,
     *                    Volume,BFaceAreax,BFaceAreay,BFaceAreaz,
     *                    FaceCentroidx,FaceCentroidy,FaceCentroidz,
     *                    BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,i1,i2,j,j1,j2,k,k1,k2
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c
      double precision, dimension(:), allocatable :: SumInvd
      double precision, dimension(:), allocatable :: SumPhiInvd
      double precision, dimension(:), allocatable :: dummyN
      double precision d,FifSfx,FifSfy,FifSfz,phif,
     *                 SumInvdFace,SumPhiInvdFace
c*********************************************************************************************
c
      dfidxT=0.
      dfidyT=0.
      dfidzT=0.
      BdfidxT=0.
      BdfidyT=0.
      BdfidzT=0.
c
      allocate(SumInvd(NumberOfNodes))
      allocate(SumPhiInvd(NumberOfNodes))
      allocate(dummyN(NumberOfNodes))
c
c--- First interpolate values to internal inodes
c
      SumInvd=0.
      SumPhiInvd=0.
c
      do i=1,NumberOfElements
        do j=1,NumbOfElementNodes(i)
c
          k=ListOfElementNodes(i,j)

          if(NodeFlag(k).eq.1) then
c
            d=dsqrt((xc(i)-x(k))**2+(yc(i)-y(k))**2+(zc(i)-z(k))**2)
            SumInvd(k)=SumInvd(k)+1./d
            SumPhiInvd(k)=SumPhiInvd(k)+FiT(i)/d
c
          endif
c
        enddo
      enddo
c
      do i=1,NumberOfElements
        do j=1,NumbOfElementNodes(i)
c
          k=ListOfElementNodes(i,j)

          if(NodeFlag(k).eq.1) then
c
            dummyN(k)=SumPhiInvd(k)/SumInvd(k)
c
          endif
c
        enddo
      enddo
c
c--- Second find for boundary nodes
c
      k2=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
c          j1=NElementBC(i,j)
c          j2=NElementBCFace(i,j)
c
          k2=k2+1
          do k=1,GlobalFaceNumberOfNodes(k2)
c          do k=1,NumberOfElementFaceNodes(j1,j2)
c
            k1=NBFaceNodes(i,j,k)
            d=dsqrt((BFaceCentroidx(i,j)-x(k1))**2+
     *                   (BFaceCentroidy(i,j)-y(k1))**2+
     *                         (BFaceCentroidz(i,j)-z(k1))**2)
            SumInvd(k1)=SumInvd(k1)+1./d
            SumPhiInvd(k1)=SumPhiInvd(k1)+BFiT(i,j)/d        
c               
          enddo
c
        enddo
      enddo
c
      k2=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
c          j1=NElementBC(i,j)
c          j2=NElementBCFace(i,j)
c
          k2=k2+1
          do k=1,GlobalFaceNumberOfNodes(k2)
c          do k=1,NumberOfElementFaceNodes(j1,j2)
c
            k1=NBFaceNodes(i,j,k)
            dummyN(k1)=SumPhiInvd(k1)/SumInvd(k1)
c               
          enddo
        enddo
      enddo
c
c--- Interpolate values to face centroid
c
      do i=1,NIFaces
c
        j=NIFaceOwner(i)
        k=NIFaceNeighbor(i) 
c
        SumInvdFace=0.
        SumPhiInvdFace=0.
        do j1=1,GlobalFaceNumberOfNodes(i)
c
          j2=NIFaceNodes(i,j1)
c
          d=dsqrt((FaceCentroidx(i)-x(j2))**2+
     *        (FaceCentroidy(i)-y(j2))**2+(FaceCentroidz(i)-z(j2))**2)
          SumInvdFace=SumInvdFace+1./d
          SumPhiInvdFace=SumPhiInvdFace+dummyN(j2)/d
c
        enddo
c
        phif=SumPhiInvdFace/SumInvdFace
c
        FifSfx=phif*FaceAreax(i)
        Fifsfy=phif*FaceAreay(i)
        Fifsfz=phif*FaceAreaz(i)
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
      deallocate(SumInvd)
      deallocate(SumPhiInvd)
      deallocate(dummyN)
c
      return
      end