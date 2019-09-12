C
C#############################################################################################
      SUBROUTINE ProcessGeometry
C#############################################################################################
c
      use Geometry1
      use Geometry3
      use Geometry4
      use Geometry5
      use User0
      use Constants1, only: big, tiny
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,j1,j2,j3,j4,k,k1,n,m,i1,i2
      double precision :: distance1,distance2,DotProduct
      double precision :: x1,y1,z1,x2,y2,z2,area1,area2
      double precision :: sx1,sy1,sz1,sx2,sy2,sz2
      double precision :: AreaBase,AreaBasex,AreaBasey,AreaBasez
      double precision :: dx,dy,dz,Productx,Producty,Productz
      double precision :: xcpyramid,ycpyramid,zcpyramid,VolumeTemp
      double precision :: centerx,centery,centerz
c*********************************************************************************************
c
      allocate(xc(NumberOfElements))
      allocate(yc(NumberOfElements))
      allocate(zc(NumberOfElements))
c
      allocate(Volume(NumberOfElements))
c
      !allocate(EdgeCentroidx(NIEdges))
      !allocate(EdgeCentroidy(NIEdges))
      !allocate(EdgeCentroidz(NIEdges))
c
      allocate(FaceCentroidx(NIFaces))
      allocate(FaceCentroidy(NIFaces))
      allocate(FaceCentroidz(NIFaces))
c
      allocate(FaceArea(NIFaces))
c
      allocate(FaceAreax(NIFaces))
      allocate(FaceAreay(NIFaces))
      allocate(FaceAreaz(NIFaces))
c
      allocate(FaceAreanx(NIFaces))
      allocate(FaceAreany(NIFaces))
      allocate(FaceAreanz(NIFaces))
c
      allocate(DistanceCF(NIFaces))
c
      allocate(DistanceCFx(NIFaces))
      allocate(DistanceCFy(NIFaces))
      allocate(DistanceCFz(NIFaces))
c
      allocate(DistanceCFux(NIFaces))
      allocate(DistanceCFuy(NIFaces))
      allocate(DistanceCFuz(NIFaces))
c
      allocate(GFactorCF(NIFaces))
c
c--- Calculate edge centroid
c
!      do k=1,NIEdges
!c
!        i=EdgeNode1(k)
!        j=EdgeNode2(k)
!c
!        EdgeCentroidx(k)=0.5*(x(i)+x(j))
!        EdgeCentroidy(k)=0.5*(y(i)+y(j))
!        EdgeCentroidz(k)=0.5*(z(i)+z(j))
!c        
!      enddo
c
c--- Calculate area and centroid of interior faces
c
      do i=1,NIFaces
c
        centerx=0.
        centery=0.
        centerz=0.
c
        do j=1,GlobalFaceNumberOfNodes(i)
c
          j1=NIFaceNodes(i,j)
c
          centerx=centerx+x(j1)     
          centery=centery+y(j1)     
          centerz=centerz+z(j1) 
c
        enddo
c          
        centerx=centerx/GlobalFaceNumberOfNodes(i)           
        centery=centery/GlobalFaceNumberOfNodes(i)           
        centerz=centerz/GlobalFaceNumberOfNodes(i)           
c
        FaceAreax(i)=0.
        FaceAreay(i)=0.
        FaceAreaz(i)=0.
        FaceArea(i)=0.
c          
        FaceCentroidx(i)=0.
        FaceCentroidy(i)=0.
        FaceCentroidz(i)=0.
c
        do j=1,GlobalFaceNumberOfNodes(i)
c
          j2=NIFaceNodes(i,j)
c
          if(j.lt.GlobalFaceNumberOfNodes(i)) then
            j3=NIFaceNodes(i,j+1)
          else
            j3=NIFaceNodes(i,1)
          endif
c          
          sx1=0.5*((y(j2)-centery)*(z(j3)-centerz)-
     *                   (y(j3)-centery)*(z(j2)-centerz))
          sy1=-0.5*((x(j2)-centerx)*(z(j3)-centerz)-
     *                   (x(j3)-centerx)*(z(j2)-centerz))
          sz1=0.5*((x(j2)-centerx)*(y(j3)-centery)-
     *                   (x(j3)-centerx)*(y(j2)-centery))
          area1=dsqrt(sx1*sx1+sy1*sy1+sz1*sz1)
c        
          x1=(centerx+x(j2)+x(j3))/3.
          y1=(centery+y(j2)+y(j3))/3.
          z1=(centerz+z(j2)+z(j3))/3.
c
          FaceAreax(i)=FaceAreax(i)+sx1
          FaceAreay(i)=FaceAreay(i)+sy1
          FaceAreaz(i)=FaceAreaz(i)+sz1
          FaceArea(i)=FaceArea(i)+area1
c          
          FaceCentroidx(i)=FaceCentroidx(i)+x1*area1
          FaceCentroidy(i)=FaceCentroidy(i)+y1*area1
          FaceCentroidz(i)=FaceCentroidz(i)+z1*area1
c
        enddo          
c
        FaceCentroidx(i)=FaceCentroidx(i)/FaceArea(i)
        FaceCentroidy(i)=FaceCentroidy(i)/FaceArea(i)
        FaceCentroidz(i)=FaceCentroidz(i)/FaceArea(i)
c
      enddo          
c
c--- Calculate area and centroid of boundary faces
c
      n=-1
      do i=1,NumberOfBCSets
        n=max(n,NBFaces(i))
      enddo
c
      j1=NumberOfBCSets
c
      allocate(BFaceCentroidx(j1,n))
      allocate(BFaceCentroidy(j1,n))
      allocate(BFaceCentroidz(j1,n))
c
      allocate(BFaceArea(j1,n))
      allocate(BFaceAreax(j1,n))
      allocate(BFaceAreay(j1,n))
      allocate(BFaceAreaz(j1,n))
      allocate(BFaceAreanx(j1,n))
      allocate(BFaceAreany(j1,n))
      allocate(BFaceAreanz(j1,n))
c
      allocate(BDistanceCF(j1,n))
      allocate(BDistanceCFx(j1,n))
      allocate(BDistanceCFy(j1,n))
      allocate(BDistanceCFz(j1,n))
      allocate(BDistanceCFux(j1,n))
      allocate(BDistanceCFuy(j1,n))
      allocate(BDistanceCFuz(j1,n))
c
      k1=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          centerx=0.
          centery=0.
          centerz=0.
          k1=k1+1
c
          do k=1,GlobalFaceNumberOfNodes(k1)
c
            j1=NBFaceNodes(i,j,k)
c
            centerx=centerx+x(j1)
            centery=centery+y(j1)
            centerz=centerz+z(j1)
c
          enddo
c          
          centerx=centerx/GlobalFaceNumberOfNodes(k1)           
          centery=centery/GlobalFaceNumberOfNodes(k1)           
          centerz=centerz/GlobalFaceNumberOfNodes(k1)           
c
          BFaceAreax(i,j)=0.
          BFaceAreay(i,j)=0.
          BFaceAreaz(i,j)=0.
          BFaceArea(i,j)=0.
c          
          BFaceCentroidx(i,j)=0.
          BFaceCentroidy(i,j)=0.
          BFaceCentroidz(i,j)=0.
c
          do k=1,GlobalFaceNumberOfNodes(k1)
c
            j2=NBFaceNodes(i,j,k)
c
            if(k.lt.GlobalFaceNumberOfNodes(k1)) then
              j3=NBFaceNodes(i,j,k+1)
            else
              j3=NBFaceNodes(i,j,1)
            endif
c          
            sx1=0.5*((y(j2)-centery)*(z(j3)-centerz)-
     *                   (y(j3)-centery)*(z(j2)-centerz))
            sy1=-0.5*((x(j2)-centerx)*(z(j3)-centerz)-
     *                   (x(j3)-centerx)*(z(j2)-centerz))
            sz1=0.5*((x(j2)-centerx)*(y(j3)-centery)-
     *                   (x(j3)-centerx)*(y(j2)-centery))
            area1=dsqrt(sx1*sx1+sy1*sy1+sz1*sz1)
c        
            x1=(centerx+x(j2)+x(j3))/3.
            y1=(centery+y(j2)+y(j3))/3.
            z1=(centerz+z(j2)+z(j3))/3.
c
            BFaceAreax(i,j)=BFaceAreax(i,j)+sx1
            BFaceAreay(i,j)=BFaceAreay(i,j)+sy1
            BFaceAreaz(i,j)=BFaceAreaz(i,j)+sz1
            BFaceArea(i,j)=BFaceArea(i,j)+area1
c          
            BFaceCentroidx(i,j)=BFaceCentroidx(i,j)+x1*area1
            BFaceCentroidy(i,j)=BFaceCentroidy(i,j)+y1*area1
            BFaceCentroidz(i,j)=BFaceCentroidz(i,j)+z1*area1
c
          enddo          
c
          BFaceCentroidx(i,j)=BFaceCentroidx(i,j)/BFaceArea(i,j)
          BFaceCentroidy(i,j)=BFaceCentroidy(i,j)/BFaceArea(i,j)
          BFaceCentroidz(i,j)=BFaceCentroidz(i,j)/BFaceArea(i,j)
c
        enddo     
c
      enddo     
c
c---- Calculate element volume and centroid
c
      do i=1,NumberOfElements
c
c---- Compute geometric center
c
        x1=0
        y1=0
        z1=0
c
        do j=1,NumbOfElementNodes(i)
c
          k=ListOfElementNodes(i,j)
          x1=x1+x(k)
          y1=y1+y(k)
          z1=z1+z(k)
c
        enddo
c
        x1=x1/NumbOfElementNodes(i)
        y1=y1/NumbOfElementNodes(i)
        z1=z1/NumbOfElementNodes(i)
c
        Volume(i)=0.
        Productx=0.
        Producty=0.
        Productz=0.
c
        do j=1,NumbOfElementFaces(i)
c
          if(NGlobalEFaces(i,j).gt.NIFaces) then          
c
            i1=BCSet(NGlobalEFaces(i,j)-NIFaces)
            i2=BCRecord(NGlobalEFaces(i,j)-NIFaces)
            AreaBase=BFaceArea(i1,i2)
            AreaBasex=BFaceAreax(i1,i2)
            AreaBasey=BFaceAreay(i1,i2)
            AreaBasez=BFaceAreaz(i1,i2)
            x2=BFaceCentroidx(i1,i2)
            y2=BFaceCentroidy(i1,i2)
            z2=BFaceCentroidz(i1,i2)
c
          else
c
            AreaBase=FaceArea(NGlobalEFaces(i,j))
            AreaBasex=FaceAreax(NGlobalEFaces(i,j))
            AreaBasey=FaceAreay(NGlobalEFaces(i,j))
            AreaBasez=FaceAreaz(NGlobalEFaces(i,j))
            x2=FaceCentroidx(NGlobalEFaces(i,j))
            y2=FaceCentroidy(NGlobalEFaces(i,j))
            z2=FaceCentroidz(NGlobalEFaces(i,j))
c
          endif
c
          dx=x2-x1
          dy=y2-y1
          dz=z2-z1
c
          xcpyramid=0.75*x2+0.25*x1
          ycpyramid=0.75*y2+0.25*y1
          zcpyramid=0.75*z2+0.25*z1
c
          VolumeTemp=dabs(AreaBasex*dx+AreaBasey*dy+AreaBasez*dz)/3.
          Volume(i)=Volume(i)+VolumeTemp
c
          Productx=Productx+xcpyramid*VolumeTemp
          Producty=Producty+ycpyramid*VolumeTemp
          Productz=Productz+zcpyramid*VolumeTemp
c
        enddo
c
        xc(i)=Productx/Volume(i)
        yc(i)=Producty/Volume(i)
        zc(i)=Productz/Volume(i)
c
      enddo

!      do i=1,NumberOfElements
!        if(volume(i).lt.0.) print*, i,volume(i)
!      enddo
!      do i=1,NumberOfElements  
!      write(90,*) i,xc(i),yc(i),zc(i)
!      enddo
c      stop
c       pause
c
c--- Calculate distance and unit vector components between owner and neighbor
c
      do k=1,NIFaces
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        DistanceCFx(k)=xc(j)-xc(i)
        DistanceCFy(k)=yc(j)-yc(i)
        DistanceCFz(k)=zc(j)-zc(i)
c
        DistanceCF(k)=dsqrt(DistanceCFx(k)**2+
     *                DistanceCFy(k)**2+DistanceCFz(k)**2)
c
        DistanceCFux(k)=DistanceCFx(k)/DistanceCF(k)
        DistanceCFuy(k)=DistanceCFy(k)/DistanceCF(k)
        DistanceCFuz(k)=DistanceCFz(k)/DistanceCF(k)
c
      enddo
c
c--- Calculate area unit vector components
c
      do k=1,NIFaces
c
        DotProduct= DistanceCFx(k)*FaceAreax(k)+
     *         DistanceCFy(k)*FaceAreay(k)+ DistanceCFz(k)*FaceAreaz(k)
c
        if(DotProduct.lt.0.) then
c
          FaceAreax(k)=-FaceAreax(k)
          FaceAreay(k)=-FaceAreay(k)
          FaceAreaz(k)=-FaceAreaz(k)
c
        endif
c
        FaceAreanx(k)=FaceAreax(k)/FaceArea(k)
        FaceAreany(k)=FaceAreay(k)/FaceArea(k)
        FaceAreanz(k)=FaceAreaz(k)/FaceArea(k)
c        
      enddo
c
c--- Calculate interpolation factors  
c
      do k=1,NIFaces
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        distance1=dsqrt((FaceCentroidx(k)-xc(i))**2+
     *         (FaceCentroidy(k)-yc(i))**2+(FaceCentroidz(k)-zc(i))**2)
        distance2=dsqrt((FaceCentroidx(k)-xc(j))**2+
     *         (FaceCentroidy(k)-yc(j))**2+(FaceCentroidz(k)-zc(j))**2)
c
        GFactorCF(k)=distance2/(distance1+distance2)
c        
      enddo
c
c--- Calculate boundary area unit vector components
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          n=NBFaceOwner(i,j)
c
          BDistanceCFx(i,j)=BFaceCentroidx(i,j)-xc(n)
          BDistanceCFy(i,j)=BFaceCentroidy(i,j)-yc(n)
          BDistanceCFz(i,j)=BFaceCentroidz(i,j)-zc(n)
c
          BDistanceCF(i,j)=dsqrt(BDistanceCFx(i,j)**2+
     *               BDistanceCFy(i,j)**2+BDistanceCFz(i,j)**2)
c
          BDistanceCFux(i,j)=BDistanceCFx(i,j)/BDistanceCF(i,j)
          BDistanceCFuy(i,j)=BDistanceCFy(i,j)/BDistanceCF(i,j)
          BDistanceCFuz(i,j)=BDistanceCFz(i,j)/BDistanceCF(i,j)
c
          DotProduct= BDistanceCFx(i,j)*BFaceAreax(i,j)+
     *                    BDistanceCFy(i,j)*BFaceAreay(i,j)+ 
     *                          BDistanceCFz(i,j)*BFaceAreaz(i,j)   
c
          if(DotProduct.lt.0.) then
c
            BFaceAreax(i,j)=-BFaceAreax(i,j)
            BFaceAreay(i,j)=-BFaceAreay(i,j)
            BFaceAreaz(i,j)=-BFaceAreaz(i,j)
c
          endif
c
          BFaceAreanx(i,j)=BFaceAreax(i,j)/BFaceArea(i,j)
          BFaceAreany(i,j)=BFaceAreay(i,j)/BFaceArea(i,j)
          BFaceAreanz(i,j)=BFaceAreaz(i,j)/BFaceArea(i,j)
c 
        enddo
      enddo
c
c--- Decompose the Surface vector into E and T along internal faces
c
!    MethodDecomposeS=1 Minimum Correction Approach
!    MethodDecomposeS=2 Orthogonal Correction Approach
!    MethodDecomposeS=3 Overrelaxed Correction Approach
c
      allocate(FaceE(NIFaces))
      allocate(FaceEx(NIFaces))
      allocate(FaceEy(NIFaces))
      allocate(FaceEz(NIFaces))
      allocate(gDiff(NIFaces))
      allocate(FaceT(NIFaces))
      allocate(FaceTx(NIFaces))
      allocate(FaceTy(NIFaces))
      allocate(FaceTz(NIFaces))
c
      if(MethodDecomposeS.eq.1) then
c
        do k=1,NIFaces
c
          DotProduct=DistanceCFux(k)*FaceAreax(k)+
     *                 DistanceCFuy(k)*FaceAreay(k)+
     *                     DistanceCFuz(k)*FaceAreaz(k)
          FaceEx(k)=DotProduct*DistanceCFux(k)
          FaceEy(k)=DotProduct*DistanceCFuy(k)
          FaceEz(k)=DotProduct*DistanceCFuz(k)
          FaceE(k)=dabs(DotProduct)
          gDiff(k)=FaceE(k)/DistanceCF(k)
c
          FaceTx(k)=FaceAreax(k)-FaceEx(k)
          FaceTy(k)=FaceAreay(k)-FaceEy(k)
          FaceTz(k)=FaceAreaz(k)-FaceEz(k)
          FaceT(k)=dsqrt(FaceTx(k)**2+FaceTy(k)**2+FaceTz(k)**2)

        enddo
c
      elseif(MethodDecomposeS.eq.2) then
c
        do k=1,NIFaces
c
          FaceEx(k)=FaceArea(k)*DistanceCFux(k)
          FaceEy(k)=FaceArea(k)*DistanceCFuy(k)
          FaceEz(k)=FaceArea(k)*DistanceCFuz(k)
          FaceE(k)=FaceArea(k)
          gDiff(k)=FaceE(k)/DistanceCF(k)
c
          FaceTx(k)=FaceAreax(k)-FaceEx(k)
          FaceTy(k)=FaceAreay(k)-FaceEy(k)
          FaceTz(k)=FaceAreaz(k)-FaceEz(k)
          FaceT(k)=dsqrt(FaceTx(k)**2+FaceTy(k)**2+FaceTz(k)**2)

        enddo
c
      elseif(MethodDecomposeS.eq.3) then
c
        do k=1,NIFaces
c
          DotProduct=DistanceCFux(k)*FaceAreax(k)+
     *                   DistanceCFuy(k)*FaceAreay(k)+
     *                           DistanceCFuz(k)*FaceAreaz(k)
          FaceEx(k)=(FaceArea(k)**2/DotProduct)*DistanceCFux(k)
          FaceEy(k)=(FaceArea(k)**2/DotProduct)*DistanceCFuy(k)
          FaceEz(k)=(FaceArea(k)**2/DotProduct)*DistanceCFuz(k)
          FaceE(k)=dsqrt(FaceEx(k)**2+FaceEy(k)**2+FaceEz(k)**2)
          gDiff(k)=FaceE(k)/DistanceCF(k)
c
          FaceTx(k)=FaceAreax(k)-FaceEx(k)
          FaceTy(k)=FaceAreay(k)-FaceEy(k)
          FaceTz(k)=FaceAreaz(k)-FaceEz(k)
          FaceT(k)=dsqrt(FaceTx(k)**2+FaceTy(k)**2+FaceTz(k)**2)

        enddo
c
      endif
c
c--- Decompose the Surface vector into E and T along boundary faces
c
      n=-1
      do i=1,NumberOfBCSets
        n=max(n,NBFaces(i))
      enddo
c
      j1=NumberOfBCSets
c
      allocate(BFaceE(j1,n))
      allocate(BFaceEx(j1,n))
      allocate(BFaceEy(j1,n))
      allocate(BFaceEz(j1,n))
      allocate(BgDiff(j1,n))
      allocate(BFaceT(j1,n))
      allocate(BFaceTx(j1,n))
      allocate(BFaceTy(j1,n))
      allocate(BFaceTz(j1,n))
c
      if(MethodDecomposeS.eq.1) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            DotProduct=BDistanceCFux(i,j)*BFaceAreax(i,j)+
     *                   BDistanceCFuy(i,j)*BFaceAreay(i,j)+
     *                     BDistanceCFuz(i,j)*BFaceAreaz(i,j)
            BFaceEx(i,j)=DotProduct*BDistanceCFux(i,j)
            BFaceEy(i,j)=DotProduct*BDistanceCFuy(i,j)
            BFaceEz(i,j)=DotProduct*BDistanceCFuz(i,j)
            BFaceE(i,j)=dabs(DotProduct)
            BgDiff(i,j)=BFaceE(i,j)/BDistanceCF(i,j)
c
            BFaceTx(i,j)=BFaceAreax(i,j)-BFaceEx(i,j)
            BFaceTy(i,j)=BFaceAreay(i,j)-BFaceEy(i,j)
            BFaceTz(i,j)=BFaceAreaz(i,j)-BFaceEz(i,j)
            BFaceT(i,j)=dsqrt(BFaceTx(i,j)**2+
     *                          BFaceTy(i,j)**2+BFaceTz(i,j)**2)

          enddo
        enddo
c
      elseif(MethodDecomposeS.eq.2) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BFaceEx(i,j)=BFaceArea(i,j)*BDistanceCFux(i,j)
            BFaceEy(i,j)=BFaceArea(i,j)*BDistanceCFuy(i,j)
            BFaceEz(i,j)=BFaceArea(i,j)*BDistanceCFuz(i,j)
            BFaceE(i,j)=BFaceArea(i,j)
            BgDiff(i,j)=BFaceE(i,j)/BDistanceCF(i,j)
c
            BFaceTx(i,j)=BFaceAreax(i,j)-BFaceEx(i,j)
            BFaceTy(i,j)=BFaceAreay(i,j)-BFaceEy(i,j)
            BFaceTz(i,j)=BFaceAreaz(i,j)-BFaceEz(i,j)
            BFaceT(i,j)=dsqrt(BFaceTx(i,j)**2+
     *                          BFaceTy(i,j)**2+BFaceTz(i,j)**2)
c
          enddo
        enddo
c
      elseif(MethodDecomposeS.eq.3) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            DotProduct=BDistanceCFux(i,j)*BFaceAreax(i,j)+
     *                     BDistanceCFuy(i,j)*BFaceAreay(i,j)+
     *                         BDistanceCFuz(i,j)*BFaceAreaz(i,j)
            BFaceEx(i,j)=(BFaceArea(i,j)**2/DotProduct)*
     *                                      BDistanceCFux(i,j)
            BFaceEy(i,j)=(BFaceArea(i,j)**2/DotProduct)*
     *                                      BDistanceCFuy(i,j)
            BFaceEz(i,j)=(BFaceArea(i,j)**2/DotProduct)*
     *                                      BDistanceCFuz(i,j)
            BFaceE(i,j)=dsqrt(BFaceEx(i,j)**2+
     *                         BFaceEy(i,j)**2+BFaceEz(i,j)**2)
            BgDiff(i,j)=BFaceE(i,j)/BDistanceCF(i,j)
c
            BFaceTx(i,j)=BFaceAreax(i,j)-BFaceEx(i,j)
            BFaceTy(i,j)=BFaceAreay(i,j)-BFaceEy(i,j)
            BFaceTz(i,j)=BFaceAreaz(i,j)-BFaceEz(i,j)
            BFaceT(i,j)=dsqrt(BFaceTx(i,j)**2+
     *                         BFaceTy(i,j)**2+BFaceTz(i,j)**2)
c
          enddo
        enddo
c
      endif
c
c---- Calculate minum and maximum quantities (Volume, area, distance,...)
c
      totalVolume=0.
      minimumVolume=big
      maximumVolume=-tiny
      do i=1,NumberOfElements
        minimumVolume=dmin1(minimumVolume,Volume(i))
        maximumVolume=dmax1(maximumVolume,Volume(i))
        totalVolume=totalVolume+Volume(i)
      enddo
c
      totalFaceArea=0.
      totalBoundaryArea=0.
      minimumFaceArea=big
      maximumFaceArea=-tiny
      do i=1,NIFaces
        minimumFaceArea=dmin1(minimumFaceArea,FaceArea(i))
        maximumFaceArea=dmax1(maximumFaceArea,FaceArea(i))
        totalFaceArea=totalFaceArea+FaceArea(i)
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
          minimumFaceArea=dmin1(minimumFaceArea,BFaceArea(i,j))
          maximumFaceArea=dmax1(maximumFaceArea,BFaceArea(i,j))
          totalFaceArea=totalFaceArea+BFaceArea(i,j)
          totalBoundaryArea=totalBoundaryArea+BFaceArea(i,j)
        enddo
      enddo
c
      minimumDistanceCF=big
      maximumDistanceCF=-tiny
      do i=1,NIFaces
        minimumDistanceCF=dmin1(minimumDistanceCF,DistanceCF(i))
        maximumDistanceCF=dmax1(maximumDistanceCF,DistanceCF(i))
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
        minimumDistanceCF=dmin1(minimumDistanceCF,BDistanceCF(i,j))
        maximumDistanceCF=dmax1(maximumDistanceCF,BDistanceCF(i,j))
        enddo
      enddo
c
      minimumX=big
      maximumX=-tiny
      minimumY=big
      maximumY=-tiny
      minimumZ=big
      maximumZ=-tiny
      do i=1,NumberOfNodes
        minimumX=dmin1(minimumX,x(i))
        maximumX=dmax1(maximumX,x(i))
        minimumY=dmin1(minimumY,y(i))
        maximumY=dmax1(maximumY,y(i))
        minimumZ=dmin1(minimumZ,z(i))
        maximumZ=dmax1(maximumZ,z(i))
      ENDDO
c
      write(11,*) 'domain extension:'
      write(11,*) minimumX," < x < ",maximumX
      write(11,*) minimumY," < y < ",maximumY
      write(11,*) minimumZ," < z < ",maximumZ
c
      write(11,*) 'minimum volume = ',minimumVolume
      write(11,*) 'maximum volume = ',maximumVolume
      write(11,*) 'total volume = ',totalVolume
c
      write(11,*) 'minimum FaceArea = ',minimumFaceArea
      write(11,*) 'maximum FaceArea = ',maximumFaceArea
      write(11,*) 'total FaceArea = ',totalFaceArea
      write(11,*) 'total Boundary Area = ',totalBoundaryArea
c
      write(11,*) 'minimum Distance between C and F= ',minimumDistanceCF
      write(11,*) 'maximum Distance between C and F= ',maximumDistanceCF
c
      Print*, 'domain extent:'
      Print*, minimumX," < x < ",maximumX
      Print*, minimumY," < y < ",maximumY
      Print*, minimumZ," < z < ",maximumZ
c
      Print*, 'minimum volume = ',minimumVolume
      Print*, 'maximum volume = ',maximumVolume
      Print*, 'total volume = ',totalVolume
c
      Print*, 'minimum FaceArea = ',minimumFaceArea
      Print*, 'maximum FaceArea = ',maximumFaceArea
      Print*, 'total FaceArea = ',totalFaceArea
      Print*, 'total Boundary Area = ',totalBoundaryArea
c
      Print*, 'minimum Distance between C and F= ',minimumDistanceCF
      Print*, 'maximum Distance between C and F= ',maximumDistanceCF
c
      return
c
      end