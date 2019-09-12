c
C#############################################################################################
c
      SUBROUTINE LocatePointInVolume(xpt,ypt,zpt,iElement)
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumbOfElementFaces,
     *                     NumbOfElementNodes
      use Geometry3, only: BCSet,BCRecord,NIFaces,NGlobalEFaces
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz,
     *                     FaceAreanx,FaceAreany,FaceAreanz,
     *                     FaceCentroidx,FaceCentroidy,FaceCentroidz,
     *                     xc,yc,zc
      use Constants1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: iElement
      integer i,j,k,i1,i2,index
      double precision xpt,ypt,zpt,xInt,yInt,zInt,x1,y1,z1,nx,ny,nz,
     *                 d,d1,d2,dot1,vx,vy,vz,denom,num,t
c*********************************************************************************************
c
      do i=1,NumberOfElements
c
        index=1
        do j=1,NumbOfElementFaces(i)
c
          k=NGlobalEFaces(i,j)
c
          if(k.gt.NIFaces) then          
c
            i1=BCSet(k-NIFaces)
            i2=BCRecord(k-NIFaces)
c
            nx=BFaceAreanx(i1,i2)
            ny=BFaceAreany(i1,i2)
            nz=BFaceAreanz(i1,i2)
c
            x1=BFaceCentroidx(i1,i2)
            y1=BFaceCentroidy(i1,i2)
            z1=BFaceCentroidz(i1,i2)
c
          else
c
            nx=FaceAreanx(k)
            ny=FaceAreany(k)
            nz=FaceAreanz(k)
c
            x1=FaceCentroidx(k)
            y1=FaceCentroidy(k)
            z1=FaceCentroidz(k)
c
          endif
c
c---- Equation of face (plane)
c
          d=-(nx*x1+ny*y1+nz*z1)
c
c---- Equation of line joining element centroid with point source
c
          vx=xpt-xc(i)
          vy=ypt-yc(i)
          vz=zpt-zc(i)
c
          denom=nx*vx+ny*vy+nz*vz
c
          if(denom.eq.0.) then
c
            cycle
c
          else
c
            num=-(nx*xpt+ny*ypt+nz*zpt+d)
c
            t=num/denom
            xInt=xpt+t*vx
            yInt=ypt+t*vy
            zInt=zpt+t*vz
c            
            d1=(xc(i)-xInt)**2+(yc(i)-yInt)**2+(zc(i)-zInt)**2
            d2=(xpt-xc(i))**2+(ypt-yc(i))**2+(zpt-zc(i))**2
c
            dot1=(xInt-xc(i))*(xpt-xc(i))+
     *                (yInt-yc(i))*(ypt-yc(i))+(zInt-zc(i))*(zpt-zc(i))
            if(dot1.gt.0.) then
c
              if(d2.gt.d1) then
c
                index=0
                exit
c
              endif              
c              
            endif
                         
          endif
c 
        enddo
c
        if(index.eq.1) then
c
          iElement=i         
          exit
c
        endif
c
      enddo
c
      return
      end