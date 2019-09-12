c
c#############################################################################################
      SUBROUTINE CorrectPresureGradientAtWalls
     *              (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
c
      use BoundaryConditions2, only:IWallSlip,IWallnoSlip,
     *                              IWallSlipOwner,IWallSlipNBFaces,
     *                              IWallSlipNumberOfBCSets,
     *                              IWallnoSlipOwner,IWallnoSlipNBFaces,
     *                              IWallnoSlipNumberOfBCSets
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,i1,i2,i3
      double precision :: nx,ny,nz
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c*********************************************************************************************
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
          BdfidxT(i2,i3)=
     *         (1-nx*nx)*dfidxT(i1)-nx*ny*dfidyT(i1)-nx*nz*dfidzT(i1)
          BdfidyT(i2,i3)=
     *         (1-ny*ny)*dfidyT(i1)-ny*nx*dfidxT(i1)-ny*nz*dfidzT(i1)
          BdfidzT(i2,i3)=
     *         (1-nz*nz)*dfidzT(i1)-nz*nx*dfidxT(i1)-nz*ny*dfidyT(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IWallnoSlip
c
          i1=IWallnoSlipOwner(i)
          i2=IWallnoSlipNumberOfBCSets(i)
          i3=IWallnoSlipNBFaces(i)
c
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
          BdfidxT(i2,i3)=
     *         (1-nx*nx)*dfidxT(i1)-nx*ny*dfidyT(i1)-nx*nz*dfidzT(i1)
          BdfidyT(i2,i3)=
     *         (1-ny*ny)*dfidyT(i1)-ny*nx*dfidxT(i1)-ny*nz*dfidzT(i1)
          BdfidzT(i2,i3)=
     *         (1-nz*nz)*dfidzT(i1)-nz*nx*dfidxT(i1)-nz*ny*dfidyT(i1)
c
        enddo
c
       return
       end
c
c#############################################################################################
      SUBROUTINE CorrectGradientAtPeriodicBoundary
     *              (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
c
      use BoundaryConditions2, only: Iperiodic,IperiodicOwner,
     *                               IperiodicNumberOfBCSets,
     *                               IperiodicNBFaces,PeriodicPair,
     *                               Icorrespondingface,a1r,b1r,c1r,
     *                               a2r,b2r,c2r,a3r,b3r,c3r,
     *                   theta,xTranslation,yTranslation,zTranslation,
     *              LRotationalPeriodicity,LTranslationalPeriodicity
      use VAriables1, only: uVelocity,vVelocity,wVelocity,BuVelocity,
     *                      BvVelocity,BwVelocity
      use Geometry4, only: xc,yc,zc,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz 
      use Geometry3, only: NBFaceOwner
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,i1,i2,i3,j1,j2,j3
      double precision :: xF1,yF1,zF1,distance1,distance2,
     *                    GFactCF,uF1,vF1,wF1
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c*********************************************************************************************
c
      if(LRotationalPeriodicity) then
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
            yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
            zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c        
            uF1=a1r(j2)*dfidxT(j1)+b1r(j2)*dfidyT(j1)+c1r(j2)*dfidzT(j1)
            vF1=a2r(j2)*dfidxT(j1)+b2r(j2)*dfidyT(j1)+c2r(j2)*dfidzT(j1)
            wF1=a3r(j2)*dfidxT(j1)+b3r(j2)*dfidyT(j1)+c3r(j2)*dfidzT(j1)
c
            BdfidxT(i2,i3)=GFactCF*dfidxT(i1)+(1.-GFactCF)*uF1
            BdfidyT(i2,i3)=GFactCF*dfidyT(i1)+(1.-GFactCF)*vF1
            BdfidzT(i2,i3)=GFactCF*dfidzT(i1)+(1.-GFactCF)*wF1
c
            BdfidxT(j2,j3)=a1r(i2)*BdfidxT(i2,i3)+
     *                  b1r(i2)*BdfidyT(i2,i3)+c1r(i2)*BdfidzT(i2,i3)
            BdfidyT(j2,j3)=a2r(i2)*BdfidxT(i2,i3)+
     *                  b2r(i2)*BdfidyT(i2,i3)+c2r(i2)*BdfidzT(i2,i3)
            BdfidzT(j2,j3)=a3r(i2)*BdfidxT(i2,i3)+
     *                  b3r(i2)*BdfidyT(i2,i3)+c3r(i2)*BdfidzT(i2,i3)
c
          endif
c
        enddo
c
      elseif(LTranslationalPeriodicity) then
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            xF1=xc(j1)+xTranslation(j2)
            yF1=yc(j1)+yTranslation(j2)
            zF1=zc(j1)+zTranslation(j2)
c
            distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
            distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            BdfidxT(i2,i3)=GFactCF*dfidxT(i1)+(1.-GFactCF)*dfidxT(j1)
            BdfidyT(i2,i3)=GFactCF*dfidyT(i1)+(1.-GFactCF)*dfidyT(j1)
            BdfidzT(i2,i3)=GFactCF*dfidzT(i1)+(1.-GFactCF)*dfidzT(j1)
c
            BdfidxT(j2,j3)=BdfidxT(i2,i3)
            BdfidyT(j2,j3)=BdfidyT(i2,i3)
c
          endif
c
        enddo
c
      endif
c
      return
      end
c
c#############################################################################################
      SUBROUTINE CorrectGradientAtOulet(dfidxT,dfidyT,dfidzT,
     *                                    BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
      use Geometry1, Only: NumberOfBCSets
      use Geometry3, Only: NBFaces,NBFaceOwner
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
      use BoundaryConditions1, only: BCType,outletTypeMomentum
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k
      double precision :: nx,ny,nz
c--------------------------------------------------------------------------
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c*********************************************************************************************
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          if(BCType(i,j).eq.'outlet') then
c
            if(outletTypeMomentum(i,j).eq.'specifiedstaticpressure'.or.
     *                  outletTypeMomentum(i,j).eq.
     *                   'specifiedaveragestaticpressure'.or.
     *           outletTypeMomentum(i,j).eq.'specifiedresistance'.or.
     *           outletTypeMomentum(i,j).eq.'fullydeveloped')then
c              
              nx=BFaceAreanx(i,j)
              ny=BFaceAreany(i,j)
              nz=BFaceAreanz(i,i)
c
              BdfidxT(i,j)=dfidxT(k)*(1.-nx*nx)-
     *                       dfidyT(k)*nx*ny-dfidzT(k)*nx*nz
              BdfidyT(i,j)=dfidyT(k)*(1.-ny*ny)-
     *                       dfidxT(k)*ny*nx-dfidzT(k)*ny*nz
              BdfidzT(i,j)=dfidzT(k)*(1.-nz*nz)-
     *                       dfidxT(k)*nz*nx-dfidyT(k)*nz*ny
c
            endif
c
          endif
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
      SUBROUTINE CorrectPressureGradientAtSymmetryBoundary
     *               (dfidxT,dfidyT,dfidzT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
      use Geometry1, Only: NumberOfBCSets
      use Geometry3, Only: NBFaces,NBFaceOwner
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz
      use BoundaryConditions1, only: BCType,outletTypeMomentum
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k
      double precision :: nx,ny,nz
c--------------------------------------------------------------------------
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c*********************************************************************************************
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          if(BCType(i,j).eq.'symmetry') then
c
            nx=BFaceAreanx(i,j)
            ny=BFaceAreany(i,j)
            nz=BFaceAreanz(i,i)
c
            BdfidxT(i,j)=dfidxT(k)*(1.-nx*nx)-
     *                       dfidyT(k)*nx*ny-dfidzT(k)*nx*nz
            BdfidyT(i,j)=dfidyT(k)*(1.-ny*ny)-
     *                       dfidxT(k)*ny*nx-dfidzT(k)*ny*nz
            BdfidzT(i,j)=dfidzT(k)*(1.-nz*nz)-
     *                       dfidxT(k)*nz*nx-dfidyT(k)*nz*ny
c
          endif
c
        enddo
      enddo
c
      return
      end