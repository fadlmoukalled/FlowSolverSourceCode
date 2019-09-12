c
c#############################################################################################
c
      SUBROUTINE UpdatePeriodic
c
c#############################################################################################
c
      use BoundaryConditions2, only: Iperiodic,IperiodicOwner,
     *                               IperiodicNumberOfBCSets,
     *                               IperiodicNBFaces,PeriodicPair,
     *                               Icorrespondingface,a1r,b1r,c1r,
     *                               a2r,b2r,c2r,a3r,b3r,c3r,
     *                               theta,LRotationalPeriodicity,
     *                               LTranslationalPeriodicity,
     *                               xTranslation,yTranslation,
     *                               zTranslation
      use VAriables1, only: uVelocity,vVelocity,wVelocity,BuVelocity,
     *                      BvVelocity,BwVelocity
      use Geometry4, only: xc,yc,zc,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz 
      use Geometry3, only: NBFaceOwner
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,j1,j2,j3
      double precision :: xF1,yF1,zF1,distance1,distance2,GFactCF,
     *                    uF1,vF1,wF1
c********************************************************************************************
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
            uF1=a1r(j2)*uVelocity(j1)+b1r(j2)*vVelocity(j1)+
     *                                         c1r(j2)*wVelocity(j1)
            vF1=a2r(j2)*uVelocity(j1)+b2r(j2)*vVelocity(j1)+
     *                                         c2r(j2)*wVelocity(j1)
            wF1=a3r(j2)*uVelocity(j1)+b3r(j2)*vVelocity(j1)+
     *                                         c3r(j2)*wVelocity(j1)
c
            BuVelocity(i2,i3)=GFactCF*uVelocity(i1)+(1.-GFactCF)*uF1
            BvVelocity(i2,i3)=GFactCF*vVelocity(i1)+(1.-GFactCF)*vF1
            BwVelocity(i2,i3)=GFactCF*wVelocity(i1)+(1.-GFactCF)*wF1
c
            BuVelocity(j2,j3)=a1r(i2)*BuVelocity(i2,i3)+
     *             b1r(i2)*BvVelocity(i2,i3)+c1r(i2)*BwVelocity(i2,i3)
            BvVelocity(j2,j3)=a2r(i2)*BuVelocity(i2,i3)+
     *             b2r(i2)*BvVelocity(i2,i3)+c2r(i2)*BwVelocity(i2,i3)
            BwVelocity(j2,j3)=a3r(i2)*BuVelocity(i2,i3)+
     *             b3r(i2)*BvVelocity(i2,i3)+c3r(i2)*BwVelocity(i2,i3)
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
     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
            GFactCF=distance2/(distance1+distance2)
c
            BuVelocity(i2,i3)=GFactCF*uVelocity(i1)+
     *                                    (1.-GFactCF)*uVelocity(j1)
            BvVelocity(i2,i3)=GFactCF*vVelocity(i1)+
     *                                    (1.-GFactCF)*vVelocity(j1)
            BwVelocity(i2,i3)=GFactCF*wVelocity(i1)+
     *                                    (1.-GFactCF)*wVelocity(j1)
c
            BuVelocity(j2,j3)=BuVelocity(i2,i3)
            BvVelocity(j2,j3)=BvVelocity(i2,i3)
            BwVelocity(j2,j3)=BwVelocity(i2,i3)
c
          endif
c
        enddo
c
      endif
c
      return
      end
