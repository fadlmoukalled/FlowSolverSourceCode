c
C#############################################################################################
      SUBROUTINE InterpolateGradientToFace
     *   (GradientInterpolationScheme,FiT,dfidxT,dfidyT,dfidzT,
     *                                  dfidxfT,dfidyfT,dfidzfT)
C#############################################################################################
      use variables1, only: mdot
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor
      use Geometry4, only: DistanceCF,DistanceCFux,DistanceCFuy,
     *                     DistanceCFuz,GFactorCF
      use PhysicalProperties1, only: Densityf,Density 
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k
      double precision gf,correction
      character*16 GradientInterpolationScheme
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:) :: dfidxfT
      double precision, dimension(:) :: dfidyfT
      double precision, dimension(:) :: dfidzfT
c*********************************************************************************************
c
      if(GradientInterpolationScheme.eq.'average') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
          gf=GFactorCF(k)
c
          dfidxfT(k)=gf*dfidxT(i)+(1.-gf)*dfidxT(j)
          dfidyfT(k)=gf*dfidyT(i)+(1.-gf)*dfidyT(j)
          dfidzfT(k)=gf*dfidzT(i)+(1.-gf)*dfidzT(j)
c
        enddo
c
      elseif(GradientInterpolationScheme.eq.'upwind') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gf=0.
          if(mdot(k).gt.0.) gf=1.          
c
          dfidxfT(k)=gf*dfidxT(i)+(1.-gf)*dfidxT(j)
          dfidyfT(k)=gf*dfidyT(i)+(1.-gf)*dfidyT(j)
          dfidzfT(k)=gf*dfidzT(i)+(1.-gf)*dfidzT(j)
c
        enddo
c
      elseif(GradientInterpolationScheme.eq.'downwind') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gf=1.
          if(mdot(k).gt.0.) gf=0.          
c
          dfidxfT(k)=gf*dfidxT(i)+(1.-gf)*dfidxT(j)
          dfidyfT(k)=gf*dfidyT(i)+(1.-gf)*dfidyT(j)
          dfidzfT(k)=gf*dfidzT(i)+(1.-gf)*dfidzT(j)
c
        enddo
c
      elseif(GradientInterpolationScheme.eq.'averagecorrected') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
          gf=GFactorCF(k)
c
          dfidxfT(k)=gf*dfidxT(i)+(1.-gf)*dfidxT(j)
          dfidyfT(k)=gf*dfidyT(i)+(1.-gf)*dfidyT(j)
          dfidzfT(k)=gf*dfidzT(i)+(1.-gf)*dfidzT(j)
c
          correction=(FiT(j)-FiT(i))/DistanceCF(k)
          correction=correction-(dfidxfT(k)*DistanceCFux(k)+
     *            dfidyfT(k)*DistanceCFuy(k)+dfidzfT(k)*DistanceCFuz(k))
c
          dfidxfT(k)=dfidxfT(k)+correction*DistanceCFux(k)
          dfidyfT(k)=dfidyfT(k)+correction*DistanceCFuy(k)
          dfidzfT(k)=dfidzfT(k)+correction*DistanceCFuz(k)
c
        enddo
c
      elseif(GradientInterpolationScheme.eq.'harmonic') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
          gf=GFactorCF(k)
c
          dfidxfT(k)=Densityf(k)*
     *           (gf*dfidxT(i)/Density(i)+(1.-gf)*dfidxT(j)/Density(j))
          dfidyfT(k)=Densityf(k)*
     *           (gf*dfidyT(i)/Density(i)+(1.-gf)*dfidyT(j)/Density(j))
          dfidzfT(k)=Densityf(k)*
     *           (gf*dfidzT(i)/Density(i)+(1.-gf)*dfidzT(j)/Density(j))
c
        enddo
c
      endif
c
      return
      end
c