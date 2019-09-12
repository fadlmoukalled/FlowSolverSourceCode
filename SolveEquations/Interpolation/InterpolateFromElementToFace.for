c
C#############################################################################################
      SUBROUTINE InterpolateElementToFace(InterpolationScheme,FiT,FiTf)
C#############################################################################################
      use variables1, only: mdot
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor
      use Geometry4, only: GFactorCF
c*********************************************************************************************      
      implicit none
c*********************************************************************************************
      integer i,j,k
      double precision gf
      character*16 InterpolationScheme
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: FiTf
c*********************************************************************************************
c
      if(InterpolationScheme.eq.'average') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
          gf=GFactorCF(k)
c
          FiTf(k)=gf*FiT(i)+(1.-gf)*FiT(j)
c
        enddo
c
      elseif(InterpolationScheme.eq.'upwind') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gf=0.
          if(mdot(k).gt.0.) gf=1.          
c
          FiTf(k)=gf*FiT(i)+(1.-gf)*FiT(j)
c
        enddo
c
      elseif(InterpolationScheme.eq.'downwind') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gf=1.
          if(mdot(k).gt.0.) gf=0.          
c
          FiTf(k)=gf*FiT(i)+(1.-gf)*FiT(j)
c
        enddo
c
      elseif(InterpolationScheme.eq.'harmonic') then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
          gf=GFactorCF(k)
c
          FiTf(k)=(FiT(i)*FiT(j))/(gf*FiT(i)+(1.-gf)*FiT(j))
c
        enddo
c
      endif
c
      return
      end
