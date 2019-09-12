c
C#############################################################################################
c
      SUBROUTINE CalculateInitialmdot
c
C#############################################################################################
c
      use User0, only: LsolveMomentum,LConvectScalar
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      mdot,Bmdot
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor,NBFaces
      use Geometry4, only: FaceAreax,FaceAreay,FaceAreaz,
     *                     BFaceAreax,BFaceAreay,BFaceAreaz,GFactorCF
      use BoundaryConditions1, only: BCType,inletTypeMomentum,
     *                               outletTypeMomentum,wallTypeMomentum
      use PhysicalProperties1, only: density,Bdensity,Densityf
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k
      double precision uVelocityf,vVelocityf,wVelocityf
c********************************************************************************************
c
c--- Initialize mdot along internal faces
c
      do k=1,NIFaces 
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        uVelocityf=GFactorCF(k)*uVelocity(i)+
     *               (1.-GFactorCF(k))*uVelocity(j)
        vVelocityf=GFactorCF(k)*vVelocity(i)+
     *               (1.-GFactorCF(k))*vVelocity(j)
        wVelocityf=GFactorCF(k)*wVelocity(i)+
     *               (1.-GFactorCF(k))*wVelocity(j)

        mdot(k)=Densityf(k)*(uVelocityf*FaceAreax(k)+
     *                   vVelocityf*FaceAreay(k)+
     *                          wVelocityf*FaceAreaz(k))
c
      enddo
c
c--- Initialize mdot along boundary faces
c
      if(LsolveMomentum) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            if(BCType(i,j).eq.'inlet') then
c
              if(inletTypeMomentum(i,j).ne.'specifiedmassflowrate') then
c
                Bmdot(i,j)=BDensity(i,j)*
     *            (BFaceAreax(i,j)*BuVelocity(i,j)+
     *                   BFaceAreay(i,j)*BvVelocity(i,j)+
     *                         BFaceAreaz(i,j)*BwVelocity(i,j))
c
              endif
c   
            elseif(BCType(i,j).eq.'outlet') then
c
              if(outletTypeMomentum(i,j).ne.'specifiedmassflowrate')then
c
                Bmdot(i,j)=BDensity(i,j)*
     *            (BFaceAreax(i,j)*BuVelocity(i,j)+
     *                   BFaceAreay(i,j)*BvVelocity(i,j)+
     *                         BFaceAreaz(i,j)*BwVelocity(i,j))
c
              endif
c
            elseif(BCType(i,j).eq.'wall') then
c
              Bmdot(i,j)=0.
c
            else
c
                Bmdot(i,j)=BDensity(i,j)*
     *            (BFaceAreax(i,j)*BuVelocity(i,j)+
     *                   BFaceAreay(i,j)*BvVelocity(i,j)+
     *                         BFaceAreaz(i,j)*BwVelocity(i,j))
c
            endif
c
          enddo
        enddo
c
      elseif(LConvectScalar) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
                Bmdot(i,j)=BDensity(i,j)*
     *            (BFaceAreax(i,j)*BuVelocity(i,j)+
     *                   BFaceAreay(i,j)*BvVelocity(i,j)+
     *                         BFaceAreaz(i,j)*BwVelocity(i,j))
c
          enddo
        enddo
c
      endif
c
	return
      end