c
C#############################################################################################
      SUBROUTINE AssembleExplicitDiffusionTerm(Variable,Gam,BGam)
C#############################################################################################
c
      use Variables1
      use Variables2
      use Geometry1
      use Geometry3
      use Geometry4
      use PhysicalProperties1
      use BoundaryConditions1, only: BCType
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k
      character*10 Variable
      double precision Gamf
      double precision FluxCflocal,FluxFflocal,FluxVflocal
      double precision, dimension(:) :: Gam
      double precision, dimension(:,:) :: BGam
c*********************************************************************************************
c
c--- Internal faces
c
      if(Variable.eq.'velx') then
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gamf=GFactorCF(k)*Gam(i)+(1.-GFactorCF(k))*Gam(j)
c
          FluxCflocal=0.
          FluxFflocal=0.
          FluxVflocal=-gamf*(uVelGradfx(k)*FaceAreax(k)+
     *         vVelGradfx(k)*FaceAreay(k)+wVelGradfx(k)*FaceAreaz(k))
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
          FluxVf(k)=FluxVf(k)+FluxVflocal
          FluxTf(k)=FluxTf(k)+FluxVflocal
c
        enddo
c
c--- Boundary faces
c
        k=NIFaces
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            k=k+1
c
c            if(BCType(i,j).ne.'wall'.and.BCType(i,j).ne.'symmetry') then
c
              gamf=Bgam(i,j)

              FluxCflocal=0.
              FluxFflocal=0.
              FluxVflocal=-gamf*(BuVelGradx(i,j)*BFaceAreax(i,j)+
     *                             BvVelGradx(i,j)*BFaceAreay(i,j)+
     *                              BwVelGradx(i,j)*BFaceAreaz(i,j))
c
              FluxCf(k)=FluxCf(k)+FluxCflocal
              FluxFf(k)=FluxFf(k)+FluxFflocal
              FluxVf(k)=FluxVf(k)+FluxVflocal
              FluxTf(k)=FluxTf(k)+FluxVflocal
c
c            endif
c
          enddo
        enddo
c
      elseif(Variable.eq.'vely') then
c
c--- Internal faces
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gamf=GFactorCF(k)*Gam(i)+(1.-GFactorCF(k))*Gam(j)
c
          FluxCflocal=0.
          FluxFflocal=0.
          FluxVflocal=-gamf*(uVelGradfy(k)*FaceAreax(k)+
     *            vVelGradfy(k)*FaceAreay(k)+wVelGradfy(k)*FaceAreaz(k))
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
          FluxVf(k)=FluxVf(k)+FluxVflocal
          FluxTf(k)=FluxTf(k)+FluxVflocal
c
        enddo
c
c--- Boundary faces
c
        k=NIFaces
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            k=k+1
c
c            if(BCType(i,j).ne.'wall'.and.BCType(i,j).ne.'symmetry') then
c
              gamf=Bgam(i,j)

              FluxCflocal=0.
              FluxFflocal=0.
              FluxVflocal=-gamf*(BuVelGrady(i,j)*BFaceAreax(i,j)+
     *                             BvVelGrady(i,j)*BFaceAreay(i,j)+
     *                                BwVelGrady(i,j)*BFaceAreaz(i,j))
c
              FluxCf(k)=FluxCf(k)+FluxCflocal
              FluxFf(k)=FluxFf(k)+FluxFflocal
              FluxVf(k)=FluxVf(k)+FluxVflocal
              FluxTf(k)=FluxTf(k)+FluxVflocal
c
c            endif
c
          enddo
        enddo
c
      elseif(Variable.eq.'velz') then
c
c--- Internal faces
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gamf=GFactorCF(k)*Gam(i)+(1.-GFactorCF(k))*Gam(j)
c
          FluxCflocal=0.
          FluxFflocal=0.
          FluxVflocal=-gamf*(uVelGradfz(k)*FaceAreax(k)+
     *            vVelGradfz(k)*FaceAreay(k)+wVelGradfz(k)*FaceAreaz(k))
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
          FluxVf(k)=FluxVf(k)+FluxVflocal
          FluxTf(k)=FluxTf(k)+FluxVflocal
c
        enddo
c
c--- Boundary faces
c
        k=NIFaces
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            k=k+1
c
c            if(BCType(i,j).ne.'wall'.and.BCType(i,j).ne.'symmetry') then
c
              gamf=Bgam(i,j)

              FluxCflocal=0.
              FluxFflocal=0.
              FluxVflocal=-gamf*(BuVelGradz(i,j)*BFaceAreax(i,j)+
     *                             BvVelGradz(i,j)*BFaceAreay(i,j)+
     *                                BwVelGradz(i,j)*BFaceAreaz(i,j))
c
              FluxCf(k)=FluxCf(k)+FluxCflocal
              FluxFf(k)=FluxFf(k)+FluxFflocal
              FluxVf(k)=FluxVf(k)+FluxVflocal
              FluxTf(k)=FluxTf(k)+FluxVflocal
c
c            endif
c
          enddo
        enddo
c      
      endif
c
      return
      end