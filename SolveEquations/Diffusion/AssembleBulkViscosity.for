c
C#############################################################################################
      SUBROUTINE AssembleBulkViscosityTerm(Variable,Gam,BGam)
C#############################################################################################
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
          FluxVflocal=-(-2./3.)*gamf*(uVelGradfx(k)+
     *                    vVelGradfy(k)+wVelGradfz(k))*FaceAreax(k)
c
          FluxCf(k)=FluxCf(k)+0.d0
          FluxFf(k)=FluxFf(k)+0.d0
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
              FluxVflocal=-(-2./3.)*gamf*(BuVelGradx(i,j)+
     *                 BvVelGrady(i,j)+BwVelGradz(i,j))*BFaceAreax(i,j)
c
              FluxCf(k)=FluxCf(k)+0.d0
              FluxFf(k)=FluxFf(k)+0.d0
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
          FluxVflocal=-(-2./3.)*gamf*(uVelGradfx(k)+
     *                     vVelGradfy(k)+wVelGradfz(k))*FaceAreay(k)
c
          FluxCf(k)=FluxCf(k)+0.d0
          FluxFf(k)=FluxFf(k)+0.d0
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
              FluxVflocal=-(-2./3.)*gamf*(BuVelGradx(i,j)+
     *              BvVelGrady(i,j)+ BwVelGradz(i,j))*BFaceAreay(i,j)
c
              FluxCf(k)=FluxCf(k)+0.d0
              FluxFf(k)=FluxFf(k)+0.d0
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
          FluxVflocal=-(-2./3.)*gamf*(uVelGradfx(k)+
     *                     vVelGradfy(k)+wVelGradfz(k))*FaceAreaz(k)
c
          FluxCf(k)=FluxCf(k)+0.d0
          FluxFf(k)=FluxFf(k)+0.d0
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
              FluxVflocal=-(-2./3.)*gamf*(BuVelGradx(i,j)+
     *              BvVelGrady(i,j)+ BwVelGradz(i,j))*BFaceAreaz(i,j)
c
              FluxCf(k)=FluxCf(k)+0.d0
              FluxFf(k)=FluxFf(k)+0.d0
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