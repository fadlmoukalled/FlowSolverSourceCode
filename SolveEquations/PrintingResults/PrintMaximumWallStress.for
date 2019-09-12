c
c#############################################################################################
c
      SUBROUTINE PrintMaximumTauWall
c
C#############################################################################################
c
      use User0, only: LSolveMomentum,LTurbulentFlow,time,LUnsteady
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BFaceArea,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz,BFaceAreax,BFaceAreay,
     *                     BFaceAreaz,BgDiffp,
     *                     BFaceTxp,BFaceTyp,BFaceTzp
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: BeDiffCoefficient
      use BoundaryConditions1, only: BoundaryType,wallTypeM
      use WallStress1
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,k
      double precision :: gDiff,dNorm,nx,ny,nz,
     *                    ShearStressX,ShearStressY,ShearStressZ,
     *                    ShearStress,areaX,areaY,areaZ,area,
     *                    uWall1,vWall1,wWall1,WallVelocity
c********************************************************************************************
      interface
c********************************************************************************************
        FUNCTION TangentialVelocity(i1)
c********************************************************************************************
          integer :: i1
          double precision :: TangentialVelocity
c********************************************************************************************
        end FUNCTION TangentialVelocity
c********************************************************************************************
        FUNCTION TangentialVelocityLaminar(i1,i2,i3)
c********************************************************************************************
          integer :: i1,i2,i3
          double precision :: TangentialVelocityLaminar
c********************************************************************************************
        end FUNCTION TangentialVelocityLaminar
c********************************************************************************************
      end interface
c********************************************************************************************
c
      if(.not.LSolveMomentum) return
c
      Variable='velx'
      call CalculateEffectiveDiffusionCoefficient(Variable)
c
c---- Print maximum shear stress along walls
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'wall'.and.
     *                  wallTypeM(i).eq.'noslip') then
c
          MaximumShearStress(i)=0.
          ShearStressX=0.
          ShearStressY=0.
          ShearStressz=0.
          ShearStress=0.
c
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            dNorm=BDistanceCFx(i,j)*BFaceAreanx(i,j)+
     *                BDistanceCFy(i,j)*BFaceAreany(i,j)+
     *                   BDistanceCFz(i,j)*BFaceAreanz(i,j)
c
            area=BFaceArea(i,j)
            areaX=BFaceAreax(i,j)
            areaY=BFaceAreay(i,j)
            areaZ=BFaceAreaz(i,j)
c
            nx=BFaceAreanx(i,j)
            ny=BFaceAreany(i,j)
            nz=BFaceAreanz(i,j)
c
            gDiff=BeDiffCoefficient(i,j)/dNorm
            ShearStressX=gDiff*(
     *               (uVelocity(k)-BuVelocity(i,j))*(1-nx**2)-
     *                        (vVelocity(k)-BvVelocity(i,j))*nx*ny-
     *                           (wVelocity(k)-BwVelocity(i,j))*nx*nz)
            ShearStressY=gDiff*(
     *               (vVelocity(k)-BvVelocity(i,j))*(1-ny**2)-
     *                        (uVelocity(k)-BuVelocity(i,j))*nx*ny-
     *                        (wVelocity(k)-BwVelocity(i,j))*ny*nz)
            ShearStressZ=gDiff*(
     *               (wVelocity(k)-BwVelocity(i,j))*(1-nz**2)-
     *                        (uVelocity(k)-BuVelocity(i,j))*nx*nz-
     *                        (vVelocity(k)-BvVelocity(i,j))*ny*nz)
c
            if(LTurbulentFlow) then
              WallVelocity=TangentialVelocity(k)
            else
              WallVelocity=TangentialVelocityLaminar(k,i,j)
            endif
            ShearStress=gDiff*WallVelocity
c
            if(ShearStress.gt.MaximumShearStress(i)) then           
c
              MaximumShearStress(i)=ShearStress
              xLocationMaxWallShearStress(i)=BFaceCentroidx(i,j)
              yLocationMaxWallShearStress(i)=BFaceCentroidy(i,j)
              zLocationMaxWallShearStress(i)=BFaceCentroidz(i,j)
c
            endif
c
          enddo
c
        endif
c        
      enddo
c
      if(.not.LUnsteady) return
      write(18,*) time
      do i=1,NumberOfBCSets
        write(18,*) xLocationMaxWallShearStress(i),
     *              yLocationMaxWallShearStress(i),
     *              zLocationMaxWallShearStress(i),
     *              MaximumShearStress(i)
      enddo
c      
      return
      end
!c
!c#############################################################################################
!c
!      SUBROUTINE calculateDeviatoricStressForPrinting
!c
!C#############################################################################################
!c
!      use User0, only: Lcompressible,LSolveTurbulenceKineticEnergy
!      use Geometry1, only: NumberOfElements,NumberOfBCSets
!      use Geometry3, only: NBFaces,NBFaceOwner
!      use Geometry4, only: BDistanceCFx,BDistanceCFy,BDistanceCFz,
!     *                     BFaceAreanx,BFaceAreany,BFaceAreanz,
!     *                     BFaceArea,BFaceCentroidx,BFaceCentroidy,
!     *                     BFaceCentroidz,BFaceAreax,BFaceAreay,
!     *                     BFaceAreaz,BgDiffp,
!     *                     BFaceTxp,BFaceTyp,BFaceTzp
!      use Variables1, only: uVelocity,vVelocity,wVelocity,
!     *                      BuVelocity,BvVelocity,BwVelocity,
!     *                      TurbulentKE, BTurbulentKE,
!     *                      uVelGradx,vVelGrady,wVelGradz,
!     *                      BuVelGradx,BvVelGrady,BwVelGradz
!      use PhysicalProperties1, only: Density,BDensity,
!     *                               eDiffCoefficient,BeDiffCoefficient
!      use Turbulence1, only: S11,S12,S13,S22,S23,S33,
!     *                       BS11,BS12,BS13,BS22,BS23,BS33,
!     *                       Tau11,Tau12,Tau13,Tau22,Tau23,Tau33,
!     *                       BTau11,BTau12,BTau13,BTau22,BTau23,BTau33 
!      use BoundaryConditions1, only: BoundaryType,wallTypeM
!      use Constants1, only: tiny,twothird
!c********************************************************************************************
!      implicit none
!c********************************************************************************************
!      character*10 Variable
!      integer :: i,j,k
!      double precision :: term1,tke1,gDiff,dNorm,nx,ny,nz,
!     *                    ShearStressX,ShearStressY,ShearStressZ,
!     *                    ShearStress,areaX,areaY,areaZ,area,
!     *                    uWall1,vWall1,wWall1,WallVelocity
!c********************************************************************************************
!      interface
!c********************************************************************************************
!        FUNCTION TangentialVelocity(i1)
!c********************************************************************************************
!          integer :: i1
!          double precision :: TangentialVelocity
!c********************************************************************************************
!        end FUNCTION TangentialVelocity
!c********************************************************************************************
!      end interface
!c********************************************************************************************
!c
!      Variable='velx'
!      call CalculateEffectiveDiffusionCoefficient(Variable)
!c
!      do i=1,NumberOfElements    
!c
!        Tau11(i)=eDiffCoefficient(i)*2.*S11(i)
!        Tau12(i)=eDiffCoefficient(i)*2.*S12(i)
!        Tau13(i)=eDiffCoefficient(i)*2.*S13(i)
!        Tau22(i)=eDiffCoefficient(i)*2.*S22(i)
!        Tau23(i)=eDiffCoefficient(i)*2.*S23(i)
!        Tau33(i)=eDiffCoefficient(i)*2.*S33(i)
!c          
!      enddo
!c
!      do i=1,NumberOfBCSets
!        do j=1,NBFaces(i)    
!c
!          BTau11(i,j)=BeDiffCoefficient(i,j)*2.*BS11(i,j)
!          BTau12(i,j)=BeDiffCoefficient(i,j)*2.*BS12(i,j)
!          BTau13(i,j)=BeDiffCoefficient(i,j)*2.*BS13(i,j)
!          BTau22(i,j)=BeDiffCoefficient(i,j)*2.*BS22(i,j)
!          BTau23(i,j)=BeDiffCoefficient(i,j)*2.*BS23(i,j)
!          BTau33(i,j)=BeDiffCoefficient(i,j)*2.*BS33(i,j)
!c          
!        enddo
!      enddo
!c
c
!      if(LSolveTurbulenceKineticEnergy) then
!c
!        if(TurbulenceModel.eq.'kklomega') then
!c
!          do i=1,NumberOfElements    
!c
!            tke1=dmax1(TurbulentKE(i)+TurbulentKL(i),0.)
!c
!            Tau11(i)=Tau11(i)-twothird*Density(i)*tke1
!            Tau22(i)=Tau22(i)-twothird*Density(i)*tke1
!            Tau33(i)=Tau33(i)-twothird*Density(i)*tke1
!c          
!          enddo
!c
!        else
!c
!          do i=1,NumberOfElements    
!c
!            tke1=dmax1(TurbulentKE(i),0.)
!c
!            Tau11(i)=Tau11(i)-twothird*Density(i)*tke1
!            Tau22(i)=Tau22(i)-twothird*Density(i)*tke1
!            Tau33(i)=Tau33(i)-twothird*Density(i)*tke1
!c          
!          enddo
!c
!        endif
!c
!      endif
!c
!      if(Lcompressible) then
!c
!        do i=1,NumberOfElements    
!c
!          term1=uVelGradx(i)+vVelGrady(i)+wVelGradz(i)
!c
!          Tau11(i)=Tau11(i)-twothird*TurbulentViscosity(i)*term1
!          Tau22(i)=Tau22(i)-twothird*TurbulentViscosity(i)*term1
!          Tau33(i)=Tau33(i)-twothird*TurbulentViscosity(i)*term1
!c          
!        enddo
!c
!      endif
!c
!
!
!
!
!      return
!      end
!c
!
!
!
!
!
!