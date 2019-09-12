c
c#############################################################################################
c
      SUBROUTINE AdjustFlow
c
c#############################################################################################
c
      use BoundaryConditions1, only: BoundaryType,outletTypeC
      use Geometry1, only: NumberOfBCSets    
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: BDistanceCFx,BDistanceCFy,BDistanceCFz
      use Variables1, only: Bmdot
      use FlowInOut1, only: mdotIn,mdotOut,mdot1,MassFlowFraction
      use Variables1, only: BuVelocity,BvVelocity,BwVelocity,
     *                      uVelocity,vVelocity,wVelocity,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz
      use Constants1, only: tiny     
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      double precision :: ratio,totalflowfraction
c********************************************************************************************
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'outlet'.and.
     *                outletTypeC(i).ne.'fullydeveloped') return
c
      enddo
c
      totalflowfraction=0.
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'outlet'.and.
     *             outletTypeC(i).eq.'fullydeveloped') then
c
          totalflowfraction=totalflowfraction+MassFlowFraction(i)
c
        endif
c
      enddo
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'outlet'.and.
     *             outletTypeC(i).eq.'fullydeveloped') then
c
          ratio=(MassFlowFraction(i)/(totalflowfraction+tiny))*
     *                               dabs(mdotIn/(mdotOut(i)+tiny))
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            Bmdot(i,j)=ratio*(Bmdot(i,j)+mdot1(i))
c            BuVelocity(i,j)=uVelocity(k)
c            BvVelocity(i,j)=vVelocity(k)
c            BwVelocity(i,j)=wVelocity(k)
c
            BuVelocity(i,j)=uVelocity(k)+
     *           (BuVelGradx(i,j)*BDistanceCFx(i,j)+
     *               BuVelGrady(i,j)*BDistanceCFy(i,j)+
     *                    BuVelGradz(i,j)*BDistanceCFz(i,j))
            BvVelocity(i,j)=vVelocity(k)+
     *           (BvVelGradx(i,j)*BDistanceCFx(i,j)+
     *               BvVelGrady(i,j)*BDistanceCFy(i,j)+
     *                    BvVelGradz(i,j)*BDistanceCFz(i,j))
            BwVelocity(i,j)=wVelocity(k)+
     *           (BwVelGradx(i,j)*BDistanceCFx(i,j)+
     *               BwVelGrady(i,j)*BDistanceCFy(i,j)+
     *                    BwVelGradz(i,j)*BDistanceCFz(i,j))
c
          enddo
c
        endif
c
      enddo
c
      return
      end
