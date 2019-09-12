c
c#############################################################################################
c
      SUBROUTINE FlowOut
c
c#############################################################################################
c
      use BoundaryConditions1, only: BoundaryType,outletTypeC
      use Geometry1, only: NumberOfBCSets    
      use Geometry3, only: NBFaces
      use Geometry4, only: BFaceAreax,BFaceAreay,BFaceAreaz
      use Variables1, only: Bmdot
      use FlowInOut1, only: mdotOut,mdot1
      use Variables1, only: BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: BDensity
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'outlet'.and.
     *                outletTypeC(i).ne.'fullydeveloped') return
c
      enddo
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'outlet'.and.
     *            outletTypeC(i).eq.'fullydeveloped') then
          mdot1(i)=0.
c
          do j=1,NBFaces(i)
c              
            Bmdot(i,j)=BDensity(i,j)*
     *                 (BuVelocity(i,j)*BFaceAreax(i,j)+
     *                         BvVelocity(i,j)*BFaceAreay(i,j)+
     *                              BwVelocity(i,j)*BFaceAreaz(i,j))
c
            if(Bmdot(i,j).lt.mdot1(i)) mdot1(i)=Bmdot(i,j)
c
          enddo
c
          mdot1(i)=dabs(mdot1(i))
c
        endif
c
      enddo
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'outlet'.and.
     *            outletTypeC(i).eq.'fullydeveloped') then
          mdotOut(i)=0.
c
          do j=1,NBFaces(i)
c
            mdotOut(i)=mdotOut(i)+mdot1(i)+Bmdot(i,j)
c
          enddo
c
        endif
c
      enddo
c
      return
      end
