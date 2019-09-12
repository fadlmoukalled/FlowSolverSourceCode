c
c#############################################################################################
c
      SUBROUTINE FlowIn
c
c#############################################################################################
c
      use BoundaryConditions1, only: BoundaryType,inletTypeC
      use Geometry1, only: NumberOfBCSets    
      use Geometry3, only: NBFaces
      use Variables1, only: Bmdot
      use FlowInOut1, only: mdotIn
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'inlet') then
c
          if(inletTypeC(i).eq.'specifiedstaticpressure'.or.
     *             inletTypeC(i).eq.'specifiedstagnationpressure') then
c              
            return
c
          endif
c
        endif
c
      enddo
c
      mdotIn=0.
c
      do i=1,NumberOfBCSets
c
        if(BoundaryType(i).eq.'inlet') then
c
          if(inletTypeC(i).eq.'specifiedvelocity'.or.
     *         inletTypeC(i).eq.'specifiedmassflowrate') then
c              
            do j=1,NBFaces(i)
c
              mdotIn=mdotIn+Bmdot(i,j)
c
            enddo
c
          endif
c
        endif
c
      enddo
c
      return
      end
