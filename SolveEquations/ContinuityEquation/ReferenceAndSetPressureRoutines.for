c
c#############################################################################################
c
      SUBROUTINE SetReferencePressure
c
c#############################################################################################
c
      use ReferencePressure1, only: LSetReferencePressure
      use User0, only: Lcompressible,LSolveMomentum
      use BoundaryConditions1, only: BCType,inletTypeMomentum,
     *                               outletTypeMomentum
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      LSetReferencePressure=.false.
      if(.not.Lcompressible) then
        if(LSolveMomentum) then
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              if(BCType(i,j).eq.'inlet') then
                if(inletTypeMomentum(i,j).eq.
     *                'specifiedstaticpressure'.or.
     *                        inletTypeMomentum(i,j).eq.
     *                        'specifiedstagnationpressure') then  
c
                   LSetReferencePressure=.true.
                   exit                    
c
                endif
              endif 
c
              if(BCType(i,j).eq.'outlet') then
                if(outletTypeMomentum(i,j).eq.
     *                'specifiedstaticpressure'.or.
     *                        outletTypeMomentum(i,j).eq.
     *                'specifiedaveragestaticpressure'.or.
     *                        outletTypeMomentum(i,j).eq.
     *                'specifiedresistance'.or.
     *                        outletTypeMomentum(i,j).eq.
     *                        'specifiedstagnationpressure') then  
c
                   LSetReferencePressure=.true.
                   exit                    
c
                endif
              endif 
c
            enddo
          enddo
        endif
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE FixPressure
c
c#############################################################################################
c
      use User0, only: FixAtLocation
      use Variables2, only: anb,bc
      use Geometry3, only:NumberofElementNeighbors 
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      i=FixAtLocation
c
      bc(i)=0
c        
      do j = 1,NumberofElementNeighbors(i)
c
        anb(i,j)=0.
c          
      enddo
c
      return
      end