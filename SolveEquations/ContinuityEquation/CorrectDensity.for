c
c#############################################################################################
c
      SUBROUTINE CorrectDensity
c
c#############################################################################################
c
      use User0, only: urfPressure
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: PressureC,BPressureC
      use PhysicalProperties1, only: Density,BDensity,drhodP,BdrhodP 
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements           
c
        Density(i)=Density(i)+urfPressure*drhodP(i)*PressureC(i)   
c      
      enddo      
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BDensity(i,j)=BDensity(i,j)+
     *             urfPressure*BdrhodP(i,j)*BPressureC(i,j)
c
        enddo
      enddo
c
      return
      end
