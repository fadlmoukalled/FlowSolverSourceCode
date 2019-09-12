c
c#############################################################################################
c
      SUBROUTINE BoundTemperature
c
c#############################################################################################
c
      use User0, only: tempmin,tempmax
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: Temperature,BTemperature
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        Temperature(i)=dmax1(tempmin,Temperature(i))
        Temperature(i)=dmin1(tempmax,Temperature(i))
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTemperature(i,j)=dmax1(tempmin,BTemperature(i,j))
          BTemperature(i,j)=dmin1(tempmax,BTemperature(i,j))
c         
        enddo
      enddo
c
      return
      end
