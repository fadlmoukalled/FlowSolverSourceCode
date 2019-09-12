c
C#############################################################################################
      SUBROUTINE LocatePointSources
C#############################################################################################
      use BoundaryConditions1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: iSource
      integer :: iElement
      double precision :: xpt,ypt,zpt
c*********************************************************************************************
c
      do iSource=1,NumberofPointSources
c
        xpt=xLocationOfPointSource(iSource)
        ypt=yLocationOfPointSource(iSource)
        zpt=zLocationOfPointSource(iSource)
c
        iElement=0
        call LocatePointInVolume(xpt,ypt,zpt,iElement)
        iElementPointSource(iSource)=iElement
c
      enddo
c
      return
      end
