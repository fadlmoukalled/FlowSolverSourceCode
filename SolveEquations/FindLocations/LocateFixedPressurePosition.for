c
C#############################################################################################
      SUBROUTINE LocateFixedPressureElement
C#############################################################################################
      use User0, only: xFixPressure,yFixPressure,
     *                 zFixPressure,FixAtLocation
c*********************************************************************************************
      implicit none
c*********************************************************************************************
c
        call LocatePointInVolume(xFixPressure,yFixPressure,
     *                              zFixPressure,FixAtLocation)
c
c         print*,xFixPressure,yFixPressure,zFixPressure,FixAtLocation
c         Pause
      return
      end
