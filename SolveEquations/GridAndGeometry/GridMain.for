c
C#############################################################################################
      SUBROUTINE Grid
C#############################################################################################
c
      use User0, only:LReadSavedGrid,LStructured,MeshType
c*********************************************************************************************
      implicit none
c*********************************************************************************************
c
      if(MeshType.eq.'neutral') then
c
        if(LReadSavedGrid) then
c
          print*,'Reading saved grid and connectivity.....Please wait'
          call ReadSavedGrid
c
        else
c
          if(LStructured) then
            print*,'Reading structured grid file...Please wait'
            call ReadStructuredGrid
          else
            print*,'Reading neutral file...Please wait'
            call ReadNeutralGrid
          endif
c
          print*,'Establishing grid connectivity...Please wait'
          call GlobalConnectivity
c
          call CalculateNumberOfNodesInGroups
c
          call SaveGrid 
c
        endif
c
      elseif(MeshType.eq.'polymesh') then
c
        if(LReadSavedGrid) then
c
          print*,'Reading saved grid and connectivity.....Please wait'
          call ReadSavedGrid
c
        else
c
          print*,'Reading polymesh file and establishing '
          print*,'grid connectivity...Please wait'
          call ReadFoamMesh
c
          call CalculateFoamNumberOfNodesInGroups
c
          call SaveGrid 
c
        endif
c
      endif
c
      print*,'Calculating geometrical quantities...Please wait'
      call ProcessGeometry
      print*,'geometry processed'
      !pause
      !stop
c
      return
      end