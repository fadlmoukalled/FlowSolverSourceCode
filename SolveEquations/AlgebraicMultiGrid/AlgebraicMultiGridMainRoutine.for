c
C#############################################################################################
c
      SUBROUTINE AlgebraicMultigrid
     *           (NF,rrF,itmax,solver,IResiduals,FResiduals,Variable1)
c
C#############################################################################################
c
      use User0, only: reAgglomorate,MultiGridCycleType, 
     *                 MaxMultiGridCycles,LTestMultiGrid,nprintMG,
     *                 nIterStartApplyingMG,MGVariable,MGType
      use MultiGrid2
      use Variables2, only:ac,anb,bc,dphi
      use Geometry1,  only:NumberOfElements,NumbOfElementNodes,
     *                     MaximumNumberofElementNodes,
     *                     ListOfElementNodes,NodeFlag,NumberOfNodes
      use Geometry3,  only:ElementNeighbor,NumberofElementNeighbors
      use MultiGrid1, only:NumberOfElementsMG,ElementNeighborMG,
     *                     NumberofElementNeighborsMG,Residuals,
     *                     acMG,bcMG,anbMG,ijbeginE,ijbeginN,dphiMG,
     *                     ijbeginNOld,NElementFacesMG,NBChildrenMaxMG
      use Geometry4, only: Volume
      use Residuals1
      use GeometricMG1
      use InteriorBoundaryNodes, only:NumberOfInteriorNodes,InteriorNode
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*35 Variable1
      character*6 solver
      integer i,j,ij,ncycle
      integer NF,itmax,iLevel,indexMGPrint
      double precision rrF,IResiduals,FResiduals
      integer, save :: indexGeometric=1
c********************************************************************************************
c
c--- set fine grid coefficients as the multigrid coefficients Level 1
c
c********************************************************************************************
c
      if(MGType.eq.'algebraic') then

        If(nIterMG.eq.nIterStartApplyingMG)  then
c
          call ResetMultigridAllocation
c
        elseif(mod(nIterMG,reAgglomorate).eq.0)  then
c
          if(Variable1.eq.MGVariable) call ResetMultigridAllocation
c
        endif
c       
      endif
c
      if(MGType.eq.'algebraic'.and.nIterMG.eq.nIterStartApplyingMG) then
c
        iLevel=1
        ijbeginE(iLevel)=0
        ijbeginN(iLevel)=0
c
        do i=1,NumberOfElements
          do j=1,NumberofElementNeighbors(i)
c
            ij=i+(j-1)*NumberOfElements
            ElementNeighborMG(ij)=ElementNeighbor(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          NumberofElementNeighborsMG(i)=NumberofElementNeighbors(i)
c
        enddo
c
        NumberOfElementsMG(iLevel)=NumberOfElements
c
        do i=1,NumberOfElements
c  
          acMG(i)=ac(i)
          bcMG(i)=bc(i)
          Residuals(i)=bc(i)
c
          do j=1,NumberofElementNeighbors(i)
            ij=i+(j-1)*NumberOfElements
c      
            anbMG(ij)=anb(i,j)
c
          enddo
c
        enddo
c
        call Agglomorate
c
        if(LTestMultiGrid) then
          indexMGPrint=1
          call PrintMultigrid(indexMGPrint)
          indexMGPrint=2
          call PrintMultigrid(indexMGPrint)
        endif
c
      elseif(MGType.eq.'algebraic'.and.mod(nIterMG,reAgglomorate)
     *                          .eq.0.and.Variable1.eq.MGVariable) then
c
        iLevel=1
        ijbeginE(iLevel)=0
        ijbeginN(iLevel)=0
c
        do i=1,NumberOfElements
          do j=1,NumberofElementNeighbors(i)
c
            ij=i+(j-1)*NumberOfElements
            ElementNeighborMG(ij)=ElementNeighbor(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          NumberofElementNeighborsMG(i)=NumberofElementNeighbors(i)
c
        enddo
c
        NumberOfElementsMG(iLevel)=NumberOfElements
c
        do i=1,NumberOfElements
c  
          acMG(i)=ac(i)
          bcMG(i)=bc(i)
          Residuals(i)=bc(i)
c
          do j=1,NumberofElementNeighbors(i)
            ij=i+(j-1)*NumberOfElements
c      
            anbMG(ij)=anb(i,j)
c
          enddo
c
        enddo
c
        call Agglomorate
c
        if(LTestMultiGrid) then
          indexMGPrint=1
          call PrintMultigrid(indexMGPrint)
          indexMGPrint=2
          call PrintMultigrid(indexMGPrint)
        endif
c
      elseif(MGType.eq.'geometricelement'.and.indexGeometric.eq.1) then
c
        call ResetMultigridAllocation
c
        iLevel=1
        ijbeginE(iLevel)=0
        ijbeginN(iLevel)=0
c
        do i=1,NumberOfElements
          do j=1,NumberofElementNeighbors(i)
c
            ij=i+(j-1)*NumberOfElements
            ElementNeighborMG(ij)=ElementNeighbor(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          NumberofElementNeighborsMG(i)=NumberofElementNeighbors(i)
c
        enddo
c
        NumberOfElementsMG(iLevel)=NumberOfElements
c
        do i=1,NumberOfElements
c  
          acMG(i)=ac(i)
          bcMG(i)=bc(i)
          Residuals(i)=bc(i)
c
          do j=1,NumberofElementNeighbors(i)
            ij=i+(j-1)*NumberOfElements
c      
            anbMG(ij)=anb(i,j)
c
          enddo
c
        enddo
c
        call Agglomorate
c
        if(LTestMultiGrid) then
          indexMGPrint=1
          call PrintMultigrid(indexMGPrint)
          indexMGPrint=2
          call PrintMultigrid(indexMGPrint)
        endif
c
        indexGeometric=0
c
      elseif(MGType.eq.'geometricnode'.and.indexGeometric.eq.1) then
c
        call ResetMultigridAllocation
c
        iLevel=1
        ijbeginE(iLevel)=0
        ijbeginN(iLevel)=0

        allocate(ijbeginNode(1))
        allocate(ijbeginENode(1))
        allocate(ijbeginNElement(1))
c
        ijbeginNode(iLevel)=0
        ijbeginENode(iLevel)=0
        ijbeginNElement(iLevel)=0
c
        do i=1,NumberOfElements
          do j=1,NumberofElementNeighbors(i)
c
            ij=i+(j-1)*NumberOfElements
            ElementNeighborMG(ij)=ElementNeighbor(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          NumberofElementNeighborsMG(i)=NumberofElementNeighbors(i)
c
        enddo
c
        NumberOfElementsMG(iLevel)=NumberOfElements
c
        allocate(NodesMG(NumberOfNodes))
        allocate(NumberOfNodesMG(1))
        allocate(FlagNodesMG(NumberOfNodes))
c
        NumberOfNodesMG(iLevel)=NumberOfNodes
        FlagNodesMG=NodeFlag
c
        do i=1,NumberOfNodes
c
          NodesMG(i)=i
c
        enddo
c
        allocate(NumbOfElementNodesMG(NumberOfElements))
        do i=1,NumberOfElements
c
          NumbOfElementNodesMG(i)=NumbOfElementNodes(i)
c
        enddo
c
        allocate(ListOfElementNodesMG
     *         (NumberOfElements*MaximumNumberofElementNodes))
        do i=1,NumberOfElements
          do j=1,NumbOfElementNodes(i)
c
            ij=ijbeginENode(iLevel)+i+(j-1)*NumberOfElements
            ListOfElementNodesMG(ij)=ListOfElementNodes(i,j)
c
          enddo
        enddo
c
        allocate(VolumeMG(NumberOfElements))
        do i=1,NumberOfElements
c
            VolumeMG(i)=Volume(i)
c
        enddo
c
        allocate(NumberOfElementsConnectedToNodeMG(1))
        allocate(ListOfElementsConnectedToNodeMG(1))
        NumberOfElementsConnectedToNodeMG=0
        ListOfElementsConnectedToNodeMG=0
c
        do i=1,NumberOfElements
c  
          acMG(i)=ac(i)
          bcMG(i)=bc(i)
          Residuals(i)=bc(i)
c
          do j=1,NumberofElementNeighbors(i)
            ij=i+(j-1)*NumberOfElements
c      
            anbMG(ij)=anb(i,j)
c
          enddo
c
        enddo
c
        call Agglomorate
c
        if(LTestMultiGrid) then
          indexMGPrint=1
          call PrintMultigrid(indexMGPrint)
          indexMGPrint=2
          call PrintMultigrid(indexMGPrint)
        endif
c
        indexGeometric=0
c
      else
c
        if(MGType.eq.'algebraic'.or.MGType.eq.'geometricelement') then
c
          do i=1,NumberOfElements
c  
            acMG(i)=ac(i)
            bcMG(i)=bc(i)
            Residuals(i)=bc(i)
c
            do j=1,NumberofElementNeighbors(i)
              ij=i+(j-1)*NumberOfElements
c      
              anbMG(ij)=anb(i,j)
c
            enddo
c
          enddo
c
          do iLevel=2,iLevelMax
c              
            call AssembleAgglomeratedLHS(iLevel)
c
          enddo
c
        elseif(MGType.eq.'geometricnode') then
c
          iLevel=1
c
          NumberOfElementsMG(iLevel)=NumberOfElements
c
          do i=1,NumberOfElements
c  
            acMG(i)=ac(i)
            bcMG(i)=bc(i)
            Residuals(i)=bc(i)
c
            do j=1,NumberofElementNeighbors(i)
              ij=i+(j-1)*NumberOfElements
c      
              anbMG(ij)=anb(i,j)
c
            enddo
c
          enddo
c
          do iLevel=2,iLevelMax
c              
            call AssembleAgglomeratedLHSGeometrical(iLevel)
c
          enddo
c
        endif
c
        if(LTestMultiGrid) then
          if(nIterMG.eq.nIterStartApplyingMG.or.
     *                      mod(nIterMG,nprintMG).eq.0)  then
c
            indexMGPrint=2
            call PrintMultigrid(indexMGPrint)
c
          endif
        endif
c
      endif
c
c--- Apply selected MutliGrid Cycle
c
      dphiMG=0.
      Residuals=0.
      ncycle=1
      if(MultiGridCycleType.eq.'vcycle') then
c
        do while(ncycle.le.MaxMultiGridCycles.and.
     *                       FResiduals.gt.rrF*IResiduals)
          call vMultiGridCycle(NF,rrF,itmax,solver)
          call CheckConvergenceMG(NF,1,FResiduals)
          ncycle=ncycle+1
        enddo
c
      elseif(MultiGridCycleType.eq.'fcycle') then
c
        do while(ncycle.le.MaxMultiGridCycles.and.
     *                       FResiduals.gt.rrF*IResiduals)
          call fMultiGridCycle(NF,rrF,itmax,solver)
          call CheckConvergenceMG(NF,1,FResiduals)
          ncycle=ncycle+1
        enddo
c
      elseif(MultiGridCycleType.eq.'wcycle') then
c
        do while(ncycle.le.MaxMultiGridCycles.and.
     *                       FResiduals.gt.rrF*IResiduals)
          call wMultiGridCycle(NF,rrF,itmax,solver)
          call CheckConvergenceMG(NF,1,FResiduals)
          ncycle=ncycle+1
        enddo
c
      endif
c
      indexMGPrint=3
      if (LTestMultiGrid) then
        if(nIterMG.eq.1.or.mod(nIterMG,nprintMG).eq.0)  then
          call PrintMultiGrid(indexMGPrint)
        endif
      endif
c
c--- Save the corrections to the original dphi array
c
      do i=1,NumberOfElements
c  
        dphi(i)=dphiMG(i)
c
      enddo
c
      return
      end
      
      
      
      
      
      
      
c
C#############################################################################################
c
      SUBROUTINE AlgebraicMultigridWallDistance
     *           (NF,rrF,itmax,solver,IResiduals,FResiduals)
c
C#############################################################################################
c
      use User0, only: reAgglomorate,MultiGridCycleType, 
     *                 MaxMultiGridCycles,LTestMultiGrid,nprintMG,
     *                 nIterStartApplyingMG,MGVariable,MGType
      use MultiGrid2
      use Variables2, only:ac,anb,bc,dphi
      use Geometry1,  only:NumberOfElements,NumbOfElementNodes,
     *                     MaximumNumberofElementNodes,
     *                     ListOfElementNodes,NodeFlag,NumberOfNodes
      use Geometry3,  only:ElementNeighbor,NumberofElementNeighbors
      use MultiGrid1, only:NumberOfElementsMG,ElementNeighborMG,
     *                     NumberofElementNeighborsMG,Residuals,
     *                     acMG,bcMG,anbMG,ijbeginE,ijbeginN,dphiMG,
     *                     ijbeginNOld,NElementFacesMG,NBChildrenMaxMG
      use Geometry4, only: Volume
      use Residuals1
      use GeometricMG1
      use InteriorBoundaryNodes, only:NumberOfInteriorNodes,InteriorNode
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*35 Variable1
      character*6 solver
      integer i,j,ij,ncycle
      integer NF,itmax,iLevel,indexMGPrint
      double precision rrF,IResiduals,FResiduals
      integer, save :: indexGeometric=1
c********************************************************************************************
c
c--- set fine grid coefficients as the multigrid coefficients Level 1
c
c********************************************************************************************
c
      if(MGType.eq.'algebraic'.and.indexGeometric.eq.1) then
c
        call ResetMultigridAllocation
c
        iLevel=1
        ijbeginE(iLevel)=0
        ijbeginN(iLevel)=0
c
        do i=1,NumberOfElements
          do j=1,NumberofElementNeighbors(i)
c
            ij=i+(j-1)*NumberOfElements
            ElementNeighborMG(ij)=ElementNeighbor(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          NumberofElementNeighborsMG(i)=NumberofElementNeighbors(i)
c
        enddo
c
        NumberOfElementsMG(iLevel)=NumberOfElements
c
        do i=1,NumberOfElements
c  
          acMG(i)=ac(i)
          bcMG(i)=bc(i)
          Residuals(i)=bc(i)
c
          do j=1,NumberofElementNeighbors(i)
            ij=i+(j-1)*NumberOfElements
c      
            anbMG(ij)=anb(i,j)
c
          enddo
c
        enddo
c
        call Agglomorate
c
        if(LTestMultiGrid) then
          indexMGPrint=1
          call PrintMultigrid(indexMGPrint)
          indexMGPrint=2
          call PrintMultigrid(indexMGPrint)
        endif
c
        indexGeometric=0
c
      elseif(MGType.eq.'geometricelement'.and.indexGeometric.eq.1) then
c
        call ResetMultigridAllocation
c
        iLevel=1
        ijbeginE(iLevel)=0
        ijbeginN(iLevel)=0
c
        do i=1,NumberOfElements
          do j=1,NumberofElementNeighbors(i)
c
            ij=i+(j-1)*NumberOfElements
            ElementNeighborMG(ij)=ElementNeighbor(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          NumberofElementNeighborsMG(i)=NumberofElementNeighbors(i)
c
        enddo
c
        NumberOfElementsMG(iLevel)=NumberOfElements
c
        do i=1,NumberOfElements
c  
          acMG(i)=ac(i)
          bcMG(i)=bc(i)
          Residuals(i)=bc(i)
c
          do j=1,NumberofElementNeighbors(i)
            ij=i+(j-1)*NumberOfElements
c      
            anbMG(ij)=anb(i,j)
c
          enddo
c
        enddo
c
        call Agglomorate
c
        if(LTestMultiGrid) then
          indexMGPrint=1
          call PrintMultigrid(indexMGPrint)
          indexMGPrint=2
          call PrintMultigrid(indexMGPrint)
        endif
c
        indexGeometric=0
c
      elseif(MGType.eq.'geometricnode'.and.indexGeometric.eq.1) then
c
        call ResetMultigridAllocation
c
        iLevel=1
        ijbeginE(iLevel)=0
        ijbeginN(iLevel)=0

        allocate(ijbeginNode(1))
        allocate(ijbeginENode(1))
        allocate(ijbeginNElement(1))
c
        ijbeginNode(iLevel)=0
        ijbeginENode(iLevel)=0
        ijbeginNElement(iLevel)=0
c
        do i=1,NumberOfElements
          do j=1,NumberofElementNeighbors(i)
c
            ij=i+(j-1)*NumberOfElements
            ElementNeighborMG(ij)=ElementNeighbor(i,j)
c
          enddo
        enddo
c
        do i=1,NumberOfElements
c
          NumberofElementNeighborsMG(i)=NumberofElementNeighbors(i)
c
        enddo
c
        NumberOfElementsMG(iLevel)=NumberOfElements
c
        allocate(NodesMG(NumberOfNodes))
        allocate(NumberOfNodesMG(1))
        allocate(FlagNodesMG(NumberOfNodes))
c
        NumberOfNodesMG(iLevel)=NumberOfNodes
        FlagNodesMG=NodeFlag
c
        do i=1,NumberOfNodes
c
          NodesMG(i)=i
c
        enddo
c
        allocate(NumbOfElementNodesMG(NumberOfElements))
        do i=1,NumberOfElements
c
          NumbOfElementNodesMG(i)=NumbOfElementNodes(i)
c
        enddo
c
        allocate(ListOfElementNodesMG
     *         (NumberOfElements*MaximumNumberofElementNodes))
        do i=1,NumberOfElements
          do j=1,NumbOfElementNodes(i)
c
            ij=ijbeginENode(iLevel)+i+(j-1)*NumberOfElements
            ListOfElementNodesMG(ij)=ListOfElementNodes(i,j)
c
          enddo
        enddo
c
        allocate(VolumeMG(NumberOfElements))
        do i=1,NumberOfElements
c
            VolumeMG(i)=Volume(i)
c
        enddo
c
        allocate(NumberOfElementsConnectedToNodeMG(1))
        allocate(ListOfElementsConnectedToNodeMG(1))
c
        do i=1,NumberOfElements
c  
          acMG(i)=ac(i)
          bcMG(i)=bc(i)
          Residuals(i)=bc(i)
c
          do j=1,NumberofElementNeighbors(i)
            ij=i+(j-1)*NumberOfElements
c      
            anbMG(ij)=anb(i,j)
c
          enddo
c
        enddo
c
        call Agglomorate
c
        if(LTestMultiGrid) then
          indexMGPrint=1
          call PrintMultigrid(indexMGPrint)
          indexMGPrint=2
          call PrintMultigrid(indexMGPrint)
        endif
c
        indexGeometric=0
c
      else
c
        if(MGType.eq.'algebraic'.or.MGType.eq.'geometricelement') then
c
          do i=1,NumberOfElements
c  
            acMG(i)=ac(i)
            bcMG(i)=bc(i)
            Residuals(i)=bc(i)
c
            do j=1,NumberofElementNeighbors(i)
              ij=i+(j-1)*NumberOfElements
c      
              anbMG(ij)=anb(i,j)
c
            enddo
c
          enddo
c
          do iLevel=2,iLevelMax
c              
            call AssembleAgglomeratedLHS(iLevel)
c
          enddo
c
        elseif(MGType.eq.'geometricnode') then
c
          iLevel=1
c
          NumberOfElementsMG(iLevel)=NumberOfElements
c
          do i=1,NumberOfElements
c  
            acMG(i)=ac(i)
            bcMG(i)=bc(i)
            Residuals(i)=bc(i)
c
            do j=1,NumberofElementNeighbors(i)
              ij=i+(j-1)*NumberOfElements
c      
              anbMG(ij)=anb(i,j)
c
            enddo
c
          enddo
c
          do iLevel=2,iLevelMax
c              
            call AssembleAgglomeratedLHSGeometrical(iLevel)
c
          enddo
c
        endif
c
        if(LTestMultiGrid) then
          if(nIterMG.eq.nIterStartApplyingMG.or.
     *                      mod(nIterMG,nprintMG).eq.0)  then
c
            indexMGPrint=2
            call PrintMultigrid(indexMGPrint)
c
          endif
        endif
c
      endif
c
c--- Apply selected MutliGrid Cycle
c
      dphiMG=0.
      Residuals=0.
      ncycle=1
      if(MultiGridCycleType.eq.'vcycle') then
c
        do while(ncycle.le.MaxMultiGridCycles.and.
     *                       FResiduals.gt.rrF*IResiduals)
          call vMultiGridCycle(NF,rrF,itmax,solver)
          call CheckConvergenceMG(NF,1,FResiduals)
          ncycle=ncycle+1
        enddo
c
      elseif(MultiGridCycleType.eq.'fcycle') then
c
        do while(ncycle.le.MaxMultiGridCycles.and.
     *                       FResiduals.gt.rrF*IResiduals)
          call fMultiGridCycle(NF,rrF,itmax,solver)
          call CheckConvergenceMG(NF,1,FResiduals)
          ncycle=ncycle+1
        enddo
c
      elseif(MultiGridCycleType.eq.'wcycle') then
c
        do while(ncycle.le.MaxMultiGridCycles.and.
     *                       FResiduals.gt.rrF*IResiduals)
          call wMultiGridCycle(NF,rrF,itmax,solver)
          call CheckConvergenceMG(NF,1,FResiduals)
          ncycle=ncycle+1
        enddo
c
      endif
c
      indexMGPrint=3
      if (LTestMultiGrid) then
        if(nIterMG.eq.1.or.mod(nIterMG,nprintMG).eq.0)  then
          call PrintMultiGrid(indexMGPrint)
        endif
      endif
c
c--- Save the corrections to the original dphi array
c
      do i=1,NumberOfElements
c  
        dphi(i)=dphiMG(i)
c
      enddo
c
      return
      end