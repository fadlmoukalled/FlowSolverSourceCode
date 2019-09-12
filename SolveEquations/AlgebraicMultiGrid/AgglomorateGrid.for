c
C#############################################################################################
c
      SUBROUTINE Agglomorate
c
C#############################################################################################
c
      use User0
      use Variables2, only:ac,anb,bc
      use Geometry1,  only:NumberOfElements
      use Geometry3,  only:NElementFaces,ElementNeighbor,
     *                     NumberofElementNeighbors
      use MultiGrid1, only:NumberOfElementsMG,ElementNeighborMG,
     *                     NumberofElementNeighborsMG,Residuals,
     *                     acMG,bcMG,anbMG,ijbeginE,ijbeginN,
     *                     ijbeginNOld,NElementFacesMG,NBChildrenMaxMG
      use MultiGrid2
c
C#############################################################################################
c
      implicit none
c********************************************************************************************
      integer iLevel
c********************************************************************************************
c
c--- Start by allocating required storage for the MultiGrid algorithm
c
      iLevel=1
      ijbeginNOld(iLevel)=ijbeginN(iLevel)
c
      NElementFacesMG(iLevel)=50 !NElementFaces
      NBChildrenMaxMG(iLevel)=50
c
      if(MGType.eq.'algebraic') then
c
        do while (iLevel.lt.maxNumberofCoarseLevels.and.
     *               NumberOfElementsMG(iLevel).Gt.2*minNumberofParents)
c
          iLevel = iLevel + 1
          call AgglomorateLevelAlgebraicly(iLevel)
          call SetMGCoefficients(iLevel)
          call assembleAgglomeratedLHS(iLevel)
c
        enddo
c
      elseif(MGType.eq.'geometricelement') then
c
        do while (iLevel.lt.maxNumberofCoarseLevels.and.
     *               NumberOfElementsMG(iLevel).Gt.2*minNumberofParents)
c
          iLevel = iLevel + 1
          call AgglomorateLevelGeometricallyElement(iLevel)
          call SetMGCoefficients(iLevel)
          call AssembleAgglomeratedLHS(iLevel)
c
        enddo
c
      elseif(MGType.eq.'geometricnode') then
c
        do while (iLevel.lt.maxNumberofCoarseLevels.and.
     *               NumberOfElementsMG(iLevel).Gt.2*minNumberofParents)
c
          iLevel = iLevel + 1
          call AgglomorateLevelGeometricallyNode(iLevel)
          call SetMGCoefficients(iLevel)
          call AssembleAgglomeratedLHSGeometrical(iLevel)
c
        enddo
c
      endif
c
      iLevelMax=iLevel
c
      return
c
      end