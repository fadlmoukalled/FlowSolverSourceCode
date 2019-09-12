c
C#############################################################################################
c
      SUBROUTINE ResetMultigridAllocation
C#############################################################################################
      use MultiGrid1
      use Geometry1, only: NumberOfElements
c********************************************************************************************
      implicit none
c********************************************************************************************
c
c--- Deallocate variables
c
      deallocate(acMG)
      deallocate(bcMG)
      deallocate(anbMG)
      deallocate(maxanb)
      deallocate(Residuals)
      deallocate(dphiMG)
c
      deallocate(NumberOfElementsMG)
      deallocate(NElementFacesMG)
      deallocate(NBChildrenMaxMG)
      deallocate(ijbeginE)
      deallocate(ijbeginN)
      deallocate(ijbeginNOld)
      deallocate(ElementNeighborMG)
      deallocate(NumberofElementNeighborsMG)
c
      deallocate(Parents)
      deallocate(Children)
      deallocate(NumberOfChildren)
      deallocate(NeighborOfParent)
      deallocate(NumberOfNeighborOfParent)
      deallocate(ElementParent)
c
c--- Reallocate variables
c
      allocate(acMG(NumberOfElements))
      allocate(bcMG(NumberOfElements))
      allocate(anbMG(NumberOfElements*60))
      allocate(maxanb(NumberOfElements))
      allocate(Residuals(NumberOfElements))
      allocate(dphiMG(NumberOfElements))
c
      allocate(NumberOfElementsMG(1))
      allocate(NElementFacesMG(1))
      allocate(NBChildrenMaxMG(1))
      allocate(ijbeginE(1))
      allocate(ijbeginN(1))
      allocate(ijbeginNOld(1))
      allocate(ElementNeighborMG(NumberOfElements*60))
      allocate(NumberofElementNeighborsMG(NumberOfElements))
c
      allocate(Parents(NumberOfElements))
      allocate(Children(NumberOfElements*60))
      allocate(NumberOfChildren(NumberOfElements))
      allocate(NeighborOfParent(NumberOfElements*60))
      allocate(NumberOfNeighborOfParent(NumberOfElements))
      allocate(ElementParent(NumberOfElements))
c
c--- Reinitialize variables
c
      acMG=0.
      bcMG=0.
      anbMG=0.
      maxanb=0.
      Residuals=0.
      dphiMG=0.
c
      NumberOfElementsMG=0
      NElementFacesMG=0
      NBChildrenMaxMG=0
      ijbeginE=0
      ijbeginN=0
      ijbeginNOld=0
      ElementNeighborMG=0
      NumberofElementNeighborsMG=0
c
      Parents=0
      Children=0
      NumberOfChildren=0
      NeighborOfParent=0
      NumberOfNeighborOfParent=0
      ElementParent=0
c
      return
      end
