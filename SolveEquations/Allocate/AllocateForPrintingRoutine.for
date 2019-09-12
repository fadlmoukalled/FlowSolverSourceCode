c
C#############################################################################################
      SUBROUTINE allocateForPrinting
C#############################################################################################
      use Turbulence1, only: S11,S12,S13,S22,S23,S33,
     *                       BS11,BS12,BS13,BS22,BS23,BS33,
     *                       StrainRate,BStrainRate,
     *                       W11,W12,W13,W22,W23,W33,
     *                       BW11,BW12,BW13,BW22,BW23,BW33,
     *                       Vorticity,BVorticity
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFacesMax
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      entry AllocateStrainRateTensor
c*********************************************************************************************
c
      allocate(S11(NumberOfElements))
      allocate(S12(NumberOfElements))
      allocate(S13(NumberOfElements))
      allocate(S22(NumberOfElements))
      allocate(S23(NumberOfElements))
      allocate(S33(NumberOfElements))
      allocate(BS11(NumberOfBCSets,NBFacesMax))
      allocate(BS12(NumberOfBCSets,NBFacesMax))
      allocate(BS13(NumberOfBCSets,NBFacesMax))
      allocate(BS22(NumberOfBCSets,NBFacesMax))
      allocate(BS23(NumberOfBCSets,NBFacesMax))
      allocate(BS33(NumberOfBCSets,NBFacesMax))
      allocate(StrainRate(NumberOfElements))
      allocate(BStrainRate(NumberOfBCSets,NBFacesMax))
c
      return
c*********************************************************************************************
      entry deAllocateStrainRateTensor
c*********************************************************************************************
c
      deallocate(S11)
      deallocate(S12)
      deallocate(S13)
      deallocate(S22)
      deallocate(S23)
      deallocate(S33)
      deallocate(BS11)
      deallocate(BS12)
      deallocate(BS13)
      deallocate(BS22)
      deallocate(BS23)
      deallocate(BS33)
      deallocate(StrainRate)
      deallocate(BStrainRate)
c
      return
c*********************************************************************************************
      entry AllocateVorticityTensor
c*********************************************************************************************
c
      allocate(W11(NumberOfElements))
      allocate(W12(NumberOfElements))
      allocate(W13(NumberOfElements))
      allocate(W22(NumberOfElements))
      allocate(W23(NumberOfElements))
      allocate(W33(NumberOfElements))
      allocate(BW11(NumberOfBCSets,NBFacesMax))
      allocate(BW12(NumberOfBCSets,NBFacesMax))
      allocate(BW13(NumberOfBCSets,NBFacesMax))
      allocate(BW22(NumberOfBCSets,NBFacesMax))
      allocate(BW23(NumberOfBCSets,NBFacesMax))
      allocate(BW33(NumberOfBCSets,NBFacesMax))
      allocate(Vorticity(NumberOfElements))
      allocate(BVorticity(NumberOfBCSets,NBFacesMax))
c
      return
c*********************************************************************************************
      entry deAllocateVorticityTensor
c*********************************************************************************************
c
      deallocate(W11)
      deallocate(W12)
      deallocate(W13)
      deallocate(W22)
      deallocate(W23)
      deallocate(W33)
      deallocate(BW11)
      deallocate(BW12)
      deallocate(BW13)
      deallocate(BW22)
      deallocate(BW23)
      deallocate(BW33)
      deallocate(Vorticity)
      deallocate(BVorticity)
c
      return
      end