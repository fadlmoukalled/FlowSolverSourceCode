c
C#############################################################################################
c
      SUBROUTINE SetMGCoefficients(iLevel)
c
C#############################################################################################
c
      use Multigrid1, only: acMG,bcMG,anbMG,dphiMG,Residuals,ijbeginE,
     *                      ijbeginN,NumberOfElementsMG,NElementFacesMG
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: iLevel,i,ij,ik
      integer :: LengthOfArray1,LengthOfArray2
      integer :: LengthOfArray3,LengthOfArray4
      double precision, dimension(:), allocatable :: anbMGt
      double precision, dimension(:), allocatable :: acMGt
      double precision, dimension(:), allocatable :: bcMGt
      double precision, dimension(:), allocatable :: Residualst
      double precision, dimension(:), allocatable :: dphiMGt
c*********************************************************************************************
c      
      LengthOfArray1=ijbeginN(iLevel)
      LengthOfArray2=ijbeginE(iLevel)
c
      allocate(anbMGt(LengthOfArray1))
      allocate(acMGt(LengthOfArray2))
      allocate(bcMGt(LengthOfArray2))
      allocate(Residualst(LengthOfArray2))
      allocate(dphiMGt(LengthOfArray2))
c
      anbMGt=anbMG
      acMGt=acMG
      bcMGt=bcMG
      Residualst=Residuals
      dphiMGt=dphiMG
c
      deallocate(anbMG)
      deallocate(acMG)
      deallocate(bcMG)
      deallocate(Residuals)
      deallocate(dphiMG)
c      
      LengthOfArray3=ijbeginN(iLevel)+
     *          NumberOfElementsMG(iLevel)*NElementFacesMG(iLevel)
      LengthOfArray4=ijbeginE(iLevel)+NumberOfElementsMG(iLevel)
c
      allocate(anbMG(LengthOfArray3))
      allocate(acMG(LengthOfArray4))
      allocate(bcMG(LengthOfArray4))
      allocate(Residuals(LengthOfArray4))
      allocate(dphiMG(LengthOfArray4))
c
      anbMG=0.
      acMG=0.
      bcMG=0.
      Residuals=0.
      dphiMG=0.
c
      anbMG(1:LengthOfArray1)=anbMGt(1:LengthOfArray1)
      acMG(1:LengthOfArray2)=acMGt(1:LengthOfArray2)
      bcMG(1:LengthOfArray2)=bcMGt(1:LengthOfArray2)
      Residuals(1:LengthOfArray2)=Residualst(1:LengthOfArray2)
      dphiMG(1:LengthOfArray2)=dphiMGt(1:LengthOfArray2)
c
      deallocate(anbMGt)
      deallocate(acMGt)
      deallocate(bcMGt)
      deallocate(Residualst)
      deallocate(dphiMGt)
c
      return
      end