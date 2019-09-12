c
C#############################################################################################
c
      SUBROUTINE AssembleAgglomeratedRHS(iLevel)
C#############################################################################################
      use User0
      use Multigrid1, only:ijbeginE,bcMG,Residuals,Parents,
     *                    NumberOfElementsMG
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,iLevel,iLevel1,ijF,ijC,iParent
c*********************************************************************************************
c
      iLevel1=iLevel-1
c
      do i=1,NumberOfElementsMG(iLevel)
c
        ijC=ijbeginE(iLevel)+i
c
        bcMG(ijC) = 0.

      enddo
c
      do i=1,NumberOfElementsMG(iLevel1)
c
        ijF=ijbeginE(iLevel1)+i
        iParent=Parents(ijF)
        ijC=ijbeginE(iLevel)+iParent
c
        bcMG(ijC) = bcMG(ijC) + Residuals(ijF)
c
      enddo
c
      return
      end
