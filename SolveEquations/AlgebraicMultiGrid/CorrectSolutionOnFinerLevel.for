c
C#############################################################################################
c
      SUBROUTINE CorrectFinerLevelSolution(iLevel)
c
C#############################################################################################
c
      use MultiGrid1, only:NumberOfElementsMG,ijbeginE,Parents,dphiMG
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer iLevel,ilevel1,i,ij,ij1,iParent
c********************************************************************************************
c
      iLevel1=iLevel-1
c
      do i=1,NumberOfElementsMG(iLevel1)
c
        ij=ijbeginE(iLevel1)+i
        iParent=Parents(ij)
        ij1=ijbeginE(iLevel)+iParent
        dphiMG(ij)=dphiMG(ij)+dphiMG(ij1)
c
      enddo
c
      return
      end
