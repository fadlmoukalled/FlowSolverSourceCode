c
C#############################################################################################
      SUBROUTINE CalculateResidualsMG(iLevel)
C#############################################################################################
      use Multigrid1, only:ijbeginE,ijbeginN,acMG,bcMG,anbMG,Residuals,
     *                    NumberofElementNeighborsMG,NumberOfElementsMG,
     *                    ElementNeighborMG,dphiMG
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer iLevel,i,j,k,ij,ij1,ij2,ij3
c*********************************************************************************************
c
	do i=1,NumberOfElementsMG(iLevel)
c
        ij=ijbeginE(iLevel)+i
        ij1=ijbeginN(iLevel)+i
        Residuals(ij)=bcMG(ij)-acMG(ij)*dphiMG(ij)
c
        do j=1,NumberofElementNeighborsMG(ij)
          ij2=ij1+(j-1)*NumberOfElementsMG(iLevel)
          k = ElementNeighborMG(ij2)
c
          if(k.gt.0) then
c
            ij3=ijbeginE(iLevel)+k
            Residuals(ij)=Residuals(ij)-anbMG(ij2)*dphiMG(ij3)
c
          endif

        enddo
c
      enddo
c
      return
      end