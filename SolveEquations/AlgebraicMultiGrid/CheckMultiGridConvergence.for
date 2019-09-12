C
C############################################################################################
      SUBROUTINE CheckConvergenceMG(NF,iLevel,AbsoluteResidual)
C############################################################################################
C
      use User0, only: LPrintMultiGridResiduals,LTestMultiGrid
      use MultiGrid1, only: acMG,anbMG,bcMG,dphiMG,ijbeginE,ijbeginN,
     *                      NumberOfElementsMG,ElementNeighborMG,
     *                      NumberofElementNeighborsMG
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,ij,ij1,ij2,NF,iLevel
      double precision AbsoluteResidual,ElementResidual
c********************************************************************************************
       AbsoluteResidual=0.
       ElementResidual=0.    
C
C---- Calculate all types of residuals	
C
	do i=1,NumberOfElementsMG(iLevel)
C
        ij=ijbeginE(iLevel)+i
        ElementResidual=bcMG(ij)-acMG(ij)*dphiMG(ij)
c
        do j = 1,NumberofElementNeighborsMG(ij)
c
          ij1=ijbeginN(iLevel)+i+(j-1)*NumberOfElementsMG(iLevel)
          k = ElementNeighborMG(ij1)
          if(k.ne.0) then
c
            ij2=ijbeginE(iLevel)+k
            ElementResidual=ElementResidual-anbMG(ij1)*dphiMG(ij2)
c
          endif

        enddo
c
        AbsoluteResidual=AbsoluteResidual+dabs(ElementResidual)
c
        enddo
        if(LPrintMultiGridResiduals) then
          print*,'Grid Level=',iLevel,'Residuals=',AbsoluteResidual
	    if(LTestMultiGrid) then
            write(14,*) 'Grid Level=',iLevel,'Residuals=',
     *                      AbsoluteResidual
          endif
        endif
c        
      return
      end