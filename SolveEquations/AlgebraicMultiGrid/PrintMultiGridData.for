c
C#############################################################################################
c
      SUBROUTINE PrintMultiGrid(indexMGPrint)
c
C#############################################################################################
c
      use MultiGrid1, only:ijbeginE,ijbeginN,NumberOfElementsMG,
     *                     NumberofElementNeighborsMG,ElementNeighborMG,
     *                     NumberOfChildren,Children,acMG,anbMG,Parents,
     *                     ElementParent,Residuals,bcMG,dphiMG
      use MultiGrid2, only: iLevelMax,nIterMG
      use User0, only: MGType
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer indexMGPrint,iLevel,nElements
      integer i,j,ij,ij1,k
      double precision sumresiduals,sumbcMG
c********************************************************************************************
c
 10   format("element",1x,"--------------NEighbors-------------------")
 20   format(i6,2x,i6,10(1x,i6))
 30   format("element",1x,"--------Parents and Children--------------")
 40   format("element",1x,"  ac   ---------anb Coefficients----------")
 50   format(i6,1P10E9.2)
 60   format(3x,"element",2x,"Parent")
 70   format(1x,i6,3x,i6)
 80   format(3x,"element",4x,"Residual",10x,"bc",10x,"dphi")
 90   format(4x,i6,3x,1PE9.2,6x,1PE9.2,4x,1PE9.2)
c
      if(indexMGPrint.eq.1) then
c
c---  Print grid connectivities
c
        do iLevel=1,iLevelMax
c        
          write(14,*) 'Connectivities for grid level  ',iLevel,
     *                ' at iteration', nIterMG
          write(14,*) 
          write(14,10) 
          write(14,*) 
          nElements=NumberOfElementsMG(iLevel)
c
          do i=1,nElements      
c
            ij=ijbeginE(iLevel)+i
            j=NumberofElementNeighborsMG(ij)
c
            ij1=ijbeginN(iLevel)+i
              
            write(14,20) i,j,
     *         (ElementNeighborMG(ij1+(k-1)*nElements),k=1,j)
          enddo
c        
        enddo
c
        if(MGType.eq.'algebraic'.or.MGType.eq.'geometricelement') then
c
          do iLevel=2,iLevelMax
c        
            write(14,*) 
            write(14,*) 'Element parent of parents for grid level',
     *                iLevel,' at iteration', nIterMG
            write(14,*) 
            write(14,60) 
            write(14,*) 
            nElements=NumberOfElementsMG(iLevel)
c
            do i=1,nElements      
c
              ij=ijbeginE(iLevel-1)+i
              j=ElementParent(ij)
              ij1=ijbeginE(iLevel-1)+j
c
              write(14,70) i,ij1
            enddo
c        
          enddo
c
        endif
c
        do iLevel=1,iLevelMax-1
c        
          write(14,*) 
          write(14,*) 'Parents and children for grid level  ',iLevel,
     *                ' at iteration', nIterMG
          write(14,*) 
          write(14,30) 
          write(14,*) 
          nElements=NumberOfElementsMG(iLevel)

          do i=1,NumberOfElementsMG(iLevel+1)      
c
            ij=ijbeginE(iLevel)+i
            j=NumberOfChildren(ij)
c
            ij1=ijbeginN(iLevel)+i
              
            write(14,20) i,j,
     *         (Children(ij1+(k-1)*nElements),k=1,j)
          enddo
c        
        enddo
c      
      elseif(indexMGPrint.eq.2) then
c      
        do iLevel=1,iLevelMax
c        
          write(14,*) 
          write(14,*) 'Main and neighbor coefficients of Level ',iLevel,
     *                ' at iteration', nIterMG
          write(14,*) 
          write(14,40) 
          write(14,*) 
          nElements=NumberOfElementsMG(iLevel)

          do i=1,nElements      
c
            ij=ijbeginE(iLevel)+i
            j=NumberofElementNeighborsMG(ij)
c
            ij1=ijbeginN(iLevel)+i
              
            write(14,50) i,acMG(ij),
     *         (anbMG(ij1+(k-1)*nElements),k=1,j)
          enddo
c        
        enddo
c      
      elseif(indexMGPrint.eq.3) then
c      
        do iLevel=1,iLevelMax
c        
          write(14,*) 
          write(14,*) 'Residuals of grid level ',iLevel,
     *                ' at iteration', nIterMG
          write(14,*) 
          write(14,80) 
          write(14,*) 
          nElements=NumberOfElementsMG(iLevel)
          sumresiduals=0.
          sumbcMG=0.
c
          do i=1,nElements      
c
            ij=ijbeginE(iLevel)+i
c
            sumbcMG=sumbcMG+bcMG(ij)
            sumresiduals=sumresiduals+Residuals(ij)
            write(14,90) i,Residuals(ij),bcMG(ij),dphiMG(ij)
          enddo
c        
          write(14,*) "sum of bc of level         ",iLevel,
     *    " = ",sumbcMG
          write(14,*) "sum of residuals of level  ",iLevel,
     *    " = ",sumresiduals
        enddo
c      
      endif
c
      return
      end
