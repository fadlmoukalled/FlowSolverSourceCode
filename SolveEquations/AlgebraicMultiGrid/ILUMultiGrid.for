c
C#############################################################################################
c
      SUBROUTINE SolveEquationsUsingILUMG(iLevel)
c
C#############################################################################################
c
      use Variables2, only: dc,rc
      use MultiGrid1, only:ijbeginE,ijbeginN,acMG,anbMG,bcMG,dphiMG,
     *                     NumberOfElementsMG,
     *                     NumberofElementNeighborsMG,ElementNeighborMG
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,i1,i2,i3,j1,j2,j3,k1,iLevel
      integer ij,ij1,ij2,ij3
      double precision res,product1,product2
c********************************************************************************************
c
c---  Start the ILU iterations by updating the residual array
c
c--------------------------------------------------------------------------------------------
c
      do i=1,NumberOfElementsMG(iLevel)
c
        ij=ijbeginE(iLevel)+i
        res = -acMG(ij)*dphiMG(ij) + bcMG(ij)
c
        do j = 1,NumberofElementNeighborsMG(ij)

          ij1=ijbeginN(iLevel)+i+(j-1)*NumberOfElementsMG(iLevel)
          k = ElementNeighborMG(ij1)
c
          if(k.ne.0) then
c
            res = res - anbMG(ij1)*dphiMG(ijbeginE(iLevel)+k)
c
          endif
c
        enddo
c
        rc(i) = res
c
      enddo
c
c--- Forward substitution
c
      do i1=1,NumberOfElementsMG(iLevel)
c
        product1 = dc(i1)*rc(i1)
c
        ij1=ijbeginE(iLevel)+i1
        i3=NumberofElementNeighborsMG(ij1)
    
        j2 = 0
c
c--- loop over neighbours of i1
c
        do while((j2+1).le.i3)
c
          j2 = j2 +1
          ij2=ijbeginN(iLevel)+i1+(j2-1)*NumberOfElementsMG(iLevel)
          j1 = ElementNeighborMG(ij2)
c
          if(j1.ne.0) then
c
c--- for all neighbour j > i do
c
            if(j1 .gt. i1 .and. j1.le. NumberOfElementsMG(iLevel)) then
c
              j3 = NumberofElementNeighborsMG(ijbeginE(iLevel)+j1) 
              i2= 0
              k = 0
c
c--- Get A[j][i]
c
              do while( ((i2+1).le.j3) .and. (k .ne. i1))
c
                i2 = i2 +1
                ij3=ijbeginN(iLevel)+j1+
     *                    (i2-1)*NumberOfElementsMG(iLevel)
                k = ElementNeighborMG(ij3)
c
              enddo
c
c--- Compute rc
c
              if(k.ne.0) then
c
                if(k .eq. i1 ) then
c
                  product2 =  anbMG(ij3)*product1 
                  rc(j1) = rc(j1) - product2
c
                else
c
                  print*,'ILU Solver Error The index for i
     * in element j is not found'
c
                endif
c
              endif
c
            endif
c
          endif   
c
        enddo
c
      enddo
c
c--- Backward substitution
c
      do i1=NumberOfElementsMG(iLevel),1,-1
c
        ij1=ijbeginE(iLevel)+i1
c
        if(i1 .lt. (NumberOfElementsMG(iLevel))) then
c
          i3=NumberofElementNeighborsMG(ij1)
          j2 = 0
c
c--- Loop over neighbours of i
c
          do while((j2+1) .le. i3)
c
            j2 = j2 +1
            ij2=ijbeginN(iLevel)+i1+
     *                    (j2-1)*NumberOfElementsMG(iLevel)
            j=ElementNeighborMG(ij2)
c
            if(j>i1)then
c
              rc(i1) = rc(i1) - anbMG(ij2)*rc(j)
c
            endif
c
          enddo
c
        endif
c
c--- Compute product D[i]*R[i]
c
        product1 = dc(i1)*rc(i1)
        rc(i1) = product1
c
c--- Update dphi
c
        dphiMG(ij1) = dphiMG(ij1) + product1
c
      enddo
c
      return
      end
