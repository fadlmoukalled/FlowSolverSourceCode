c
C#############################################################################################
c
      SUBROUTINE ILUFactorizationMG(iLevel)
c
C#############################################################################################
c
      use Variables2, only: dc,rc
      use MultiGrid1, only:ijbeginE,ijbeginN,acMG,anbMG,bcMG,dphiMG,
     *                     NumberOfElementsMG,
     *                     NumberofElementNeighborsMG,ElementNeighborMG
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,i1,i2,i3,i4,j2,j3,j4,iLevel
      integer ij,ij1,ij4,ij5
c********************************************************************************************
c
c---Start the ILU Factorization by setting the dc and rc arrays to zero
c
      dc=0.  
      rc=0.
c
      do i=1,NumberOfElementsMG(iLevel)
c
        ij=ijbeginE(iLevel)+i
        dc(i)=acMG(ij)
c
      enddo
c
      do i1=1,NumberOfElementsMG(iLevel)
c
        ij=ijbeginE(iLevel)+i1
        dc(i1) = 1.0/dc(i1)
        rc(i1) = bcMG(ij)
c    
        i3 =NumberofElementNeighborsMG(ij)
    
        if(i1.ne.NumberOfElementsMG(iLevel)-1) then
c
c--- loop over neighbours of iElement
c
          j2 = 1
          do while(j2.le.i3)
c
            ij1=ijbeginN(iLevel)+i1+(j2-1)*NumberOfElementsMG(iLevel)
            j4 = ElementNeighborMG(ij1)
c
c--- for all neighbour j > i do
c
            if((j4.gt.i1) .and. (j4.le.NumberOfElementsMG(iLevel))) then
              ij4=ijbeginE(iLevel)+j4
                j3 = NumberofElementNeighborsMG(ij4)
                i2= 0
                i4 = -1
c
c--- find i4 index to get A[j][i4]
c
              do while((i2.le.j3) .and. (i4 .ne. i1))
                i2 = i2 + 1
                ij5=ijbeginN(iLevel)+j4+
     *                     (i2-1)*NumberOfElementsMG(iLevel)
                i4 =  ElementNeighborMG(ij5) 
              enddo
c
c--- Compute A[j][i]*D[i]*A[i][j]
c
              if(i4.eq. i1) then
c
                dc(j4) = dc(j4) - anbMG(ij5)*dc(i1)*anbMG(ij1)
c
              else
c
                print*,'the index for i in j is not found'
c
              endif
            endif
c
            j2 = j2 + 1
c
          enddo
        endif
      enddo
c
      return
      end
