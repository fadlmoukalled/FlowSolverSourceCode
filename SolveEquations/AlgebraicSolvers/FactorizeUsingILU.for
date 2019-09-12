c
C#############################################################################################
c
      SUBROUTINE ILUFactorization
c
C#############################################################################################
c
      use Variables2, only: ac,anb,bc,dc,rc
      use Geometry1,  only: NumberOfElements
      use Geometry3,  only: NBFacesMax,ElementNeighbor,
     *                      NumberofElementNeighbors
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,i1,i2,i3,i4,j1,j2,j3,j4
c********************************************************************************************
c
c---Start the ILU Factorization by setting the dc and rc arrays to zero
c
      dc=0.  
      rc=0.
c
      do i=1,NumberOfElements
c
        dc(i)=ac(i)
c
      enddo
c
      do i1=1,NumberOfElements
c
        dc(i1) = 1.0/dc(i1)
        rc(i1) = bc(i1)
c    
        i3 =NumberofElementNeighbors(i1)
    
        if(i1.ne.NumberOfElements-1) then
c
c--- loop over neighbours of iElement
c
          j2 = 1
          do while(j2.le.i3)
            j4 = ElementNeighbor(i1,j2)
c
c--- for all neighbour j > i do
c
            if((j4.gt.i1) .and. (j4.le.numberOfElements)) then
                j3 = NumberofElementNeighbors(j4)
                i2= 0
                i4 = -1
c
c--- find i4 index to get A[j][i4]
c
              do while((i2.le.j3) .and. (i4 .ne. i1))
                i2 = i2 + 1
                i4 =  ElementNeighbor(j4,i2) 
              enddo
c
c--- Compute A[j][i]*D[i]*A[i][j]
c
              if(i4.eq. i1) then
c
                dc(j4) = dc(j4) - anb(j4,i2)*dc(i1)*anb(i1,j2)
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
