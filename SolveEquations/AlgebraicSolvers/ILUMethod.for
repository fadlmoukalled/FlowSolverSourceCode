c
C#############################################################################################
c
      SUBROUTINE SolveEquationsUsingILU
c
C#############################################################################################
c
      use Variables2, only: ac,anb,bc,dc,rc,dphi
      use Geometry1,  only: NumberOfElements
      use Geometry3,  only: ElementNeighbor,NumberofElementNeighbors
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,i1,i2,i3,j1,j2,j3,k1
      double precision res,product1,product2
c********************************************************************************************
c
c---  Start the ILU iterations by updating the residual array
c
c--------------------------------------------------------------------------------------------
c
      do i=1,NumberOfElements
c
        res = -ac(i)*dphi(i) + bc(i)
c
        do j = 1,NumberofElementNeighbors(i)

          k = ElementNeighbor(i,j)
c
          if(k.ne.0) then
c
            res = res - anb(i,j)*dphi(k)
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
      do i1=1,NumberOfElements
c
        product1 = dc(i1)*rc(i1)
c
        i3=NumberofElementNeighbors(i1)
    
        j2 = 0
c
c--- loop over neighbours of i
c
        do while((j2+1).le.i3)
c
          j2 = j2 +1
          j1 = ElementNeighbor(i1,j2)
c
          if(j1.ne.0) then
c
c--- for all neighbour j > i do
c
            if(j1 .gt. i1 .and. j1.le. NumberOfElements) then
c
              j3 = NumberofElementNeighbors(j1) 
              i2= 0
              k = 0
c
c--- Get A[j][i]
c
              do while( ((i2+1).le.j3) .and. (k .ne. i1))
c
                i2 = i2 +1
                k = ElementNeighbor(j1,i2)
c
              enddo
c
c--- Compute rc
c
              if(k.ne.0) then
c
                if(k .eq. i1 ) then
c
                  product2 =  anb(j1,i2)*product1 
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
      do i1=NumberOfElements,1,-1
c
        if(i1 .lt. (NumberOfElements)) then
c
          i3=NumberofElementNeighbors(i1)
          j2 = 0
c
c--- Loop over neighbours of i
c
          do while((j2+1) .le. i3)
c
            j2 = j2 +1
             j=ElementNeighbor(i1,j2)
c
            if(j>i1)then
c
              rc(i1) = rc(i1) - anb(i1,j2)*rc(j)
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
        dphi(i1) = dphi(i1) + product1
c
      enddo
c
      return
      end