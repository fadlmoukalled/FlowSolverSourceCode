c
C#############################################################################################
      SUBROUTINE sort(arr)
C#############################################################################################
      use swapData, only:swapDPvariables,maskedSwapDPvariables
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      interface
c*********************************************************************************************
        SUBROUTINE nrerror(string)
          character(LEN=*), intent(in) :: string
        end SUBROUTINE nrerror
c*********************************************************************************************
      end interface
c
c*********************************************************************************************
      double precision, dimension(:), intent(inout) :: arr
      integer, parameter :: NN=15, NSTACK=50
c
c*********************************************************************************************
c
c      Sorts an array arr into ascending numerical order using the Quicksort algorithm. arr is
c      replaced on output by its sorted rearrangement.
c      Parameters: NN is the size of subarrays sorted by straight insertion and NSTACK is the
c      required auxiliary storage.
c
c*********************************************************************************************
c
      double precision :: a
      integer :: n,k,i,j,jstack,l,r
      integer, dimension(nstack) :: istack
c
c*********************************************************************************************
c
      n=size(arr)
      jstack=0
      l=1
      r=n
c
      do
        if (r-l < NN) then                   
          do j=l+1,r
c
            a=arr(j)
c
            do i=j-1,l,-1
c
              if (arr(i) <= a) exit
              arr(i+1)=arr(i)
c
            enddo
c
            arr(i+1)=a
c
          enddo
c
          if (jstack == 0) RETURN
          r=istack(jstack)                          
          l=istack(jstack-1)                        
          jstack=jstack-2
c
        else                                         
c
          k=(l+r)/2                                   
          call swapDPvariables(arr(k),arr(l+1))                  
          call maskedSwapDPvariables(arr(l),arr(r),arr(l)>arr(r))
          call maskedSwapDPvariables(arr(l+1),arr(r),arr(l+1)>arr(r))
          call maskedSwapDPvariables(arr(l),arr(l+1),arr(l)>arr(l+1))
          i=l+1                                       
          j=r
          a=arr(l+1)                                  
c
          do                                          
            do                                        
c
              i=i+1
              if (arr(i) >= a) exit
c
            enddo
            do                                       
c
              j=j-1
              if (arr(j) <= a) exit
c
            enddo
c
            if (j < i) exit                             
            call swapDPvariables(arr(i),arr(j))                 
c
          enddo
c
          arr(l+1)=arr(j)                              
          arr(j)=a
          jstack=jstack+2
c
          if (jstack > NSTACK) call nrerror("sort:NSTACK too small")
          if (r-i+1 >= j-l) then
c
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
c
          else
c
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
c
          endif
        endif
      enddo
c
      return
      end