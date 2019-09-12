c
C#############################################################################################
      SUBROUTINE CalculateLastrField
C#############################################################################################
      use User0, only: NumberOfrFieldsToSolve,LSolverField
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use VolumeOfFluid1, only: rField,BrField
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k,n,indexrField
      double precision :: sum
c********************************************************************************************
c
      indexrField=0
      do i=1,NumberOfrFieldsToSolve
        if(.not.LSolverField(i)) then
          indexrField=1
          k=i
        endif
      enddo    
c
      if(indexrField.eq.1) then
c
        do i=1,NumberOfElements
          rField(i,k)=1.
        enddo
c
        do j=1,NumberOfrFieldsToSolve
          if(LSolverField(j)) then
            do i=1,NumberOfElements
              rField(i,k)=rField(i,k)-rField(i,j)
            enddo    
          endif
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
            BrField(i,j,k)=1.
          enddo
        enddo
c
        do n=1,NumberOfrFieldsToSolve
          if(LSolverField(n)) then
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
                BrField(i,j,k)=BrField(i,j,k)-BrField(i,j,n)
              enddo
            enddo
          endif
        enddo
c
      else
c          
c--- Bound the r fields          
c
        do i=1,NumberOfElements
          sum=0.
          do j=1,k
            sum=sum+rField(i,j)
          enddo    
          do j=1,k
            rField(i,j)=rField(i,j)/dmax1(sum,tiny)
          enddo    
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
            sum=0.
            do n=1,k
              sum=sum+BrField(i,j,n)
            enddo
            do n=1,k
              BrField(i,j,n)=BrField(i,j,n)/dmax1(sum,tiny)
            enddo
          enddo
        enddo
c          
      endif
c      
      return
      end