c
C#############################################################################################
c
      SUBROUTINE CrankNicolsonrField
c
C#############################################################################################
c
      use User0, only: NumberOfrFieldsToSolve,LSolverField 
      use VolumeOfFluid1, only: rField,rFieldOld,BrField,BrFieldOld
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,j,irField
c
c********************************************************************************************
c
      do irField=1,NumberOfrFieldsToSolve
c
        if(LSolverField(irField)) then
c        
          do i=1,NumberOfElements
c
            rField(i,irField)=2.*rField(i,irField)-rFieldOld(i,irField)
c
          enddo       
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BrField(i,j,irField)=2.*BrField(i,j,irField)-
     *                                      BrFieldOld(i,j,irField)
c
            enddo
          enddo        
c        
        endif
c
      enddo
c
      return
      end