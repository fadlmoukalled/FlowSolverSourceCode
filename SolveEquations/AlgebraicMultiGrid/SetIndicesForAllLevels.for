c
C#############################################################################################
c
      SUBROUTINE SetLevelIndices(iLevel)
c
C#############################################################################################
      use MultiGrid1, only:ijbeginE,ijbeginN,NumberOfElementsMG,
     *                     NElementFacesMG,NBChildrenMaxMG
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer, dimension(:), allocatable :: ijbegint
      integer :: i,iLevel
c*********************************************************************************************
c
      allocate(ijbegint(iLevel))
c
c--- Elements
c
      do i=1,iLevel-1
        ijbegint(i)=ijbeginE(i)
      enddo
c
      ijbegint(iLevel)=ijbegint(iLevel-1)+
     *               NumberOfElementsMG(iLevel-1)
c
      deallocate(ijbeginE)
      allocate(ijbeginE(iLevel))
c
      ijbeginE=ijbegint
c
c--- Neighbors
c
      do i=1,iLevel-1
        ijbegint(i)=ijbeginN(i)
      enddo
c
      ijbegint(iLevel)=ijbegint(iLevel-1)+
     *           NumberOfElementsMG(iLevel-1)*NElementFacesMG(iLevel-1)
c      ijbegint(iLevel)=ijbegint(iLevel-1)+
c     *           NumberOfElementsMG(iLevel-1)*NBChildrenMaxMG(iLevel-1)
c
      deallocate(ijbeginN)
      allocate(ijbeginN(iLevel))
c
      ijbeginN=ijbegint
c
      deallocate(ijbegint)
c      
      Return
      end
