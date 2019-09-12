c
C#############################################################################################
c
      SUBROUTINE BoundTurbulentKE
c
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TurbulentKE,BTurbulentKE
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TurbulentKE(i)=dmax1(TurbulentKE(i),tiny)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentKE(i,j)=dmax1(BTurbulentKE(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE BoundTurbulentED
c
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TurbulentED,BTurbulentED
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TurbulentED(i)=dmax1(TurbulentED(i),tiny)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentED(i,j)=dmax1(BTurbulentED(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE BoundTurbulentOmega
c
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TurbulentOmega,BTurbulentOmega
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TurbulentOmega(i)=dmax1(TurbulentOmega(i),tiny)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentOmega(i,j)=dmax1(BTurbulentOmega(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE BoundTReTheta
c
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TReTheta,BTReTheta
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TReTheta(i)=dmax1(TReTheta(i),tiny)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTReTheta(i,j)=dmax1(BTReTheta(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE BoundTGamma
c
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TGamma,BTGamma
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TGamma(i)=dmin1(TGamma(i),1.)
        TGamma(i)=dmax1(TGamma(i),0.)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTGamma(i,j)=dmin1(BTGamma(i,j),1.)
          BTGamma(i,j)=dmax1(BTGamma(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE BoundTfRelaxation
c
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TfRelaxation,BTfRelaxation
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TfRelaxation(i)=dmin1(TfRelaxation(i),1.d6)
        TfRelaxation(i)=dmax1(TfRelaxation(i),0.)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTfRelaxation(i,j)=dmin1(BTfRelaxation(i,j),1.d6)
          BTfRelaxation(i,j)=dmax1(BTfRelaxation(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE BoundTurbulentV2
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TurbulentV2,BTurbulentV2
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TurbulentV2(i)=dmax1(TurbulentV2(i),tiny)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentV2(i,j)=dmax1(BTurbulentV2(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE BoundTurbulentZeta
c
C#############################################################################################
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TurbulentZeta,BTurbulentZeta
      use Constants1, only: twothird
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TurbulentZeta(i)=dmin1(TurbulentZeta(i),2.)
        TurbulentZeta(i)=dmax1(TurbulentZeta(i),0.)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentZeta(i,j)=dmin1(BTurbulentZeta(i,j),2.)
          BTurbulentZeta(i,j)=dmax1(BTurbulentZeta(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE BoundTurbulentKL
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: TurbulentKL,BTurbulentKL
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        TurbulentKL(i)=dmax1(TurbulentKL(i),tiny)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BTurbulentKL(i,j)=dmax1(BTurbulentKL(i,j),0.)
c         
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE BoundModifiedED
c
c#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Variables1, only: ModifiedED,BModifiedED
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        ModifiedED(i)=dmax1(ModifiedED(i),tiny)
c         
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BModifiedED(i,j)=dmax1(BModifiedED(i,j),0.)
c         
        enddo
      enddo
c
      return
      end