c
c#############################################################################################
c
      SUBROUTINE CalculateStrainRateTensor
c
C#############################################################################################
c
      use Variables1, only: uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz
      use Turbulence1, only:S11,S12,S13,S22,S23,S33,
     *                      BS11,BS12,BS13,BS22,BS23,BS33,
     *                      StrainRate,BStrainRate
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c 
      do i=1,NumberOfElements    
c
        S11(i)=uVelGradx(i)
        S12(i)=0.5*(uVelGrady(i)+vVelGradx(i))
        S13(i)=0.5*(uVelGradz(i)+wVelGradx(i))
        S22(i)=vVelGrady(i)
        S23(i)=0.5*(vVelGradz(i)+wVelGrady(i))
        S33(i)=wVelGradz(i)
c          
        StrainRate(i)=dsqrt(2.*
     *         (S11(i)*S11(i)+S12(i)*S12(i)+S13(i)*S13(i)+
     *          S12(i)*S12(i)+S22(i)*S22(i)+S23(i)*S23(i)+
     *          S13(i)*S13(i)+S23(i)*S23(i)+S33(i)*S33(i)))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BS11(i,j)=BuVelGradx(i,j)
          BS12(i,j)=0.5*(BuVelGrady(i,j)+BvVelGradx(i,j))
          BS13(i,j)=0.5*(BuVelGradz(i,j)+BwVelGradx(i,j))
          BS22(i,j)=BvVelGrady(i,j)
          BS23(i,j)=0.5*(BvVelGradz(i,j)+BwVelGrady(i,j))
          BS33(i,j)=BwVelGradz(i,j)
c          
          BStrainRate(i,j)=dsqrt(2.*
     *       (BS11(i,j)*BS11(i,j)+BS12(i,j)*BS12(i,j)+
     *          BS13(i,j)*BS13(i,j)+BS12(i,j)*BS12(i,j)+
     *             BS22(i,j)*BS22(i,j)+BS23(i,j)*BS23(i,j)+
     *              BS13(i,j)*BS13(i,j)+BS23(i,j)*BS23(i,j)+
     *                   BS33(i,j)*BS33(i,j)))
c
        enddo
      enddo
c
      return
      end