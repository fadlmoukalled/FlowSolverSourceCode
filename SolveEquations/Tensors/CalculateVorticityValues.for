c
c#############################################################################################
c
      SUBROUTINE CalculateVorticityTensor
c
C#############################################################################################
c
      use Variables1, only: uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz
      use Turbulence1, only:W11,W12,W13,W22,W23,W33,
     *                      BW11,BW12,BW13,BW22,BW23,BW33,
     *                      Vorticity,BVorticity
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use User0, only: LRotation
      use Coriolis1, only: AngularVelocityX,AngularVelocityY,
     *                     AngularVelocityZ

c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j
c********************************************************************************************
c 
      do i=1,NumberOfElements    
c
        W11(i)=0.
        W12(i)=0.5*(uVelGrady(i)-vVelGradx(i))
        W13(i)=0.5*(uVelGradz(i)-wVelGradx(i))
        W22(i)=0.
        W23(i)=0.5*(vVelGradz(i)-wVelGrady(i))
        W33(i)=0.
c          
        if(LRotation) then
c
          W12(i)=W12(i)-AngularVelocityZ
          W13(i)=W13(i)+AngularVelocityY
          W23(i)=W23(i)-AngularVelocityX
c
        endif     
c
        Vorticity(i)=dsqrt(2.*
     *         (W11(i)*W11(i)+W12(i)*W12(i)+W13(i)*W13(i)+
     *          (-W12(i))*(-W12(i))+W22(i)*W22(i)+W23(i)*W23(i)+
     *          (-W13(i))*(-W13(i))+(-W23(i))*(-W23(i))+W33(i)*W33(i)))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BW11(i,j)=0.
          BW12(i,j)=0.5*(BuVelGrady(i,j)-BvVelGradx(i,j))
          BW13(i,j)=0.5*(BuVelGradz(i,j)-BwVelGradx(i,j))
          BW22(i,j)=0.
          BW23(i,j)=0.5*(BvVelGradz(i,j)-BwVelGrady(i,j))
          BW33(i,j)=0.
c          
          if(LRotation) then
c
            BW12(i,j)=BW12(i,j)-AngularVelocityZ
            BW13(i,j)=BW13(i,j)+AngularVelocityY
            BW23(i,j)=BW23(i,j)-AngularVelocityX
c
          endif     
c
          BVorticity(i,j)=dsqrt(2.*
     *       (BW11(i,j)*BW11(i,j)+BW12(i,j)*BW12(i,j)+
     *          BW13(i,j)*BW13(i,j)+(-BW12(i,j))*(-BW12(i,j))+
     *             BW22(i,j)*BW22(i,j)+BW23(i,j)*BW23(i,j)+
     *              (-BW13(i,j))*(-BW13(i,j))+(-BW23(i,j))*(-BW23(i,j))+
     *                   BW33(i,j)*BW33(i,j)))
c
        enddo
      enddo
c
      return
      end
