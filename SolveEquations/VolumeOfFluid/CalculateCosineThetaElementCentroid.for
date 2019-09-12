c
c#############################################################################################
c
      SUBROUTINE CalculateCosineThetaElement
c
c#############################################################################################
      use Geometry1, only: NumberOfElements
      use VolumeOfFluid1, only: cosThetaE
      use TransferrField1, only: rFieldGradxT,rFieldGradyT,
     *                           rFieldGradzT
      use Variables1, only: uVelocity,vVelocity,wVelocity
      use Constants1, only: tiny 
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i
      double precision :: rFieldGrad,VelMagnitude
c********************************************************************************************

      do i=1,NumberOfElements
c
        rFieldGrad=dsqrt(rFieldGradxT(i)*rFieldGradxT(i)+
     *                   rFieldGradyT(i)*rFieldGradyT(i)+
     *                   rFieldGradzT(i)*rFieldGradzT(i))
        VelMagnitude=dsqrt(uVelocity(i)*uVelocity(i)+
     *                     vVelocity(i)*vVelocity(i)+
     *                     wVelocity(i)*wVelocity(i))
        cosThetaE(i)=dabs(uVelocity(i)*rFieldGradxT(i)+
     *                    vVelocity(i)*rFieldGradyT(i)+
     *                    wVelocity(i)*rFieldGradzT(i))/
     *                          (rFieldGrad*VelMagnitude+tiny)
        if(cosThetaE(i).gt.1.) cosThetaE(i)=1.
c
      enddo
c
      return
      end
