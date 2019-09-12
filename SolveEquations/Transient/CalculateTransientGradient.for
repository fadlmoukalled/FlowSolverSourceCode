c
C#############################################################################################
c
      SUBROUTINE ComputeTransientGradient(phi,phiOld,phiOldOld,dphidt)
c
C#############################################################################################
c
      use User0, only: TransientScheme,dt
      use Geometry1, only:NumberOfElements
      use Transient1, only:dtOld,dtOldOld
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: factor1,factor2,factor3,dtCN
      double precision, dimension(:) :: phi,phiOld,phiOldOld,dphidt
c********************************************************************************************
c
      if(TransientScheme.eq.'euler') then
c
        do i=1,NumberOfElements
c
          dphidt(i)=(phi(i)-phiOld(i))/dt
c
        enddo
c
      elseif(TransientScheme.eq.'cranknicolson1') then
c
        dtCN=dt/2.
c
        do i=1,NumberOfElements
c
          dphidt(i)=(phi(i)-phiOld(i))/dtCN
c
        enddo
c
      elseif(TransientScheme.eq.'cranknicolson2') then
c
          factor1=(dtOld/dt)/(dt+dtOld)
          factor2=1./(dt+dtOld)-(dtOldOld/dt)/(dtOld+dtOldOld)
          factor3=-(dtOld/dt)/(dtOld+dtOldOld)
c
        do i=1,NumberOfElements
c
          dphidt(i)=factor1*phi(i)+factor2*phiOld(i)+
     *                                    factor3*phiOldOld(i)
c
        enddo
c
      elseif(TransientScheme.eq.'adamsmoulton') then
c
        factor1=1./dt+1./(dt+dtOld)
        factor3=(dtOld/dt)/(dtOld+dtOldOld)
        factor2=factor1+factor3
c
        do i=1,NumberOfElements
c
          dphidt(i)=factor1*phi(i)-factor2*phiOld(i)+
     *                                    factor3*phiOldOld(i)
c
        enddo
c
      endif
c
      return
      end