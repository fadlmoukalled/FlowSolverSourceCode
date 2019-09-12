C
C############################################################################################
      SUBROUTINE CalculateResiduals(NF,FiT)
C############################################################################################
C
      use User0
      use Residuals1
      use constants1
      use Geometry1, only: NumberOfElements
      use Geometry3, only: ElementNeighbor,NumberofElementNeighbors
      use Variables2, only: ac,anb,bc,dphi
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,NF
      double precision AbsoluteResidual,RMSResiduals,acphic
      double precision, dimension(:) :: FiT
c********************************************************************************************
      ResorAbs(NF)=0.
      ResorMax(NF)=-1.
      ResorScaled(NF)=-1
      if(nf.eq.4) ResorScaled(NF)=0.
c
      AbsoluteResidual=0.
      RMSResiduals=0.
      acphic=-1.
C
C---- Start calculating the maximum ac*phic
C
	do i=1,NumberOfElements
C
        acphic=dmax1(acphic,dabs(ac(i)*FiT(i)))
C
	enddo
C
C---- Calculate all types of residuals	
C
	do i=1,NumberOfElements
C
        AbsoluteResidual=bc(i)-ac(i)*dphi(i)
c
        do j = 1,NumberofElementNeighbors(i)
c
          k = ElementNeighbor(i,j)
          if(k.ne.0) then
c
            AbsoluteResidual=AbsoluteResidual-anb(i,j)*dphi(k)
c
          endif

        enddo
C
        AbsoluteResidual=dabs(AbsoluteResidual)
        ResorAbs(NF)=ResorAbs(NF)+AbsoluteResidual
        ResorMax(NF)=dmax1(ResorMax(NF),AbsoluteResidual)
        RMSResiduals=RMSResiduals+AbsoluteResidual*AbsoluteResidual
c
        if(NF.ne.4) then
c
          ResorScaled(NF)=dmax1(ResorScaled(NF),
     *                           AbsoluteResidual/(acphic+Tiny))
c
        endif
C
	enddo
C	
      ResorRMS(NF)=DSQRT(RMSResiduals/NumberOfElements)
C
      return
      end
