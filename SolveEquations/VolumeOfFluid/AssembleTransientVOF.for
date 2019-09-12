c
C#############################################################################################
c
      SUBROUTINE AssembleTransientTermrField
     *                       (Variable,FiT,FiTold,FiToldold)
c
C#############################################################################################
c
      use User0, only: dt,TransientSchemerField
      use Transient1
      use Variables3
      use Geometry1, only: NumberOfElements
      use Geometry4, only: Volume
      use MultiGrid2, only: nIter
      use VolumeOfFluid1, only: cosThetaE
c********************************************************************************************
      implicit none      
c********************************************************************************************
      integer i,j,k
      character*10 Variable
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: FiTold
      double precision, dimension(:) :: FiToldold
      double precision :: FluxCElocal,FluxCEoldlocal,FluxCEoldoldlocal
      double precision :: dtCN
      double precision :: factor1,factor2,factor3,factor
      double precision :: phiRright,phiRleft
c********************************************************************************************
c
      if(TransientSchemerField.eq.'euler') then
c
        do i=1,NumberOfElements
c
          FluxCElocal=Volume(i)/dt
          FluxCEoldlocal=-Volume(i)/dt
c
          FluxCE(i)=FluxCE(i)+FluxCElocal
          FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
          FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
        enddo
c
      elseif(TransientSchemerField.eq.'cranknicolson1') then
c
        dtCN=dt/2.
c
        do i=1,NumberOfElements
c
          FluxCElocal=Volume(i)/dtCN
          FluxCEoldlocal=-Volume(i)/dtCN
c
          FluxCE(i)=FluxCE(i)+FluxCElocal
          FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
c
          FluxTE(i)=FluxTE(i)+
     *                FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)
c
        enddo
c
      elseif(TransientSchemerField.eq.'cranknicolson2') then
c
        factor1=(dtOld/dt)/(dt+dtOld)
        factor2=1./(dt+dtOld)-(dtOldOld/dt)/(dtOld+dtOldOld)
        factor3=-(dtOld/dt)/(dtOld+dtOldOld)
c
        do i=1,NumberOfElements
c
          FluxCElocal=factor1*Volume(i)
          FluxCEoldlocal=factor2*Volume(i)
          FluxCEoldoldlocal=factor3*Volume(i)
c
          FluxCE(i)=FluxCE(i)+FluxCElocal
          FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
          FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
          FluxTE(i)=FluxTE(i)+
     *            FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
        enddo
c
      elseif(TransientSchemerField.eq.'adamsmoulton') then
c
        factor1=1./dt+1./(dt+dtOld)
        factor3=(dtOld/dt)/(dtOld+dtOldOld)
        factor2=factor1+factor3
c
        do i=1,NumberOfElements
c
          FluxCElocal=factor1*Volume(i)
          FluxCEoldlocal=-factor2*Volume(i)
          FluxCEoldoldlocal=factor3*Volume(i)
c
          FluxCE(i)=FluxCE(i)+FluxCElocal
          FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
          FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
          FluxTE(i)=FluxTE(i)+
     *            FluxCElocal*FiT(i)+FluxCEoldlocal*FiTold(i)+
     *                                FluxCEoldoldlocal*FiToldold(i)
c
        enddo
c
      elseif(TransientSchemerField.eq.'tics1.75') then
c
        call CalculateCosineThetaElement
c
        do i=1,NumberOfElements
c
          FluxCElocal=Volume(i)/dt
          FluxCEoldlocal=0.
          FluxCEoldoldlocal=0.
c
          factor=cosThetaE(i)*cosThetaE(i)*cosThetaE(i)*cosThetaE(i)
c
          phiRright=
     *       dmax1(dmin1(1.75*FiT(i)-0.75*FiTold(i),1.),0.)*factor+
     *       dmax1(dmin1(1.5*FiT(i)-0.5*FiTold(i),1.),0.)*(1.-factor)
          phiRleft=
     *    dmax1(dmin1(1.75*FiTold(i)-0.75*FiToldold(i),1.),0.)*factor+
     *    dmax1(dmin1(1.5*FiTold(i)-0.5*FiToldold(i),1.),0.)*(1.-factor)
c                
          FluxCE(i)=FluxCE(i)+FluxCElocal
          FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
          FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
          FluxTE(i)=FluxTE(i)+FluxCElocal*(phiRright-phiRleft)
c
        enddo
c
      elseif(TransientSchemerField.eq.'tics2.5') then
c
        call CalculateCosineThetaElement
c
        do i=1,NumberOfElements
c
          FluxCElocal=Volume(i)/dt
          FluxCEoldlocal=0.
          FluxCEoldoldlocal=0.
c
          factor=cosThetaE(i)*cosThetaE(i)*cosThetaE(i)*cosThetaE(i)
c
          phiRright=
     *       dmax1(dmin1(2.5*FiT(i)-1.5*FiTold(i),1.),0.)*factor+
     *       dmax1(dmin1(1.5*FiT(i)-0.5*FiTold(i),1.),0.)*(1.-factor)
          phiRleft=
     *    dmax1(dmin1(2.5*FiTold(i)-1.5*FiToldold(i),1.),0.)*factor+
     *    dmax1(dmin1(1.5*FiTold(i)-0.5*FiToldold(i),1.),0.)*(1.-factor)
c                
          FluxCE(i)=FluxCE(i)+FluxCElocal
          FluxCEold(i)=FluxCEold(i)+FluxCEoldlocal
          FluxCEoldold(i)=FluxCEoldold(i)+FluxCEoldoldlocal
c
          FluxTE(i)=FluxTE(i)+FluxCElocal*(phiRright-phiRleft)
c
        enddo
c
      endif
c
      return
      end