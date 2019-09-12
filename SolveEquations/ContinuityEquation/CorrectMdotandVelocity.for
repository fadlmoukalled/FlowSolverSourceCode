c
c#############################################################################################
c
      SUBROUTINE Correctmdot
c
c#############################################################################################
c
      use User0, only: LCompressible,LFreeSurfaceFlow  
      use Variables1, only: PressureC,PressCGradx,PressCGrady,
     *                      PressCGradz,
     *                      PressCGradfx,PressCGradfy,PressCGradfz,mdot,
     *                      Du2Velocity,Dv2Velocity,Dw2Velocity
      use Variables2, only: FluxCf,FluxFf
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor
      use PhysicalProperties1, only: Densityf,drhodp     
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*16 interpolation
      integer :: i,j,k
      double precision :: rhof
c********************************************************************************************
c
      if(LCompressible) then
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          rhof=Densityf(k)
          mdot(k)=mdot(k)+FluxCf(k)*PressureC(i)+FluxFf(k)*PressureC(j)+
     *            dmax1(mdot(k)/rhof,0.)*drhodP(i)*PressureC(i)-
     *              dmax1(-mdot(k)/rhof,0.)*drhodP(j)*PressureC(j)
c
        enddo
c
      else
c
        if(LFreeSurfaceFlow) then
c
          do k=1,NIFaces
c
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            rhof=Densityf(k)
            mdot(k)=mdot(k)+
     *          rhof*FluxCf(k)*PressureC(i)+rhof*FluxFf(k)*PressureC(j)
c
          enddo
c
        else
c
          do k=1,NIFaces
c
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            mdot(k)=mdot(k)+
     *                FluxCf(k)*PressureC(i)+FluxFf(k)*PressureC(j)
c
          enddo
c
        endif
c
      endif
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CorrectBoundarymdot
c
c#############################################################################################
c
      use User0, only: Lcompressible,MethodDecomposeS,LFreeSurfaceFlow
      use Variables1, only: PressureC,BPressureC,Bmdot,
     *                      Du2Velocity,Dv2Velocity,Dw2Velocity,
     *                      BuVelocity,BvVelocity,BwVelocity,BPressure,
     *                      PressCGradx,PressCGrady,PressCGradz,
     *                      uVelocity,vVelocity,wVelocity,Pressure,
     *                      xVeldirection,yVeldirection,zVeldirection
      use Variables2, only: FluxCf,FluxFf
      use BoundaryConditions2
      use Geometry3, only: NIFaces,NBFaces,NBFaceOwner
      use PhysicalProperties1, only: BDensity,BSpecificHeat,Bdrhodp
      use Geometry4, only: BFaceAreax,BFaceAreay,BFaceAreaz,
     *                     BDistanceCFux,BDistanceCFuy,BDistanceCFuz,
     *                     BDistanceCF,xc,yc,zc,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz,
     *                     BFaceAreax,BFaceAreay,BFaceAreaz
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j,j1,j2,j3
      double precision :: rhof,sfx,sfy,sfz,ex,ey,ez,cf,DuSf,DvSf,DwSf,
     *                    eDuSf,eDvSf,eDwSf,DuEf,DvEf,DwEf,geoDiff,
     *                    Magnitude,velocity2,
     *                    TotalTemp,dpdub,c1,xF1,yF1,zF1,
     *                    distance1,distance2,GFactCF,
     *                    DistCFx,DistCFy,DistCFz,
     *                    DistCF,DistCFux,DistCFuy,DistCFuz,
     *                    pF1x,pF1y,pF1z,
     *                    uF1,vF1,wF1,Ub1,Ub2,term1,term2,term3,term4,
     *                    correction,velMagnitude
c********************************************************************************************
c
      if(.not.Lcompressible) then
c
c----------------------------------------------------------------------
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c        
            rhof=1.
            if(LFreeSurfaceFlow) rhof=BDensity(i2,i3)
c
            BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
            BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
            BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
            Bmdot(i2,i3)=Bmdot(i2,i3)+rhof*FluxCf(i4)*PressureC(i1)+
     *                              rhof*FluxFf(i4)*PressureC(j1)
c
            BuVelocity(j2,j3)=a1r(i2)*BuVelocity(i2,i3)+
     *           b1r(i2)* BvVelocity(i2,i3)+c1r(i2)* BwVelocity(i2,i3)
            BvVelocity(j2,j3)=a2r(i2)*BuVelocity(i2,i3)+
     *           b2r(i2)* BvVelocity(i2,i3)+c2r(i2)* BwVelocity(i2,i3)
            BwVelocity(j2,j3)=a3r(i2)*BuVelocity(i2,i3)+
     *           b3r(i2)* BvVelocity(i2,i3)+c3r(i2)* BwVelocity(i2,i3)
c
            Bmdot(j2,j3)=-Bmdot(i2,i3)
c
          endif
c
        enddo
c
c----------------------------------------------------------------------
c
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i1=IinletSpecifiedStaticPressureOwner(i)
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
          geoDiff=FluxFf(i4)
          if(LFreeSurfaceFlow) geoDiff=rhof*FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
c
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          velMagnitude=Bmdot(i2,i3)/
     *       (rhof*(xVeldirection(i2,i3)*sfx+
     *            yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
c          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
c     *                          Du2Velocity(i1)*PressCGradx(i1)
c          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
c     *                          Dv2Velocity(i1)*PressCGrady(i1)
c          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
c     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationPressure
c
          i1=IinletSpecifiedStagnationPressureOwner(i)
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
          geoDiff=FluxFf(i4)
          if(LFreeSurfaceFlow) geoDiff=rhof*FluxFf(i4)
c
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
c
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          velMagnitude=Bmdot(i2,i3)/
     *       (rhof*(xVeldirection(i2,i3)*sfx+
     *            yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
c          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
c     *                          Du2Velocity(i1)*PressCGradx(i1)
c          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
c     *                          Dv2Velocity(i1)*PressCGrady(i1)
c          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
c     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo

c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedStaticPressure
c
          i1=IoutletSpecifiedStaticPressureOwner(i)
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
          geoDiff=FluxFf(i4)
          if(LFreeSurfaceFlow) geoDiff=rhof*FluxFf(i4)
c
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i1=IoutletSpecifiedAverageStaticPressureOwner(i)
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
          geoDiff=FluxFf(i4)
          if(LFreeSurfaceFlow) geoDiff=rhof*FluxFf(i4)
c
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedResistance
c
          i1=IoutletSpecifiedResistanceOwner(i)
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          rhof=BDensity(i2,i3)
          geoDiff=FluxFf(i4)
          if(LFreeSurfaceFlow) geoDiff=rhof*FluxFf(i4)
c
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
c
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
c
      elseif(Lcompressible) then
c
c----------------------------------------------------------------------
c
        do i=1,IinletSpecifiedVelocity
c
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
c
          Bmdot(i2,i3)=Bmdot(i2,i3)+Bmdot(i2,i3)*BdrhodP(i2,i3)*
     *                            BPressureC(i2,i3)/Bdensity(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
c
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStaticPressure
c
          i1=IinletSpecifiedStaticPressureOwner(i)
          i2=IinletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStaticPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
c
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          rhof=BDensity(i2,i3)
c
          velMagnitude=Bmdot(i2,i3)/
     *       (rhof*(xVeldirection(i2,i3)*sfx+
     *            yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
c          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
c     *                          Du2Velocity(i1)*PressCGradx(i1)
c          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
c     *                          Dv2Velocity(i1)*PressCGrady(i1)
c          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
c     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedStagnationPressure
c
          i1=IinletSpecifiedStagnationPressureOwner(i)
          i2=IinletSpecifiedStagnationPressureNumberOfBCSets(i)
          i3=IinletSpecifiedStagnationPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
c
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          rhof=BDensity(i2,i3)
c
          velMagnitude=Bmdot(i2,i3)/
     *       (rhof*(xVeldirection(i2,i3)*sfx+
     *            yVeldirection(i2,i3)*sfy+zVeldirection(i2,i3)*sfz))
	    BuVelocity(i2,i3)=velMagnitude*xVeldirection(i2,i3)
	    BvVelocity(i2,i3)=velMagnitude*yVeldirection(i2,i3)
	    BwVelocity(i2,i3)=velMagnitude*zVeldirection(i2,i3)
c
c          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
c     *                          Du2Velocity(i1)*PressCGradx(i1)
c          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
c     *                          Dv2Velocity(i1)*PressCGrady(i1)
c          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
c     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IpressureFarField
c
          i1=IpressureFarFieldOwner(i)
          i2=IpressureFarFieldNumberOfBCSets(i)
          i3=IpressureFarFieldNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletTransmissive
c
          i1=IoutletTransmissiveOwner(i)
          i2=IoutletTransmissiveNumberOfBCSets(i)
          i3=IoutletTransmissiveNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Ioutletsupersonic
c
          i1=IoutletsupersonicOwner(i)
          i2=IoutletsupersonicNumberOfBCSets(i)
          i3=IoutletsupersonicNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
c
          Bmdot(i2,i3)=Bmdot(i2,i3)+Bmdot(i2,i3)*BdrhodP(i2,i3)*
     *                            BPressureC(i2,i3)/Bdensity(i2,i3)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedStaticPressure
c
          i1=IoutletSpecifiedStaticPressureOwner(i)
          i2=IoutletSpecifiedStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedStaticPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedAverageStaticPressure
c
          i1=IoutletSpecifiedAverageStaticPressureOwner(i)
          i2=IoutletSpecifiedAverageStaticPressureNumberOfBCSets(i)
          i3=IoutletSpecifiedAverageStaticPressureNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedResistance
c
          i1=IoutletSpecifiedResistanceOwner(i)
          i2=IoutletSpecifiedResistanceNumberOfBCSets(i)
          i3=IoutletSpecifiedResistanceNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          geoDiff=FluxFf(i4)
          Bmdot(i2,i3)=Bmdot(i2,i3)-
     *               geoDiff*(BPressureC(i2,i3)-PressureC(i1))
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
c
          BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
          BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
          BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
c
        enddo
c----------------------------------------------------------------------
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
c
            i4=NIFaces
c
            do j=1,i2-1
c
              i4=i4+NBFaces(j)
c
            enddo
c
            i4=i4+i3
c
            j3=Icorrespondingface(i2,i3)
            j1=NBFaceOwner(j2,j3)
c
            BuVelocity(i2,i3)=BuVelocity(i2,i3)-
     *                          Du2Velocity(i1)*PressCGradx(i1)
            BvVelocity(i2,i3)=BvVelocity(i2,i3)-
     *                          Dv2Velocity(i1)*PressCGrady(i1)
            BwVelocity(i2,i3)=BwVelocity(i2,i3)-
     *                          Dw2Velocity(i1)*PressCGradz(i1)
            Bmdot(i2,i3)=Bmdot(i2,i3)+FluxCf(i4)*PressureC(i1)+
     *                              FluxFf(i4)*PressureC(j1)
c
            BuVelocity(j2,j3)=a1r(i2)*BuVelocity(i2,i3)+
     *           b1r(i2)* BvVelocity(i2,i3)+c1r(i2)* BwVelocity(i2,i3)
            BvVelocity(j2,j3)=a2r(i2)*BuVelocity(i2,i3)+
     *           b2r(i2)* BvVelocity(i2,i3)+c2r(i2)* BwVelocity(i2,i3)
            BwVelocity(j2,j3)=a3r(i2)*BuVelocity(i2,i3)+
     *           b3r(i2)* BvVelocity(i2,i3)+c3r(i2)* BwVelocity(i2,i3)
c
            Bmdot(j2,j3)=-Bmdot(i2,i3)
c
          endif
c
        enddo
c----------------------------------------------------------------------
      endif
c----------------------------------------------------------------------
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CorrectVelocity
c
c#############################################################################################
c
      use User0, only: MethodCalcGradientContinuity,
     *                 nIterGradientContinuity,LimitGradientContinuity,
     *                 LimitGradientContinuityMethod
      use Geometry1, only: NumberOfElements
      use Variables1, only: PressureC,PressCGradx,PressCGrady,
     *                      PressCGradz,BPressureC,
     *                      BPressCGradx,BPressCGrady,BPressCGradz,
     *                      uVelocity,vVelocity,wVelocity,
     *                      Du2Velocity,Dv2Velocity,Dw2Velocity
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------
          character*10 Variable
          integer :: MethodCalcGradient,nIterGradientPhi
          logical :: LimitGradient
          integer :: LimitGradientMethod
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c       
      Variable='pressc'
      call Gradient(Variable,MethodCalcGradientContinuity,
     *         PressureC,PressCGradx,PressCGrady,PressCGradz,
     *         BPressureC,BPressCGradx,BPressCGrady,BPressCGradz,
     *         nIterGradientContinuity,LimitGradientContinuity,
     *         LimitGradientContinuityMethod)


      do i=1,NumberOfElements           
c
        uVelocity(i)=uVelocity(i)-Du2Velocity(i)*PressCGradx(i) 
        vVelocity(i)=vVelocity(i)-Dv2Velocity(i)*PressCGrady(i) 
        wVelocity(i)=wVelocity(i)-Dw2Velocity(i)*PressCGradz(i) 
c      
      enddo      
c
      return
      end