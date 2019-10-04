c
c#############################################################################################
c
      SUBROUTINE AssemblePressureCorrection
c
c#############################################################################################
c
      use User0, only: MethodDecomposeS,LUnsteady,dt,Lcompressible
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces,NIFaceOwner,NIFaceNeighbor,
     *                     NBFaceOwner
      use Geometry4, only: FaceAreax,FaceAreay,FaceAreaz,
     *                     Volume,xc,yc,zc,
     *                     DistanceCFux,DistanceCFuy,DistanceCFuz,
     *                     DistanceCF,BFaceAreax,BFaceAreay,BFaceAreaz,
     *                     BDistanceCFux,BDistanceCFuy,BDistanceCFuz,
     *                     BDistanceCF,BFaceArea,
     *                     BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
      use Variables1, only: Du2Velocity,Dv2Velocity,Dw2Velocity,mdot,
     *                      Bmdot,BuVelocity,BvVelocity,BwVelocity,
     *                      BTemperature,BPressure,uVelocity,vVelocity,
     *                      wVelocity,Temperature,Pressure
      use Variables2, only: FluxCf,FluxFf,FluxVf,FluxTf
      use Variables3, only: FluxCE,FluxTE
      use BoundaryConditions2
      use PhysicalProperties1, only:Density,DensityOld,DensityOldOld,
     *                              drhodp,Bdrhodp,Densityf,BDensity,
     *                              BSpecificHeat,RGas
      use Constants1, only: tiny
      use ArteryResistance1, only: ArteryResistance,geoDiffB,geoDiffSum
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j,k,iBFace,j1,j2,j3
      Character*16 Interpolation
      double precision :: rhof,sfx,sfy,sfz,ex,ey,ez,cf,DuSf,DvSf,DwSf,
     *                    DuEf,DvEf,DwEf,DuTf,DvTf,DwTf,geoDiff,
     *                    FluxCfLocal,FluxFfLocal,Magnitude,eDuSf,
     *                    eDvSf,eDwSf,TotalTemperature,dpdub,c1,diff,
     *                    velocity,velocity2
c
      double precision :: nx,ny,nz,vdotInfinity,vdotin,cinfinity,cin,
     *                    vdotb,cb,Tb,rhoInfinity,sbEntropy,rhob,
     *                    vel2,xdir,ydir,zdir,FluxCfLocal1,ratio,uF1,
     *                    vF1,Ub1,Ub2,term1,term2,term3,term4
      double precision :: distance1,distance2,GFactCF,
     *                    DistCFx,DistCFy,DistCFz,xF1,yF1,zF1,DuF1,DvF1,
     *                    DwF1,DistCF,DistCFux,DistCFuy,DistCFuz,
     *                    dotproduct
c
      double precision, save, dimension(:), allocatable :: Du2f
      double precision, save, dimension(:), allocatable :: Dv2f
      double precision, save, dimension(:), allocatable :: Dw2f
      double precision, save, dimension(:), allocatable :: drhodpf
      double precision, save, dimension(:), allocatable :: drhodt
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE ComputeTransientGradient(phi,phiOld,phiOldOld,dphidt)
c--------------------------------------------------------------
          double precision, dimension(:) :: phi,phiOld,phiOldOld,dphidt
c--------------------------------------------------------------
        end SUBROUTINE ComputeTransientGradient
c--------------------------------------------------------------
       SUBROUTINE InterpolateElementToFace(InterpolationScheme,FiT,FiTf)
c--------------------------------------------------------------
          character*16 InterpolationScheme
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: FiTf
        end SUBROUTINE InterpolateElementToFace
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      allocate(Du2f(NIFaces))
      allocate(Dv2f(NIFaces))
      allocate(Dw2f(NIFaces))

      interpolation='average'
      call InterpolateElementToFace(interpolation,Du2Velocity,Du2f)
      call InterpolateElementToFace(interpolation,Dv2Velocity,Dv2f)
      call InterpolateElementToFace(interpolation,Dw2Velocity,Dw2f)

      do k=1,NIFaces
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        rhof=Densityf(k)
        sfx=FaceAreax(k)
        sfy=FaceAreay(k)
        sfz=FaceAreaz(k)
        ex=DistanceCFux(k)
        ey=DistanceCFuy(k)
        ez=DistanceCFuz(k)
        cf=DistanceCF(k)
c
c--- difference in pressure gradients term
c
        if(MethodDecomposeS.eq.1) then
c
          DuSf=Du2f(k)*sfx
          DvSf=Dv2f(k)*sfy
          DwSf=Dw2f(k)*sfz
c
          dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
          DuEf=dotproduct*ex
          DvEf=dotproduct*ey
          DwEf=dotproduct*ez
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
        elseif(MethodDecomposeS.eq.2) then
c
          DuSf=Du2f(k)*sfx
          DvSf=Dv2f(k)*sfy
          DwSf=Dw2f(k)*sfz
c          
          Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
          DuEf=ex*Magnitude
          DvEf=ey*Magnitude
          DwEf=ez*Magnitude
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
        elseif(MethodDecomposeS.eq.3) then
c
          DuSf=Du2f(k)*sfx
          DvSf=Dv2f(k)*sfy
          DwSf=Dw2f(k)*sfz
          Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
          eDuSf=DuSf/Magnitude
          eDvSf=DvSf/Magnitude
          eDwSf=DwSf/Magnitude
c
          dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
          DuEf=Magnitude*ex/dotproduct
          DvEf=Magnitude*ey/dotproduct
          DwEf=Magnitude*ez/dotproduct
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
        endif
c
        FluxCfLocal=rhof*geoDiff
        FluxFfLocal=-rhof*geoDiff
c
        FluxCf(k)= FluxCfLocal
        FluxFf(k)= FluxFfLocal
        FluxTf(k)= mdot(k)
c
      enddo
c
      iBFace=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          iBFace=iBFace+1
c
          FluxTf(iBFace)= Bmdot(i,j)
c
        enddo
      enddo
c
      if(Lcompressible) then
c
        allocate(drhodpf(NIFaces))
        call InterpolateElementToFace(interpolation,drhodp,drhodpf)
c
        do k=1,NIFaces
c        
          FluxCflocal=(drhodpf(k)/Densityf(k))*dmax1(mdot(k),0.)
          FluxFflocal=-(drhodpf(k)/Densityf(k))*dmax1(-mdot(k),0.)
c
          FluxCf(k)=FluxCf(k)+FluxCflocal
          FluxFf(k)=FluxFf(k)+FluxFflocal
c
        enddo
c
        iBFace=NIFaces
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            iBFace=iBFace+1
c
            FluxCflocal=
     *          (Bdrhodp(i,j)/BDensity(i,j))*dmax1(Bmdot(i,j),0.)
c
            FluxCf(iBFace)=FluxCf(iBFace)+FluxCflocal
c
          enddo
        enddo
c
        if(Lunsteady) then
c
          allocate(drhodt(NumberOfElements))
          call ComputeTransientGradient
     *            (Density,DensityOld,DensityOldOld,drhodt)
c
          do i=1,NumberOfElements
c
            FluxCE(i)=FluxCE(i)+drhodp(i)*Volume(i)/dt
            FluxTE(i)=FluxTE(i)+drhodt(i)*Volume(i)

c
          enddo
c
          deallocate(drhodt)
c      
        endif
c
        deallocate(drhodpf)
c      
      endif
c
c----------------------------------------------------------------------
c--- Modify along inlet and oulet boundaries for incompressible and compressible flows
c----------------------------------------------------------------------
c
      if(.not.Lcompressible) then
c
c----------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i1=IoutletspecifiedVelocityOwner(i)
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxCfLocal
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxCfLocal
c
        enddo
c----------------------------------------------------------------------
        geoDiffSum=0.
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          geoDiffB(i2,i3)=geoDiff
          geoDiffSum(i2)=geoDiffSum(i2)+geoDiff
c
        enddo
c
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
c
          FluxCfLocal=rhof*geoDiffB(i2,i3)
          FluxFfLocal=FluxCfLocal
c
          FluxCfLocal=FluxCfLocal*
     *      (1+(geoDiffSum(i2)-geoDiffB(i2,i3))*ArteryResistance(i2))/
     *          (1.+ArteryResistance(i2)*geoDiffSum(i2))
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxFfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletFullyDeveloped
c
          i1=IoutletFullyDevelopedOwner(i)
          i2=IoutletFullyDevelopedNumberOfBCSets(i)
          i3=IoutletFullyDevelopedNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i1=IinletSpecifiedVelocityOwner(i)
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxCfLocal
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=FluxCfLocal
c
          velocity2=BuVelocity(i2,i3)**2+
     *                BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2
c
          FluxCfLocal=Bmdot(i2,i3)*FluxCfLocal/
     *          (Bmdot(i2,i3)-FluxCfLocal*BDensity(i2,i3)*velocity2)
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxFfLocal
c
        enddo
c----------------------------------------------------------------------
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          if(LRotationalPeriodicity) then
c
            xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
            yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
            zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
          elseif(LTranslationalPeriodicity) then
c
            xF1=xc(j1)+xTranslation(j2)
            yF1=yc(j1)+yTranslation(j2)
            zF1=zc(j1)+zTranslation(j2)
c
          endif
c
          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
c
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          ex=DistCFux
          ey=DistCFuy
          ez=DistCFuz
          cf=DistCF
c
          DuF1=a1r(j2)*Du2Velocity(j1)+b1r(j2)*Dv2Velocity(j1)+
     *                                        c1r(j2)*Dw2Velocity(j1)
          DvF1=a2r(j2)*Du2Velocity(j1)+b2r(j2)*Dv2Velocity(j1)+
     *                                        c2r(j2)*Dw2Velocity(j1)
          DwF1=a3r(j2)*Du2Velocity(j1)+b3r(j2)*Dv2Velocity(j1)+
     *                                        c3r(j2)*Dw2Velocity(j1)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*DuF1
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*DvF1
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*DwF1
            DwSf=DwSf*sfz
c
            DuEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ex
            DvEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ey
            DwEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*DuF1
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*DvF1
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*DwF1
            DwSf=DwSf*sfz
c          
            DuEf=ex*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DvEf=ey*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DwEf=ez*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*Du2Velocity(j1)
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*Dv2Velocity(j1)
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*Dw2Velocity(j1)
            DwSf=DwSf*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            DuEf=Magnitude*ex/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
            DvEf=Magnitude*ey/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
            DwEf=Magnitude*ez/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
          endif
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxFfLocal
c
        enddo
c
c----------------------------------------------------------------------
c
      elseif(Lcompressible) then
c
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
          FluxCf(i4)=FluxCf(i4)
          FluxFf(i4)=FluxCf(i4)
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletspecifiedVelocity
c
          i1=IoutletspecifiedVelocityOwner(i)
          i2=IoutletspecifiedVelocityNumberOfBCSets(i)
          i3=IoutletspecifiedVelocityNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4)
          FluxFf(i4)=0.
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxCfLocal
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxCfLocal
c
        enddo
c----------------------------------------------------------------------
        geoDiffSum=0.
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          geoDiffB(i2,i3)=geoDiff
          geoDiffSum(i2)=geoDiffSum(i2)+geoDiff
c
        enddo
c
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
c
          FluxCfLocal=rhof*geoDiffB(i2,i3)
          FluxFfLocal=FluxCfLocal
c
          FluxCfLocal=FluxCfLocal*
     *      (1+(geoDiffSum(i2)-geoDiffB(i2,i3))*ArteryResistance(i2))/
     *          (1.+ArteryResistance(i2)*geoDiffSum(i2))
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxFfLocal
c
        enddo
c----------------------------------------------------------------------
        do i=1,IoutletSpecifiedMassFlowRate
c
          i1=IoutletSpecifiedMassFlowRateOwner(i)
          i2=IoutletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IoutletSpecifiedMassFlowRateNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
c
        enddo
c----------------------------------------------------------------------
c----------------------------------------------------------------------
        do i=1,Iinletsupersonic
c
          i1=IinletsupersonicOwner(i)
          i2=IinletsupersonicNumberOfBCSets(i)
          i3=IinletsupersonicNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedVelocity
c
          i1=IinletSpecifiedVelocityOwner(i)
          i2=IinletSpecifiedVelocityNumberOfBCSets(i)
          i3=IinletSpecifiedVelocityNBFaces(i)
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
          FluxCf(i4)=FluxCf(i4) 
          FluxFf(i4)=0.
c
        enddo
c----------------------------------------------------------------------
        do i=1,IinletSpecifiedMassFlowRate
c
          i1=IinletSpecifiedMassFlowRateOwner(i)
          i2=IinletSpecifiedMassFlowRateNumberOfBCSets(i)
          i3=IinletSpecifiedMassFlowRateNBFaces(i)
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
          FluxCf(i4)=0.
          FluxFf(i4)=0.
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxCfLocal
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
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal=rhof*geoDiff
c
          velocity2=BuVelocity(i2,i3)**2+
     *                BvVelocity(i2,i3)**2+BwVelocity(i2,i3)**2
          TotalTemperature=BTemperature(i2,i3)+
     *                     velocity2/(2.*BSpecificHeat(i2,i3))
          dpdub=-BPressure(i2,i3)*Bmdot(i2,i3)/
     *                       (rhof*TotalTemperature*Rgas)
          c1=dpdub*FluxCfLocal/(rhof+dpdub*FluxCfLocal)
c
          diff=Bdrhodp(i2,i3)*Bmdot(i2,i3)/BDensity(i2,i3)
c
          FluxCf(i4)=FluxCfLocal*(1.-c1)+diff*c1
          FluxFf(i4)=FluxCfLocal
c
        enddo
c
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
c          FluxCfLocal1=rhof*geoDiff/BFaceArea(i2,i3)
c
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
          ratio=BSpecificHeat(i2,i3)/(BSpecificHeat(i2,i3)-RGas)
c
          vdotInfinity=uVelocityFarField(i2)*nx+
     *             vVelocityFarField(i2)*ny+wVelocityFarField(i2)*nz
          vdotin=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
          cinfinity=dsqrt(ratio*RGas*TemperatureFarField(i2))
          cin=dsqrt(ratio*RGas*Temperature(i1))
c
          vdotb=0.5*(vdotInfinity+vdotin)-(cinfinity-cin)/(ratio-1.)
          cb=-(vdotInfinity-vdotin)*(ratio-1.)/4.+0.5*(cinfinity+cin)
          Tb=(cb*cb)/(ratio*RGas)
          rhoInfinity=pressureFarField(i2)/
     *                            (RGas*TemperatureFarField(i2))
c
          if(Bmdot(i2,i3).lt.0.) then
c
            sbEntropy=cinfinity**2/(ratio*(rhoInfinity**(ratio-1)))
c
          else
c
            sbEntropy=cin**2/(ratio*(Density(i1)**(ratio-1)))
c
          endif           
c
          rhob=(cb*cb/(ratio*sbEntropy))**(1./(ratio-1.))
c
          FluxCfLocal=rhob*(1.+vdotb/cb)*BFaceArea(i2,i3)*
     *     geoDiff/(BFaceArea(i2,i3)+geoDiff*BPressure(i2,i3)/cb)
c          FluxCfLocal=rhob*(1.+vdotb/cb)*BFaceArea(i2,i3)*
c     *     FluxCfLocal1/(rhob+FluxCfLocal1*BPressure(i2,i3)/cb)
c
c          FluxCfLocal=rhob*(1.+vdotb/cb)*BFaceArea(i2,i3)*
c     *         (FluxCfLocal1/2.+cin/(2.*ratio*Pressure(i1)))/
c     *         (1.+0.5*FluxCfLocal1*ratio*BPressure(i2,i3)/cb)
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=rhof*geoDiff
c
        enddo
c
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
          ex=BDistanceCFux(i2,i3)
          ey=BDistanceCFuy(i2,i3)
          ez=BDistanceCFuz(i2,i3)
          cf=BDistanceCF(i2,i3)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c
            dotproduct=ex*DuSf+ey*DvSf+ez*DwSf
            DuEf=dotproduct*ex
            DvEf=dotproduct*ey
            DwEf=dotproduct*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf   
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
c          
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DuEf=ex*Magnitude
            DvEf=ey*Magnitude
            DwEf=ez*Magnitude
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf  
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=Du2Velocity(i1)*sfx
            DvSf=Dv2Velocity(i1)*sfy
            DwSf=Dw2Velocity(i1)*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            dotproduct=ex*eDuSf+ey*eDvSf+ez*eDwSf
            DuEf=Magnitude*ex/dotproduct
            DvEf=Magnitude*ey/dotproduct
            DwEf=Magnitude*ez/dotproduct
c          =
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf
c
          endif
c
          FluxCfLocal1=rhof*geoDiff/BFaceArea(i2,i3)
c
          nx=BFaceAreanx(i2,i3)
          ny=BFaceAreany(i2,i3)
          nz=BFaceAreanz(i2,i3)
          ratio=BSpecificHeat(i2,i3)/(BSpecificHeat(i2,i3)-RGas)
c
          vdotin=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
          cin=dsqrt(ratio*RGas*Temperature(i1))
c
          vdotb=vdotin
          cb=cin
          Tb=(cb*cb)/(ratio*RGas)
c
          sbEntropy=cin**2/(ratio*(Density(i1)**(ratio-1)))
c
          rhob=(cb*cb/(ratio*sbEntropy))**(1./(ratio-1.))
c
          FluxCfLocal=rhob*(1.+vdotb/cb)*BFaceArea(i2,i3)*
     *         FluxCfLocal1/(1.+FluxCfLocal1*BPressure(i2,i3)/cb)
c          FluxCfLocal=rhob*FluxCfLocal1*(1.+vdotb/cb)*BFaceArea(i2,i3)/
c     *                         (1.+FluxCfLocal1*cb/ratio)
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=rhof*geoDiff
c
        enddo
c
c----------------------------------------------------------------------
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
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
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
c
          if(LRotationalPeriodicity) then
c
            xF1=a1r(j2)*xc(j1)+b1r(j2)*yc(j1)+c1r(j2)*zc(j1)
            yF1=a2r(j2)*xc(j1)+b2r(j2)*yc(j1)+c2r(j2)*zc(j1)
            zF1=a3r(j2)*xc(j1)+b3r(j2)*yc(j1)+c3r(j2)*zc(j1)
c
          elseif(LTranslationalPeriodicity) then
c
            xF1=xc(j1)+xTranslation(j2)
            yF1=yc(j1)+yTranslation(j2)
            zF1=zc(j1)+zTranslation(j2)
c
          endif
c
          distance1=dsqrt((BFaceCentroidx(i2,i3)-xc(i1))**2+
     *                         (BFaceCentroidy(i2,i3)-yc(i1))**2+
     *                            (BFaceCentroidz(i2,i3)-zc(i1))**2)
          distance2=dsqrt((BFaceCentroidx(i2,i3)-xF1)**2+
     *                         (BFaceCentroidy(i2,i3)-yF1)**2+
     *                            (BFaceCentroidz(i2,i3)-zF1)**2)
c
          GFactCF=distance2/(distance1+distance2)
c
          DistCFx=xF1-xc(i1)
          DistCFy=yF1-yc(i1)
          DistCFz=zF1-zc(i1)
c
          DistCF=dsqrt(DistCFx**2+DistCFy**2+DistCFz**2)
c
          DistCFux=DistCFx/DistCF
          DistCFuy=DistCFy/DistCF
          DistCFuz=DistCFz/DistCF
c
          ex=DistCFux
          ey=DistCFuy
          ez=DistCFuz
          cf=DistCF
c
          DuF1=a1r(j2)*Du2Velocity(j1)+b1r(j2)*Dv2Velocity(j1)+
     *                                        c1r(j2)*Dw2Velocity(j1)
          DvF1=a2r(j2)*Du2Velocity(j1)+b2r(j2)*Dv2Velocity(j1)+
     *                                        c2r(j2)*Dw2Velocity(j1)
          DwF1=a3r(j2)*Du2Velocity(j1)+b3r(j2)*Dv2Velocity(j1)+
     *                                        c3r(j2)*Dw2Velocity(j1)
c
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*DuF1
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*DvF1
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*DwF1
            DwSf=DwSf*sfz
c
            DuEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ex
            DvEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ey
            DwEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ez
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
          elseif(MethodDecomposeS.eq.2) then
c
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*DuF1
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*DvF1
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*DwF1
            DwSf=DwSf*sfz
c          
            DuEf=ex*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DvEf=ey*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
            DwEf=ez*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
          elseif(MethodDecomposeS.eq.3) then
c
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*Du2Velocity(j1)
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*Dv2Velocity(j1)
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*Dw2Velocity(j1)
            DwSf=DwSf*sfz
            Magnitude=dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)    
            eDuSf=DuSf/Magnitude
            eDvSf=DvSf/Magnitude
            eDwSf=DwSf/Magnitude
c
            DuEf=Magnitude*ex/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
            DvEf=Magnitude*ey/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
            DwEf=Magnitude*ez/(ex*eDuSf+ey*eDvSf+ez*eDwSf)
c          
            geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
          endif
c
          FluxCfLocal=rhof*geoDiff
          FluxFfLocal=-rhof*geoDiff
c
          FluxCf(i4)=FluxCfLocal
          FluxFf(i4)=FluxFfLocal
c
        enddo
c
      endif
c
      deallocate(Du2f)
      deallocate(Dv2f)
      deallocate(Dw2f)
c
      return
      end