c
C#############################################################################################
c
      SUBROUTINE AssembleViscousWorkTerm
c
C#############################################################################################
      use User0, only: LTurbulentFlow,LCompressible,
     *                 LSolveTurbulenceKineticEnergy
      use Variables1, only: uVelGradfx,uVelGradfy,uVelGradfz,
     *                      vVelGradfx,vVelGradfy,vVelGradfz,
     *                      wVelGradfx,wVelGradfy,wVelGradfz,
     *                      uVelGradx,uVelGrady,uVelGradz,
     *                      vVelGradx,vVelGrady,vVelGradz,
     *                      wVelGradx,wVelGrady,wVelGradz,
     *                      BuVelGradx,BuVelGrady,BuVelGradz,
     *                      BvVelGradx,BvVelGrady,BvVelGradz,
     *                      BwVelGradx,BwVelGrady,BwVelGradz,
     *                      uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      TurbulentKE,BTurbulentKE
      use Variables3, only: FluxTE
      use Geometry1, only: NumberOfElements
      use Geometry3, only: NIFaceOwner,NIFaceNeighbor,NIFaces,NBFaces,
     *                     NBFaceOwner
      use Geometry4, only: FaceAreax,FaceAreay,FaceAreaz,GFactorCF,
     *                     BFaceAreax,BFaceAreay,BFaceAreaz,xc,yc,zc,
     *                     BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
      use PhysicalProperties1, only:Viscosity,TurbulentViscosity,
     *                              Density,Densityf,BViscosity,
     *                              BTurbulentViscosity,BDensity
      use Constants1, only: twothird
      use BoundaryConditions2
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k,i1,i2,i3,i4,j1,j2,j3,j4,k1
      double precision :: FluxTElocal,viscf,tkef,nx,ny,nz,dNorm,
     *                    Tau11,Tau12,Tau13,Tau22,Tau23,Tau33,
     *                    uVelf,vVelf,wVelf,term1,term2,term3,
     *                    uWall1,vWall1,wWall1,dotproduct,
     *                    uVel1,vVel1,wVel1,xF1,yF1,zF1,
     *                    distance1,distance2,GFactCF,densf
c********************************************************************************************
      interface
c********************************************************************************************
       SUBROUTINE InterpolateElementToFace(InterpolationScheme,FiT,FiTf)
c********************************************************************************************
          character*16 InterpolationScheme
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: FiTf
        end SUBROUTINE InterpolateElementToFace
c********************************************************************************************
      endinterface
c********************************************************************************************************
c
c--- Interior faces
c
      do k=1,NIFaces
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        if(LTurbulentFlow) then
c
          viscf=(TurbulentViscosity(i)+Viscosity(i))*GFactorCF(k)+
     *           (TurbulentViscosity(j)+Viscosity(j))*(1.-GFactorCF(k))
c
        else
c
          viscf=Viscosity(i)*GFactorCF(k)+Viscosity(j)*(1.-GFactorCF(k))
c
        endif
c
        Tau11=2.*viscf*uVelGradfx(k)
        Tau12=viscf*(uVelGradfy(k)+vVelGradfx(k))
        Tau13=viscf*(uVelGradfz(k)+wVelGradfx(k))
        Tau22=2.*viscf*vVelGradfy(k)
        Tau23=viscf*(vVelGradfz(k)+wVelGradfy(k))
        Tau33=2.*viscf*wVelGradfz(k)
c
        if(Lcompressible) then
c
          term1=uVelGradfx(k)+vVelGradfy(k)+wVelGradfz(k)
c
          Tau11=Tau11-twothird*viscf*term1
          Tau22=Tau22-twothird*viscf*term1
          Tau33=Tau33-twothird*viscf*term1
c          
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          tkef=dmax1(TurbulentKE(i)*GFactorCF(k)+
     *                    TurbulentKE(j)*(1.-GFactorCF(k)),0.)
c
          Tau11=Tau11-twothird*Densityf(k)*tkef
          Tau22=Tau22-twothird*Densityf(k)*tkef
          Tau33=Tau33-twothird*Densityf(k)*tkef
c
        endif
c
        uVelf=uVelocity(i)*GFactorCF(k)+uVelocity(j)*(1.-GFactorCF(k))
        vVelf=vVelocity(i)*GFactorCF(k)+vVelocity(j)*(1.-GFactorCF(k))
        wVelf=wVelocity(i)*GFactorCF(k)+wVelocity(j)*(1.-GFactorCF(k))
c
        term1=(Tau11*uVelf+Tau12*vVelf+Tau13*wVelf)*FaceAreax(k)
        term2=(Tau12*uVelf+Tau22*vVelf+Tau23*wVelf)*FaceAreay(k)
        term3=(Tau13*uVelf+Tau23*vVelf+Tau33*wVelf)*FaceAreaz(k)
c
        FluxTElocal=term1+term2+term3
c
        FluxTE(i)= FluxTE(i)-FluxTElocal
        FluxTE(j)= FluxTE(j)+FluxTElocal
c
      enddo
c----------------------------------------------------------------------
c
c--- Boundary faces
c
c----------------------------------------------------------------------
      do i=1,IWallnoSlip
c
        i1=IWallnoSlipOwner(i)
        i2=IWallnoSlipNumberOfBCSets(i)
        i3=IWallnoSlipNBFaces(i)
c
        viscf=BViscosity(i2,i3)
        if(LTurbulentFlow) viscf=BTurbulentViscosity(i2,i3)
c
        nx=BFaceAreanx(i2,i3)
        ny=BFaceAreany(i2,i3)
        nz=BFaceAreanz(i2,i3)
        dNorm=BDistanceCFx(i2,i3)*nx+
     *               BDistanceCFy(i2,i3)*ny+BDistanceCFz(i2,i3)*nz
c
        uWall1=BuVelocity(i2,i3)-uVelocity(i1)
        vWall1=BvVelocity(i2,i3)-vVelocity(i1)
        wWall1=BwVelocity(i2,i3)-wVelocity(i1)
c
        dotproduct=uWall1*nx+vWall1*ny+wWall1*nz
c
        uWall1=uWall1-dotproduct*nx
        vWall1=vWall1-dotproduct*ny
        wWall1=wWall1-dotproduct*nz
c
        Tau11=-(Bviscosity(i2,i3)/dNorm)*uWall1
        Tau12=-(Bviscosity(i2,i3)/dNorm)*vWall1
        Tau13=-(Bviscosity(i2,i3)/dNorm)*wWall1
c          
        term1=BuVelocity(i2,i3)*Tau11*BFaceAreax(i2,i3)
        term2=BvVelocity(i2,i3)*Tau12*BFaceAreay(i2,i3)
        term3=BwVelocity(i2,i3)*Tau13*BFaceAreaz(i2,i3)
c
        FluxTElocal=term1+term2+term3
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo
c----------------------------------------------------------------------
      do i=1,Iinletsupersonic
c
        i1=IinletsupersonicOwner(i)
        i2=IinletsupersonicNumberOfBCSets(i)
        i3=IinletsupersonicNBFaces(i)
c
        viscf=BViscosity(i2,i3)
        if(LTurbulentFlow) viscf=BTurbulentViscosity(i2,i3)
c
        Tau11=2.*viscf*BuVelGradx(i2,i3)
        Tau12=viscf*(BuVelGrady(i2,i3)+BvVelGradx(i2,i3))
        Tau13=viscf*(BuVelGradz(i2,i3)+BwVelGradx(i2,i3))
        Tau22=2.*viscf*BvVelGrady(i2,i3)
        Tau23=viscf*(BvVelGradz(i2,i3)+BwVelGrady(i2,i3))
        Tau33=2.*viscf*BwVelGradz(i2,i3)
c
        if(Lcompressible) then
c
          term1=BuVelGradx(i2,i3)+BvVelGrady(i2,i3)+BwVelGradz(i2,i3)
c
          Tau11=Tau11-twothird*viscf*term1
          Tau22=Tau22-twothird*viscf*term1
          Tau33=Tau33-twothird*viscf*term1
c          
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          tkef=dmax1(BTurbulentKE(i2,i3),0.)
c
          Tau11=Tau11-twothird*BDensity(i2,i3)*tkef
          Tau22=Tau22-twothird*BDensity(i2,i3)*tkef
          Tau33=Tau33-twothird*BDensity(i2,i3)*tkef
c
        endif
c
        uVelf=BuVelocity(i2,i3)
        vVelf=BvVelocity(i2,i3)
        wVelf=BwVelocity(i2,i3)
c
        term1=(Tau11*uVelf+Tau12*vVelf+Tau13*wVelf)*BFaceAreax(i2,i3)
        term2=(Tau12*uVelf+Tau22*vVelf+Tau23*wVelf)*BFaceAreay(i2,i3)
        term3=(Tau13*uVelf+Tau23*vVelf+Tau33*wVelf)*BFaceAreaz(i2,i3)
c
        FluxTElocal=term1+term2+term3
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo
c----------------------------------------------------------------------
      do i=1,IinletSpecifiedStaticTemperature
c
        i1=IinletSpecifiedStaticTemperatureOwner(i)
        i2=IinletSpecifiedStaticTemperatureNumberOfBCSets(i)
        i3=IinletSpecifiedStaticTemperatureNBFaces(i)
c
        viscf=BViscosity(i2,i3)
        if(LTurbulentFlow) viscf=BTurbulentViscosity(i2,i3)
c
        Tau11=2.*viscf*BuVelGradx(i2,i3)
        Tau12=viscf*(BuVelGrady(i2,i3)+BvVelGradx(i2,i3))
        Tau13=viscf*(BuVelGradz(i2,i3)+BwVelGradx(i2,i3))
        Tau22=2.*viscf*BvVelGrady(i2,i3)
        Tau23=viscf*(BvVelGradz(i2,i3)+BwVelGrady(i2,i3))
        Tau33=2.*viscf*BwVelGradz(i2,i3)
c
        if(Lcompressible) then
c
          term1=BuVelGradx(i2,i3)+BvVelGrady(i2,i3)+BwVelGradz(i2,i3)
c
          Tau11=Tau11-twothird*viscf*term1
          Tau22=Tau22-twothird*viscf*term1
          Tau33=Tau33-twothird*viscf*term1
c          
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          tkef=dmax1(BTurbulentKE(i2,i3),0.)
c
          Tau11=Tau11-twothird*BDensity(i2,i3)*tkef
          Tau22=Tau22-twothird*BDensity(i2,i3)*tkef
          Tau33=Tau33-twothird*BDensity(i2,i3)*tkef
c
        endif
c
        uVelf=BuVelocity(i2,i3)
        vVelf=BvVelocity(i2,i3)
        wVelf=BwVelocity(i2,i3)
c
        term1=(Tau11*uVelf+Tau12*vVelf+Tau13*wVelf)*BFaceAreax(i2,i3)
        term2=(Tau12*uVelf+Tau22*vVelf+Tau23*wVelf)*BFaceAreay(i2,i3)
        term3=(Tau13*uVelf+Tau23*vVelf+Tau33*wVelf)*BFaceAreaz(i2,i3)
c
        FluxTElocal=term1+term2+term3
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo
c----------------------------------------------------------------------
      do i=1,IinletSpecifiedStagnationTemperature
c
        i1=IinletSpecifiedStagnationTemperatureOwner(i)
        i2=IinletSpecifiedStagnationTemperatureNumberOfBCSets(i)
        i3=IinletSpecifiedStagnationTemperatureNBFaces(i)
c
        viscf=BViscosity(i2,i3)
        if(LTurbulentFlow) viscf=BTurbulentViscosity(i2,i3)
c
        Tau11=2.*viscf*BuVelGradx(i2,i3)
        Tau12=viscf*(BuVelGrady(i2,i3)+BvVelGradx(i2,i3))
        Tau13=viscf*(BuVelGradz(i2,i3)+BwVelGradx(i2,i3))
        Tau22=2.*viscf*BvVelGrady(i2,i3)
        Tau23=viscf*(BvVelGradz(i2,i3)+BwVelGrady(i2,i3))
        Tau33=2.*viscf*BwVelGradz(i2,i3)
c
        if(Lcompressible) then
c
          term1=BuVelGradx(i2,i3)+BvVelGrady(i2,i3)+BwVelGradz(i2,i3)
c
          Tau11=Tau11-twothird*viscf*term1
          Tau22=Tau22-twothird*viscf*term1
          Tau33=Tau33-twothird*viscf*term1
c          
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          tkef=dmax1(BTurbulentKE(i2,i3),0.)
c
          Tau11=Tau11-twothird*BDensity(i2,i3)*tkef
          Tau22=Tau22-twothird*BDensity(i2,i3)*tkef
          Tau33=Tau33-twothird*BDensity(i2,i3)*tkef
c
        endif
c
        uVelf=BuVelocity(i2,i3)
        vVelf=BvVelocity(i2,i3)
        wVelf=BwVelocity(i2,i3)
c
        term1=(Tau11*uVelf+Tau12*vVelf+Tau13*wVelf)*BFaceAreax(i2,i3)
        term2=(Tau12*uVelf+Tau22*vVelf+Tau23*wVelf)*BFaceAreay(i2,i3)
        term3=(Tau13*uVelf+Tau23*vVelf+Tau33*wVelf)*BFaceAreaz(i2,i3)
c
        FluxTElocal=term1+term2+term3
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo

c----------------------------------------------------------------------
      do i=1,Ioutletsupersonic
c
        i1=IoutletsupersonicOwner(i)
        i2=IoutletsupersonicNumberOfBCSets(i)
        i3=IoutletsupersonicNBFaces(i)
c
        FluxTElocal=0.
c
        FluxTE(i1)= FluxTE(i1)+FluxTElocal
c
      enddo

c----------------------------------------------------------------------
      do i=1,IoutletFullyDevelopedEnergy
c
        i1=IoutletFullyDevelopedEnergyOwner(i)
        i2=IoutletFullyDevelopedEnergyNumberOfBCSets(i)
        i3=IoutletFullyDevelopedEnergyNBFaces(i)
c
        FluxTElocal=0.
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo
c----------------------------------------------------------------------
      do i=1,Isymmetry
c
        i1=IsymmetryOwner(i)
        i2=IsymmetryNumberOfBCSets(i)
        i3=IsymmetryNBFaces(i)
c
        viscf=BViscosity(i2,i3)
        if(LTurbulentFlow) viscf=BTurbulentViscosity(i2,i3)
c
        dNorm=BDistanceCFx(i2,i3)*BFaceAreanx(i2,i3)+
     *                BDistanceCFy(i2,i3)*BFaceAreany(i2,i3)+
     *                BDistanceCFz(i2,i3)*BFaceAreanz(i2,i3)
c
        nx=BFaceAreanx(i2,i3)
        ny=BFaceAreany(i2,i3)
        nz=BFaceAreanz(i2,i3)
c
        dotproduct=uVelocity(i1)*nx+vVelocity(i1)*ny+wVelocity(i1)*nz
        uVel1=dotproduct*nx
        vVel1=dotproduct*ny
        wVel1=dotproduct*nz
c
        Tau11=-2.*viscf*uVel1/dNorm        
        Tau12=-2.*viscf*vVel1/dNorm        
        Tau13=-2.*viscf*wVel1/dNorm        
c
        uVelf=BuVelocity(i2,i3)
        vVelf=BvVelocity(i2,i3)
        wVelf=BwVelocity(i2,i3)
c
        term1=Tau11*uVelf*BFaceAreax(i2,i3)
        term2=Tau12*vVelf*BFaceAreay(i2,i3)
        term3=Tau13*wVelf*BFaceAreaz(i2,i3)
c
        FluxTElocal=term1+term2+term3
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo
c----------------------------------------------------------------------
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
        do k1=1,j2-1
c
          j4=j4+NBFaces(k1)
c
        enddo
c
        j4=j4+j3
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
     *                              (BFaceCentroidz(i2,i3)-zF1)**2)
c
        GFactCF=distance2/(distance1+distance2)
c
        if(LTurbulentFlow) then
c
          viscf=(TurbulentViscosity(i1)+Viscosity(i1))*GFactCF+
     *           (TurbulentViscosity(j1)+Viscosity(j1))*(1.-GFactCF)
c
        else
c
          viscf=Viscosity(i1)*GFactCF+Viscosity(j1)*(1.-GFactCF)
c
        endif
c
        Tau11=2.*viscf*(uVelGradx(i1)*GFactCF+
     *                           uVelGradx(j1)*(1.-GFactCF))
        Tau12=viscf*(
     *        uVelGrady(i1)*GFactCF+uVelGrady(j1)*(1.-GFactCF)+
     *            vVelGradx(i1)*GFactCF+vVelGradx(j1)*(1.-GFactCF))
        Tau13=viscf*(
     *        uVelGradz(i1)*GFactCF+uVelGradz(j1)*(1.-GFactCF)+
     *            wVelGradx(i1)*GFactCF+wVelGradx(j1)*(1.-GFactCF))
        Tau22=2.*viscf*(vVelGrady(i1)*GFactCF+
     *                           vVelGrady(j1)*(1.-GFactCF))
        Tau23=viscf*(
     *        vVelGradz(i1)*GFactCF+vVelGradz(j1)*(1.-GFactCF)+
     *            wVelGrady(i1)*GFactCF+wVelGrady(j1)*(1.-GFactCF))
        Tau33=2.*viscf*(wVelGradz(i1)*GFactCF+
     *                           wVelGradz(j1)*(1.-GFactCF))
c
        if(Lcompressible) then
c
          term1=uVelGradx(i1)*GFactCF+uVelGradx(j1)*(1.-GFactCF)+
     *            vVelGrady(i1)*GFactCF+vVelGrady(j1)*(1.-GFactCF)+
     *              wVelGradz(i1)*GFactCF+wVelGradz(j1)*(1.-GFactCF)
c
          Tau11=Tau11-twothird*viscf*term1
          Tau22=Tau22-twothird*viscf*term1
          Tau33=Tau33-twothird*viscf*term1
c          
        endif
c
        if(LSolveTurbulenceKineticEnergy) then
c
          tkef=dmax1(TurbulentKE(i1)*GFactCF+
     *                    TurbulentKE(j1)*(1.-GFactCF),0.)
c
          densf=Density(i1)*GFactCF+Density(j1)*(1.-GFactCF)
          Tau11=Tau11-twothird*densf*tkef
          Tau22=Tau22-twothird*densf*tkef
          Tau33=Tau33-twothird*densf*tkef
c
        endif
c
        uVelf=uVelocity(i1)*GFactCF+uVelocity(j1)*(1.-GFactCF)
        vVelf=vVelocity(i1)*GFactCF+vVelocity(j1)*(1.-GFactCF)
        wVelf=wVelocity(i1)*GFactCF+wVelocity(j1)*(1.-GFactCF)
c
        term1=(Tau11*uVelf+Tau12*vVelf+Tau13*wVelf)*BFaceAreax(i2,i3)
        term2=(Tau12*uVelf+Tau22*vVelf+Tau23*wVelf)*BFaceAreay(i2,i3)
        term3=(Tau13*uVelf+Tau23*vVelf+Tau33*wVelf)*BFaceAreaz(i2,i3)
c
        FluxTElocal=term1+term2+term3
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo
c----------------------------------------------------------------------
      do i=1,Iaxis
c
        i1=IaxisOwner(i)
        i2=IaxisNumberOfBCSets(i)
        i3=IaxisNBFaces(i)
c
        FluxTElocal=0.
c
        FluxTE(i1)= FluxTE(i1)-FluxTElocal
c
      enddo
c
      return
      end