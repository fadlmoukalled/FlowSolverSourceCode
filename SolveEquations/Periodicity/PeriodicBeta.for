c
c#############################################################################################
c
      SUBROUTINE updateperiodicBeta
c
c#############################################################################################
c
      use User0, only: MethodDecomposeS
      use BoundaryConditions2, only: Iperiodic,IperiodicOwner,
     *                               IperiodicNumberOfBCSets,
     *                               IperiodicNBFaces,PeriodicPair,
     *                               Icorrespondingface,a1r,b1r,c1r,
     *                               a2r,b2r,c2r,a3r,b3r,c3r,
     *                               xTranslation,yTranslation,
     *                               zTranslation,xTranslationUse,
     *                               yTranslationUse,zTranslationUse,
     *                               periodicBeta,periodicMdot,
     *                               relaxBeta,BetaIterations,
     *                               nIterCorrectBeta
      use PhysicalProperties1, only: BDensity,Densityf
      use Variables1, only: uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity,
     *                      mdot,Bmdot,Du2Velocity,Dv2Velocity,
     *                      Dw2Velocity
      use Geometry4, only: xc,yc,zc,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz,BFaceAreax,BFaceAreay,
     *                     BFaceAreaz,FaceAreax,FaceAreay,FaceAreaz,
     *                     GFactorCF,DistanceCFux,DistanceCFuy,
     *                     DistanceCFuz,DistanceCF
      use Geometry3, only: NBFaceOwner,NIFaces,NBFaces,NIFaceOwner,
     *                     NIFaceNeighbor
      use Geometry1, only: NumberOfElements
      use MultiGrid2, only: nIter
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j,j1,j2,j3,k
      double precision :: xF1,yF1,zF1,distance1,distance2,GFactCF,
     *                    rhof,sfx,sfy,sfz,DistCFx,DistCFy,DistCFz,
     *                    DistCF,DistCFux,DistCFuy,DistCFuz,
     *                    ex,ey,ez,cf,Du1f,Dv1f,Dw1f,DuSf,DvSf,
     *                    DwSf,DuEf,DvEf,DwEf,geoDiff,Magnitude,
     *                    eDuSf,eDvSf,eDwSf,mdotin,rhogeodiff,dpdxprime,
     *                    gf,LTranslation,betax,betay,betaz
c
c********************************************************************************************
c
      if(mod(nIter,nIterCorrectBeta).ne.0) return
c
      dpdxprime=0.
c
      do k=1,BetaIterations
c
        mdotin=0.
        rhogeodiff=0.
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          j2=PeriodicPair(i2)         
c
          if(j2.gt.i2) then
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
            xF1=xc(j1)+xTranslation(j2)
            yF1=yc(j1)+yTranslation(j2)
            zF1=zc(j1)+zTranslation(j2)
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
            rhof=BDensity(i2,i3)
            sfx=BFaceAreax(i2,i3)
            sfy=BFaceAreay(i2,i3)
            sfz=BFaceAreaz(i2,i3)
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
c--- difference in pressure gradients term
c
            if(MethodDecomposeS.eq.1) then
c
              DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*Du2Velocity(j1)
              DuSf=DuSf*sfx
              DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*Dv2Velocity(j1)
              DvSf=DvSf*sfy
              DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*Dw2Velocity(j1)
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
              DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*Du2Velocity(j1)
              DuSf=DuSf*sfx
              DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*Dv2Velocity(j1)
              DvSf=DvSf*sfy
              DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*Dw2Velocity(j1)
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
            mdotin=mdotin+Bmdot(i2,i3)
            rhogeodiff=rhogeodiff+rhof*geoDiff
c
          endif
c
        enddo
c
        dpdxprime=(periodicMdot-mdotin)/rhogeodiff
        periodicBeta=periodicBeta+relaxBeta*dpdxprime
c
      enddo
c
      do i=1,Iperiodic
c
        i1=IperiodicOwner(i)
        i2=IperiodicNumberOfBCSets(i)
        i3=IperiodicNBFaces(i)
        i4=NIFaces
c
        j2=PeriodicPair(i2)         
c
        if(j2.gt.i2) then
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
          BuVelocity(i2,i3)=BuVelocity(i2,i3)+
     *                                    Du2Velocity(i1)*dpdxprime
          BvVelocity(i2,i3)=BvVelocity(i2,i3)+
     *                                    Dv2Velocity(i1)*dpdxprime
          BwVelocity(i2,i3)=BwVelocity(i2,i3)+
     *                                    Dw2Velocity(i1)*dpdxprime
          BuVelocity(j2,j3)=BuVelocity(i2,i3)
          BvVelocity(j2,j3)=BvVelocity(i2,i3)
          BwVelocity(j2,j3)=BwVelocity(i2,i3)
c
          xF1=xc(j1)+xTranslation(j2)
          yF1=yc(j1)+yTranslation(j2)
          zF1=zc(j1)+zTranslation(j2)
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
          rhof=BDensity(i2,i3)
          sfx=BFaceAreax(i2,i3)
          sfy=BFaceAreay(i2,i3)
          sfz=BFaceAreaz(i2,i3)
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
c--- difference in pressure gradients term
c
          if(MethodDecomposeS.eq.1) then
c
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*Du2Velocity(j1)
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*Dv2Velocity(j1)
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*Dw2Velocity(j1)
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
            DuSf=GFactCF*Du2Velocity(i1)+(1.-GFactCF)*Du2Velocity(j1)
            DuSf=DuSf*sfx
            DvSf=GFactCF*Dv2Velocity(i1)+(1.-GFactCF)*Dv2Velocity(j1)
            DvSf=DvSf*sfy
            DwSf=GFactCF*Dw2Velocity(i1)+(1.-GFactCF)*Dw2Velocity(j1)
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
          Bmdot(i2,i3)=Bmdot(i2,i3)+rhof*geoDiff*dpdxprime
          Bmdot(j2,j3)=-Bmdot(i2,i3)
c
        endif
c
      enddo
c
      LTranslation=dsqrt(xTranslationUse**2+
     *                 yTranslationUse**2+zTranslationUse**2)
c
      betax=xTranslationUse/LTranslation
      betay=yTranslationUse/LTranslation
      betaz=zTranslationUse/LTranslation
c
      do k=1,NumberOfElements
c
        uVelocity(k)=uVelocity(k)+Du2Velocity(k)*dpdxprime*betax  
        vVelocity(k)=vVelocity(k)+Dv2Velocity(k)*dpdxprime*betay  
        wVelocity(k)=wVelocity(k)+Dw2Velocity(k)*dpdxprime*betaz  
c
      enddo
c
      do k=1,NIFaces
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
        gf=GFactorCF(k)
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
        Du1f=gf*Du2Velocity(i)+(1.-gf)*Du2Velocity(j)
        Dv1f=gf*Dv2Velocity(i)+(1.-gf)*Dv2Velocity(j)
        Dw1f=gf*Dw2Velocity(i)+(1.-gf)*Dw2Velocity(j)
c
c--- difference in pressure gradients term
c
        if(MethodDecomposeS.eq.1) then
c
          DuSf=Du1f*sfx
          DvSf=Dv1f*sfy
          DwSf=Dw1f*sfz
c
          DuEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ex
          DvEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ey
          DwEf=(ex*DuSf+ey*DvSf+ez*DwSf)*ez
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
        elseif(MethodDecomposeS.eq.2) then
c
          DuSf=Du1f*sfx
          DvSf=Dv1f*sfy
          DwSf=Dw1f*sfz
c          
          DuEf=ex*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
          DvEf=ey*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
          DwEf=ez*dsqrt(DuSf*DuSf+DvSf*DvSf+DwSf*DwSf)
c          
          geoDiff=dsqrt(DuEf*DuEf+DvEf*DvEf+DwEf*DwEf)/cf         
c
        elseif(MethodDecomposeS.eq.3) then
c
          DuSf=Du1f*sfx
          DvSf=Dv1f*sfy
          DwSf=Dw1f*sfz
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
        mdot(k)=mdot(k)+rhof*geoDiff*dpdxprime
c
      enddo
c
      return
      end