c
C#############################################################################################
      SUBROUTINE calculateSprime(Variable)
C#############################################################################################
c
      use PhysicalProperties1, only: Conductivity11,Conductivity12,
     *                               Conductivity13,Conductivity22,
     *                               Conductivity23,Conductivity33,
     *                               BConductivity11,BConductivity12,
     *                               BConductivity13,BConductivity22,
     *                               BConductivity23,BConductivity33,
     *                 DiffusionCoefficient11,DiffusionCoefficient12,
     *                 DiffusionCoefficient13,DiffusionCoefficient22,
     *                 DiffusionCoefficient23,DiffusionCoefficient33,
     *                 BDiffusionCoefficient11,BDiffusionCoefficient12,
     *                 BDiffusionCoefficient13,BDiffusionCoefficient22,
     *                 BDiffusionCoefficient23,BDiffusionCoefficient33,
     *                 InterpolationSchemeGamaScalar
      use User0, only: MethodDecomposeSprime,LSolveLambdaELEEquation,
     *                 InterpolationSchemeGamaEnergy
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces,NIFaceOwner,NIFaceNeighbor
      use Geometry4, only: GFactorCF,FaceAreax,FaceAreay,BFaceAreax,
     *                     BFaceAreay,DistanceCFux,DistanceCFuy,
     *                     BDistanceCFux,BDistanceCFuy,DistanceCF,
     *                     BDistanceCF,FaceAreanx,FaceAreany,
     *                     BFaceAreanx,BFaceAreany,
     *                     FaceAreap,FaceAreaxp,FaceAreayp,BFaceAreap,
     *                     BFaceAreaxp,BFaceAreayp,FaceEp,FaceExp,
     *                     FaceEyp,gDiffp,FaceTp,FaceTxp,FaceTyp,
     *                     BFaceEp,BFaceExp,BFaceEyp,BgDiffp,BFaceTp,
     *                     BFaceTxp,BFaceTyp,BDistanceCFuz,
     *                     FaceAreaz,FaceAreazp,BFaceAreaz,BFaceAreazp,
     *                     DistanceCFuz,FaceEzp,FaceTzp,BFaceEzp,
     *                     BFaceTzp,FaceAreanz,BFaceAreanz
      use Variables1, only: mdot
      use Scalar2
      use Constants1, only: tiny
      use MultiGrid2, only: nIter
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k,i5
      double precision :: gamf11,gamf12,gamf13,gamf21,gamf22,gamf23,
     *                    gamf31,gamf32,gamf33,DotProduct,
     *                    DotProduct1,DotProduct2,Ratio,gf
      character*10 Variable
c*********************************************************************************************
c
      if(LSolveLambdaELEEquation.and.nIter.gt.1.and.
     *                                  Variable.eq.'lambda') return
      if(Variable.eq.'temp'.or.Variable.eq.'lambda') then 
c
        if(InterpolationSchemeGamaEnergy.eq.'average') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=GFactorCF(k)
c
            gamf11=gf*Conductivity11(i)+(1.-gf)*Conductivity11(j)
            gamf12=gf*Conductivity12(i)+(1.-gf)*Conductivity12(j)
            gamf13=gf*Conductivity13(i)+(1.-gf)*Conductivity13(j)
            gamf21=gamf12
            gamf22=gf*Conductivity22(i)+(1.-gf)*Conductivity22(j)
            gamf23=gf*Conductivity23(i)+(1.-gf)*Conductivity23(j)
            gamf31=gamf13
            gamf32=gamf23
            gamf33=gf*Conductivity33(i)+(1.-gf)*Conductivity33(j)
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *              gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *              gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *              gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                   FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        elseif(InterpolationSchemeGamaEnergy.eq.'upwind') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=0.
            if(mdot(k).gt.0.) gf=1.          
c
            gamf11=gf*Conductivity11(i)+(1.-gf)*Conductivity11(j)
            gamf12=gf*Conductivity12(i)+(1.-gf)*Conductivity12(j)
            gamf13=gf*Conductivity13(i)+(1.-gf)*Conductivity13(j)
            gamf21=gamf12
            gamf22=gf*Conductivity22(i)+(1.-gf)*Conductivity22(j)
            gamf23=gf*Conductivity23(i)+(1.-gf)*Conductivity23(j)
            gamf31=gamf13
            gamf32=gamf23
            gamf33=gf*Conductivity33(i)+(1.-gf)*Conductivity33(j)
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *              gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *              gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *              gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                   FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        elseif(InterpolationSchemeGamaEnergy.eq.'downwind') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=1.
            if(mdot(k).gt.0.) gf=0.          
c
            gamf11=gf*Conductivity11(i)+(1.-gf)*Conductivity11(j)
            gamf12=gf*Conductivity12(i)+(1.-gf)*Conductivity12(j)
            gamf13=gf*Conductivity13(i)+(1.-gf)*Conductivity13(j)
            gamf21=gamf12
            gamf22=gf*Conductivity22(i)+(1.-gf)*Conductivity22(j)
            gamf23=gf*Conductivity23(i)+(1.-gf)*Conductivity23(j)
            gamf31=gamf13
            gamf32=gamf23
            gamf33=gf*Conductivity33(i)+(1.-gf)*Conductivity33(j)
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *              gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *              gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *              gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                   FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        elseif(InterpolationSchemeGamaEnergy.eq.'harmonic') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=GFactorCF(k)
c
            gamf11=Conductivity11(i)*Conductivity11(j)
            if(gamf11.ne.0.) gamf11=gamf11/
     *             (gf*Conductivity11(i)+(1.-gf)*Conductivity11(j))
            gamf12=Conductivity12(i)*Conductivity12(j)
            if(gamf12.ne.0.) gamf12=gamf12/
     *             (gf*Conductivity12(i)+(1.-gf)*Conductivity12(j))
            gamf13=Conductivity13(i)*Conductivity13(j)
            if(gamf13.ne.0.) gamf13=gamf13/
     *             (gf*Conductivity13(i)+(1.-gf)*Conductivity13(j))
            gamf21=gamf12
            gamf22=Conductivity22(i)*Conductivity22(j)
            if(gamf22.ne.0.) gamf22=gamf22/
     *             (gf*Conductivity22(i)+(1.-gf)*Conductivity22(j))
            gamf23=Conductivity23(i)*Conductivity23(j)
            if(gamf23.ne.0.) gamf23=gamf23/
     *             (gf*Conductivity23(i)+(1.-gf)*Conductivity23(j))
            gamf31=gamf13
            gamf32=gamf23
            gamf33=Conductivity33(i)*Conductivity33(j)
            if(gamf33.ne.0.) gamf33=gamf33/
     *             (gf*Conductivity33(i)+(1.-gf)*Conductivity33(j))
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *              gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *              gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *              gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                   FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        endif
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BFaceAreaxp(i,j)=BConductivity11(i,j)*BFaceAreax(i,j)+
     *                          BConductivity12(i,j)*BFaceAreay(i,j)+
     *                              BConductivity13(i,j)*BFaceAreaz(i,j)
            BFaceAreayp(i,j)=BConductivity12(i,j)*BFaceAreax(i,j)+
     *                          BConductivity22(i,j)*BFaceAreay(i,j)+
     *                              BConductivity23(i,j)*BFaceAreaz(i,j)
            BFaceAreazp(i,j)=BConductivity13(i,j)*BFaceAreax(i,j)+
     *                          BConductivity23(i,j)*BFaceAreay(i,j)+
     *                              BConductivity33(i,j)*BFaceAreaz(i,j)
            BFaceAreap(i,j)=dsqrt(BFaceAreaxp(i,j)**2+
     *                          BFaceAreayp(i,j)**2+BFaceAreazp(i,j)**2)
c
          enddo
        enddo
c
      else
c      
        i5=iScalarVariable
c
        if(InterpolationSchemeGamaScalar(i5).eq.'average') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=GFactorCF(k)
c
            gamf11=gf*DiffusionCoefficient11(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient11(j,i5)
            gamf12=gf*DiffusionCoefficient12(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient12(j,i5)
            gamf13=gf*DiffusionCoefficient13(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient13(j,i5)
            gamf21=gamf12
            gamf22=gf*DiffusionCoefficient22(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient22(j,i5)
            gamf23=gf*DiffusionCoefficient23(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient23(j,i5)
            gamf31=gamf13
            gamf32=gamf23
            gamf33=gf*DiffusionCoefficient33(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient33(j,i5)
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *                          gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *                          gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *                          gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                               FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        elseif(InterpolationSchemeGamaScalar(i5).eq.'upwind') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=0.
            if(mdot(k).gt.0.) gf=1.
c
            gamf11=gf*DiffusionCoefficient11(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient11(j,i5)
            gamf12=gf*DiffusionCoefficient12(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient12(j,i5)
            gamf13=gf*DiffusionCoefficient13(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient13(j,i5)
            gamf21=gamf12
            gamf22=gf*DiffusionCoefficient22(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient22(j,i5)
            gamf23=gf*DiffusionCoefficient23(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient23(j,i5)
            gamf31=gamf13
            gamf32=gamf23
            gamf33=gf*DiffusionCoefficient33(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient33(j,i5)
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *                          gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *                          gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *                          gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                               FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        elseif(InterpolationSchemeGamaScalar(i5).eq.'downwind') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=1.
            if(mdot(k).gt.0.) gf=0.
c
            gamf11=gf*DiffusionCoefficient11(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient11(j,i5)
            gamf12=gf*DiffusionCoefficient12(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient12(j,i5)
            gamf13=gf*DiffusionCoefficient13(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient13(j,i5)
            gamf21=gamf12
            gamf22=gf*DiffusionCoefficient22(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient22(j,i5)
            gamf23=gf*DiffusionCoefficient23(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient23(j,i5)
            gamf31=gamf13
            gamf32=gamf23
            gamf33=gf*DiffusionCoefficient33(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient33(j,i5)
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *                          gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *                          gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *                          gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                               FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        elseif(InterpolationSchemeGamaScalar(i5).eq.'harmonic') then
c
          do k=1,NIFaces
c        
            i=NIFaceOwner(k)
            j=NIFaceNeighbor(k)
c
            gf=GFactorCF(k)
c
            gamf11=DiffusionCoefficient11(i,i5)*
     *                         DiffusionCoefficient11(i,i5)
            if(gamf11.ne.0.) 
     *            gamf11=gamf11/(gf*DiffusionCoefficient11(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient11(j,i5))
            gamf12=DiffusionCoefficient12(i,i5)*
     *                         DiffusionCoefficient12(i,i5)
            if(gamf12.ne.0.) 
     *            gamf12=gamf12/(gf*DiffusionCoefficient12(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient12(j,i5))
            gamf13=DiffusionCoefficient13(i,i5)*
     *                         DiffusionCoefficient13(i,i5)
            if(gamf13.ne.0.) 
     *            gamf13=gamf13/(gf*DiffusionCoefficient13(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient13(j,i5))
            gamf21=gamf12
            gamf22=DiffusionCoefficient22(i,i5)*
     *                         DiffusionCoefficient22(i,i5)
            if(gamf22.ne.0.) 
     *            gamf22=gamf22/(gf*DiffusionCoefficient22(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient22(j,i5))
            gamf23=DiffusionCoefficient23(i,i5)*
     *                         DiffusionCoefficient23(i,i5)
            if(gamf23.ne.0.) 
     *            gamf23=gamf23/(gf*DiffusionCoefficient23(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient23(j,i5))
            gamf31=gamf13
            gamf32=gamf23
            gamf33=DiffusionCoefficient33(i,i5)*
     *                         DiffusionCoefficient33(i,i5)
            if(gamf33.ne.0.) 
     *            gamf33=gamf33/(gf*DiffusionCoefficient33(i,i5)+
     *                   (1.-gf)*DiffusionCoefficient33(j,i5))
c            
            FaceAreaxp(k)=gamf11*FaceAreax(k)+
     *                          gamf21*FaceAreay(k)+gamf31*FaceAreaz(k)
            FaceAreayp(k)=gamf12*FaceAreax(k)+
     *                          gamf22*FaceAreay(k)+gamf32*FaceAreaz(k)
            FaceAreazp(k)=gamf13*FaceAreax(k)+
     *                          gamf23*FaceAreay(k)+gamf33*FaceAreaz(k)
            FaceAreap(k)=dsqrt(FaceAreaxp(k)**2+
     *                               FaceAreayp(k)**2+FaceAreazp(k)**2)
c
          enddo
c
        endif
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BFaceAreaxp(i,j)=
     *            BDiffusionCoefficient11(i,j,i5)*BFaceAreax(i,j)+
     *               BDiffusionCoefficient12(i,j,i5)*BFaceAreay(i,j)+
     *                  BDiffusionCoefficient13(i,j,i5)*BFaceAreaz(i,j)
            BFaceAreayp(i,j)=
     *            BDiffusionCoefficient12(i,j,i5)*BFaceAreax(i,j)+
     *               BDiffusionCoefficient22(i,j,i5)*BFaceAreay(i,j)+
     *                  BDiffusionCoefficient23(i,j,i5)*BFaceAreaz(i,j)
            BFaceAreazp(i,j)=
     *            BDiffusionCoefficient13(i,j,i5)*BFaceAreax(i,j)+
     *               BDiffusionCoefficient23(i,j,i5)*BFaceAreay(i,j)+
     *                  BDiffusionCoefficient33(i,j,i5)*BFaceAreaz(i,j)
            BFaceAreap(i,j)=dsqrt(BFaceAreaxp(i,j)**2+
     *                         BFaceAreayp(i,j)**2+BFaceAreazp(i,j)**2)
c
          enddo
        enddo
c      
      endif
c
c--- Decompose S' into E and T along internal faces
c
      if(MethodDecomposeSprime.eq.1) then
c
        do k=1,NIFaces
c
          DotProduct=DistanceCFux(k)*FaceAreaxp(k)+
     *      DistanceCFuy(k)*FaceAreayp(k)+DistanceCFuz(k)*FaceAreazp(k)
          FaceExp(k)=DotProduct*DistanceCFux(k)
          FaceEyp(k)=DotProduct*DistanceCFuy(k)
          FaceEzp(k)=DotProduct*DistanceCFuz(k)
          FaceEp(k)=dabs(DotProduct)
          gDiffp(k)=FaceEp(k)/DistanceCF(k)
c
          FaceTxp(k)=FaceAreaxp(k)-FaceExp(k)
          FaceTyp(k)=FaceAreayp(k)-FaceEyp(k)
          FaceTzp(k)=FaceAreazp(k)-FaceEzp(k)
          FaceTp(k)=dsqrt(FaceTxp(k)**2+FaceTyp(k)**2+FaceTzp(k)**2)
c
        enddo
c
      elseif(MethodDecomposeSprime.eq.2) then
c
        do k=1,NIFaces
c
          FaceExp(k)=FaceAreap(k)*DistanceCFux(k)
          FaceEyp(k)=FaceAreap(k)*DistanceCFuy(k)
          FaceEzp(k)=FaceAreap(k)*DistanceCFuz(k)
          FaceEp(k)=FaceAreap(k)
          gDiffp(k)=FaceEp(k)/DistanceCF(k)
c
          FaceTxp(k)=FaceAreaxp(k)-FaceExp(k)
          FaceTyp(k)=FaceAreayp(k)-FaceEyp(k)
          FaceTzp(k)=FaceAreazp(k)-FaceEzp(k)
          FaceTp(k)=dsqrt(FaceTxp(k)**2+FaceTyp(k)**2+FaceTzp(k)**2)
c
        enddo
c
      elseif(MethodDecomposeSprime.eq.3) then
c
        do k=1,NIFaces
c
          DotProduct=DistanceCFux(k)*FaceAreaxp(k)+
     *       DistanceCFuy(k)*FaceAreayp(k)+DistanceCFuz(k)*FaceAreazp(k)
          FaceExp(k)=(FaceAreap(k)**2/DotProduct)*DistanceCFux(k)
          FaceEyp(k)=(FaceAreap(k)**2/DotProduct)*DistanceCFuy(k)
          FaceEzp(k)=(FaceAreap(k)**2/DotProduct)*DistanceCFuz(k)
          FaceEp(k)=dsqrt(FaceExp(k)**2+FaceEyp(k)**2+FaceEzp(k)**2)
          gDiffp(k)=FaceEp(k)/DistanceCF(k)
c
          FaceTxp(k)=FaceAreaxp(k)-FaceExp(k)
          FaceTyp(k)=FaceAreayp(k)-FaceEyp(k)
          FaceTzp(k)=FaceAreazp(k)-FaceEzp(k)
          FaceTp(k)=dsqrt(FaceTxp(k)**2+FaceTyp(k)**2+FaceTzp(k)**2)
c
        enddo
c
      elseif(MethodDecomposeSprime.eq.4) then
c
        do k=1,NIFaces
c
          dotproduct1=FaceAreaxp(k)*FaceAreanx(k)+
     *                     FaceAreayp(k)*FaceAreany(k)+
     *                             FaceAreazp(k)*FaceAreanz(k)
          dotproduct2=DistanceCFux(k)*FaceAreanx(k)+
     *                     DistanceCFuy(k)*FaceAreany(k)+
     *                             DistanceCFuz(k)*FaceAreanz(k)
c
          ratio=dotproduct1/dotproduct2
          FaceExp(k)=ratio*DistanceCFux(k)
          FaceEyp(k)=ratio*DistanceCFuy(k)
          FaceEzp(k)=ratio*DistanceCFuz(k)
          FaceEp(k)=dabs(ratio)
c
          FaceTxp(k)=FaceAreaxp(k)-FaceExp(k)
          FaceTyp(k)=FaceAreayp(k)-FaceEyp(k)
          FaceTzp(k)=FaceAreazp(k)-FaceEzp(k)
          FaceTp(k)=dsqrt(FaceTxp(k)**2+FaceTyp(k)**2+FaceTzp(k)**2)
c
          gDiffp(k)=FaceEp(k)/DistanceCF(k)
c
        enddo
c
      endif
c
c--- Decompose S' into E and T along boundary faces
c
      if(MethodDecomposeSprime.eq.1) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            DotProduct=BDistanceCFux(i,j)*BFaceAreaxp(i,j)+
     *                   BDistanceCFuy(i,j)*BFaceAreayp(i,j)+
     *                     BDistanceCFuz(i,j)*BFaceAreazp(i,j)
            BFaceExp(i,j)=DotProduct*BDistanceCFux(i,j)
            BFaceEyp(i,j)=DotProduct*BDistanceCFuy(i,j)
            BFaceEzp(i,j)=DotProduct*BDistanceCFuz(i,j)
            BFaceEp(i,j)=dabs(DotProduct)
            BgDiffp(i,j)=BFaceEp(i,j)/BDistanceCF(i,j)
c
            BFaceTxp(i,j)=BFaceAreaxp(i,j)-BFaceExp(i,j)
            BFaceTyp(i,j)=BFaceAreayp(i,j)-BFaceEyp(i,j)
            BFaceTzp(i,j)=BFaceAreazp(i,j)-BFaceEzp(i,j)
            BFaceTp(i,j)=dsqrt(BFaceTxp(i,j)**2+
     *                BFaceTyp(i,j)**2+BFaceTzp(i,j)**2)
c
          enddo
        enddo
c
      elseif(MethodDecomposeSprime.eq.2) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            BFaceExp(i,j)=BFaceAreap(i,j)*BDistanceCFux(i,j)
            BFaceEyp(i,j)=BFaceAreap(i,j)*BDistanceCFuy(i,j)
            BFaceEzp(i,j)=BFaceAreap(i,j)*BDistanceCFuz(i,j)
            BFaceEp(i,j)=BFaceAreap(i,j)
            BgDiffp(i,j)=BFaceEp(i,j)/BDistanceCF(i,j)
c
            BFaceTxp(i,j)=BFaceAreaxp(i,j)-BFaceExp(i,j)
            BFaceTyp(i,j)=BFaceAreayp(i,j)-BFaceEyp(i,j)
            BFaceTzp(i,j)=BFaceAreazp(i,j)-BFaceEzp(i,j)
            BFaceTp(i,j)=dsqrt(BFaceTxp(i,j)**2+
     *                 BFaceTyp(i,j)**2+BFaceTzp(i,j)**2)
c
          enddo
        enddo
c
      elseif(MethodDecomposeSprime.eq.3) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            DotProduct=BDistanceCFux(i,j)*BFaceAreaxp(i,j)+
     *                   BDistanceCFuy(i,j)*BFaceAreayp(i,j)+
     *                     BDistanceCFuz(i,j)*BFaceAreazp(i,j)
            BFaceExp(i,j)=(BFaceAreap(i,j)**2/DotProduct)*
     *                                      BDistanceCFux(i,j)
            BFaceEyp(i,j)=(BFaceAreap(i,j)**2/DotProduct)*
     *                                      BDistanceCFuy(i,j)
            BFaceEzp(i,j)=(BFaceAreap(i,j)**2/DotProduct)*
     *                                      BDistanceCFuz(i,j)
            BFaceEp(i,j)=dsqrt(BFaceExp(i,j)**2+
     *                       BFaceEyp(i,j)**2+BFaceEzp(i,j)**2)
            BgDiffp(i,j)=BFaceEp(i,j)/BDistanceCF(i,j)
c
            BFaceTxp(i,j)=BFaceAreaxp(i,j)-BFaceExp(i,j)
            BFaceTyp(i,j)=BFaceAreayp(i,j)-BFaceEyp(i,j)
            BFaceTzp(i,j)=BFaceAreazp(i,j)-BFaceEzp(i,j)
            BFaceTp(i,j)=dsqrt(BFaceTxp(i,j)**2+
     *                       BFaceTyp(i,j)**2+BFaceTzp(i,j)**2)
c
          enddo
        enddo
c
      elseif(MethodDecomposeSprime.eq.4) then
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            dotproduct1=BFaceAreaxp(i,j)*BFaceAreanx(i,j)+
     *                     BFaceAreayp(i,j)*BFaceAreany(i,j)+
     *                             BFaceAreazp(i,j)*BFaceAreanz(i,j)
            dotproduct2=BDistanceCFux(i,j)*BFaceAreanx(i,j)+
     *                     BDistanceCFuy(i,j)*BFaceAreany(i,j)+
     *                             BDistanceCFuz(i,j)*BFaceAreanz(i,j)
c
            ratio=dotproduct1/dotproduct2
            BFaceExp(i,j)=ratio*BDistanceCFux(i,j)
            BFaceEyp(i,j)=ratio*BDistanceCFuy(i,j)
            BFaceEzp(i,j)=ratio*BDistanceCFuz(i,j)
            BFaceEp(i,j)=dabs(ratio)
c
            BFaceTxp(i,j)=BFaceAreaxp(i,j)-BFaceExp(i,j)
            BFaceTyp(i,j)=BFaceAreayp(i,j)-BFaceEyp(i,j)
            BFaceTzp(i,j)=BFaceAreazp(i,j)-BFaceEzp(i,j)
            BFaceTp(i,j)=dsqrt(BFaceTxp(i,j)**2+
     *                           BFaceTyp(i,j)**2+BFaceTzp(i,j)**2)
c
            BgDiffp(i,j)=BFaceEp(i,j)/BDistanceCF(i,j)
c
          enddo
        enddo
c
      endif
c
      return
      end