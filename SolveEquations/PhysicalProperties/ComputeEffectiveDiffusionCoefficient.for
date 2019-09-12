c
c#############################################################################################
c
      SUBROUTINE CalculateEffectiveDiffusionCoefficient(Variable)
c
C#############################################################################################
c
      use User0, only: LTurbulentFlow,TurbulenceModel,WallTreatment,
     *                 LNegativeSpalartAllmaras
      use PhysicalProperties1
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Scalar2
      use Turbulence1
      use Variables1, only: ModifiedED,BModifiedED
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,k
      double precision :: sigTKElocal,sigTEDlocal,sigR
c********************************************************************************************
c
      if(LTurbulentFlow) then
c
        if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.
     *                          Variable.eq.'velz') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=Viscosity(i)+TurbulentViscosity(i) 
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                             BTurbulentViscosity(i,j) 
c
            enddo
          enddo
c
          if(WallTreatment.eq.'wallfunctions') then
c
            call MomentumWallFunctions
c
          endif
c
        elseif(Variable.eq.'temp') then
c
          if(TurbulenceModel.eq.'kklomega') then
c          
            call CalculateKKLWAlfaTheta         
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Conductivity(i)+
     *               SpecificHeat(i)*Density(i)*AlfaTheta(i)
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BConductivity(i,j)+
     *            BSpecificHeat(i,j)*BDensity(i,j)*BAlfaTheta(i,j)
c
              enddo
            enddo
c
          else
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Conductivity(i)+
     *           SpecificHeat(i)*TurbulentViscosity(i)/SigT
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BConductivity(i,j)+
     *            BSpecificHeat(i,j)*BTurbulentViscosity(i,j)/SigT
c
              enddo
            enddo
c
            if(WallTreatment.eq.'wallfunctions') then
c
              call EnergyWallFunctions
c
            endif
c
          endif
c
        elseif(Variable.eq.'tke') then
c
          if(TurbulenceModel.eq.'komegasst'.or.
     *                TurbulenceModel.eq.'komegabsl'.or.
     *                   TurbulenceModel.eq.'sstgama'.or.
     *                   TurbulenceModel.eq.'sstgamaretheta') then
c
            do i=1,NumberOfElements
c
              sigTKElocal=F1factor(i)*sigTKE1+(1.-F1factor(i))*sigTKE2
c
              eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)*sigTKElocal 
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                sigTKElocal=BF1factor(i,j)*sigTKE1+
     *                                   (1.-BF1factor(i,j))*sigTKE2
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                          BTurbulentViscosity(i,j)*sigTKElocal
c
              enddo
            enddo
c
            if(TurbulenceModel.eq.'sstgamaretheta') 
     *                                   call CalculateTGammaEff
c
          elseif(TurbulenceModel.eq.'kklmodel') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)*sigTKE 
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                            BTurbulentViscosity(i,j)*sigTKE
c
              enddo
            enddo
c
          elseif(TurbulenceModel.eq.'kklomega') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=
     *               Viscosity(i)+Density(i)*AlfaT(i)/sigTKE 
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=
     *              BViscosity(i,j)+BDensity(i,j)*BAlfaT(i,j)/sigTKE
c
              enddo
            enddo
c
          elseif(TurbulenceModel.eq.'kepsilonrng') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=sigTKERNG(i)*
     *                    (Viscosity(i)+TurbulentViscosity(i))
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BsigTKERNG(i,j)*
     *                (BViscosity(i,j)+BTurbulentViscosity(i,j))
c
              enddo
            enddo
c
          else
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)/sigTKE 
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                            BTurbulentViscosity(i,j)/sigTKE
c
              enddo
            enddo
c
          endif
c
        elseif(Variable.eq.'ted') then
c
          if(TurbulenceModel.eq.'kepsilonrng') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=sigTEDRNG(i)*
     *                    (Viscosity(i)+TurbulentViscosity(i))
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BsigTEDRNG(i,j)*
     *                (BViscosity(i,j)+BTurbulentViscosity(i,j))
c
              enddo
            enddo
c
          else 
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)/sigTED 
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                             BTurbulentViscosity(i,j)/sigTED  
c
              enddo
            enddo
c
          endif
c
        elseif(Variable.eq.'tomega') then
c
          if(TurbulenceModel.eq.'komegasst'.or.
     *                TurbulenceModel.eq.'komegabsl'.or.
     *                   TurbulenceModel.eq.'sstgama'.or.
     *                   TurbulenceModel.eq.'sstgamaretheta') then
c
            do i=1,NumberOfElements
c
              sigTEDlocal=F1factor(i)*sigTED1+(1.-F1factor(i))*sigTED2
c
              eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)*sigTEDlocal
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                sigTEDlocal=BF1factor(i,j)*sigTED1+
     *                             (1.-BF1factor(i,j))*sigTED2
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                             BTurbulentViscosity(i,j)*sigTEDlocal
c
              enddo
            enddo
c
          elseif(TurbulenceModel.eq.'kklomega') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=
     *               Viscosity(i)+Density(i)*AlfaT(i)/sigTED
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=
     *               BViscosity(i,j)+BDensity(i,j)*BAlfaT(i,j)/sigTED
c
              enddo
            enddo
c
          else
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)/sigTED 
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                             BTurbulentViscosity(i,j)/sigTED
c
              enddo
            enddo
c
          endif
c
        elseif(Variable.eq.'tkl') then
c
          if(TurbulenceModel.eq.'kklmodel') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)*sigKL
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                             BTurbulentViscosity(i,j)*sigKL
c
              enddo
            enddo
c
          elseif(TurbulenceModel.eq.'kklomega') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=Viscosity(i)
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=BViscosity(i,j)
c
              enddo
            enddo
c
          endif  
c
        elseif(Variable.eq.'tgamma') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)/sigTGamma
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                             BTurbulentViscosity(i,j)/sigTGamma
c
            enddo
          enddo
c
        elseif(Variable.eq.'tretheta') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=(Viscosity(i)+
     *                            TurbulentViscosity(i))*sigTReTheta
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=(BViscosity(i,j)+
     *                             BTurbulentViscosity(i,j))*sigTReTheta
c
            enddo
          enddo
c
        elseif(Variable.eq.'frelax') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=LScale(i)*LScale(i)
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BLScale(i,j)*BLScale(i,j)
c
            enddo
          enddo
c
        elseif(Variable.eq.'tv2') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=Viscosity(i)+
     *                            TurbulentViscosity(i)/sigTKE 
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BViscosity(i,j)+
     *                            BTurbulentViscosity(i,j)/sigTKE
c
            enddo
          enddo
c
        elseif(Variable.eq.'tzeta') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=TurbulentViscosity(i)/sigTKE 
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BTurbulentViscosity(i,j)/sigTKE
c
            enddo
          enddo
c
        elseif(Variable.eq.'med') then
c
          if(TurbulenceModel.eq.'spalartallmaras') then
c
            if(.not.LNegativeSpalartAllmaras) then
c
              do i=1,NumberOfElements
c
                eDiffCoefficient(i)=(Viscosity(i)+
     *                            Density(i)*ModifiedED(i))/sigMED
c
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
c
                  BeDiffCoefficient(i,j)=(BViscosity(i,j)+
     *                         BDensity(i,j)*BModifiedED(i,j))/sigMED
c
                enddo
              enddo
c
            elseif(LNegativeSpalartAllmaras) then
c
              do i=1,NumberOfElements
c
                if(ModifiedED(i).ge.0.) then
c
                  eDiffCoefficient(i)=(Viscosity(i)+
     *                            Density(i)*ModifiedED(i))/sigMED
c
                else
c
                  eDiffCoefficient(i)=(Viscosity(i)+
     *                fnCoefficient(i)*Density(i)*ModifiedED(i))/sigMED
c
                endif
c
              enddo
c
              do i=1,NumberOfBCSets
                do j=1,NBFaces(i)
c
                  if(BModifiedED(i,j).ge.0.) then
c
                    BeDiffCoefficient(i,j)=(BViscosity(i,j)+
     *                         BDensity(i,j)*BModifiedED(i,j))/sigMED
c
                  else
c
                    BeDiffCoefficient(i,j)=(BViscosity(i,j)+
     *                         BfnCoefficient(i,j)*BDensity(i,j)*
     *                                       BModifiedED(i,j))/sigMED
c
                  endif
c
                enddo
              enddo
c
            endif
c
          elseif(TurbulenceModel.eq.'wrayagarwal') then
c
            do i=1,NumberOfElements
c
              sigR=f1WA(i)*(sigKW-sigKE)+sigKE
              eDiffCoefficient(i)=
     *                   Viscosity(i)+Density(i)*sigR*ModifiedED(i)
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                sigR=Bf1WA(i,j)*(sigKW-sigKE)+sigKE
                BeDiffCoefficient(i,j)=
     *              BViscosity(i,j)+BDensity(i,j)*sigR*BModifiedED(i,j)
c
              enddo
            enddo
c
          elseif(TurbulenceModel.eq.'kepsilonrt') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=
     *                   Viscosity(i)+TurbulentViscosity(i)/sigRT
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=
     *              BViscosity(i,j)+BTurbulentViscosity(i,j)/sigRT
c
              enddo
            enddo
c
          elseif(TurbulenceModel.eq.'nut92') then
c
            do i=1,NumberOfElements
c
              eDiffCoefficient(i)=
     *                   Viscosity(i)+Cnut0*TurbulentViscosity(i)
c
            enddo
c
            do i=1,NumberOfBCSets
              do j=1,NBFaces(i)
c
                BeDiffCoefficient(i,j)=
     *              BViscosity(i,j)+Cnut0*BTurbulentViscosity(i,j)
c
              enddo
            enddo
c
          endif
c
        else
c
          k=iScalarVariable
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=DiffusionCoefficient(i,k)+
     *        SpecificHeatScalar(i,k)*TurbulentViscosity(i)/SigScalar(k)
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BDiffusionCoefficient(i,j,k)+
     *                         BSpecificHeatScalar(i,j,k)*
     *                             TurbulentViscosity(i)/SigScalar(k)
c
            enddo
          enddo
c
          if(WallTreatment.eq.'wallfunctions') then
c
            call ScalarWallFunctions
c
          endif
c
        endif
c
      else
c
        if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.
     *                            Variable.eq.'velz') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=Viscosity(i)
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BViscosity(i,j)
c
            enddo
          enddo
c
        elseif(Variable.eq.'temp') then
c
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=Conductivity(i)
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BConductivity(i,j)
c
            enddo
          enddo
c
        else
c
          k=iScalarVariable
          do i=1,NumberOfElements
c
            eDiffCoefficient(i)=DiffusionCoefficient(i,k)
c
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              BeDiffCoefficient(i,j)=BDiffusionCoefficient(i,j,k)
c
            enddo
          enddo
c
        endif
c
      endif
c
      return
      end
