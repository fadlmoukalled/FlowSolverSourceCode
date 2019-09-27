c
C#############################################################################################
c
      SUBROUTINE HRSchemesCorrections(Variable,ConvectionScheme,
     *               Bleed,HRFramework,FiT,dfidxT,dfidyT,dfidzT)
c
c#############################################################################################
c
      use geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor
      use Variables1, only: mdot
      use Variables2, only: FluxTf
      use Scalar2
      use PhysicalProperties1, only: SpecificHeat,SpecificHeatScalar
      use Geometry4, only: GFactorCF,DistanceCFx,DistanceCFy,DistanceCFz
      use VolumeOfFluid2, only: i,j,k
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*20 ConvectionScheme
      character*10 Variable
      character*4 HRFramework
      double precision Bleed
      double precision, save :: phic,phid,phiu,phiTeldaC,diff,phiHR,
     *                          cpf,gf,rf,psirf,phiTeldaf,dcfx,dcfy,dcfz
      integer, save :: iUpwind
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
c********************************************************************************************
c      
      if(ConvectionScheme.eq.'upwind') return
c
c--- Implement HR schemes (NVF and TVD formulations)
c
      do k=1,NIFaces
c        
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        if(mdot(k).gt.0.) then
          
          iUpwind=i
          phic=FiT(i)
          phid=FiT(j)
          dcfx=DistanceCFx(k)
          dcfy=DistanceCFy(k)
          dcfz=DistanceCFz(k)
c
        else
c
          iUpwind=j
          phic=FiT(j)
          phid=FiT(i)
          dcfx=-DistanceCFx(k)
          dcfy=-DistanceCFy(k)
          dcfz=-DistanceCFz(k)
c
        endif
c          
        phiu=phid-2.*(dfidxT(iUpwind)*dcfx+
     *                dfidyT(iUpwind)*dcfy+dfidzT(iUpwind)*dcfz)
c
c--- NVF implementation
c
        if(HRFramework.eq.'nvf') then
c
          diff=phid-phiu
c
            if(dabs(diff).lt.1.e-9) then 
c
              phiTeldaC=0.
c
            else
c
              phiTeldaC=(phic-phiu)/diff
c
          endif
c
          call NVFLimiter(ConvectionScheme,phiTeldaC,phiTeldaf)
c
          phiHR=phiTeldaf*diff+phiu
          if(ConvectionScheme.ne.'gamma') 
     *                       phiHR=(1.-Bleed)*phiHR+Bleed*phic
c
        elseif(HRFramework.eq.'tvd') then
c
c--- TVD Implementation
c
          diff=phid-phic
c
          if(dabs(diff).lt.1.e-9) then 
c
            rf=0.
c
          else
c
            rf=(phic-phiu)/diff
c
          endif
c
          call TVDLimiter(ConvectionScheme,rf,psirf)
c
          phiHR=phic+0.5*psirf*diff
          phiHR=(1.-Bleed)*phiHR+Bleed*phic
c
        endif
c
          cpf=1.
          gf=GFactorCF(k)
c
          if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.
     *       Variable.eq.'velz'.or.Variable.eq.'htotal'.or.
     *        Variable.eq.'tgamma'.or.Variable.eq.'tretheta'.or.
     *          Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *            Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *               Variable.eq.'tkl'.or.Variable.eq.'rfield'.or.
     *               Variable.eq.'tv2'.or.Variable.eq.'tzeta') then
            cpf=1.
c
          elseif(Variable.eq.'temp') then
c
            cpf=gf*SpecificHeat(i)+(1.-gf)*SpecificHeat(j)
c
          else
c
            cpf=gf*SpecificHeatScalar(i,iScalarVariable)+
     *            (1.-gf)*SpecificHeatScalar(j,iScalarVariable)
c
          endif
c
          FluxTf(k)=FluxTf(k)+mdot(k)*(phiHR-phic)*cpf                
c
      enddo            
c
      return
      end