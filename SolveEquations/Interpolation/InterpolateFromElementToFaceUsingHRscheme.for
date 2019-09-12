c
C#############################################################################################
      SUBROUTINE InterpolateToFaceUsingHRscheme(ConvectionScheme,
     *               Bleed,NVF,TVD,FiT,dfidxT,dfidyT,dfidzT,FiTf)
c#############################################################################################
      use geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor
      use Variables1, only: mdot
      use Geometry4, only: GFactorCF,DistanceCFx,DistanceCFy,DistanceCFz
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      character*20 ConvectionScheme
      character*10 Variable
      logical NVF,TVD
      double precision :: Bleed
      double precision, save:: phic,phid,phiu,phiTeldaC,diff,phiHR,gf,
     *                 rf,psirf,phiTeldaf,dcfx,dcfy,dcfz
      integer i,j,k,iUpwind
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:) :: FiTf
c*********************************************************************************************
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
     *                   dfidyT(iUpwind)*dcfy+dfidzT(iUpwind)*dcfz)
c
c--- NVF implementation
c
        if(NVF) then
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
          phiHR=(1.-Bleed)*phiHR+Bleed*phic
c
        elseif(TVD) then
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
        else
c
          phiHR=phic
c
        endif
c
        FiTf(k)=phiHR
c
      enddo            
c
      return
      end
