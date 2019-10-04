c
c#############################################################################################
c
      SUBROUTINE NVFLimiter(ConvectionScheme,phiTeldaC,phiTeldaf)
c
c#############################################################################################
      use User0, only: dt,BettaConvection
      use VolumeOfFluid2, only: i,j,k
      use VolumeOfFluid1, only: cosThetaF
      use Geometry4, only: Volume,GFactorCF
      use Variables1, only: mdot
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*20 ConvectionScheme
      double precision :: phiTeldaC,phiTeldaf
      double precision :: phiSuperB,phiSMART,factor,phiHRIC,Volumef,
     *                    courantNumber
c********************************************************************************************
c     NVF implementation of convection schemes
c********************************************************************************************
c     upwind    (bounded)
c     downwind  (not bounded)
c     smart     (bounded)
c     minmod    (bounded)
c     boundedcd (bounded)
c     cd        (unbounded)
c     osher     (bounded)
c     stoic     (bounded)
c     muscl     (bounded)
c     superbee  (bounded)
c     gamma     (bounded)
c     stacs     (bounded)
c     hric      (bounded)
c     cicsam    (bounded)
c--------------------------------------------------------------------------------------------
c
      if(ConvectionScheme.eq.'upwind') then
c
        phiTeldaf=phiTeldaC
        return
c
      elseif(ConvectionScheme.eq.'smart') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.ge.0..AND.phiTeldaC.le.1.) 
     *    phiTeldaf=dmin1(3.*phiTeldaC,0.375+0.75*phiTeldaC,
     *            (2.+phiTeldaC)/3.)
        return
c
      elseif(ConvectionScheme.eq.'minmod') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.ge.0..AND.phiTeldaC.le.1.) 
     *   phiTeldaf=dmin1(1.5*phiTeldaC,0.5*phiTeldaC+0.5)
c
        return
c
      elseif(ConvectionScheme.eq.'boundedcd') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.ge.0..AND.phiTeldaC.le.1.)
     *                       phiTeldaf=0.5*phiTeldaC+0.5
        return
c
      elseif(ConvectionScheme.eq.'cd') then
c
        phiTeldaf=0.5*phiTeldaC+0.5
        return
c
      elseif(ConvectionScheme.eq.'osher') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.ge.0..AND.phiTeldaC.le.1.)
     *     phiTeldaf=dmin1(1.5*phiTeldaC,1.)
        return
c
      elseif(ConvectionScheme.eq.'stoic') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.ge.0..AND.phiTeldaC.le.1.) 
     *    phiTeldaf=dmin1(0.5*phiTeldaC+0.5,0.375+0.75*phiTeldaC,
     *            (2.+phiTeldaC)/3.)
        return
c
      elseif(ConvectionScheme.eq.'muscl') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.ge.0..AND.phiTeldaC.le.1.) 
     *    phiTeldaf=dmin1(2.*phiTeldaC,0.25+phiTeldaC,1.)
        return
c
      elseif(ConvectionScheme.eq.'superbee') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.ge.0..AND.phiTeldaC.le.1.) 
     *    phiTeldaf=dmin1(2.*phiTeldaC,0.5*phiTeldaC+0.5,
     *                                            1.5*phiTeldaC,1.)
        return
c
      elseif(ConvectionScheme.eq.'downwind') then
c
        phiTeldaf=1.
c
        return
c
      elseif(ConvectionScheme.eq.'gamma') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.gt.0..and.phiTeldaC.lt.BettaConvection) then
          factor=0.5/BettaConvection
          phiTeldaf=(-factor*phiTeldaC+(1.+factor))*phiTeldaC
        elseif(phiTeldaC.ge.BettaConvection.and.phiTeldaC.lt.1.) then
          phiTeldaf=0.5+0.5*phiTeldaC
        endif
c
        return
c
      elseif(ConvectionScheme.eq.'stacs') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.gt.0..and.phiTeldaC.lt.1.) 
     *               phiTeldaf=dmin1(2.*phiTeldaC,0.85+0.15*phiTeldaC)
        phiSuperB=phiTeldaf
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.gt.0..and.phiTeldaC.lt.1.) 
     *         phiTeldaf=dmin1(dmin1(3.*phiTeldaC,
     *                  0.375+0.75*phiTeldaC),0.9+0.1*phiTeldaC)
        phiSMART=phiTeldaf
c
        factor=cosThetaF(k)*cosThetaF(k)*cosThetaF(k)*cosThetaF(k)
c
        phiTeldaf=phiSuperB*factor+phiSMART*(1.-factor)
c
        return
c
      elseif(ConvectionScheme.eq.'hric') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.gt.0..and.phiTeldaC.lt.1.) 
     *               phiTeldaf=dmin1(2.*phiTeldaC,0.85+0.15*phiTeldaC)
        phiSuperB=phiTeldaf
c
        factor=dsqrt(cosThetaF(k))
        phiHRIC=phiSuperB*factor+phiTeldaC*(1.-factor)
c
        Volumef=GFactorCF(k)*Volume(i)+(1.-GFactorCF(k))*Volume(j)
        courantNumber=dabs(mdot(k)*dt/Volumef)
c
        if(courantNumber.lt.0.3) then
c
          phiTeldaf=phiHRIC
c
        elseif(courantNumber.lt.0.7) then
c
          phiTeldaf=phiTeldaC+
     *              (phiHRIC-phiTeldaC)*(0.7-courantNumber)/0.4
c
        else
c
          phiTeldaf=phiTeldaC
c
        endif
c
        return
c
      elseif(ConvectionScheme.eq.'cicsam') then
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.gt.0..and.phiTeldaC.lt.1.) 
     *               phiTeldaf=dmin1(2.*phiTeldaC,0.85+0.15*phiTeldaC)
        phiSuperB=phiTeldaf
c
        phiTeldaf=phiTeldaC
        if(phiTeldaC.gt.0..and.phiTeldaC.lt.1.) 
     *         phiTeldaf=dmin1(dmin1(3.*phiTeldaC,
     *                  0.375+0.75*phiTeldaC),0.9+0.1*phiTeldaC)
        phiSMART=phiTeldaf
c
        factor=cosThetaF(k)
        Volumef=GFactorCF(k)*Volume(i)+(1.-GFactorCF(k))*Volume(j)
        courantNumber=dabs(mdot(k)*dt/Volumef)
c
        phiTeldaf=(courantNumber*phiTeldaC+
     *           (1.-courantNumber)*phiSuperB)*factor+
     *                 (courantNumber*phiTeldaC+(1.-courantNumber)*
     *                                           phiSMART)*(1.-factor)
c
        return
c        
      else
c
        print*, 'convection scheme ',ConvectionScheme,' not implemented'
        print*, 'program has stoped'
        stop
c
      endif
c
      return
      end