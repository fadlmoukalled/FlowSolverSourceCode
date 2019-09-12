c
C#############################################################################################
c
      SUBROUTINE TVDLimiter(ConvectionScheme,rf,psirf)
c
C#############################################################################################
      use User0, only: dt
      use VolumeOfFluid2, only: i,j,k
      use VolumeOfFluid1, only: cosThetaF
      use Geometry4, only: Volume,GFactorCF
      use Variables1, only: mdot
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*20 ConvectionScheme
      double precision rf,psirf,factor1,factor2,factor
      double precision psiSuperB,psiSMART,Volumef,courantNumber
c********************************************************************************************
c     TVD implementation of convection schemes
c********************************************************************************************
c     upwind   (bounded)
c     downwind (not bounded)
c     fromm    (not bounded)
c     sou      (not bounded)
c     cd       (not bounded)
c     quick    (not bounded)
c     hquick   (harmonic based on quick) 
c     superbee (bounded)
c     minmod   (bounded)
c     osher    (bounded)
c     vanleer  (bounded)
c     muscl    (bounded)
c     smart    (bounded)
c     umist    (bounded)
c     charm    (bounded)
c--------------------------------------------------------------------------------------------
c
      if(ConvectionScheme.eq.'upwind') then
c
        psirf=0.
        return
c
      elseif(ConvectionScheme.eq.'downwind') then
c
        psirf=2.
        return
c
      elseif(ConvectionScheme.eq.'fromm') then
c
        psirf=0.5*(1.+rf)
        return
c
      elseif(ConvectionScheme.eq.'sou') then
c
        psirf=rf
        return
c
      elseif(ConvectionScheme.eq.'cd') then
c
        psirf=1.
        return
c
      elseif(ConvectionScheme.eq.'quick') then
c
        psirf=(3.+rf)/4.
        return
c
      elseif(ConvectionScheme.eq.'superbee') then
c
        psirf=dmax1(0.,dmin1(1.,2.*rf),dmin1(2.,rf))
        return
c
      elseif(ConvectionScheme.eq.'minmod') then
c
        psirf=dmax1(0.,dmin1(1.,rf))
        return
c
      elseif(ConvectionScheme.eq.'osher') then
c
        psirf=dmax1(0.,dmin1(2.,rf))
        return
c
      elseif(ConvectionScheme.eq.'vanleer') then
c
        psirf=(rf+dabs(rf))/(1.+dabs(rf))
        return
c
      elseif(ConvectionScheme.eq.'muscl') then
c
        psirf=dmax1(0.,dmin1(2.*rf,0.5*(1.+rf),2.))
        return
c
      elseif(ConvectionScheme.eq.'smart') then
c
        psirf=dmax1(0.,dmin1(2.*rf,0.75*rf+0.25,4.))
        return
c
      elseif(ConvectionScheme.eq.'hquick') then
c
        psirf=2.*(rf+dabs(rf))/(rf+3.)
        return
c
      elseif(ConvectionScheme.eq.'umist') then
c
        psirf=dmax1(0.,dmin1(2.*rf,0.25+0.75*rf,0.75+0.25*rf,2.))
        return
c
      elseif(ConvectionScheme.eq.'charm') then
c
        if(rf.le.0.) then
          psirf=0.
        else
          psirf=rf*(3.*rf+1.)/((rf+1.)**2)
        endif
        return
c
      elseif(ConvectionScheme.eq.'stacs') then
c
        psirf=dmax1(0.,dmin1(1.,2.*rf),dmin1(2.,rf))
        psiSuperB=psirf
c        
        psirf=dmax1(0.,dmin1(2.*rf,0.5*(1.+rf),2.))
        psiSMART=psirf
c
        factor=cosThetaF(k)*cosThetaF(k)*cosThetaF(k)*cosThetaF(k)
c
        psirf=psiSuperB*factor+psiSMART*(1.-factor)
c
        return
c
      elseif(ConvectionScheme.eq.'hric') then
c
        psirf=dmax1(0.,dmin1(1.,2.*rf),dmin1(2.,rf))
        psiSuperB=psirf
c
        Volumef=GFactorCF(k)*Volume(i)+(1.-GFactorCF(k))*Volume(j)
        courantNumber=dabs(mdot(k)*dt/Volumef)
c
        if(courantNumber.lt.0.3) then
c
          factor=dsqrt(cosThetaF(k))
c
        elseif(courantNumber.lt.0.7) then
c
          factor=dsqrt(cosThetaF(k))*(0.7-courantNumber)/0.4
c
        else
c
          factor=0.
c
        endif
c
        psirf=psiSuperB*factor
c
        return
c
      elseif(ConvectionScheme.eq.'cicsam') then
c
        psirf=dmax1(0.,dmin1(1.,2.*rf),dmin1(2.,rf))
        psiSuperB=psirf
c        
        psirf=dmax1(0.,dmin1(2.*rf,0.5*(1.+rf),2.))
        psiSMART=psirf
c
        Volumef=GFactorCF(k)*Volume(i)+(1.-GFactorCF(k))*Volume(j)
        courantNumber=dabs(mdot(k)*dt/Volumef)
c
        factor1=1.-courantNumber
        factor2=cosThetaF(k)
c
        psirf=factor1*(factor2*psiSuperB+(1.-factor2)*psiSMART)
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