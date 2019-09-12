c
c#############################################################################################
c
      SUBROUTINE CalculateCompressibilityCorrection
c
C#############################################################################################
c
      use User0, only: CompressibilityCorrectionMethod,
     *                 Zeman1,Zeman2
      use Geometry1, only: NumberOfElements
      use Turbulence1, only: XiStarFmt
      use PhysicalProperties1, only: SpecificHeat,RGas
      use Variables1, only: TurbulentKE,Temperature
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i
      double precision :: term1,term2,XiStar,Fmt,Mt,Mt0
c********************************************************************************************
c
      if(CompressibilityCorrectionMethod.eq.'sarkar') then
c
        XiStar=1.
c
        do i=1,NumberOfElements    
c
          term1=SpecificHeat(i)/(SpecificHeat(i)-RGas)
          term2=term1*RGas*Temperature(i)
          Fmt=2.*dmax1(TurbulentKE(i),0.)/dmax1(term2,tiny)
          XiStarFmt(i)=XiStar*Fmt
c
        enddo
c
      elseif(CompressibilityCorrectionMethod.eq.'zeman') then
c
        XiStar=0.75d0
c
        do i=1,NumberOfElements    
c
          term1=SpecificHeat(i)/(SpecificHeat(i)-RGas)
          Mt0=Zeman1*dsqrt(2./(dmax1(term1+1,tiny)))
          term2=term1*RGas*Temperature(i)
          Mt=dsqrt(2.*dmax1(TurbulentKE(i),0.)/dmax1(term2,tiny))
          XiStarFmt(i)=0.
c
          if(Mt.ge.Mt0) then
c
            term2=0.5*(term1+1.)*(Mt-Mt0)*(Mt-Mt0)/(Zeman2*Zeman2)
            Fmt=1.-dexp(-term2)
            XiStarFmt(i)=XiStar*Fmt
c
          endif
c
        enddo
c
      elseif(CompressibilityCorrectionMethod.eq.'wilcox') then
c
        XiStar=2.d0
        Mt0=0.0625d0   !this is mt0*mt0
c
        do i=1,NumberOfElements    
c
          term1=SpecificHeat(i)/(SpecificHeat(i)-RGas)
          term2=term1*RGas*Temperature(i)
          Mt=2.*dmax1(TurbulentKE(i),0.)/dmax1(term2,tiny)
          XiStarFmt(i)=0.
c
          if(Mt.ge.Mt0) then
c
            Fmt=Mt-Mt0
            XiStarFmt(i)=XiStar*Fmt
c
          endif
c
        enddo
c
      elseif(CompressibilityCorrectionMethod.eq.'nicoetala') then
c
        XiStar=1.
        Mt0=0.1
c
        do i=1,NumberOfElements    
c
          term1=SpecificHeat(i)/(SpecificHeat(i)-RGas)
          term2=term1*RGas*Temperature(i)
          Mt=2.*dmax1(TurbulentKE(i),0.)/dmax1(term2,tiny)
          Fmt=Mt
          XiStarFmt(i)=XiStar*Fmt
c
          if(Mt.gt.Mt0) then
c
            Fmt=0.1
            XiStarFmt(i)=XiStar*Fmt
c
          endif
c
        enddo
c
      elseif(CompressibilityCorrectionMethod.eq.'nicoetalb') then
c
        XiStar=1.
        Mt0=0.1
c
        do i=1,NumberOfElements    
c
          term1=SpecificHeat(i)/(SpecificHeat(i)-RGas)
          term2=term1*RGas*Temperature(i)
          Mt=2.*dmax1(TurbulentKE(i),0.)/dmax1(term2,tiny)
          Fmt=10*Mt*Mt
          XiStarFmt(i)=XiStar*Fmt
c
          if(Mt.gt.Mt0) then
c
            Fmt=0.1
            XiStarFmt(i)=XiStar*Fmt
c
          endif
c
        enddo
c
      endif
c
      return
      end
