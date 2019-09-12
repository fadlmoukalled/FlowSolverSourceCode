c
C#############################################################################################
      FUNCTION Flength1(x)
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision :: x
      double precision :: FLength1
c*********************************************************************************************
c
      if(x.lt.400.) then
        FLength1=39.8189-(119.270d-4+132.567d-6*x)*x
      elseif(x.ge.400..and.x.lt.596.) then 
        FLength1=263.404+(-123.939d-2+(194.548d-5-101.695d-8*x)*x)*x
      elseif(x.ge.596..and.x.lt.1200.) then 
        FLength1=0.5-3.d-4*(x-596.)
      elseif(x.ge.1200.) then 
        FLength1=0.3188
      endif
c
      end FUNCTION FLength1
c
C#############################################################################################
      FUNCTION ReThetaC(x)
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision :: x
      double precision :: ReThetaC
c*********************************************************************************************
c
      if(x.le.1870.) then
        ReThetaC=-396.035d-2+(10120.656d-4+
     *        (-868.23d-6+(696.506d-9-174.105d-12*x)*x)*x)*x
      elseif(x.gt.1870.) then 
        ReThetaC=x-(593.11+0.482*(x-1870.))
      endif
c
      end FUNCTION ReThetaC
c
C#############################################################################################
      FUNCTION ReThetaTeqPrime(x,y,z)
!     x=Tu    y=LambdaTheta  z:dLambdadTheta
C#############################################################################################
      use Constants1, only: tiny
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision :: x,y,z
      double precision :: ReThetaTeqPrime
c*********************************************************************************************
      interface
c*********************************************************************************************
        FUNCTION FPrimeLambdaTheta(x,y,z)
          double precision :: x,y,z
          double precision :: FPrimeLambdaTheta
        end FUNCTION FPrimeLambdaTheta
      end interface
c*********************************************************************************************
c
      if(x.le.1.3) then
        ReThetaTeqPrime=(1173.51-589.428*x+
     *              0.2196/dmax1(x*x,tiny))*FPrimeLambdaTheta(x,y,z)
      elseif(x.gt.1.3) then 
        ReThetaTeqPrime=331.5*((x-0.5658)**(-0.671))*
     *                                      FPrimeLambdaTheta(x,y,z)
      endif
c
      end FUNCTION ReThetaTeqPrime
c
C#############################################################################################
      FUNCTION ReThetaTeq(x,y)
!     x=Tu    y=LambdaTheta
C#############################################################################################
      use Constants1, only: tiny
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision :: x,y
      double precision :: ReThetaTeq
c*********************************************************************************************
      interface
c*********************************************************************************************
        FUNCTION FLambdaTheta(x,y)
          double precision :: x,y
          double precision :: FLambdaTheta
        end FUNCTION FLambdaTheta
      end interface
c*********************************************************************************************
c
      if(x.le.1.3) then
        ReThetaTeq=
     *    (1173.51-589.428*x+0.2196/dmax1(x*x,tiny))*FLambdaTheta(x,y)
      elseif(x.gt.1.3) then 
        ReThetaTeq=331.5*((x-0.5658)**(-0.671))*FLambdaTheta(x,y)
      endif
c
      ReThetaTeq=dmax1(ReThetaTeq,20.) 
c
      end FUNCTION ReThetaTeq
c
C#############################################################################################
      FUNCTION FPrimeLambdaTheta(x,y,z)
!     x=Tu    y=LambdaTheta   z:dLambdadTheta
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision :: x,y,z
      double precision :: FPrimeLambdaTheta
c*********************************************************************************************
c
      if(y.le.0.) then
        FPrimeLambdaTheta=
     *           z*(12.986+(247.32+1217.067*y)*y)*dexp(-(x/1.5)**1.5)
      elseif(y.gt.0.) then 
        FPrimeLambdaTheta=9.625*z*dexp(-35.*y)*dexp(-(x/0.5))
      endif
c
      end FUNCTION FPrimeLambdaTheta
c
C#############################################################################################
      FUNCTION FLambdaTheta(x,y)
!     x=Tu    y=LambdaTheta
C#############################################################################################
      implicit none
c*********************************************************************************************
      double precision :: x,y
      double precision :: FLambdaTheta
c*********************************************************************************************
c
      if(y.le.0.) then
        FLambdaTheta=
     *        1.+((12.986+(123.66+405.689*y)*y)*y)*dexp(-(x/1.5)**1.5)
      elseif(y.gt.0.) then 
        FLambdaTheta=1.+0.275*(1.-dexp(-35.*y))*dexp(-(x/0.5))
      endif
c
      end FUNCTION FLambdaTheta