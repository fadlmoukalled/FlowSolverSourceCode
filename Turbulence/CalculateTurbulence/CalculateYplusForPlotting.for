c
c#############################################################################################
c
      SUBROUTINE CalculateYplusPlotting
c
c#############################################################################################
      use User0, only: urff2Coefficient,urffmuCoefficient,
     *                 TurbulenceModel
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner,NBFacesmax
      use Geometry4, only: BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     xc,yc,zc,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz
      use BoundaryConditionsTurbulence2, 
     *      only: IwallTurbulence,IWallTurbulenceOwner,
     *            IWallTurbulenceNumberOfBCSets,IWallTurbulenceNBFaces
      use Turbulence1, only: fmuCoefficient,f2Coefficient,LTKE,LTED,
     *                       BfmuCoefficient,ReT,BReT,ustar,cmu25
      use Tecplot1, only: yplusPlot,uplusPlot,ByplusPlot,BuplusPlot
      use WallDistance1, only: iTau,jTau,BiTau,BjTau,WallDistance,
     *                         BWallDistance
      use Variables1, only: TurbulentKE,TurbulentED,BTurbulentKE,
     *                      BTurbulentED,uVelocity,vVelocity,wVelocity,
     *                      BuVelocity,BvVelocity,BwVelocity
      use PhysicalProperties1, only: Density,Viscosity,
     *                               BDensity,BViscosity,
     *                               BeDiffCoefficient,eDiffCoefficient
      use Constants1, only: tiny,big
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,k,i1,i2,i3
      double precision :: distance,nx,ny,nz,dNorm,uWall1,vWall1,wWall1,
     *                    WallVelocity,TauWall,uTau,dplus,f2Old,
     *                    fmuOld,dotproduct,BWallDistance1
c********************************************************************************************
      interface
c********************************************************************************************
        FUNCTION TangentialVelocity(i1)
c********************************************************************************************
          integer :: i1
          double precision :: TangentialVelocity
c********************************************************************************************
        end FUNCTION TangentialVelocity
c********************************************************************************************
        FUNCTION BTangentialVelocity(i1,i2)
c********************************************************************************************
          integer :: i1,i2
          double precision :: BTangentialVelocity
c********************************************************************************************
        end FUNCTION BTangentialVelocity
c********************************************************************************************
      end interface
c********************************************************************************************
c
      allocate(yplusPlot(NumberOfElements))
      allocate(uplusPlot(NumberOfElements))
      allocate(ByplusPlot(NumberOfBCSets,NBFacesMax))
      allocate(BuplusPlot(NumberOfBCSets,NBFacesMax))
c
      variable='velx'
      call CalculateEffectiveDiffusionCoefficient(Variable)
c
      do i=1,NumberOfElements
c
        i2=iTau(i)
        i3=jTau(i)
        i1=NBFaceOwner(i2,i3)
c
        dNorm=WallDistance(i1)
c
        WallVelocity=TangentialVelocity(i1)
        TauWall=BeDiffCoefficient(i2,i3)*WallVelocity/dNorm
c
        uTau=dsqrt(TauWall/BDensity(i2,i3))
        yplusPlot(i)=WallDistance(i)*
     *               BDensity(i2,i3)*uTau/BViscosity(i2,i3)
c
        WallVelocity=TangentialVelocity(i)
        uplusPlot(i)=WallVelocity/dmax1(uTau,tiny)
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          i2=BiTau(i,j)
          i3=BjTau(i,j)
          i1=NBFaceOwner(i2,i3)
c
          dNorm=WallDistance(i1)
c
          WallVelocity=TangentialVelocity(i1)
          TauWall=BeDiffCoefficient(i2,i3)*WallVelocity/dNorm
c
          uTau=dsqrt(TauWall/BDensity(i2,i3))
          ByplusPlot(i,j)=BWallDistance(i,j)*BDensity(i2,i3)*
     *                                   uTau/BViscosity(i2,i3)
c
          WallVelocity=BTangentialVelocity(i,j)
          BuplusPlot(i,j)=WallVelocity/dmax1(uTau,tiny)
c
        enddo
      enddo
c
      return
      end
c
c#############################################################################################
c
      SUBROUTINE CalculateYplus1(Pr1,Prp1,yplT,sig)
c
C#############################################################################################
c
      use Turbulence1, only: sigT,cappa,cc,ctrans
      use Constants1, only: big,tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: n0,i,index
      double precision :: Pr1,Prp1,yplT,sig
      double precision :: tol,a1,a2,d
      data n0,tol/200,1.e-5/
c********************************************************************************************
c
      a1=cappa*Pr1/sig
      a2=cappa*(cc+Prp1)
c      
      i=1
      index=0
      if(yplT.eq.0.) yplT=ctrans
c      
      do while(i.le.n0.and.index.eq.0)
c
        d=(dlog(yplT)-a1*yplT+a2)/(1./(yplT+tiny)-a1)
        yplT=yplT-d
c
        if(dabs(d).lt.tol) then
          index=1
        endif
c      
        i=i+1
        if(i.eq.n0.and.index.eq.0)  
     *         print*, 'Thermal yplus iterations diverged'              
c
      enddo
      if(i.ge.n0) print*,i,dabs(d)
c      
      return
      end          
c
c#############################################################################################
c
      SUBROUTINE Newton(ut,yplus1,dNorm,rho,visc)
c
C#############################################################################################
c
      use turbulence1, only: cappa,cc
      use Constants1, only: tiny
c********************************************************************************************
      implicit none
c********************************************************************************************
      double precision :: ut,yplus1,dNorm,rho,visc
      double precision :: term1,term2,term3,termp1,termp2,termp3,termp4,
     *                    termp5,dx
      double precision, save :: tolerance=1.e-5
      integer, save :: nIterNewton=100
      integer :: i
c********************************************************************************************
      yplus1=dmax1(yplus1,tiny)
      do i=1,nIterNewton
c
        term1=(ut/yplus1)**4+(ut/((dlog(yplus1))/cappa+cc))**4
        term2=(dNorm*rho/visc)*(term1**0.25)
        term3=yplus1-term2
c      
        termp1=0.25*(dNorm*rho/visc)*(term1**(-0.75))
        termp2=(ut**4)
        termp3=4.*termp2*(-4./(yplus1**5))
        termp4=termp2*(-4.*(1./(cappa*yplus1))/
     *          (((1./cappa)*dlog(yplus1)+cc)**5))
        termp5=1.-termp1*(termp3+termp4)
c
        dx=term3/termp5
        yplus1=yplus1-dx
        if(dabs(dx).lt.tolerance) return 
c       
      enddo
c
      print*,'number of yplus iterations exceeded'
      print*,'residuals=',dabs(dx)
c
      return
      end