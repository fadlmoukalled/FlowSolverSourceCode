c
c#############################################################################################
c
      SUBROUTINE MinusTwoThirdsRhoK(Variable)
c
C#############################################################################################
c
      use User0 
      use Turbulence1
      use Variables1, only: TurbulentKE,BTurbulentKE,
     *                      TurbulentKL,BTurbulentKL
      use Geometry4, only: Volume
      use PhysicalProperties1, only: Density, BDensity
      use Variables3, only: FluxTE
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Constants1, only: TwoThird
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer:: i,j
      character*10 Variable
      character*10 Variable2
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE Gradient(Variable,MethodCalcGradient,
     *         FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT,
     *         nIterGradientPhi,LimitGradient,LimitGradientMethod)
c--------------------------------------------------------------------------------
          character*10 Variable
          integer :: MethodCalcGradient,nIterGradientPhi
          logical :: LimitGradient
          integer :: LimitGradientMethod
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------------
        end SUBROUTINE Gradient
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c      
      if(.not.LSolveTurbulenceKineticEnergy) return
c
      if(Variable.eq.'velx') then
c
        if(TurbulenceModel.eq.'kklomega') then
c
          do i=1,NumberOfElements    
c
            rhok(i)=TwoThird*Density(i)*(TurbulentKE(i)+TurbulentKL(i))
c        
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Brhok(i,j)=TwoThird*BDensity(i,j)*
     *                           (BTurbulentKE(i,j)+BTurbulentKL(i,j))
c
            enddo
          enddo
c
        else
c
          do i=1,NumberOfElements    
c
            rhok(i)=TwoThird*Density(i)*TurbulentKE(i)
c        
          enddo
c
          do i=1,NumberOfBCSets
            do j=1,NBFaces(i)
c
              Brhok(i,j)=TwoThird*BDensity(i,j)*BTurbulentKE(i,j)
c
            enddo
          enddo
c
        endif
c
        Variable2='rhok'
        call Gradient(Variable2,MethodCalcGradientMomentum,
     *         rhok,drhokdx,drhokdy,drhokdz,Brhok,Bdrhokdx,Bdrhokdy,
     *            Bdrhokdz,nIterGradientMomentum, LimitGradientMomentum,
     *               LimitGradientMomentumMethod)
c
      endif
c
      if(Variable.eq.'velx') then
c
        do i=1,NumberOfElements    
c
          FluxTE(i)=FluxTE(i)+drhokdx(i)*Volume(i)
c        
        enddo
c
      elseif(Variable.eq.'vely') then
c
        do i=1,NumberOfElements    
c
          FluxTE(i)=FluxTE(i)+drhokdy(i)*Volume(i)
c        
        enddo
c
      elseif(Variable.eq.'velz') then
c
        do i=1,NumberOfElements    
c
          FluxTE(i)=FluxTE(i)+drhokdz(i)*Volume(i)
c        
        enddo
c
      endif
c
      return
      end
