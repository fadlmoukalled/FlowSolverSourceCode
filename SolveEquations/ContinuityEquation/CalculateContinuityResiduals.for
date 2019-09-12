c
c#############################################################################################
c
      SUBROUTINE CalculateResidualsContinuity
c
c#############################################################################################
c
      use User0, only: LUnsteady,LFreeSurfaceFlow
      use Residuals1
      use BoundaryConditions1, only: BCType
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NIFaces,NBFaces,
     *                     NIFaceOwner,NIFaceNeighbor,NBFaceOwner
      use Geometry4, only: Volume
      use Variables1, only: mdot,Bmdot
      use PhysicalProperties1, only: Density,DensityOld,DensityOldOld,
     *                               Densityf,BDensity
c********************************************************************************************
      implicit none
c********************************************************************************************
      double precision :: MaximumImbalance,rmsImbalance,AbsImbalance
      double precision :: mdotin
      integer :: i,j,k
c
      double precision, save, dimension(:), allocatable :: massImbalance
      double precision, save, dimension(:), allocatable :: drhodt
c********************************************************************************************
      interface
c********************************************************************************************
        SUBROUTINE ComputeTransientGradient(phi,phiOld,phiOldOld,dphidt)
c--------------------------------------------------------------
          double precision, dimension(:) :: phi,phiOld,phiOldOld,dphidt
c--------------------------------------------------------------
        end SUBROUTINE ComputeTransientGradient
c--------------------------------------------------------------
      end interface
c--------------------------------------------------------------
c
      if(LFreeSurfaceFlow) then
c
        allocate(massImbalance(NumberOfElements))
        massImbalance=0.
c      
        ResorMax(5)=0.
        ResorAbs(5)=0.
        ResorRMS(5)=0.
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          massImbalance(i)=massImbalance(i)+mdot(k)/Densityf(k)    
          massImbalance(j)=massImbalance(j)-mdot(k)/Densityf(k)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            massImbalance(k)=massImbalance(k)+Bmdot(i,j)/BDensity(i,j)
c
          enddo  
        enddo  
c
        MaximumImbalance=0.
        OverallImbalance=0.
        AbsImbalance=0.
        rmsImbalance=0.
c
        do i=1,NumberOfElements
c
          OverallImbalance=OverallImbalance+massImbalance(i)          
          AbsImbalance=AbsImbalance+dabs(massImbalance(i))          
          MaximumImbalance=
     *                dmax1(MaximumImbalance,dabs(massImbalance(i)))
          if(MaximumImbalance.eq.dabs(massImbalance(i))) iMaxImbalance=i
          rmsImbalance=rmsImbalance+massImbalance(i)*massImbalance(i)
c
        enddo
c
        ResorMax(5)=MaximumImbalance
        ResorAbs(5)=AbsImbalance
        ResorRMS(5)=dsqrt(rmsImbalance/NumberOfElements)      
c
        mdotIn=0.
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            if(BCType(i,j).eq.'inlet') then
c              
              mdotIn=mdotIn+Bmdot(i,j)/BDensity(i,j)
c
            endif
c
            if(BCType(i,j).eq.'pressurefarfield') then
c              
              if(Bmdot(i,j).lt.0.) then
c
                mdotIn=mdotIn+Bmdot(i,j)/BDensity(i,j)
c
              endif
c
            endif
c
            if(BCType(i,j).eq.'periodic') then
c              
              if(Bmdot(i,j).lt.0.) then
c
                mdotIn=mdotIn+Bmdot(i,j)/BDensity(i,j)
c
              endif
c
            endif
c
          enddo
        enddo
c
        if(mdotIn.eq.0.) then
          ResorScaled(5)=0.
        else
          ResorScaled(5)=ResorMax(5)/dabs(mdotIn)
        endif
c
        deallocate(massImbalance)
c
      else
c
        allocate(massImbalance(NumberOfElements))
        massImbalance=0.
c      
        ResorMax(5)=0.
        ResorAbs(5)=0.
        ResorRMS(5)=0.
c
        do k=1,NIFaces
c
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          massImbalance(i)=massImbalance(i)+mdot(k)    
          massImbalance(j)=massImbalance(j)-mdot(k)
c
        enddo
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            k=NBFaceOwner(i,j)
c
            massImbalance(k)=massImbalance(k)+Bmdot(i,j)    
c
          enddo  
        enddo  
c
        if(LUnsteady) then
c
          allocate(drhodt(NumberOfElements))
          call ComputeTransientGradient
     *            (Density,DensityOld,DensityOldOld,drhodt)
c
          do i=1,NumberOfElements
c
            massImbalance(i)=massImbalance(i)+drhodt(i)*Volume(i) 
c
          enddo
c
          deallocate(drhodt)
c
        endif
c
        MaximumImbalance=0.
        OverallImbalance=0.
        AbsImbalance=0.
        rmsImbalance=0.
c
        do i=1,NumberOfElements
c
          OverallImbalance=OverallImbalance+massImbalance(i)          
          AbsImbalance=AbsImbalance+dabs(massImbalance(i))          
          MaximumImbalance=
     *                dmax1(MaximumImbalance,dabs(massImbalance(i)))
          if(MaximumImbalance.eq.dabs(massImbalance(i))) iMaxImbalance=i
          rmsImbalance=rmsImbalance+massImbalance(i)*massImbalance(i)
c
        enddo
c
        ResorMax(5)=MaximumImbalance
        ResorAbs(5)=AbsImbalance
        ResorRMS(5)=dsqrt(rmsImbalance/NumberOfElements)      
c
        mdotIn=0.
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c
            if(BCType(i,j).eq.'inlet') then
c              
              mdotIn=mdotIn+Bmdot(i,j)
c
            endif
c
            if(BCType(i,j).eq.'pressurefarfield') then
c              
              if(Bmdot(i,j).lt.0.) then
c
                mdotIn=mdotIn+Bmdot(i,j)
c
              endif
c
            endif
c
            if(BCType(i,j).eq.'periodic') then
c              
              if(Bmdot(i,j).lt.0.) then
c
                mdotIn=mdotIn+Bmdot(i,j)
c
              endif
c
            endif
c
          enddo
        enddo
c
        if(mdotIn.eq.0.) then
          ResorScaled(5)=0.
        else
          ResorScaled(5)=ResorMax(5)/dabs(mdotIn)
        endif
c
        deallocate(massImbalance)
c
      endif
c
      return
      end