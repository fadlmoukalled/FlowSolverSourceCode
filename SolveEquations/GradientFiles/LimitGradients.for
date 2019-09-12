c
c#############################################################################################
      SUBROUTINE LimitGradientVariable(FiT,dfidxT,dfidyT,dfidzT,
     *                 BFiT,BdfidxT,BdfidyT,BdfidzT,LimitGradientMethod)
c#############################################################################################
      implicit none
c*********************************************************************************************
      integer :: LimitGradientMethod
c
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c*********************************************************************************************
      interface
c*********************************************************************************************
        SUBROUTINE LimitGradientMethod1(FiT,dfidxT,dfidyT,dfidzT,
     *                                  BFiT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE LimitGradientMethod1
c--------------------------------------------------------------------------
        SUBROUTINE LimitGradientMethod2(FiT,dfidxT,dfidyT,dfidzT,
     *                                  BFiT,BdfidxT,BdfidyT,BdfidzT)
c--------------------------------------------------------------------------
          double precision, dimension(:) :: FiT
          double precision, dimension(:) :: dfidxT
          double precision, dimension(:) :: dfidyT
          double precision, dimension(:) :: dfidzT
          double precision, dimension(:,:) :: BFiT
          double precision, dimension(:,:) :: BdfidxT
          double precision, dimension(:,:) :: BdfidyT
          double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
        end SUBROUTINE LimitGradientMethod2
c*********************************************************************************************
      end interface
c*********************************************************************************************
c
      if(LimitGradientMethod.eq.1) then
c
        call LimitGradientMethod1(FiT,dfidxT,dfidyT,dfidzT,
     *                            BFiT,BdfidxT,BdfidyT,BdfidzT)
c
      elseif(LimitGradientMethod.eq.2) then
c
        call LimitGradientMethod2(FiT,dfidxT,dfidyT,dfidzT,
     *                            BFiT,BdfidxT,BdfidyT,BdfidzT)
c
      endif
c
      return
      end
c
c#############################################################################################
      SUBROUTINE LimitGradientMethod1(FiT,dfidxT,dfidyT,dfidzT,
     *                                BFiT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
      use constants1, only: big,tiny
      use Geometry1, only: NumberOfElements,NumberOfBCSets,
     *                     NumbOfElementFaces
      use Geometry3, only: NBFaces,NGlobalEFaces,ElementNeighbor,
     *                     NBFaceOwner,NumberofElementNeighbors
      use Geometry4, only: FaceCentroidx,FaceCentroidy,FaceCentroidz,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz,
     *                     xc,yc,zc
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k,kE,kf
      double precision :: rk,deltaMinMax,phif
c
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c
      double precision, dimension(:), allocatable :: phiminL
      double precision, dimension(:), allocatable :: phimaxL
      double precision, dimension(:), allocatable :: deltaMin
      double precision, dimension(:), allocatable :: deltaMax
      double precision, dimension(:), allocatable :: cellLimiter
c
      data rk/1./
c
c*********************************************************************************************
c
      allocate(phiminL(NumberOfElements))
      allocate(phimaxL(NumberOfElements))
      allocate(deltaMin(NumberOfElements))
      allocate(deltaMax(NumberOfElements))
      allocate(cellLimiter(NumberOfElements))
c
      phiminL= big
      phimaxL=-big
      cellLimiter=1.
c
      do i=1,NumberOfElements
c
        phiminL(i)=dmin1(phiminL(i),FiT(i))
        phimaxL(i)=dmax1(phimaxL(i),FiT(i))
c
      enddo
c
c--- Interior faces
c
      do i=1,NumberOfElements
        do j=1,NumberofElementNeighbors(i)
c
          kE=ElementNeighbor(i,j)
c
          if(kE.ne.0) then
c
            phiminL(i)=dmin1(phiminL(i),FiT(KE))
            phimaxL(i)=dmax1(phimaxL(i),FiT(KE))
c
          endif
c
        enddo
      enddo
c
c--- Boundary faces
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          kE=NBFaceOwner(i,j)
c
          phiminL(kE)=dmin1(phiminL(kE),BFiT(i,j))
          phimaxL(kE)=dmax1(phimaxL(kE),BFiT(i,j))
c
        enddo
      enddo
c
c--- All elements
c
      do i=1,NumberOfElements
c
        deltaMax(i)=phimaxL(i)-FiT(i)
        deltaMin(i)=phiminL(i)-FiT(i)
c
      enddo
c
      if(rk.gt.0.and.rk.lt.1.) then
c
        do i=1,NumberOfElements
c
          deltaMinMax=(1./rk-1.)*(phimaxL(i)-phiminL(i))
          deltaMax(i)=deltaMax(i)+deltaMinMax
          deltaMin(i)=deltaMin(i)-deltaMinMax
c
        enddo
c
      endif
c
c--- Interior faces
c
      do i=1,NumberOfElements
        do j=1,NumberofElementNeighbors(i)
c
          kE=ElementNeighbor(i,j)
          kf=NGlobalEFaces(i,j)
c
          if(kE.ne.0) then
c
            phif=dfidxT(i)*(FaceCentroidx(kf)-xc(i))+
     *               dfidyT(i)*(FaceCentroidy(kf)-yc(i))+
     *                     dfidzT(i)*(FaceCentroidz(kf)-zc(i))
c
            if(phif.gt.(deltaMax(i)+tiny)) then
c
              cellLimiter(i)=dmin1(cellLimiter(i),
     *                         deltaMax(i)/(phif+tiny))
c
            elseif(phif.lt.(deltaMin(i)-tiny)) then
c
              cellLimiter(i)=dmin1(cellLimiter(i),
     *                         deltaMin(i)/(phif+tiny))
c
            endif
c
          endif
c
        enddo
      enddo
c
c--- Boundary faces
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          kE=NBFaceOwner(i,j)
c
            phif=BdfidxT(i,j)*(BFaceCentroidx(i,j)-xc(kE))+
     *               BdfidyT(i,j)*(BFaceCentroidy(i,j)-yc(kE))+
     *                    BdfidzT(i,j)*(BFaceCentroidz(i,j)-zc(kE))
c
            if(phif.gt.(deltaMax(kE)+tiny)) then
c
              cellLimiter(kE)=dmin1(cellLimiter(kE),
     *                         deltaMax(kE)/(phif+tiny))
c
            elseif(phif.lt.(deltaMin(kE)-tiny)) then
c
              cellLimiter(kE)=dmin1(cellLimiter(kE),
     *                         deltaMin(kE)/(phif+tiny))
c
            endif
c
        enddo
      enddo
c
c--- Limit the gradient at cell centers
c
      do i=1,NumberOfElements
c      
        dfidxT(i)=cellLimiter(i)*dfidxT(i)     
        dfidyT(i)=cellLimiter(i)*dfidyT(i)     
        dfidzT(i)=cellLimiter(i)*dfidzT(i)     
c      
      enddo
c
c--- Update Gradients along boundaries
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          BdfidxT(i,j)=dfidxT(k)
          BdfidyT(i,j)=dfidyT(k)
          BdfidzT(i,j)=dfidzT(k)
c
        enddo
      enddo
c
      deallocate(phiminL)
      deallocate(phimaxL)
      deallocate(deltaMax)
      deallocate(deltaMin)
      deallocate(cellLimiter)
c
      return
      end
c
c#############################################################################################
      SUBROUTINE LimitGradientMethod2(FiT,dfidxT,dfidyT,dfidzT,
     *                                BFiT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################c
      use constants1, only: epsilon,big,tiny
      use Geometry1, only: NumberOfElements,NumberOfBCSets,
     *                     NumbOfElementFaces
      use Geometry3, only: NBFaces,NGlobalEFaces,ElementNeighbor,
     *                     NBFaceOwner,NumberofElementNeighbors
      use Geometry4, only: FaceCentroidx,FaceCentroidy,FaceCentroidz,
     *                     BFaceCentroidx,BFaceCentroidy,BFaceCentroidz,
     *                     xc,yc,zc
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j,k,kE,kf
      double precision :: phiMin,phiMax,phif,eps2,ztemp
      double precision :: deltaP,deltaM,DeltaP2,DeltaM2,part1,part2
c
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c--------------------------------------------------------------------------
      double precision, dimension(:), allocatable :: phiminL
      double precision, dimension(:), allocatable :: phimaxL
      double precision, dimension(:), allocatable :: cellLimiter
c*********************************************************************************************
      interface
c*********************************************************************************************
        SUBROUTINE PhiFieldMinMax(FiT,BFiT,phiMin,phiMax)
c--------------------------------------------------------------------------
          double precision :: phiMin,phiMax
          double precision, dimension(:) :: FiT
          double precision, dimension(:,:) :: BFiT
c--------------------------------------------------------------------------
        end SUBROUTINE PhiFieldMinMax
c*********************************************************************************************
      end interface
c*********************************************************************************************
c
      call PhiFieldMinMax(FiT,BFiT,phiMin,phiMax)
c
      eps2=(epsilon*(phiMax-phiMin))**2
c
      allocate(phiminL(NumberOfElements))
      allocate(phimaxL(NumberOfElements))
      allocate(cellLimiter(NumberOfElements))
c
      phiminL= big
      phimaxL=-big
      cellLimiter=1.e30
c
      do i=1,NumberOfElements
c
        phiminL(i)=dmin1(phiminL(i),FiT(i))
        phimaxL(i)=dmax1(phimaxL(i),FiT(i))
c
      enddo
c
c--- Interior faces
c
      do i=1,NumberOfElements
        do j=1,NumberofElementNeighbors(i)
c
          kE=ElementNeighbor(i,j)
c
          if(kE.ne.0) then
c
            phiminL(i)=dmin1(phiminL(i),FiT(KE))
            phimaxL(i)=dmax1(phimaxL(i),FiT(KE))
c
          endif
c
        enddo
      enddo
c
c--- Boundary faces
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          kE=NBFaceOwner(i,j)
c
          phiminL(kE)=dmin1(phiminL(kE),BFiT(i,j))
          phimaxL(kE)=dmax1(phimaxL(kE),BFiT(i,j))
c
        enddo
      enddo
c
c--- Interior faces
c
      do i=1,NumberOfElements
        do j=1,NumberofElementNeighbors(i)
c
          kE=ElementNeighbor(i,j)
          kf=NGlobalEFaces(i,j)
c
          if(kE.ne.0) then
c
            deltaM=dfidxT(i)*(FaceCentroidx(kf)-xc(i))+
     *               dfidyT(i)*(FaceCentroidy(kf)-yc(i))+
     *                     dfidzT(i)*(FaceCentroidz(kf)-zc(i))
c
            if(deltaM.ge.0.) then
c
              deltaP=phimaxL(i)-FiT(i)
c
            else
c
              deltaP=phiminL(i)-FiT(i)
c
            endif
c
            deltaP2=deltaP*deltaP
            deltaM2=deltaM*deltaM
            part1=(deltaP2+eps2)*deltaM+2.*deltaM2*deltaP
            part2=deltaP2+2.*deltaM2+deltaP*deltaM+eps2
            Ztemp=1.
            phif=part1/(dsign(ztemp,deltaM)*(dabs(deltaM)+tiny)*Part2)
c
            cellLimiter(i)=dmin1(1.,phif,cellLimiter(i))
c
          endif
c
        enddo
      enddo
c
c--- Boundary faces
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          kE=NBFaceOwner(i,j)
c
            deltaM=BdfidxT(i,j)*(BFaceCentroidx(i,j)-xc(kE))+
     *               BdfidyT(i,j)*(BFaceCentroidy(i,j)-yc(kE))+
     *                  BdfidzT(i,j)*(BFaceCentroidz(i,j)-zc(kE))
c
            if(deltaM.ge.0.) then
c
              deltaP=phimaxL(kE)-BFiT(i,j)
c
            else
c
              deltaP=phiminL(kE)-BFiT(i,j)
c
            endif
c
            deltaP2=deltaP*deltaP
            deltaM2=deltaM*deltaM
            part1=(deltaP2+eps2)*deltaM+2.*deltaM2*deltaP
            part2=deltaP2+2.*deltaM2+deltaP*deltaM+eps2
            Ztemp=1.
            phif=part1/(dsign(ztemp,deltaM)*(dabs(deltaM)+tiny)*Part2)
c
            cellLimiter(kE)=dmin1(1.,phif,cellLimiter(kE))
c
        enddo
      enddo
c
c--- Limit the gradient at cell centers
c
      do i=1,NumberOfElements
c      
        dfidxT(i)=cellLimiter(i)*dfidxT(i)     
        dfidyT(i)=cellLimiter(i)*dfidyT(i)     
        dfidzT(i)=cellLimiter(i)*dfidzT(i)     
c      
      enddo
c
c--- Update Gradients along boundaries
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          BdfidxT(i,j)=dfidxT(k)
          BdfidyT(i,j)=dfidyT(k)
          BdfidzT(i,j)=dfidzT(k)
c
        enddo
      enddo
c
      deallocate(phiminL)
      deallocate(phimaxL)
      deallocate(cellLimiter)
c
      return
      end
c
c#############################################################################################
      SUBROUTINE PhiFieldMinMax(FiT,BFiT,phiMin,phiMax)
c#############################################################################################
      use constants1, only: big
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,j
      double precision :: phiMin,phiMax
      double precision, dimension(:) :: FiT
      double precision, dimension(:,:) :: BFiT
c*********************************************************************************************
c
      PhiMin=big
      phiMax=-big
c
      do i=1,NumberOfElements
c
        phiMin=dmin1(phiMin,FiT(i))
        phiMax=dmax1(phiMax,FiT(i))
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          phiMin=dmin1(phiMin,BFiT(i,j))
          phiMax=dmax1(phiMax,BFiT(i,j))
c
        enddo
      enddo
c
      return
      end