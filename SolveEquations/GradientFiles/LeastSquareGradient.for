c
c#############################################################################################
      SUBROUTINE LeastSquares
     *         (FiT,dfidxT,dfidyT,dfidzT,BFiT,BdfidxT,BdfidyT,BdfidzT)
c#############################################################################################
      use Geometry1, only:NumberOfNodes,ListOfElementNodes,
     *                    NumberOfElements,NumbOfElementFaces,
     *                    x,y,z,NumberOfBCSets,NBCDataRecords
      use Geometry3, only:NIFaceOwner,NIFaceNeighbor,
     *                    NIFaces,NBFaces,NBFaceOwner
      use Geometry4, only:xc,yc,zc,Volume,BFaceCentroidx,
     *                    BFaceCentroidy,BFaceCentroidz
      use User0,     only: InvDistancePower
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      double precision :: w1,w2,w3,w4,w5,w6,w7,w8,w9,
     *                    denom,denom1,denom2,denom3
      integer :: i,j,k
      double precision, dimension(:) :: FiT
      double precision, dimension(:) :: dfidxT
      double precision, dimension(:) :: dfidyT
      double precision, dimension(:) :: dfidzT
      double precision, dimension(:,:) :: BFiT
      double precision, dimension(:,:) :: BdfidxT
      double precision, dimension(:,:) :: BdfidyT
      double precision, dimension(:,:) :: BdfidzT
c********************************************************************************************
      double precision :: deltaPhi,deltax,deltay,deltaz,weight
      double precision, dimension(:), allocatable :: a11
      double precision, dimension(:), allocatable :: a22
      double precision, dimension(:), allocatable :: a33
      double precision, dimension(:), allocatable :: a12
      double precision, dimension(:), allocatable :: a13
      double precision, dimension(:), allocatable :: a23
      double precision, dimension(:), allocatable :: b1
      double precision, dimension(:), allocatable :: b2
      double precision, dimension(:), allocatable :: b3
c*********************************************************************************************
c
      dfidxT=0.
      dfidyT=0.
      dfidzT=0.
      BdfidxT=0.
      BdfidyT=0.
      BdfidzT=0.
c
      allocate(a11(NumberOfElements))
      allocate(a22(NumberOfElements))
      allocate(a33(NumberOfElements))
      allocate(a12(NumberOfElements))
      allocate(a13(NumberOfElements))
      allocate(a23(NumberOfElements))
      allocate(b1(NumberOfElements))
      allocate(b2(NumberOfElements))
      allocate(b3(NumberOfElements))
c
      a11=0.
      a22=0.
      a33=0.
      a12=0.
      a13=0.
      a23=0.
      b1=0.
      b2=0.
      b3=0.
c
      do k=1,NIFaces
c        
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
c
        deltaPhi=Fit(j)-Fit(i)
        deltax=xc(j)-xc(i)
        deltay=yc(j)-yc(i)
        deltaz=zc(j)-zc(i)
c
        weight=1./dsqrt(deltax*deltax+deltay*deltay+deltaz*deltaz)
        weight=weight**InvDistancePower
c
        w1=weight*deltax*deltax
        w2=weight*deltax*deltay
        w3=weight*deltax*deltaz
        w4=weight*deltay*deltay
        w5=weight*deltay*deltaz
        w6=weight*deltaz*deltaz
        w7=weight*deltax*deltaPhi
        w8=weight*deltay*deltaPhi
        w9=weight*deltaz*deltaPhi
c
        a11(i)=a11(i)+w1
        a12(i)=a12(i)+w2
        a13(i)=a13(i)+w3
        a22(i)=a22(i)+w4
        a23(i)=a23(i)+w5
        a33(i)=a33(i)+w6
        b1(i)=b1(i)+w7
        b2(i)=b2(i)+w8
        b3(i)=b3(i)+w9
c
        a11(j)=a11(j)+w1
        a12(j)=a12(j)+w2
        a13(j)=a13(j)+w3
        a22(j)=a22(j)+w4
        a23(j)=a23(j)+w5
        a33(j)=a33(j)+w6
        b1(j)=b1(j)+w7
        b2(j)=b2(j)+w8
        b3(j)=b3(j)+w9
c
      enddo
c
c--- Boundary faces
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          k=NBFaceOwner(i,j)
c
          deltaPhi=BFiT(i,j)-Fit(k)
          deltax=BFaceCentroidx(i,j)-xc(k)
          deltay=BFaceCentroidy(i,j)-yc(k)
          deltaz=BFaceCentroidz(i,j)-zc(k)
c
          weight=1./dsqrt(deltax*deltax+deltay*deltay+deltaz*deltaz)
          weight=weight**InvDistancePower
c
          a11(k)=a11(k)+weight*deltax*deltax
          a12(k)=a12(k)+weight*deltax*deltay
          a13(k)=a13(k)+weight*deltax*deltaz
          a22(k)=a22(k)+weight*deltay*deltay
          a23(k)=a23(k)+weight*deltay*deltaz
          a33(k)=a33(k)+weight*deltaz*deltaz
          b1(k)=b1(k)+weight*deltax*deltaPhi
          b2(k)=b2(k)+weight*deltay*deltaPhi
          b3(k)=b3(k)+weight*deltaz*deltaPhi
c
        enddo
      enddo
c
      do i=1,NumberOfElements
c
        denom=a11(i)*(a22(i)*a33(i)-a23(i)*a23(i))-
     *        a12(i)*(a12(i)*a33(i)-a13(i)*a23(i))+
     *        a13(i)*(a12(i)*a23(i)-a13(i)*a22(i))
        denom1=b1(i)*(a22(i)*a33(i)-a23(i)*a23(i))-
     *        a12(i)*(b2(i)*a33(i)-b3(i)*a23(i))+
     *        a13(i)*(b2(i)*a23(i)-b3(i)*a22(i))
        denom2=a11(i)*(b2(i)*a33(i)-b3(i)*a23(i))-
     *        b1(i)*(a12(i)*a33(i)-a13(i)*a23(i))+
     *        a13(i)*(a12(i)*b3(i)-a13(i)*b2(i))
        denom3=a11(i)*(a22(i)*b3(i)-b2(i)*a23(i))-
     *        a12(i)*(a12(i)*b3(i)-a13(i)*b2(i))+
     *        b1(i)*(a12(i)*a23(i)-a13(i)*a22(i))


        dfidxT(i)=denom1/denom
        dfidyT(i)=denom2/denom
        dfidzT(i)=denom3/denom
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
      deallocate(a11)
      deallocate(a22)
      deallocate(a33)
      deallocate(a12)
      deallocate(a13)
      deallocate(a23)
      deallocate(b1)
      deallocate(b2)
      deallocate(b3)
c
      return
      end