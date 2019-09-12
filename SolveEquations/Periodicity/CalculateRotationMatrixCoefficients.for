c
C#############################################################################################
      SUBROUTINE RotationMatrixCoefficients
C#############################################################################################
      use User0
      use Geometry1, only: NumberOfBCSets,NumbOfElementNodes
      use Geometry3, only: NIFaces,NBFaces,NBFacesMax,NBFaceOwner,
     *                     NBFaceNeighbor,ElementNeighbor,NGlobalEFaces
      use Geometry4, only: BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
      use BoundaryConditions1, only: BoundaryType
      use BoundaryConditions2, only:PeriodicPair,theta,a1r,b1r,c1r,
     *                              a2r,b2r,c2r,a3r,b3r,c3r,
     *                              a1Axis,a2Axis,a3Axis,
     *                              Icorrespondingface,Iperiodic,
     *             LRotationalPeriodicity,LTranslationalPeriodicity,
     *                      xTranslation,yTranslation,zTranslation,
     *                 xTranslationUse,yTranslationUse,zTranslationUse,
     *                             LPeriodicImplicit,IperiodicOwner,
     *                      IperiodicNumberOfBCSets,IperiodicNBFaces,
     *                      ElementPeriodicFace
      use Constants1, only: big,pi,epsilonPeriodic
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer :: i,i1,i2,i3,i4,i4desired,j,j1,j2,j3,k
      double precision :: thetaR,xF1,yF1,zF1,distance1,distance2,
     *                    distance3,difference,minDiff,
     *                    difference1,difference2,difference3
c*********************************************************************************************
c
c---- Translation periodicity
c
c     R=I
c
c---- rotation through an angle theta around the axis (a1,a2,a3)
c     in 2D the axis is (0,0,a3)
c
c       | A1  B1  C1 |
c     R=| A2  B2  C2 |=I+(W^2)(1-cos(theta))+Wsin(theta)
c       | A3  B3  C3 |
c        
c    W is the skew symmetric tensor given by        
c       | 0  -a3  a2 |
c     W=| a3  0  -a1 |
c       |-a2  a1   0 |
c        
c     A1=1-(a2^2+a3^2)(1-cos(theta))   
c     A2=a1a2(1-cos(theta))+a3sin(theta)
c     A3=a1a3(1-cos(theta))-a2sin(theta)
c     B1=a1a2(1-cos(theta))-a3sin(theta) 
c     B2=1-(a1^2+a3^2)(1-cos(theta)) 
c     B3=a2a3(1-cos(theta))+a1sin(theta) 
c     C1=a1a3(1-cos(theta))+a2sin(theta) 
c     C2=a2a3(1-cos(theta))-a1sin(theta)
c     C3=1-(a1^2+a2^2)(1-cos(theta)) 
c        
c---- Tranformation of velocity : V(F1)=RxV(C2)
c
C********************************************************************************************
c
c
c--- In two dimensional situations a1 and a2 are zero ==> A3=B3=C1=C2=0 C3=1 
c    for every periodic face the user has to specify the rotation vector
c    and rotation angle (positive in the clockwise direction)
c    for every face the user has to specify the corresponding transformed face     
c
      allocate(Icorrespondingface(NumberOfBCSets,NBFacesMax))
c
      if(LRotationalPeriodicity) then
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'periodic') then
c
            thetaR=theta(i)*pi/180.
c
            a1r(i)=1.-(a2Axis(i)**2+a3Axis(i)**2)*(1.-dcos(thetaR))
            a2r(i)=a1Axis(i)*a2Axis(i)*(1.-dcos(thetaR))+
     *                                       a3Axis(i)*dsin(thetaR)
            a3r(i)=a1Axis(i)*a3Axis(i)*(1.-dcos(thetaR))-
     *                                       a2Axis(i)*dsin(thetaR)
            b1r(i)=a1Axis(i)*a2Axis(i)*(1.-dcos(thetaR))-
     *                                       a3Axis(i)*dsin(thetaR)
            b2r(i)=1.-(a1Axis(i)**2+a3Axis(i)**2)*(1.-dcos(thetaR))
            b3r(i)=a2Axis(i)*a3Axis(i)*(1.-dcos(thetaR))+
     *                                       a1Axis(i)*dsin(thetaR)
            c1r(i)=a1Axis(i)*a3Axis(i)*(1.-dcos(thetaR))+
     *                                       a2Axis(i)*dsin(thetaR)
            c2r(i)=a2Axis(i)*a3Axis(i)*(1.-dcos(thetaR))-
     *                                       a1Axis(i)*dsin(thetaR)
            c3r(i)=1.-(a1Axis(i)**2+a2Axis(i)**2)*(1.-dcos(thetaR))
c
          endif
c
        enddo
c
c--- select one set of the periodic boundary couples and establish face correspondence 
c    and the positive direction of rotation
c
        do i1=1,NumberOfBCSets
c
          if(BoundaryType(i1).eq.'periodic') then
c
            i2=PeriodicPair(i1)
c
            if(i2.gt.i1) then
c
              thetaR=theta(i1)*pi/180.
c
              do i3=1,NBFaces(i1)
c
                xF1=a1r(i1)*BFaceCentroidx(i1,i3)+
     *                           b1r(i1)*BFaceCentroidy(i1,i3)+
     *                                c1r(i1)*BFaceCentroidz(i1,i3)
                yF1=a2r(i1)*BFaceCentroidx(i1,i3)+
     *                           b2r(i1)*BFaceCentroidy(i1,i3)+
     *                                c2r(i1)*BFaceCentroidz(i1,i3)
                zF1=a3r(i1)*BFaceCentroidx(i1,i3)+
     *                           b3r(i1)*BFaceCentroidy(i1,i3)+
     *                                c3r(i1)*BFaceCentroidz(i1,i3)
c                distance1=xF1**2+yF1**2+zF1**2
c                minDiff=big         
c
                do i4=1,NBFaces(i2)
c
c                  distance2=BFaceCentroidx(i2,i4)**2+
c     *                             BFaceCentroidy(i2,i4)**2+
c     *                                BFaceCentroidz(i2,i4)**2
c                  difference=dabs(distance2-distance1)                 
                  difference1=dabs(BFaceCentroidx(i2,i4)-xF1)
                  difference2=dabs(BFaceCentroidy(i2,i4)-yF1)
                  difference3=dabs(BFaceCentroidz(i2,i4)-zF1)
c
c                  if(difference.lt.minDiff) then
                  if(difference1.lt.epsilonPeriodic.and.
     *               difference2.lt.epsilonPeriodic.and.
     *               difference3.lt.epsilonPeriodic) then
c
c                    minDiff=difference
                    i4desired=i4
c
                  endif
c
                enddo            
c              
                if(LPeriodicImplicit) then
c
                  NBFaceNeighbor(i1,i3)=NBFaceOwner(i2,i4desired)
                  NBFaceNeighbor(i2,i4desired)=NBFaceOwner(i1,i3)
c
                endif
c
                Icorrespondingface(i1,i3)=i4desired
                Icorrespondingface(i2,i4desired)=i3
c
              enddo          
c          
            endif
c          
          endif
c          
        enddo
c
      elseif(LTranslationalPeriodicity) then
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'periodic') then
c
            a1r(i)=1.
            a2r(i)=0.
            a3r(i)=0.
            b1r(i)=0.
            b2r(i)=1.
            b3r(i)=0.
            c1r(i)=0.
            c2r(i)=0.
            c3r(i)=1.
c
          endif
c
        enddo
c
        do i1=1,NumberOfBCSets
c
          if(BoundaryType(i1).eq.'periodic') then
c
            i2=PeriodicPair(i1)
            if(i2.gt.i1) then
c
              do i3=1,NBFaces(i1)
c
                xF1=BFaceCentroidx(i1,i3)
                yF1=BFaceCentroidy(i1,i3)
                zF1=BFaceCentroidz(i1,i3)
                minDiff=big         
c
                do i4=1,NBFaces(i2)
c
                  distance1=dabs(dabs(BFaceCentroidx(i2,i4)-xF1)-
     *                                          dabs(xTranslation(i1)))
                  distance2=dabs(dabs(BFaceCentroidy(i2,i4)-yF1)-
     *                                          dabs(yTranslation(i1)))
                  distance3=dabs(dabs(BFaceCentroidz(i2,i4)-zF1)-
     *                                          dabs(zTranslation(i1)))
                  difference=dmax1(distance1,distance2,distance3)
c
                  if(difference.lt.minDiff) then
c
                    minDiff=difference
                    i4desired=i4
c
                  endif
c
                enddo            
c              
                if(LPeriodicImplicit) then
c
                  NBFaceNeighbor(i1,i3)=NBFaceOwner(i2,i4desired)
                  NBFaceNeighbor(i2,i4desired)=NBFaceOwner(i1,i3)
c
                endif
c
                Icorrespondingface(i1,i3)=i4desired
                Icorrespondingface(i2,i4desired)=i3
c
              enddo          
c          
            endif
c          
          endif
c          
        enddo
c
        do i1=1,NumberOfBCSets
c
          if(BoundaryType(i1).eq.'periodic') then
c
            xTranslationUse=dabs(xTranslation(i1))         
            yTranslationUse=dabs(yTranslation(i1))         
            zTranslationUse=dabs(zTranslation(i1))         
            exit
c
          endif
c
        enddo
c
      endif
c              
      if(LPeriodicImplicit) then
c
        allocate(ElementPeriodicFace(Iperiodic))

        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          do j=1,i2-1
c
            i4=i4+NBFaces(j)
c
          enddo
c
          i4=i4+i3
c
          j2=PeriodicPair(i2)         
          j3=Icorrespondingface(i2,i3)
          j1=NBFaceOwner(j2,j3)
c
          do j=1,NumbOfElementNodes(i1)
c
            k=NGlobalEFaces(i1,j)
c
            if(k.eq.i4) then
c
              ElementPeriodicFace(i)=j
              ElementNeighbor(i1,j)=j1
              exit
c
            endif
c
          enddo
c
        enddo
c
      endif
c
      return
      end