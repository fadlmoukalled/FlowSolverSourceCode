c
C#############################################################################################
      SUBROUTINE SetBoundaryConditions
C#############################################################################################
      use Geometry1, only : NumberOfBCSets
      use Geometry3, only : NBFaces
      use BoundaryConditions1
      use BoundaryConditionsScalar1
      use BoundaryConditionsrField1
      use User0
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k
c*********************************************************************************************
c
      if(LSolveLambdaELEEquation) then
c          
        do i=1,NumberOfBCSets
c
            if(BoundaryType(i).eq.'closed') then
               BoundaryType(i)='wall' 
               wallTypeL(i)='vonneumann'
            else
               BoundaryType(i)='wall' 
               wallTypeL(i)='dirichlet'
            endif
c
        enddo
c          
      endif
c      
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          BCType(i,j)=BoundaryType(i)
          wallTypeMomentum(i,j)=wallTypeM(i)
          wallTypeContinuity(i,j)=wallTypeC(i)
          wallTypeEnergy(i,j)=wallTypeE(i)
          wallTypeLambda(i,j)=wallTypeL(i)
c
          outletTypeMomentum(i,j)=outletTypeM(i)
          outletTypeContinuity(i,j)=outletTypeC(i)
          outletTypeEnergy(i,j)=outletTypeE(i)
c
          inletTypeMomentum(i,j)=inletTypeM(i)
          inletTypeContinuity(i,j)=inletTypeC(i)
          inletTypeEnergy(i,j)=inletTypeE(i)
c
          do k=1,NumberOfScalarsToSolve
c
            wallTypeScalar(i,j,k)=wallTypeS(i,k)
            outletTypeScalar(i,j,k)=outletTypeS(i,k)
            inletTypeScalar(i,j,k)=inletTypeS(i,k)
c
          enddo
c
          do k=1,NumberOfrFieldsToSolve
c
            wallTyperField(i,j,k)=wallTypeR(i,k)
            outletTyperField(i,j,k)=outletTypeR(i,k)
            inletTyperField(i,j,k)=inletTypeR(i,k)
c
          enddo
c
        enddo
      enddo
c
      return
      end