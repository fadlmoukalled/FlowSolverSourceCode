c
C#############################################################################################
c
      SUBROUTINE AssembleGlobalMatrixElementFluxes
c
C#############################################################################################
c
      use Geometry1, only:NumberOfElements
      use Variables2, only: ac,acold,acoldold,bc
      use Variables3, only: FluxCE,FluxTE,FluxCEold,FluxCEoldold
      use User0, only: LUnsteady
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,i1,j1,iBFace
c********************************************************************************************
c
c--- Add element flux contribution to elements
c
      do i=1,NumberOfElements
c      
        ac(i)=ac(i)+FluxCE(i)
        bc(i)=bc(i)-FluxTE(i)
c
      enddo
c
      if(LUnsteady) then
c
        do i=1,NumberOfElements
c      
          acold(i)=acold(i)+FluxCEold(i)
          acoldold(i)=acoldold(i)+FluxCEoldold(i)
c
        enddo
c
      endif
c
      return
      end
c
C#############################################################################################
c
      SUBROUTINE AssembleGlobalMatrixFaceFluxes
c
C#############################################################################################
c
      use Geometry1, only:NumberOfBCSets,NumberOfElements
      use Variables2
      use Geometry3
      use BoundaryConditions2, only:Iperiodic,IperiodicOwner,
     *                              IperiodicNumberOfBCSets,
     *                              LPeriodicImplicit,IperiodicNBFaces,
     *                              ElementPeriodicFace
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer i,j,k,i1,i2,i3,i4,j1,iBFace
c********************************************************************************************

      do k=1,NIFaces
c
        i=NIFaceOwner(k)
        j=NIFaceNeighbor(k)
        i1=iOwnerNeighbor(k)
        j1=iNeighborOwner(k)
c
c--- Assemble fluxes for owner cell
c
        ac(i) = ac(i) + FluxCf(k)
        anb(i,i1) = anb(i,i1) + FluxFf(k)
        bc(i)= bc(i) - FluxTf(k)
c 
c--- Assemble fluxes for neighbour cell
c   
        ac(j)  = ac(j)  - FluxFf(k)
        anb(j,j1) = anb(j,j1) - FluxCf(k)
        bc(j)   = bc(j) + FluxTf(k)
c
      enddo
c
c--- assemble fluxes of boundary faces
c
      iBFace=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c      
          k=NBFaceOwner(i,j)
          iBFace=iBFace+1
c
c--- Assemble fluxes for owner cell
c
          ac(k) = ac(k) + FluxCf(iBFace);
          bc(k) = bc(k) - FluxTf(iBFace);
c
        enddo
      enddo
c
c--- Add contribution to coefficients for periodic boundaries
c
      if(LPeriodicImplicit) then
c
        do i=1,Iperiodic
c
          i1=IperiodicOwner(i)
          i2=IperiodicNumberOfBCSets(i)
          i3=IperiodicNBFaces(i)
          i4=NIFaces
c
          do j1=1,i2-1
c
            i4=i4+NBFaces(j1)
c
          enddo
c
          i4=i4+i3
c
          j1 = ElementPeriodicFace(i)
          anb(i1,j1) = anb(i1,j1) + FluxFf(i4)
c
        enddo
c
      endif
c      
      return
      end