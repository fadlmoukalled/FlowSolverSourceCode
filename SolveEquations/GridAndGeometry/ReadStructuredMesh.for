c
C#############################################################################################
      SUBROUTINE ReadStructuredGrid
C#############################################################################################
c
      use User0, only: GridScalex,GridScaley,GridScalez
      use Geometry1
      use Geometry2
      use Structured1
c*********************************************************************************************
      implicit none
c*********************************************************************************************
      integer i,j,k,m,l1,m1,n1,ij,i1,j1,k1
c*********************************************************************************************
c
      read(1,*) NumberOfElementGroups
      read(1,*) l1,m1,n1
      NumberOfNodes=l1*m1*n1
      NumberOfElements=(l1-1)*(m1-1)*(n1-1)
c      NumberOfElementGroups=1
      NumberOfBCSets=6
      NumberOfCoordDir=3
      NumberOfVelComp=3
c
      allocate(xtemp(l1,m1,n1))
      allocate(ytemp(l1,m1,n1))
      allocate(ztemp(l1,m1,n1))
c
      read(1,*) (((xtemp(i,j,k),i=1,l1),j=1,m1),k=1,n1),
     +          (((ytemp(i,j,k),i=1,l1),j=1,m1),k=1,n1),
     +          (((ztemp(i,j,k),i=1,l1),j=1,m1),k=1,n1)
c
c---- Nodal Coordinates
c
      allocate(x(NumberOfNodes))
      allocate(y(NumberOfNodes))
      allocate(z(NumberOfNodes))
c
      do k=1,n1
        do j=1,m1
          do i=1,l1
c
            x((k-1)*l1*m1+(j-1)*l1+i)=xtemp(i,j,k)     
            y((k-1)*l1*m1+(j-1)*l1+i)=ytemp(i,j,k)     
            z((k-1)*l1*m1+(j-1)*l1+i)=ztemp(i,j,k)     
c
          enddo
        enddo
      enddo
c
      deallocate(xtemp)
      deallocate(ytemp)
      deallocate(ztemp)
c
      x=x*GridScalex
      y=y*GridScaley
      z=z*GridScalez
c
c---- Element/Cell Connectivity
c
      n=8
      allocate(NTypeGeometry(NumberOfElements))
      allocate(NumbOfElementNodes(NumberOfElements))
      allocate(ListOfElementNodes(NumberOfElements,n))
      allocate(NumbOfElementFaces(NumberOfElements))
c      allocate(NumbOfElementEdges(NumberOfElements))
c
      NTypeGeometry=4
      NumbOfElementNodes=8
      NumbOfElementFaces=6
c      NumbOfElementEdges=12
      k1=0
      do k=1,n1-1
        do j=1,m1-1
          do i=1,l1-1
c
            k1=k1+1
            ij=(k-1)*l1*m1+(j-1)*l1+i
            ListOfElementNodes(k1,1)=ij+l1
            ListOfElementNodes(k1,2)=ij+l1+1
            ListOfElementNodes(k1,3)=ij+l1+l1*m1
            ListOfElementNodes(k1,4)=ij+l1+l1*m1+1
            ListOfElementNodes(k1,5)=ij
            ListOfElementNodes(k1,6)=ij+1
            ListOfElementNodes(k1,7)=ij+l1*m1
            ListOfElementNodes(k1,8)=ij+l1*m1+1
c
          enddo
        enddo
      enddo
c
c---- Element Group Information
c
      allocate(NumberOfGroupElements(NumberOfElementGroups))
      allocate(MaterialGroupType(NumberOfElementGroups))
      allocate(NumberOfGroupFlags(NumberOfElementGroups))
      allocate(GroupName(NumberOfElementGroups))
      allocate(NElementsInGroup
     *          (NumberOfElementGroups,NumberOfElements))
c
      ElementGroup='group1'
      NGroup='Group1'
      NumberOfGroup=1
      NElementRead='nelem'
      NumberOfGroupElements(NumberOfGroup)=NumberOfElements
      material='air'
      MaterialGroupType(NumberOfGroup)=5
      NFlagsRead='read'
      NumberOfGroupFlags(NumberOfGroup)=1
      GroupName(NumberOfGroup)='group'
      allocate(NGroupFlags(1))
      NGroupFlags(1)=1
      do i=1,NumberOfElements
        NElementsInGroup(1,i)=i
      enddo
c
c---- Boundary Conditions
c
      allocate(BoundaryName(NumberOfBCSets))
      allocate(NBCDataType(NumberOfBCSets))
      allocate(NBCDataRecords(NumberOfBCSets))
c
      allocate(NodeBCTemp(NumberOfBCSets,NumberOfElements))
      allocate(NElementBCTemp(NumberOfBCSets,NumberOfElements))
      allocate(NElementBCTypeTemp(NumberOfBCSets,NumberOfElements))
      allocate(NElementBCFaceTemp(NumberOfBCSets,NumberOfElements))
c
      do i=1,NumberOfBCSets
c
        BoundaryConditions='BConditions'      
        BoundaryName(i)='i'
        NBCDataType(i)=1
        if(i.eq.1.or.i.eq.3) then
          NBCDataRecords(i)=(l1-1)*(m1-1)
        elseif(i.eq.2.or.i.eq.4) then
          NBCDataRecords(i)=(m1-1)*(n1-1)
        elseif(i.eq.5.or.i.eq.6) then
          NBCDataRecords(i)=(l1-1)*(n1-1)
        endif
c
        NBCValuesPerRecord=0
        NBCCode=6
c



        if(i.eq.1) then 
c
          k1=1
          j=0
          do i1=1,l1-1
            do j1=1,m1-1
            ij=(k1-1)*(l1-1)*(m1-1)+(j1-1)*(l1-1)+i1
            j=j+1
            NElementBCTemp(i,j)=ij
            NElementBCTypeTemp(i,j)=3
            NElementBCFaceTemp(i,j)=1
          enddo
        enddo
c
        elseif(i.eq.2) then 
c
          i1=l1-1
          j=0
          do j1=1,m1-1
            do k1=1,n1-1
            ij=(k1-1)*(l1-1)*(m1-1)+(j1-1)*(l1-1)+i1
            j=j+1
            NElementBCTemp(i,j)=ij
            NElementBCTypeTemp(i,j)=3
            NElementBCFaceTemp(i,j)=2
          enddo
        enddo
c
        elseif(i.eq.3) then 
c
          k1=n1-1
          j=0
          do i1=1,l1-1
            do j1=1,m1-1
            ij=(k1-1)*(l1-1)*(m1-1)+(j1-1)*(l1-1)+i1
            j=j+1
            NElementBCTemp(i,j)=ij
            NElementBCTypeTemp(i,j)=3
            NElementBCFaceTemp(i,j)=3
          enddo
        enddo
c
        elseif(i.eq.4) then 
c
          i1=1
          j=0
          do j1=1,m1-1
            do k1=1,n1-1
            ij=(k1-1)*(l1-1)*(m1-1)+(j1-1)*(l1-1)+i1
            j=j+1
            NElementBCTemp(i,j)=ij
            NElementBCTypeTemp(i,j)=3
            NElementBCFaceTemp(i,j)=4
          enddo
        enddo
c
        elseif(i.eq.5) then 
c
          j1=m1-1
          j=0
          do i1=1,l1-1
            do k1=1,n1-1
            ij=(k1-1)*(l1-1)*(m1-1)+(j1-1)*(l1-1)+i1
            j=j+1
            NElementBCTemp(i,j)=ij
            NElementBCTypeTemp(i,j)=3
            NElementBCFaceTemp(i,j)=5
          enddo
        enddo
c
        elseif(i.eq.6) then 
c
          j1=1
          j=0
          do i1=1,l1-1
            do k1=1,n1-1
            ij=(k1-1)*(l1-1)*(m1-1)+(j1-1)*(l1-1)+i1
            j=j+1
            NElementBCTemp(i,j)=ij
            NElementBCTypeTemp(i,j)=3
            NElementBCFaceTemp(i,j)=6
          enddo
        enddo
c
        endif
c        
      enddo
c
      ntemp=-1
      do i=1,NumberOfBCSets
c
        ntemp=max(ntemp,NBCDataRecords(i))
c
      enddo
c
      allocate(NodeBC(NumberOfBCSets,ntemp))
      allocate(NElementBC(NumberOfBCSets,ntemp))
      allocate(NElementBCType(NumberOfBCSets,ntemp))
      allocate(NElementBCFace(NumberOfBCSets,ntemp))
c
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
c
          NodeBC(i,j)=NodeBCTemp(i,j)
          NElementBC(i,j)=NElementBCTemp(i,j)
          NElementBCType(i,j)=NElementBCTypeTemp(i,j)
          NElementBCFace(i,j)=NElementBCFaceTemp(i,j)
c
        enddo
      enddo
c
      deallocate(NodeBCTemp)
      deallocate(NElementBCTemp)
      deallocate(NElementBCTypeTemp)
      deallocate(NElementBCFaceTemp)
c
      close(1)
c
      return
      end