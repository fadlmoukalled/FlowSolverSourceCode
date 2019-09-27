!************************************************************************************************
!      Module to read polymesh grid, establish its connectivity, and print results to paraview
!************************************************************************************************
     MODULE ReadpolyMesh   
!************************************************************************************************
     type Face
       integer nPoints, neighbour, owner
       integer,allocatable :: PointsID(:)
     end type Face
!
     type Boundary
       integer :: nFaces, startFace
       character (len = 32) :: Boundary_name
     end type Boundary
!
     type Cell
       integer :: nFaces, nPoints
       integer,allocatable :: FacesID(:), PointsID(:)
     end type Cell

     integer :: nPoints, nFaces, nNeighbour, nBoundaries, nCells
     double precision,allocatable :: points(:,:)
     type(face), allocatable :: Faces(:)
     type(Boundary), allocatable :: Boundaries(:)
     type(Cell), allocatable :: Cells(:)

     !
     contains
!
!************************************************************************************************
     SUBROUTINE Readit
!************************************************************************************************
     use User0, only: PolyMeshDirectory,GridScalex,GridScaley,GridScalez
     use Geometry1, only: NumberOfNodes,NumberOfElements,NumberOfElementGroups,NumberOfBCSets,&
                          x,y,z,NumbOfElementFaces,NumbOfElementNodes,NumberOfGroupElements,&
                          ListOfElementNodes,NTypeGeometry,MaterialGroupType,NElementsInGroup,&
                          NBCDataType,NBCDataRecords,NodeBC,NElementBC,NElementBCType,&
                          NElementBCFace,NodeFlag,MaximumNumberofElementNodes
     use Geometry2, only: ListOfElementNodesTemp,NumberOfGroupFlags,GroupName,NElementsInGroupTemp,&
                             BoundaryName,NGroupFlags
     use Geometry3, only: NFacesTotal,NIFaces,NBFaces,NIFaceNodes,NumberOfElementFaceNodes,&
                          NBFacesMax,NIFaceOwner,NIFaceNeighbor,NBFaceOwner,NBFaceNeighbor,&
                          NGlobalEFaces,GlobalFaceNumberOfNodes,NBFaceNodes,&
                          LocalElementFaceNodes,BCSet,BCRecord,NBFacesTotal,iOwnerNeighbor,&
                          iNeighborOwner,ElementNeighbor,NumberofElementNeighbors,NElementFaces,&
                           MaximumNumberOfFaceNodes
     use InteriorBoundaryNodes, only: InteriorNodeTemp,BoundaryNodeTemp,InteriorNode,&
                                         BoundaryNode,NumberOfInteriorNodes,NumberOfBoundaryNodes
!************************************************************************************************
     implicit none
!************************************************************************************************
     integer :: i,i1,j,k,k1,k2,k3
     logical :: add
     integer :: IOstatus
     integer :: iOwner,iNeighbor,ntemp,n,iEdge,indexEdge,iElement
     character (len = 70) :: x9(10)
     character (len = 70) :: L0, L1 
!       
!----Reading points
!
     open(unit = 7, file = trim(PolyMeshDirectory)//"/points")
!
     i=0
     do 
       L0 = L1
       read (7,"(a)",iostat=IOstatus) L1
       if (IOstatus/=0) exit 
        if((index(L1,"("))>0) then
          read(L0,*) nPoints
          exit
        end if
     end do
!
     allocate(points(nPoints,3))
!
     do i=1,nPoints
!
       read (7,*) x9(1), x9(2), x9(3)
       read( x9(1)(2:20), * ) points(i,1)
       read( x9(2)(1:20), * ) points(i,2)
       read( x9(3)(1:(len(trim(x9(3)))-1)), * ) points(i,3)
       points(i,1)=points(i,1)*GridScalex
       points(i,2)=points(i,2)*GridScaley
       points(i,3)=points(i,3)*GridScalez
     end do
!
     close(7)
!        
!----Reading Faces
!
     open(unit = 7, file = trim(PolyMeshDirectory)//"/faces")
!        
     i=0
     L1 = ""
     do 
       L0 = L1
       read (7,"(a)",iostat=IOstatus) L1
       if (IOstatus/=0) exit 
       if((index(L1,"("))>0) then
         read(L0,*) nFaces
         exit
       end if
     end do
!        
     allocate(Faces(nFaces))
!        
     do i=1,nFaces
       read( 7, "(A)" ) L1
       read( L1(1:1), * ) faces(i)%nPoints
!
       allocate(faces(i)%PointsID(faces(i)%nPoints))
!           
       read(L1,*) x9(1:faces(i)%nPoints)
!
!----Reading Face Points
!
       read( x9(1)(3:20), * ) faces(i)%PointsID(1)
!
       do j=2,(faces(i)%nPoints-1)
         read( x9(j)(1:20), * ) faces(i)%PointsID(j)
       end do
!
       read( x9(faces(i)%nPoints)(1:(len(trim(x9(faces(i)%nPoints)))-1)), * ) faces(i)%PointsID(faces(i)%nPoints)
       faces(i)%PointsID(:) = faces(i)%PointsID(:)+1
!
     end do
!
     close(7)
!        
!----Reading owners
!
     open(unit = 7, file = trim(PolyMeshDirectory)//"/owner")
!
     do 
       read (7,"(a)",iostat=IOstatus) L1
       if (IOstatus/=0) exit 
       if((index(L1,"("))>0) exit
     end do
!        
     do i=1,nFaces
       read( 7, * ) faces(i)%owner
       faces(i)%owner = faces(i)%owner+1
     end do
!
     close(7)
!        
!----Reading neighbours
!
     open(unit = 7, file = trim(PolyMeshDirectory)//"/neighbour")
!     
     L1 = ""
     do 
       L0 = L1
       read (7,"(a)",iostat=IOstatus) L1
       if (IOstatus/=0) exit 
       if((index(L1,"("))>0) then
         read(L0,*) nNeighbour
         exit
        end if
     end do
!
     faces(:)%neighbour = -1
!        
     do i=1,nNeighbour
       read( 7, * ) faces(i)%neighbour
       faces(i)%neighbour = faces(i)%neighbour+1
     end do
!
     close(7)
!        
!----Reading Boundaries
!
     open(unit = 7, file = trim(PolyMeshDirectory)//"/boundary")
!        
     L1 = ""
     do 
       L0 = L1
       read (7,"(a)",iostat=IOstatus) L1
       if (IOstatus/=0) exit 
       if((index(L1,"("))>0) then
         read(L0,*) nBoundaries
         exit
       end if
     end do
!        
     allocate(Boundaries(nBoundaries))
!        
     i=1
     do 
       L0 = L1
       read (7,"(a)",iostat=IOstatus) L1
       if (IOstatus/=0) exit
!             
       if((index(L1,"{"))>0) then
         read(L0,*) Boundaries(i)%Boundary_name
       end if
!            
       if((index(L1,"nFaces"))>0) then
         read( L1(1:(len(trim(L1))-1)), "(A)" ) x9(1)
         read(x9(1),*) x9(2), Boundaries(i)%nFaces
       end if
!           
       if((index(L1,"startFace"))>0) then
         read( L1(1:(len(trim(L1))-1)), "(A)" ) x9(1)
         read(x9(1),*) x9(2), Boundaries(i)%startFace
         Boundaries(i)%StartFace = Boundaries(i)%StartFace+1
         i=i+1
       end if
!            
     end do
!
     close(7)
! 
!----Establish Global Connectivity        
!        
     NumberOfNodes=nPoints     
     allocate(x(NumberOfNodes)) 
     allocate(y(NumberOfNodes)) 
     allocate(z(NumberOfNodes)) 
!        
     do i=1,NumberOfNodes
       x(i)=points(i,1)
       y(i)=points(i,2)
       z(i)=points(i,3)
     enddo
!
     NFacesTotal=nFaces
     NumberOfBCSets=nBoundaries
     NIFaces=NFacesTotal
     NBFaces=0
!
!---- Extract number of elements
!
     NumberOfElements=-1
     do i=1,NFacesTotal
       NumberOfElements=max(NumberOfElements,faces(i)%owner,faces(i)%neighbour)
     enddo
!
     allocate(NBFaces(NumberOfBCSets))
!
     NBFacesMax=-1
     do i=1,NumberOfBCSets
       NBFaces(i)=Boundaries(i)%nFaces
       NBFacesMax=max(NBFacesMax,NBFaces(i))
       NIFaces=NIFaces-NBFaces(i)
     enddo
!
     allocate(NIFaceOwner(NIFaces))
     allocate(NIFaceNeighbor(NIFaces))
     allocate(NBFaceOwner(NumberOfBCSets,NBFacesMax))
     allocate(NBFaceNeighbor(NumberOfBCSets,NBFacesMax))
!       
     do i=1,NIFaces
       NIFaceOwner(i)=faces(i)%owner
       NIFaceNeighbor(i)=faces(i)%neighbour
     enddo
!       
     k=NIFaces
     do i=1,NumberOfBCSets
       do j=1,NBFaces(i)
         k=k+1
         NBFaceOwner(i,j)=faces(k)%owner
         NBFaceNeighbor(i,j)=-1
       enddo
     enddo
!
     allocate(NumbOfElementFaces(NumberOfElements))
!
     NumbOfElementFaces=0
     do i=1,NIFaces
       iOwner=faces(i)%owner
       iNeighbor=faces(i)%neighbour
       NumbOfElementFaces(iOwner)=NumbOfElementFaces(iOwner)+1
       NumbOfElementFaces(iNeighbor)=NumbOfElementFaces(iNeighbor)+1
     enddo
     do i=NIFaces+1,NFacesTotal
       iOwner=faces(i)%owner
       NumbOfElementFaces(iOwner)=NumbOfElementFaces(iOwner)+1
     enddo
!
     NElementFaces=-1
     do i=1,NumberOfElements
       NElementFaces=max(NElementFaces,NumbOfElementFaces(i))   
     enddo   
!
     allocate(NGlobalEFaces(NumberOfElements,NElementFaces))
!
     NumbOfElementFaces=0
     do i=1,NIFaces
       iOwner=faces(i)%owner
       iNeighbor=faces(i)%neighbour
       NumbOfElementFaces(iOwner)=NumbOfElementFaces(iOwner)+1
       NGlobalEFaces(iOwner,NumbOfElementFaces(iOwner))=i
       NumbOfElementFaces(iNeighbor)=NumbOfElementFaces(iNeighbor)+1
       NGlobalEFaces(iNeighbor,NumbOfElementFaces(iNeighbor))=i
     enddo
!
     do i=NIFaces+1,NFacesTotal
       iOwner=faces(i)%owner
       NumbOfElementFaces(iOwner)=NumbOfElementFaces(iOwner)+1
       NGlobalEFaces(iOwner,NumbOfElementFaces(iOwner))=i
     enddo
!
     nCells=NumberOfElements
     allocate(cells(NumberOfElements))
!            
     do i=1,NumberOfElements   
       Cells(i)%nFaces=NumbOfElementFaces(i)
       allocate(Cells(i)%FacesID(NumbOfElementFaces(i)))
       do j=1,NumbOfElementFaces(i)
         Cells(i)%FacesID(j)=NGlobalEFaces(i,j)
       enddo
     enddo
!
     MaximumNumberOfFaceNodes=0
     do i=1,NFacesTotal
       MaximumNumberOfFaceNodes=max(MaximumNumberOfFaceNodes,faces(i)%nPoints)
     enddo
!
     allocate(NIFaceNodes(NIFaces,MaximumNumberOfFaceNodes))
     allocate(GlobalFaceNumberOfNodes(NFacesTotal))
     allocate(NumberOfElementFaceNodes(NumberOfElements,NElementFaces))
     allocate(NBFaceNodes(NumberOfBCSets,NBFacesMax,MaximumNumberOfFaceNodes))
! 
     do i=1,NIFaces
       GlobalFaceNumberOfNodes(i)= faces(i)%nPoints
       do j=1,GlobalFaceNumberOfNodes(i)
         NIFaceNodes(i,j)=faces(i)%PointsID(j)
       enddo
     enddo
!
     k=NIFaces
     do i=1,NumberOfBCSets
       do j=1,NBFaces(i)
         k=k+1
         GlobalFaceNumberOfNodes(k)= faces(k)%nPoints
         do k1=1,faces(k)%nPoints
           NBFaceNodes(i,j,k1)= faces(k)%PointsID(k1)
         enddo
       enddo
     enddo
!
!---- Set the maximum possible number of element nodes
!
     n=200
     allocate(ListOfElementNodesTemp(NumberOfElements,n))
     allocate(NumbOfElementNodes(NumberOfElements))
     ListOfElementNodesTemp=-1
     NumbOfElementNodes=0
!
     do i=1,NumberOfElements
       k1=0
       do j=1,NumbOfElementFaces(i)
         k2=NGlobalEFaces(i,j)
         do k=1,faces(k2)%nPoints
           if(j.eq.1) then
             k1=k1+1
             NumbOfElementNodes(i)=NumbOfElementNodes(i)+1
             ListOfElementNodesTemp(i,k1)=faces(k2)%PointsID(k)
           elseif(k1.lt.n) then
             k1=k1+1
             NumbOfElementNodes(i)=NumbOfElementNodes(i)+1
             ListOfElementNodesTemp(i,k1)=faces(k2)%PointsID(k)
             i1loop: do i1=1,k1-1
               if(faces(k2)%PointsID(k).eq.ListOfElementNodesTemp(i,i1))  then
                 ListOfElementNodesTemp(i,k1)=-1
                 NumbOfElementNodes(i)=NumbOfElementNodes(i)-1
                 k1=k1-1
                 exit i1loop
               endif
             enddo i1loop
           endif
         enddo
       enddo
     enddo
!
     MaximumNumberofElementNodes=-1
     do i=1,NumberOfElements
       MaximumNumberofElementNodes=max(MaximumNumberofElementNodes,NumbOfElementNodes(i))
     enddo      
     allocate(ListOfElementNodes(NumberOfElements,MaximumNumberofElementNodes))
     ListOfElementNodes=-1
     do i=1,NumberOfElements
       do j=1,NumbOfElementNodes(i)
         ListOfElementNodes(i,j)= ListOfElementNodesTemp(i,j)
       enddo      
     enddo    
!
     deallocate(ListOfElementNodesTemp)
!
     do i=1,NumberOfElements
       Cells(i)%nPoints=NumbOfElementNodes(i)
       allocate(Cells(i)%PointsID(NumbOfElementNodes(i)))
       do j=1,NumbOfElementNodes(i)
         Cells(i)%PointsID(j)=ListOfElementNodes(i,j)
       enddo
     enddo
!
!--- Sort nodes for elements
!
     allocate(NTypeGeometry(NumberOfElements))
     NTypeGeometry=1
!
!--- Set element geometry type (not used with polymesh)
!
     do i=1,NumberOfElements
       if(NumbOfElementNodes(i).eq.8) then
         NTypeGeometry(i)=4
       elseif(NumbOfElementNodes(i).eq.6) then
         NTypeGeometry(i)=5
       elseif(NumbOfElementNodes(i).eq.4) then
         NTypeGeometry(i)=6
       elseif(NumbOfElementNodes(i).eq.5) then
         NTypeGeometry(i)=7
       endif
     enddo
!
     do i=1,NumberOfElements
!
       k=NGlobalEFaces(i,1)
       k1=0
       do j=1,faces(k)%nPoints
         ListOfElementNodes(i,j)=faces(k)%PointsID(j)
         k1=k1+1
       enddo
!            
       jloop: do j=2,NumbOfElementFaces(i)
         k=NGlobalEFaces(i,j)
!
         do k2=1,faces(k)%nPoints
           add=.true.
           k3loop: do k3=1,k1
             if(ListOfElementNodes(i,k3).eq.faces(k)%PointsID(k2)) then
               add=.false.
               exit k3loop
             endif
           enddo k3loop
           if(add) then
             k1=k1+1
             ListOfElementNodes(i,k1)=faces(k)%PointsID(k2)
           endif
           if(k1.eq.NumbOfElementNodes(i)) exit jloop
         enddo
!
       enddo jloop
!
     enddo
!      
     allocate(LocalElementFaceNodes(NumberOfElements,NElementFaces,MaximumNumberOfFaceNodes))
!
     do i=1,NumberOfElements
       do j=1,NumbOfElementFaces(i)
!           
         k=NGlobalEFaces(i,j)
         NumberOfElementFaceNodes(i,j)=GlobalFaceNumberOfNodes(k)
!
         do k1=1,GlobalFaceNumberOfNodes(k)
           LocalElementFaceNodes(i,j,k1)=faces(k)%PointsID(k1)
         enddo
!
       enddo
     enddo
!
!---- Element Group Information
!
     NumberOfElementGroups=1
!
     allocate(NumberOfGroupElements(NumberOfElementGroups))
     allocate(MaterialGroupType(NumberOfElementGroups))
     allocate(NumberOfGroupFlags(NumberOfElementGroups))
     allocate(GroupName(NumberOfElementGroups))
     allocate(NElementsInGroup(NumberOfElementGroups,NumberOfElements))
     allocate(NGroupFlags(1))
!
     NumberOfGroupElements(1)=NumberOfElements
     NumberOfGroupFlags(1)=1
     NGroupFlags=1
!
     do i=1,NumberOfElementGroups
       do j=1,NumberOfGroupElements(i)      
         NElementsInGroup(i,j)=j
       enddo
     enddo
!      
     allocate(BoundaryName(NumberOfBCSets))
     allocate(NBCDataType(NumberOfBCSets))
     allocate(NBCDataRecords(NumberOfBCSets))
!
     do i=1,NumberOfBCSets
       NBCDataRecords(i)=NBFaces(i)
       BoundaryName(i)=Boundaries(i)%Boundary_name
     enddo
!      
     ntemp=-1
     do i=1,NumberOfBCSets
       ntemp=max(ntemp,NBCDataRecords(i))
     enddo
!
     allocate(NodeBC(NumberOfBCSets,ntemp))
     allocate(NElementBC(NumberOfBCSets,ntemp))
     allocate(NElementBCType(NumberOfBCSets,ntemp))
     allocate(NElementBCFace(NumberOfBCSets,ntemp))
!
!---- Process boundary faces (Assumes BC type 1, i.e. element based)
!
     i1=0
     do i=1,NumberOfBCSets
       i1=i1+NBCDataRecords(i)
     enddo
!
     allocate(BCSet(i1))
     allocate(BCRecord(i1))
!
     NBFacesTotal=0
!
     do i=1,NumberOfBCSets
       do j=1,NBCDataRecords(i)
!
         NBFacesTotal=NBFacesTotal+1
         BCSet(NBFacesTotal)=i
         BCRecord(NBFacesTotal)=j
!
       enddo
     enddo
!
!--- Locate the element neighbors, OwnerNeighbor, and NeighborOwner connectivity
!
     allocate(iOwnerNeighbor(NIFaces))
     allocate(iNeighborOwner(NIFaces))
     allocate(ElementNeighbor(NumberOfElements,NElementFaces))
!
     iOwnerNeighbor=0
     iNeighborOwner=0
     ElementNeighbor=0
!
     do i=1,NIFaces
!
       iElement=NIFaceOwner(i)
       iNeighbor=NIFaceNeighbor(i)
!
       do j=1,NumbOfElementFaces(iElement)
!
         k=NGlobalEFaces(iElement,j)
!
         if(k.eq.i) then
!
           iOwnerNeighbor(i)=j
           ElementNeighbor(iElement,j)=iNeighbor
           exit
!
         endif
!
       enddo
!
       do j=1,NumbOfElementFaces(iNeighbor)
!
         k=NGlobalEFaces(iNeighbor,j)
!
         if(k.eq.i) then
!
           iNeighborOwner(i)=j 
           ElementNeighbor(iNeighbor,j)=iElement
           exit
!
         endif
!
       enddo
!
     enddo
!
     allocate(NumberofElementNeighbors(NumberOfElements))
!
     do i=1,NumberOfElements
!
       NumberofElementNeighbors(i)=NumbOfElementFaces(i)
!
     enddo
!
!--- Create flags for nodes (interior node ---> flag=1,boundary node ---> flag=0) 
!
     k=NIFaces
     do i=1,NumberOfBCSets
       do j=1,NBFaces(i)
         k=k+1
         NElementBC(i,j)=faces(k)%owner
         NElementBCFace(i,j)=k
       enddo
     enddo
!
     allocate(NodeFlag(NumberOfNodes))
     NodeFlag=1
!      
     k=NIFaces
     do i=1,NumberOfBCSets
       do j=1,NBFaces(i)
         k=k+1
         do k1=1,faces(k)%nPoints
           NodeFlag(NBFaceNodes(i,j,k1))= 0
         enddo
       enddo
     enddo
!
     allocate(InteriorNodeTemp(NumberOfNodes))
     allocate(BoundaryNodeTemp(NumberOfNodes))
!
     j=0
     k=0   
     do i=1,NumberOfNodes
!
       if(NodeFlag(i).eq.1) then
!
         j=j+1
         InteriorNodeTemp(j)=i     
!
       elseif(NodeFlag(i).eq.0) then
!
         k=k+1
         BoundaryNodeTemp(k)=i
!
       endif
!
     enddo
!
     NumberOfInteriorNodes=j
     NumberOfBoundaryNodes=k
!
     allocate(InteriorNode(NumberOfInteriorNodes))
     allocate(BoundaryNode(NumberOfBoundaryNodes))
!
     InteriorNode(1:NumberOfInteriorNodes)=InteriorNodeTemp(1:NumberOfInteriorNodes)
     BoundaryNode(1:NumberOfBoundaryNodes)=BoundaryNodeTemp(1:NumberOfBoundaryNodes)
!
     deallocate(InteriorNodeTemp)
     deallocate(BoundaryNodeTemp)
!
     return
     end SUBROUTINE Readit
!************************************************************************************************
     SUBROUTINE digit_to_ch ( digit, ch )
!************************************************************************************************
!
!--- DIGIT_TO_CH returns the character representation of a decimal digit.
!    Discussion: Instead of CHAR, we now use the ACHAR function, which
!                guarantees the ASCII collating sequence.
!    Example:
!
!    DIGIT   CH 
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!    Licensing: This code is distributed under the GNU LGPL license. 
!    Modified:04 August 1999
!    Author: John Burkardt
!
!    Parameters:
!      Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!      Output, character CH, the corresponding character.
!************************************************************************************************
     implicit none
!************************************************************************************************
     character ch
     integer ( kind = 4 ) digit
!************************************************************************************************
!
     if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
     else
       ch = '*'
     end if
!
     return
     end SUBROUTINE digit_to_ch
!************************************************************************************************
     SUBROUTINE get_unit ( iunit )
!************************************************************************************************
!
!    GET_UNIT returns a free FORTRAN unit number.
!
!    Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!    Licensing:This code is distributed under the GNU LGPL license. 
!    Modified:26 October 2008
!    Author:John Burkardt
!    Parameters:
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!************************************************************************************************
     implicit none
!************************************************************************************************
     integer ( kind = 4 ) i
     integer ( kind = 4 ) ios
     integer ( kind = 4 ) iunit
     logical lopen
!************************************************************************************************
!
     iunit = 0
!
     do i = 1, 99
       if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
         inquire ( unit = i, opened = lopen, iostat = ios )
         if ( ios == 0 ) then
           if ( .not. lopen ) then
             iunit = i
             return
           end if
         end if
       end if
     end do
!
     return
     end SUBROUTINE get_unit
!************************************************************************************************
     SUBROUTINE i4_to_s_left ( i4, s )
!************************************************************************************************
!
!    I4_TO_S_LEFT converts an I4 to a left-justified string.
!    Discussion: An I4 is an integer ( kind = 4 ).
!    Example: Assume that S is 6 characters long:
!
!         I4  S
!          1  1
!         -1  -1
!          0  0
!       1952  1952
!     123456  123456
!    1234567  ******  <-- Not enough room!
!
!    Licensing: This code is distributed under the GNU LGPL license. 
!    Modified:28 July 2000
!    Author:John Burkardt
!
!    Parameters:
!      Input, integer ( kind = 4 ) I4, an integer to be converted.
!      Output, character ( len = * ) S, the representation of the integer.
!      The integer will be left-justified.  If there is not enough space,
!      the string will be filled with stars.
!************************************************************************************************
     implicit none
!************************************************************************************************
     character c
     integer ( kind = 4 ) i
     integer ( kind = 4 ) i4
     integer ( kind = 4 ) idig
     integer ( kind = 4 ) ihi
     integer ( kind = 4 ) ilo
     integer ( kind = 4 ) ipos
     integer ( kind = 4 ) ival
     character ( len = * ) s
!************************************************************************************************
!
     s = ' '
!
     ilo = 1
     ihi = len ( s )
!
     if ( ihi <= 0 ) then
       return
     end if
!
!--- Make a copy of the integer.
!
     ival = i4
!
!--- Handle the negative sign.
!
     if ( ival < 0 ) then
       if ( ihi <= 1 ) then
         s(1:1) = '*'
         return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
     end if
!
!--- The absolute value of the integer goes into S(ILO:IHI).
!
     ipos = ihi
!
!--- Find the last digit of IVAL, strip it off, and stick it into the string.
!
     do
       idig = mod ( ival, 10 )
       ival = ival / 10
       if ( ipos < ilo ) then
         do i = 1, ihi
           s(i:i) = '*'
         end do
         return
       end if
       call digit_to_ch ( idig, c )
       s(ipos:ipos) = c
       ipos = ipos - 1
       if ( ival == 0 ) then
         exit
       end if
     end do
!
!--- Shift the string to the left.
!
     s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
     s(ilo+ihi-ipos:ihi) = ' '
 
     return
     end SUBROUTINE i4_to_s_left
!
!************************************************************************************************
     SUBROUTINE makeVTU(output_unit)
!************************************************************************************************
     use User0, only:LprintParaviewNodeBased,LprintParaviewCellBased
!************************************************************************************************
     implicit none
!************************************************************************************************
     character (len = 20) s_nPoints, s_nCells
     integer :: Cellsize, i, j, output_unit, offsett
     double precision, allocatable :: CellID(:), NodeID(:)
!************************************************************************************************
!
     print*,'writing Paraview file ....please wait'
     Cellsize = 0
     do i=1,nCells
       Cellsize = Cellsize+Cells(i)%nPoints+1
     end do
!   
     call i4_to_s_left(nPoints,s_nPoints)
     call i4_to_s_left(nCells,s_nCells)
!
!--- Writing Headers
!
     write ( output_unit,'(a)')  & 
      '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//NEW_LINE('A')// &
      '<UnstructuredGrid>'//NEW_LINE('A')//'<Piece NumberOfPoints="'// trim (s_nPoints) // &
      '" NumberOfCells="' // trim (s_nCells) // '">'
!
!--- Writing Points
!
     write ( output_unit,'(a)') '<Points>'//NEW_LINE('A')// &
      '<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'
     do i=1,nPoints
       write ( output_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) Points(i,:)
     end do
     write ( output_unit,'(a)') '</DataArray>'//NEW_LINE('A')//'</Points>'
!
!--- writing Cells
!
     write ( output_unit,'(a)') '<Cells>'
!
!--- Connectivities
!
     write ( output_unit,'(a)') '<DataArray type="Int64" Name="connectivity" format="ascii">'
     do i=1,nCells
       write ( output_unit,*) (Cells(i)%PointsID(:)-1)
     end do
     write ( output_unit,'(a)') '</DataArray>'
!
!--- offsets
!
     write ( output_unit,'(a)') '<DataArray type="Int64" Name="offsets" format="ascii">'
         j = 0
         do i=1,nCells
            j = j + Cells(i)%nPoints
            write ( output_unit,'(i8)', advance = 'no') j    
         end do
     write ( output_unit,*)
     write ( output_unit,'(a)') '</DataArray>'
!
!--- types
!
     write ( output_unit,'(a)') '<DataArray type="Int64" Name="types" format="ascii">'
     do i=1,nCells
       write ( output_unit, '(a)', advance = 'no') '42 '
     end do
     write ( output_unit,*)
     write ( output_unit,'(a)') '</DataArray>'
!
!--- faces  
!
     write ( output_unit,'(a)') '<DataArray type="Int64" Name="faces" format="ascii">'
     do i=1,nCells
       write ( output_unit, *) Cells(i)%nFaces
       do j=1,Cells(i)%nFaces
         write ( output_unit,*) Faces(Cells(i)%FacesID(j))%nPoints, (Faces(Cells(i)%FacesID(j))%PointsID(:)-1)
       end do  
     end do
     write ( output_unit,'(a)') '</DataArray>'
!
!--- faces offsets
!
     write ( output_unit,'(a)') '<DataArray type="Int64" Name="faceoffsets" format="ascii">'
     offsett = 0
     do i=1,nCells
       do j = 1,Cells(i)%nFaces
         offsett = offsett + (Faces(Cells(i)%FacesID(j))%nPoints+1)
       end do
        offsett = offsett + 1
        write ( output_unit,'(i10)', advance = 'no') offsett
     end do
     write ( output_unit,*)
     write ( output_unit,'(a)') '</DataArray>'
     write ( output_unit,'(a)') '</Cells>' !printing Cells finishes here
!
     allocate(CellID(nCells))
     do i=1,nCells
       CellID(i)=i-1
     enddo
!
     allocate(NodeID(nPoints))
     do i=1,nPoints
       NodeID(i)=i-1
     enddo
!
!--- Printing Variables as point data
!
     write ( output_unit,'(a)') '<PointData Scalars="scalars">'
!
!--- Put here all point data that you want to print
!
     call printScalar(NodeID, 'NodeID', output_unit, 0)
     if(LprintParaviewNodeBased) call PrintVariablesBasedOnNodes(output_unit)
     write ( output_unit,'(a)') '</PointData>'
!   
!--- Printing Variables as cell data
!
     write ( output_unit,'(a)') '<CellData Scalars="scalars">'
!
!--- Put here all cell data that you want to print
!
     call printScalar(CellID, 'CellID', output_unit, 1)
     if(LprintParaviewCellBased) call PrintVariablesBasedOnCells(output_unit)
     write ( output_unit,'(a)') '</CellData>'
!
!--- vtu file closure
!
     write ( output_unit,'(a)') '</Piece>'//NEW_LINE('A')//'</UnstructuredGrid>' &
            //NEW_LINE('A')//'</VTKFile>'
!     
     end SUBROUTINE makeVTU     
!************************************************************************************************
     SUBROUTINE printScalar(Variable, NameV, VTU_unit, ArrayType)
!************************************************************************************************
     implicit none
!************************************************************************************************
     double precision, allocatable :: Variable(:)
     character(len=*), intent(in) :: NameV
     integer :: VTU_unit, i, ArrayType
!************************************************************************************************
!
!--- ArrayType = 0 >> point data
!--- ArrayType = 1 >> Cell data
!
     write ( VTU_unit,'(a)') '<DataArray type="Float64" Name="' &
             //trim (NameV)//'" Format="ascii">'
!    
     if (ArrayType==0)  then    
       do i=1,nPoints
         write(VTU_unit,*) Variable(i)
       end do
     else if (ArrayType == 1) then
       do i=1,nCells
         write(VTU_unit,*) Variable(i)
       end do
     else
       print*, "Error ArrayType = 0 >> point data; ArrayType = 1 >> Cell data"
     end if
!
     write ( VTU_unit,'(a)')
     write ( VTU_unit,'(a)') '</DataArray>'
!
     end SUBROUTINE printScalar
!************************************************************************************************
     SUBROUTINE printVector(Variable, DimV, NameV, VTU_unit, ArrayType)
!************************************************************************************************
     implicit none
!************************************************************************************************
     double precision, allocatable :: Variable(:,:)
     character(len=*), intent(in) :: NameV
     integer :: VTU_unit, i, DimV, ArrayType
     character (len = 8) :: s_DimV
!************************************************************************************************
!
     call i4_to_s_left(DimV,s_DimV)
!    
     write ( VTU_unit,'(a)') '<DataArray type="Float64" Name="'// &
          NameV//'" NumberOfComponents="'//trim(s_DimV)//'" Format="ascii">'
!    
     if (ArrayType==0)  then    
       do i=1,nPoints
         write(VTU_unit,*) Variable(i,:)
       end do
     else if (ArrayType == 1) then
       do i=1,nCells
         write(VTU_unit,*) Variable(i,:)
       end do
     else
       print*, "Error ArrayType = 0 >> point data; ArrayType = 1 >> Cell data"
     end if
!    
     write ( VTU_unit,'(a)')
     write ( VTU_unit,'(a)') '</DataArray>'
!    
     end SUBROUTINE printVector     
!************************************************************************************************
SUBROUTINE InitializeV
!
!************************************************************************************************
!
     use User0, only: PolyMeshDirectory,nInterPoints
     use Variables1, only: uVelocity,vVelocity,wVelocity
     use Geometry4, only: xc,yc,zc !,BFaceCentroidx,BFaceCentroidy,BFaceCentroidz
     use constants1, only: pi
     
     implicit none
     integer :: i,j
!
    type observedVelocities
       double precision :: xc, yc, zc !location x y z
       double precision :: uVelocity, vVelocity
    end type observedVelocities
    
    Type Interpolation
        !cell indexes of vertically interpolated elements
        integer, allocatable :: CellID(:)
        !Z coordinate of vertically interpolated elements
        double precision, allocatable :: Zc(:)
        !number of cells in this column cells
        integer :: nCells
    end Type Interpolation
    
    integer :: IOstatus, k, l, InterCellsID(nInterPoints)     
    double precision :: magnitude, angle, deltaC1, deltaC2, zc2(nCells), &
            Xmin, Xmax, Ymin, Ymax, rSqrd, Nu, Nv, D, zArray(nCells), &
            p = (0.4) !p is the power law constant        
    double precision, allocatable :: Difference2D(:)
    integer :: CellsArray(nCells)
    LOGICAL, allocatable :: mk(:) !used as mask while searching
    integer , allocatable :: ascendingID(:)
    integer :: nObservations
    type(observedVelocities), allocatable :: observedVelocity(:)
    type(Interpolation), allocatable :: VerInter(:) !for each observation point we have verInter set
    
    zc2 = zc
    
    open(unit = 7, file = trim(PolyMeshDirectory)//"/observations.txt")
    !counting the number of observations
    i = 0
    do 
      read (7,*,iostat=IOstatus) 
      if (IOstatus/=0) exit 
      i = i+1
    end do
    nObservations = i
    allocate(observedVelocity(nObservations))
    allocate(VerInter(nObservations))
    close(7)
    
    !Reading Observations
    open(unit = 7, file = trim(PolyMeshDirectory)//"/observations.txt")
    do i=1,nObservations
       read (7,*) observedVelocity(i)%xc, observedVelocity(i)%yc, &
                  observedVelocity(i)%zc, magnitude, angle
       angle = angle*pi/180.
       observedVelocity(i)%uVelocity = sin(angle)*magnitude
       observedVelocity(i)%vVelocity = cos(angle)*magnitude
    end do
    
    close(7)
    
    !Changing Altitudes to Hights
    call AltitudeToHight(xc,yc,zc,nCells,PolyMeshDirectory)
    
    allocate(Difference2D(nCells))
    allocate(mk(nCells))

    !Interpolating vertically
    do i=1,nObservations
        deltaC1 = infinity()
        deltaC2 = infinity()
        CellsArray = 0
        zArray = 0
        mk = .true.
        !Find the limits in x and y direction for intepolating according to cell size
        Difference2D = sqrt(abs(observedVelocity(i)%xc-xc)**2 &
                          + abs(observedVelocity(i)%yc-yc)**2)
        !l is the index of the closest cell to the measurment in xy directions
        l = minloc(Difference2D, 1, mk)
        Xmin = minVal(Points(Cells(l)%pointsID(:),1))
        Xmax = maxVal(Points(Cells(l)%pointsID(:),1))
        Ymin = minVal(Points(Cells(l)%pointsID(:),2))
        Ymax = maxVal(Points(Cells(l)%pointsID(:),2))
        deltaC1 = sqrt((Xmin-Xmax)**2+(Ymin-Ymax)**2)
        j = 0
        !Interpolating vertically
        do k = 1,nCells
            Difference2D = sqrt(abs(observedVelocity(i)%xc-xc)**2 &
                              + abs(observedVelocity(i)%yc-yc)**2)
            !l is the index of the closest cell to the measurment
             l = minloc(Difference2D,1,mk)
             mk(l) = .FALSE.
             deltaC2 = sqrt((observedVelocity(i)%xc-xc(l))**2 &
                          + (observedVelocity(i)%yc-yc(l))**2)
             if (deltaC2.gt.(1.2*deltaC1)) exit
             if (.not.((xc(l)<Xmin).OR.(xc(l)>Xmax).OR.(yc(l)<Ymin).OR.(yc(l)>Ymax))) then
                 !use the power law to assign velocity on index l
                 uVelocity(l) = observedVelocity(i)%uVelocity*(zc(l)/observedVelocity(i)%zc)**p
                 vVelocity(l) = observedVelocity(i)%vVelocity*(zc(l)/observedVelocity(i)%zc)**p
                 j = j+1
                 CellsArray(j) = l  !saving the index of the cell 
                 zArray(j) = zc2(l)  !saving the z of the cell
             end if
        end do
        allocate(VerInter(i)%CellID(j))
        allocate(VerInter(i)%Zc(j))
        VerInter(i)%CellID(:) = CellsArray(1:j)
        VerInter(i)%Zc(:) = zArray(1:j)
        VerInter(i)%nCells = j
    end do
    deallocate(mk)
    deallocate(Difference2D)
    zc = zc2
    !interpolating horizontally
    !allocate mk and ascendingID with same size of array that I'm searching
    allocate(mk(nObservations)) 
    allocate(ascendingID(nObservations))
    allocate(Difference2D(nObservations))
    do i=1,nCells
        mk = .true. 
        !finding the first nInterPoints closest observations
        Difference2D = sqrt(abs(observedVelocity(:)%xc-xc(i))**2 &
                          + abs(observedVelocity(:)%yc-yc(i))**2)
        do j = 1, nObservations
           !ascendingID is array of observations ID's sorted from nearest
           !e.g. u of nearest observation is observedVelocity(ascendingID(1))%uVelocity
           ascendingID(j) =  MINLOC(Difference2D, 1,mk)     
           mk(ascendingID(j)) = .FALSE.
        end do
        !search nearest arrays for cell with the nearest z
        do j=1,nInterPoints
           l = minloc(abs(VerInter(ascendingID(j))%zc-zc(i)),1) !l is the index of the cell in the column cells
           InterCellsID(j) = VerInter(ascendingID(j))%cellID(l)
        end do
        !interpolating
        Nu = 0 !numerator of u equation
        Nv = 0 !numerator of v equation
        D = 0  !denomenator
        do j=1,nInterPoints
           rSqrd = ((xc(InterCellsID(j))-xc(i))**2+(yc(InterCellsID(j))-yc(i))**2)
           if (rSqrd ==0 ) goto 20
           Nu = Nu + (uVelocity(InterCellsID(j))/rSqrd)
           Nv = Nv + (vVelocity(InterCellsID(j))/rSqrd)
           D = D + (1/rSqrd)  
!           write(90,*) InterCellsID(j), uVelocity(InterCellsID(j)) , Nu
        end do
        uVelocity(i) = Nu/D  
        vVelocity(i) = Nv/D
!        write(90,*) uVelocity(i), vVelocity(i)
        20 continue
    end do

     end SUBROUTINE InitializeV
!************************************************************************************************
     real FUNCTION infinity() 
!************************************************************************************************
     implicit none 
     real :: x 
     x = huge(1.) 
     infinity = x + x 
     end FUNCTION infinity     
!     
     end MODULE ReadpolyMesh
!     
!************************************************************************************************
     SUBROUTINE PrintVariablesBasedOnCells(VTK_unit)
!************************************************************************************************
     use ReadpolyMesh
     use User0, only: LSolveMomentum,LSolveContinuity,LSolveTurbulenceKineticEnergy,&
                     LSolveModifiedED,LSolveTurbulenceDissipationRate,LSolveTurbulenceV2Equation,&
                     LSolveTurbulenceZetaEquation,LSolveTurbulencefRelaxationEquation,&
                     LSolveTurbulenceSpecificDissipationRate,LSolveTurbulenceGammaEquation,&              
                     LSolveTurbulenceReynoldsThetaEquation,LSolveTurbulentKL,LSolveEnergy,&
                     TurbulenceModel,LTurbulentFlow,LBuoyancy,LtestTurbulenceModel,EnergyEquation,&
                     NumberOfrFieldsToSolve,LSolverField,NumberOfScalarsToSolve,LSolveScalar,&
                     LCompressible,LFreeSurfaceFlow,rFieldName,ScalarName,LUnsteady,LPrintGradients,&
                     LSolveLambdaELEEquation
     use Geometry1, only: NumberOfElements,NumbOfElementFaces,NumbOfElementNodes,ListOfElementNodes
     use Geometry3, only: NGlobalEFaces,NFacesTotal
     use Variables1, only: uVelocity,vVelocity,WVelocity,Pressure,TurbulentKE,ModifiedED,TurbulentED,&
                          TurbulentV2,TurbulentZeta,TfRelaxation,TurbulentOmega,TGamma,TGammaEff,&
                          TReTheta,TurbulentKL,TurbulenceProduction,TurbulenceProductionB,Temperature,&
                          Htotal,MachNumber,uVelGradx,uVelGrady,uVelGradz,vVelGradx,vVelGrady,vVelGradz,&
                          wVelGradx,wVelGrady,wVelGradz,PressGradx,PressGrady,PressGradz,&
                          TKEGradx,TKEGrady,TKEGradz,TEDGradx,TEDGrady,TEDGradz,TurbulentV2Gradx,&
                          TurbulentV2Grady,TurbulentV2Gradz,TurbulentZetaGradx,TurbulentZetaGrady,&
                          TurbulentZetaGradz,TfRelaxationGradx,TfRelaxationGrady,TfRelaxationGradz,&
                          TOmegaGradx,TOmegaGrady,TOmegaGradz,TGammaGradx,TGammaGrady,TGammaGradz,&
                          TReThetaGradx,TReThetaGrady,TReThetaGradz,TurbulentKLGradx,TurbulentKLGrady,&
                          TurbulentKLGradz,ModifiedEDGradx,ModifiedEDGrady,ModifiedEDGradz,TempGradx,TempGrady,&
                          TempGradz,HtotalGradx,HtotalGrady,HtotalGradz,LambdaELE,LambdaELEGradx,&
                          LambdaELEGrady,LambdaELEGradz,InitialVelDivergence,FinalVelDivergence
     use PhysicalProperties1, only: Density,Viscosity,TurbulentViscosity,SpecificHeat,ReferenceTemperature,&
                                   eDiffCoefficient
     use Turbulence1, only: TurbulentViscosityTl,TurbulentViscosityTs,ProductionKT,ProductionKL,ReT,&
                            StrainRate,Vorticity
     use ReferenceValues1, only: UInfinity
     use constants1, only: twothird,tiny
     use WallDistance1, only: WallDistance
     use VolumeOfFluid1, only: rField,rFieldGradx,rFieldGrady,rFieldGradz
     use Scalar1, only: Scalar,ScalarGradx,ScalarGrady,ScalarGradz
     use Tecplot1, only: yplusPlot,uplusPlot,ByplusPlot,BuplusPlot
!************************************************************************************************
     implicit none 
!************************************************************************************************
     character*10 :: Variable1
     integer :: iScalar,irField
     integer :: VTK_unit
     character*10 :: part1
     character*5 :: part2,part3,part4
     character*15 :: name1,name2,name3
     integer ::   namesize
     double precision, allocatable :: Momentum(:,:)
     double precision, allocatable :: ScalarToPrint(:)
1    format(1x,',"',A15,'" ')
!********************************************************************************************
!
     if(LUnsteady) LtestTurbulenceModel=.false.   !in order not to modify normal distance to wall
!    
     allocate(ScalarToPrint(NumberOfElements))
!
     if(LsolveMomentum.or.LSolveLambdaELEEquation) then    
       allocate(Momentum(NumberOfElements,3))
       Momentum(:,1)=uVelocity(:)
       Momentum(:,2)=vVelocity(:)
       Momentum(:,3)=WVelocity(:)
       call printVector(Momentum, 3, 'Velocity(m/s)', VTK_unit,1)
       deallocate(Momentum)
       ScalarToPrint=dsqrt(uVelocity*uVelocity+vVelocity*vVelocity+wVelocity*wVelocity)
       call printScalar(ScalarToPrint,'Velocity-Magnitude(m/s)', VTK_unit,1)
     endif
     if(LSolveLambdaELEEquation) then    
       call printScalar(LambdaELE,'Lambda-Eulerian-Lagrangian-Equation', VTK_unit,1)
       call printScalar(InitialVelDivergence,'Initial-velocity-divergence(1/s)', VTK_unit,1)
       call printScalar(FinalVelDivergence,'Final-velocity-divergence(1/s)', VTK_unit,1)
     endif
     if(LSolveContinuity) then
       call printScalar(Pressure,'Pressure(Pa)', VTK_unit,1)
     endif
     if(LSolveTurbulenceKineticEnergy) then
       call printScalar(TurbulentKE,'Turbulent-kinetic-energy(m2/s2)', VTK_unit,1)
     endif
     if(LSolveModifiedED) then
       call printScalar(ModifiedED,'Modified-eddy-diffusivity(m2/s)', VTK_unit,1)
     endif
     if(LSolveTurbulenceDissipationRate) then
       call printScalar(TurbulentED,'Turbulent-dissipation-rate(m2/s3)', VTK_unit,1)
     endif
     if(LSolveTurbulenceV2Equation) then
       call printScalar(TurbulentV2,'Turbulent-v2(m2/s2)', VTK_unit,1)
     endif
     if(LSolveTurbulenceZetaEquation) then
       call printScalar(TurbulentZeta,'Turbulent-zeta', VTK_unit,1)
     endif
     if(LSolveTurbulencefRelaxationEquation) then
       call printScalar(TfRelaxation,'f-relaxation(m2/s)', VTK_unit,1)
     endif
     if(LSolveTurbulenceSpecificDissipationRate) then
       call printScalar(TurbulentOmega,'Turbulent-specific-dissipation-rate(1/s)', VTK_unit,1)
     endif
     if(LSolveTurbulenceGammaEquation) then
       call printScalar(TGamma,'Intermittency(Turbulent-gama)', VTK_unit,1)
     endif
     if(TurbulenceModel.eq.'sstgamaretheta') then
       call printScalar(TGammaEff,'Intermittency-effective', VTK_unit,1)
     endif
     if(LSolveTurbulenceReynoldsThetaEquation) then
       call printScalar(TReTheta,'Momentum-thickness(turbulent-Re-Theta)', VTK_unit,1)
     endif
     if(LSolveTurbulentKL) then
       call printScalar(TurbulentKL,'Laminar-kinetic-energy(m2/s2)', VTK_unit,1)
     endif
     if(LTurbulentFlow.and.TurbulenceModel.eq.'kklomega') then
       ScalarToPrint=TurbulentKE+TurbulentKL
       call printScalar(ScalarToPrint,'Total-fluctuation-energy(m2/s2)', VTK_unit,1)
     endif
     if(LSolveMomentum) then
       call printScalar(Viscosity,'Laminar-viscosity(Pa.s)', VTK_unit,1)
     endif
     if(LTurbulentFlow) then
       call printScalar(TurbulentViscosity,'Turbulent-viscosity(Pa.s)', VTK_unit,1)
     endif
     if(LTurbulentFlow.and.TurbulenceModel.eq.'kklomega') then
       ScalarToPrint=Density*TurbulentViscosityTl
       call printScalar(ScalarToPrint,'Turbulent-viscosity-(large-scale)(Pa.s)', VTK_unit,1)
       ScalarToPrint=Density*TurbulentViscosityTs
       call printScalar(ScalarToPrint,'Turbulent-viscosity(small-scale)(Pa.s)', VTK_unit,1)
     endif
     if(LTurbulentFlow) then
       ScalarToPrint=Viscosity+TurbulentViscosity
       call printScalar(ScalarToPrint,'Effective-viscosity(Pa.s)', VTK_unit,1)
       ScalarToPrint=TurbulentViscosity/Viscosity
       call printScalar(ScalarToPrint,'Turbulent-viscosity-ratio', VTK_unit,1)
     endif
     if(LSolveTurbulenceKineticEnergy) then
       if(TurbulenceModel.eq.'kklomega') then         
         ScalarToPrint=100.*dsqrt(twothird*(TurbulentKE+TurbulentKL))/dmax1(Uinfinity,tiny)
       else
         ScalarToPrint=100.*dsqrt(twothird*TurbulentKE)/dmax1(Uinfinity,tiny)
       endif
       call printScalar(ScalarToPrint,'Turbulent-intensity(%)', VTK_unit,1)
     endif
     if(LTurbulentFlow) then
       if(TurbulenceModel.ne.'spalartallmaras'.and.TurbulenceModel.ne.'wrayagarwal') then
         if(TurbulenceModel.eq.'kklomega') then
           ScalarToPrint=Density*ProductionKT   ! Production of k
           call printScalar(ScalarToPrint,'Production-of-k(Kg/(m.s3))', VTK_unit,1)
           ScalarToPrint=Density*ProductionKL   ! Production of laminar k
           call printScalar(ScalarToPrint,'Production-of-laminar-k(Kg/(m.s3))', VTK_unit,1)
         elseif(TurbulenceModel.eq.'sstgamaretheta') then
           call printScalar(TurbulenceProduction,'Production-of-k(Kg/(m.s3))', VTK_unit,1)
         else       
           call printScalar(TurbulenceProduction,'Production-of-k(Kg/(m.s3))', VTK_unit,1)
         endif
       endif
       if(LBuoyancy) then 
           call printScalar(TurbulenceProductionB,'Production-of-k-by-buoyancy(Kg/(m.s3))', VTK_unit,1)
       endif
     endif
     if(LTurbulentFlow.and.TurbulenceModel.ne.'wrayagarwal'.and.TurbulenceModel.ne.'spalartallmaras') then 
       call CalculateTurbulentReynoldsNumber
       call printScalar(ReT,'Turbulent-Reynolds-number', VTK_unit,1)
     endif
     if(LTurbulentFlow.and.LtestTurbulenceModel) then
       call CalculateYplusPlotting
       call printScalar(yplusPlot,'yplus', VTK_unit,1)
       call printScalar(uplusPlot,'uplus', VTK_unit,1)
       deallocate(yplusPlot)
       deallocate(uplusPlot)
       deallocate(ByplusPlot)
       deallocate(BuplusPlot)

     endif
     if(LTurbulentFlow) then
       call printScalar(WallDistance,'Normal-distance-to-wall(m)', VTK_unit,1)
     endif
     if(LSolveEnergy) then
       call printScalar(Temperature,'Temperature(K)', VTK_unit,1)
       ScalarToPrint=Temperature+0.5*(uVelocity*uVelocity+vVelocity*vVelocity+wVelocity*wVelocity)/SpecificHeat   
       call printScalar(ScalarToPrint,'Total-Temperature(K)', VTK_unit,1)
       ScalarToPrint=SpecificHeat*(Temperature-ReferenceTemperature)
       call printScalar(ScalarToPrint,'Static-enthalpy(J/Kg)', VTK_unit,1)
     endif
     if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
       call printScalar(Htotal,'Total-enthalpy(J/Kg)', VTK_unit,1)
     endif
     if(LSolveEnergy.and.LTurbulentFlow) then
       Variable1='temp'
       call CalculateEffectiveDiffusionCoefficient(Variable1)
       call printScalar(eDiffCoefficient,'Effective-thermal-conductivity(W/(m.K))', VTK_unit,1)
     endif
!     
      do irField=1,NumberOfrFieldsToSolve
        if(LSolverField(irField)) then
          ScalarToPrint(:)=rField(:,irField)  !rfield
          call printScalar(ScalarToPrint,rFieldName(irField), VTK_unit,1)
        endif
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LSolveScalar(iScalar)) then
          ScalarToPrint(:)=Scalar(:,iScalar)  !scalar
          call printScalar(ScalarToPrint,ScalarName(iScalar), VTK_unit,1)
        endif
      enddo
      if(LCompressible) then
        call printScalar(Density,'Density(Kg/m3)', VTK_unit,1)
      endif
      if(.not.LCompressible) then
        if(LFreeSurfaceFlow.and.LSolveMomentum) then
          call printScalar(Density,'Density(Kg/m3)', VTK_unit,1)
        endif
      endif
      if(LCompressible) then
        call printScalar(MachNumber,'Mach-Number', VTK_unit,1)
      endif
!
      if(LPrintGradients) then
!
        if(LSolveMomentum) then
          call printScalar(uVelGradx,'uVelGradx(1/s)', VTK_unit,1)
          call printScalar(uVelGrady,'uVelGrady(1/s)', VTK_unit,1)
          call printScalar(uVelGradz,'uVelGradz(1/s)', VTK_unit,1)
          call printScalar(vVelGradx,'vVelGradx(1/s)', VTK_unit,1)
          call printScalar(vVelGrady,'vVelGrady(1/s)', VTK_unit,1)
          call printScalar(vVelGradz,'vVelGradz(1/s)', VTK_unit,1)
          call printScalar(wVelGradx,'wVelGradx(1/s)', VTK_unit,1)
          call printScalar(wVelGrady,'wVelGrady(1/s)', VTK_unit,1)
          call printScalar(wVelGradz,'wVelGradz(1/s)', VTK_unit,1)
          if(.not.LTurbulentFlow) then
            call AllocateStrainRateTensor
            call CalculateStrainRateTensor
          endif
          call printScalar(StrainRate,'Strain-rate(1/s)', VTK_unit,1)
          if(.not.LTurbulentFlow) call deAllocateStrainRateTensor
          if(.not.LTurbulentFlow) then 
            call AllocateVorticityTensor
            call CalculateVorticityTensor
          endif
          call printScalar(Vorticity,'Vorticity(1/s)', VTK_unit,1)
          if(.not.LTurbulentFlow) call deAllocateVorticityTensor
        endif
        if(LSolveLambdaELEEquation) then
          call printScalar(uVelGradx,'uVelGradx(1/s)', VTK_unit,1)
          call printScalar(uVelGrady,'uVelGrady(1/s)', VTK_unit,1)
          call printScalar(uVelGradz,'uVelGradz(1/s)', VTK_unit,1)
          call printScalar(vVelGradx,'vVelGradx(1/s)', VTK_unit,1)
          call printScalar(vVelGrady,'vVelGrady(1/s)', VTK_unit,1)
          call printScalar(vVelGradz,'vVelGradz(1/s)', VTK_unit,1)
          call printScalar(wVelGradx,'wVelGradx(1/s)', VTK_unit,1)
          call printScalar(wVelGrady,'wVelGrady(1/s)', VTK_unit,1)
          call printScalar(wVelGradz,'wVelGradz(1/s)', VTK_unit,1)
          call printScalar(LambdaELEGradx,'LambdaGradx', VTK_unit,1)
          call printScalar(LambdaELEGrady,'LambdaGrady', VTK_unit,1)
          call printScalar(LambdaELEGradz,'LambdaGradz', VTK_unit,1)
        endif
        if(LSolveContinuity) then
          call printScalar(PressGradx,'PressureGradx(Pa/m)', VTK_unit,1)
          call printScalar(PressGrady,'PressureGrady(Pa/m)', VTK_unit,1)
          call printScalar(PressGradz,'PressureGradz(Pa/m)', VTK_unit,1)
        endif
        if(LSolveTurbulenceKineticEnergy) then
          call printScalar(TKEGradx,'TKEGradx(m/s2)', VTK_unit,1)
          call printScalar(TKEGrady,'TKEGrady(m/s2)', VTK_unit,1)
          call printScalar(TKEGradz,'TKEGradz(m/s2)', VTK_unit,1)
        endif
        if(LSolveTurbulenceDissipationRate) then
          call printScalar(TEDGradx,'TEDGradx(m/s2)', VTK_unit,1)
          call printScalar(TEDGrady,'TEDGrady(m/s2)', VTK_unit,1)
          call printScalar(TEDGradz,'TEDGradz(m/s2)', VTK_unit,1)
        endif
        if(LSolveTurbulenceV2Equation) then
          call printScalar(TurbulentV2Gradx,'v2Gradx(m/s2)', VTK_unit,1)
          call printScalar(TurbulentV2Grady,'v2Grady(m/s2)', VTK_unit,1)
          call printScalar(TurbulentV2Gradz,'v2Gradz(m/s2)', VTK_unit,1)
        endif
        if(LSolveTurbulenceZetaEquation) then
          call printScalar(TurbulentZetaGradx,'ZetaGradx(m/s2)', VTK_unit,1)
          call printScalar(TurbulentZetaGrady,'ZetaGrady(m/s2)', VTK_unit,1)
          call printScalar(TurbulentZetaGrady,'ZetaGradz(m/s2)', VTK_unit,1)
        endif
        if(LSolveTurbulencefRelaxationEquation) then
          call printScalar(TfRelaxationGradx,'fGradx(m/s2)', VTK_unit,1)
          call printScalar(TfRelaxationGrady,'fGrady(m/s2)', VTK_unit,1)
          call printScalar(TfRelaxationGradz,'fGradz(m/s2)', VTK_unit,1)
        endif
        if(LSolveTurbulenceSpecificDissipationRate) then
          call printScalar(TOmegaGradx,'TOmegaGradx(1/(m.s))', VTK_unit,1)
          call printScalar(TOmegaGrady,'TOmegaGrady(1/(m.s))', VTK_unit,1)
          call printScalar(TOmegaGradz,'TOmegaGradz(1/(m.s))', VTK_unit,1)
        endif
        if(LSolveTurbulenceGammaEquation) then
          call printScalar(TGammaGradx,'TGammaGradx(1/m)', VTK_unit,1)
          call printScalar(TGammaGrady,'TGammaGrady(1/m)', VTK_unit,1)
          call printScalar(TGammaGradz,'TGammaGradz(1/m)', VTK_unit,1)
        endif
        if(LSolveTurbulenceReynoldsThetaEquation) then
          call printScalar(TReThetaGradx,'TReThetaGradx(1/m)', VTK_unit,1)
          call printScalar(TReThetaGrady,'TReThetaGrady(1/m)', VTK_unit,1)
          call printScalar(TReThetaGradz,'TReThetaGradz(1/m)', VTK_unit,1)
        endif
        if(LSolveTurbulentKL) then
          call printScalar(TurbulentKLGradx,'TKLGradx(m/s2)', VTK_unit,1)
          call printScalar(TurbulentKLGrady,'TKLGrady(m/s2)', VTK_unit,1)
          call printScalar(TurbulentKLGradz,'TKLGradz(m/s2)', VTK_unit,1)
        endif
        if(LSolveModifiedED) then
          call printScalar(ModifiedEDGradx,'MEDGradx(m/s)', VTK_unit,1)
          call printScalar(ModifiedEDGrady,'MEDGrady(m/s)', VTK_unit,1)
          call printScalar(ModifiedEDGradz,'MEDGradz(m/s)', VTK_unit,1)
        endif
        if(LSolveEnergy) then
          call printScalar(TempGradx,'TempGradx(K/m)', VTK_unit,1)
          call printScalar(TempGrady,'TempGrady(K/m)', VTK_unit,1)
          call printScalar(TempGradz,'TempGradz(K/m)', VTK_unit,1)
        endif
        if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
          call printScalar(HtotalGradx,'HTotalGradx(J/(Kg.m))', VTK_unit,1)
          call printScalar(HtotalGrady,'HTotalGrady(J/(Kg.m))', VTK_unit,1)
          call printScalar(HtotalGradz,'HTotalGradz(J/(Kg.m))', VTK_unit,1)
        endif
        do irField=1,NumberOfrFieldsToSolve
          if(LSolverField(irField)) then
            part1=trim(rFieldName(irField))
            part2="Gradx"
            part3="Grady"
            part4="Gradz"
            namesize=len(part1)+5
            name1=trim(part1)//part2
            name2=trim(part1)//part3
            name3=trim(part1)//part4
            write(12,1,advance='no') name1(1:namesize)
            write(12,1,advance='no') name2(1:namesize)
            write(12,1,advance='no') name3(1:namesize)
            ScalarToPrint=rFieldGradx(:,irField) 
            call printScalar(ScalarToPrint,name1, VTK_unit,1)
            ScalarToPrint=rFieldGrady(:,irField)
            call printScalar(ScalarToPrint,name2, VTK_unit,1)
            ScalarToPrint=rFieldGradz(:,irField)
            call printScalar(ScalarToPrint,name3, VTK_unit,1)
          endif
        enddo
        do iScalar=1,NumberOfScalarsToSolve
          if(LSolveScalar(iScalar)) then
            part1=trim(ScalarName(iScalar))
            part2="Gradx"
            part3="Grady"
            part4="Gradz"
            namesize=len(part1)+5
            name1=trim(part1)//part2
            name2=trim(part1)//part3
            name3=trim(part1)//part4
            write(12,1,advance='no') name1(1:namesize)
            write(12,1,advance='no') name2(1:namesize)
            write(12,1,advance='no') name3(1:namesize)
            ScalarToPrint=ScalarGradx(:,iScalar) 
            call printScalar(ScalarToPrint,name1, VTK_unit,1)
            ScalarToPrint=ScalarGrady(:,iScalar) 
            call printScalar(ScalarToPrint,name2, VTK_unit,1)
            ScalarToPrint=ScalarGradz(:,iScalar) 
            call printScalar(ScalarToPrint,name3, VTK_unit,1)
          endif
        enddo
      endif
!      
      deallocate(ScalarToPrint)
!    
      return
    end SUBROUTINE PrintVariablesBasedOnCells
!     
!************************************************************************************************
     SUBROUTINE PrintVariablesBasedOnNodes(VTK_unit)
!************************************************************************************************
     use ReadpolyMesh
     use User0, only: LSolveMomentum,LSolveContinuity,LSolveTurbulenceKineticEnergy,&
                      LSolveModifiedED,LSolveTurbulenceDissipationRate,LSolveTurbulenceV2Equation,&
                      LSolveTurbulenceZetaEquation,LSolveTurbulencefRelaxationEquation,&
                      LSolveTurbulenceSpecificDissipationRate,LSolveTurbulenceGammaEquation,&              
                      LSolveTurbulenceReynoldsThetaEquation,LSolveTurbulentKL,LSolveEnergy,&
                      TurbulenceModel,LTurbulentFlow,LBuoyancy,LtestTurbulenceModel,EnergyEquation,&
                      NumberOfrFieldsToSolve,LSolverField,NumberOfScalarsToSolve,LSolveScalar,&
                      LCompressible,LFreeSurfaceFlow,rFieldName,ScalarName,LUnsteady,LPrintGradients,&
                     LSolveLambdaELEEquation
     use Geometry1, only: NumberOfElements,NumbOfElementFaces,NumbOfElementNodes,ListOfElementNodes,&
                          NumberOfNodes,NumberOfBCSets
     use Geometry3, only: NGlobalEFaces,NFacesTotal,NBFacesMax,NBFaces,NBFaceOwner
     use Variables1, only: uVelocity,vVelocity,wVelocity,Pressure,TurbulentKE,ModifiedED,TurbulentED,&
                           TurbulentV2,TurbulentZeta,TfRelaxation,TurbulentOmega,TGamma,TGammaEff,&
                           TReTheta,TurbulentKL,TurbulenceProduction,TurbulenceProductionB,Temperature,&
                           Htotal,MachNumber,BuVelocity,BvVelocity,BwVelocity,BPressure,BTurbulentKE,BModifiedED,&
                           BTurbulentED,BTurbulentV2,BTurbulentZeta,BTfRelaxation,BTurbulentOmega,BTGamma,&
                           BTReTheta,BTurbulentKL,BTurbulenceProduction,BTurbulenceProductionB,BTemperature,&
                           BHtotal,BMachNumber,uVelGradx,uVelGrady,uVelGradz,vVelGradx,vVelGrady,vVelGradz,&
                           wVelGradx,wVelGrady,wVelGradz,PressGradx,PressGrady,PressGradz,&
                           TKEGradx,TKEGrady,TKEGradz,TEDGradx,TEDGrady,TEDGradz,TurbulentV2Gradx,&
                           TurbulentV2Grady,TurbulentV2Gradz,TurbulentZetaGradx,TurbulentZetaGrady,&
                           TurbulentZetaGradz,TfRelaxationGradx,TfRelaxationGrady,TfRelaxationGradz,&
                           TOmegaGradx,TOmegaGrady,TOmegaGradz,TGammaGradx,TGammaGrady,TGammaGradz,&
                           TReThetaGradx,TReThetaGrady,TReThetaGradz,TurbulentKLGradx,TurbulentKLGrady,&
                           TurbulentKLGradz,ModifiedEDGradx,ModifiedEDGrady,ModifiedEDGradz,TempGradx,TempGrady,&
                           TempGradz,HtotalGradx,HtotalGrady,HtotalGradz,BuVelGradx,BuVelGrady,BuVelGradz,&
                           BvVelGradx,BvVelGrady,BvVelGradz,BwVelGradx,BwVelGrady,BwVelGradz,BPressGradx,&
                           BPressGrady,BPressGradz,BTKEGradx,BTKEGrady,BTKEGradz,BTEDGradx,BTEDGrady,BTEDGradz,&
                           BTurbulentV2Gradx,BTurbulentV2Grady,BTurbulentV2Gradz,BTurbulentZetaGradx,&
                           BTurbulentZetaGrady,BTurbulentZetaGradz,BTfRelaxationGradx,BTfRelaxationGrady,&
                           BTfRelaxationGradz,BTOmegaGradx,BTOmegaGrady,BTOmegaGradz,BTGammaGradx,BTGammaGrady,&
                           BTGammaGradz,BTReThetaGradx,BTReThetaGrady,BTReThetaGradz,BTurbulentKLGradx,BTurbulentKLGrady,&
                           BTurbulentKLGradz,BModifiedEDGradx,BModifiedEDGrady,BModifiedEDGradz,BTempGradx,BTempGrady,&
                           BTempGradz,BHtotalGradx,BHtotalGrady,BHtotalGradz,LambdaELE,BLambdaELE,LambdaELEGradx,&
                           LambdaELEGrady,LambdaELEGradz,BLambdaELEGradx,BLambdaELEGrady,BLambdaELEGradz,&
                           InitialVelDivergence,FinalVelDivergence
     use PhysicalProperties1, only: Density,Viscosity,TurbulentViscosity,SpecificHeat,ReferenceTemperature,&
                                    eDiffCoefficient,BDensity,BViscosity,BTurbulentViscosity,BSpecificHeat,&
                                    BeDiffCoefficient
     use Turbulence1, only: TurbulentViscosityTl,TurbulentViscosityTs,ProductionKT,ProductionKL,ReT,&
                            BTurbulentViscosityTl,BTurbulentViscosityTs,BReT,StrainRate,Vorticity,BStrainRate,BVorticity
     use ReferenceValues1, only: UInfinity
     use constants1, only: twothird,tiny
     use WallDistance1, only: WallDistance,BWallDistance
     use VolumeOfFluid1, only: rField,BrField,rFieldGradx,rFieldGrady,rFieldGradz,BrFieldGradx,BrFieldGrady,BrFieldGradz
     use Scalar1, only: Scalar,BScalar,ScalarGradx,ScalarGrady,ScalarGradz,BScalarGradx,BScalarGrady,BScalarGradz
     use Tecplot1, only: dummyE,BdummyE,dummyN,yplusPlot,ByplusPlot,uplusPlot,BuplusPlot,&
                         SumInvd,SumPhiInvd
!
!************************************************************************************************
     implicit none 
!************************************************************************************************
     character*10 :: Variable1
     integer :: i,j,k,iScalar,irField
     integer :: VTK_unit
     character*10 :: part1
     character*5 :: part2,part3,part4
     character*15 :: name1,name2,name3
     integer ::   namesize
     double precision, allocatable :: Momentum(:,:)
  1  format(1x,',"',A15,'" ')
!************************************************************************************************
!
     if(LUnsteady) LtestTurbulenceModel=.false.   !in order not to modify normal distance to wall
!
     allocate(dummyE(NumberOfElements))
     allocate(BdummyE(NumberOfBCSets,NBFacesMax))
     allocate(dummyN(NumberOfNodes))
     allocate(SumInvd(NumberOfNodes))
     allocate(SumPhiInvd(NumberOfNodes))
!
     if(LsolveMomentum.or.LSolveLambdaELEEquation) then    
       allocate(Momentum(NumberOfNodes,3))
       dummyE=uVelocity
       BdummyE=BuVelocity
       call InterpolateToNodesFromElements
       Momentum(:,1)=dummyN(:)
       dummyE=vVelocity
       BdummyE=BvVelocity
       call InterpolateToNodesFromElements
       Momentum(:,2)=dummyN(:)
       dummyE=wVelocity
       BdummyE=BwVelocity
       call InterpolateToNodesFromElements
       Momentum(:,3)=dummyN(:)
       call printVector(Momentum, 3, 'Velocity(m/s)', VTK_unit,0)
       deallocate(Momentum)
       dummyE=dsqrt(uVelocity*uVelocity+vVelocity*vVelocity+wVelocity*wVelocity)
       BdummyE=dsqrt(BuVelocity*BuVelocity+BvVelocity*BvVelocity+BwVelocity*BwVelocity)
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Velocity-Magnitude(m/s)', VTK_unit,0)
     endif
     if(LSolveLambdaELEEquation) then    
       dummyE=LambdaELE
       BdummyE=BLambdaELE
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Lambda-Eulerian-Lagrangian-Equation', VTK_unit,0)
       dummyE=InitialVelDivergence
       do i=1,NumberOfBCSets
         do j=1,NBFaces(i)
           k=NBFaceOwner(i,j)
           BdummyE(i,j)=InitialVelDivergence(k)
         enddo
       enddo
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Initial-velocity-divergence(1/s)', VTK_unit,0)
       dummyE=FinalVelDivergence
       do i=1,NumberOfBCSets
         do j=1,NBFaces(i)
           k=NBFaceOwner(i,j)
           BdummyE(i,j)=FinalVelDivergence(k)
         enddo
       enddo
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Final-velocity-divergence(1/s)', VTK_unit,0)
     endif
     if(LSolveContinuity) then
       dummyE=Pressure
       BdummyE=BPressure
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Pressure(Pa)', VTK_unit,0)
     endif
     if(LSolveTurbulenceKineticEnergy) then
       dummyE=TurbulentKE
       BdummyE=BTurbulentKE
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-kinetic-energy(m2/s2)', VTK_unit,0)
     endif
     if(LSolveModifiedED) then
       dummyE=ModifiedED
       BdummyE=BModifiedED
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Modified-eddy-diffusivity(m2/s)', VTK_unit,0)
     endif
     if(LSolveTurbulenceDissipationRate) then
       dummyE=TurbulentED
       BdummyE=BTurbulentED
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-dissipation-rate(m2/s3)', VTK_unit,0)
     endif
     if(LSolveTurbulenceV2Equation) then
       dummyE=TurbulentV2
       BdummyE=BTurbulentV2
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-v2(m2/s2)', VTK_unit,0)
     endif
     if(LSolveTurbulenceZetaEquation) then
       dummyE=TurbulentZeta
       BdummyE=BTurbulentZeta
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-zeta', VTK_unit,0)
     endif
     if(LSolveTurbulencefRelaxationEquation) then
       dummyE=TfRelaxation
       BdummyE=BTfRelaxation
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'f-relaxation(m2/s)', VTK_unit,0)
     endif
     if(LSolveTurbulenceSpecificDissipationRate) then
       dummyE=TurbulentOmega
       BdummyE=BTurbulentOmega
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-specific-dissipation-rate(1/s)', VTK_unit,0)
     endif
     if(LSolveTurbulenceGammaEquation) then
       dummyE=TGamma
       BdummyE=BTGamma
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Intermittency(Turbulent-gama)', VTK_unit,0)
     endif
     if(TurbulenceModel.eq.'sstgamaretheta') then
       dummyE=TGammaEff
       do i=1,NumberOfBCSets
         do j=1,NBFaces(i)
           k=NBFaceOwner(i,j)
           BdummyE(i,j)=TGammaEff(k)
         enddo
       enddo
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Intermittency-effective', VTK_unit,0)
     endif
     if(LSolveTurbulenceReynoldsThetaEquation) then
       dummyE=TReTheta
       BdummyE=BTReTheta
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Momentum-thickness(turbulent-Re-Theta)', VTK_unit,0)
     endif
     if(LSolveTurbulentKL) then
       dummyE=TurbulentKL
       BdummyE=BTurbulentKL
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Laminar-kinetic-energy(m2/s2)', VTK_unit,0)
     endif
     if(LTurbulentFlow.and.TurbulenceModel.eq.'kklomega') then
       dummyE=TurbulentKE+TurbulentKL
       BdummyE=BTurbulentKE+BTurbulentKL
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Total-fluctuation-energy(m2/s2)', VTK_unit,0)
     endif
     if(LSolveMomentum) then
       dummyE=Viscosity
       BdummyE=BViscosity
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Laminar-viscosity(Pa.s)', VTK_unit,0)
     endif
     if(LTurbulentFlow) then
       dummyE=TurbulentViscosity
       BdummyE=BTurbulentViscosity
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-viscosity(Pa.s)', VTK_unit,0)
     endif
     if(LTurbulentFlow.and.TurbulenceModel.eq.'kklomega') then
       dummyE=Density*TurbulentViscosityTl
       BdummyE=BDensity*BTurbulentViscosityTl
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-viscosity-(large-scale)(Pa.s)', VTK_unit,0)
       dummyE=Density*TurbulentViscosityTs
       BdummyE=BDensity*BTurbulentViscosityTs
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-viscosity(small-scale)(Pa.s)', VTK_unit,0)
     endif
     if(LTurbulentFlow) then
       dummyE=Viscosity+TurbulentViscosity
       BdummyE=BViscosity+BTurbulentViscosity
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Effective-viscosity(Pa.s)', VTK_unit,0)
       dummyE=TurbulentViscosity/Viscosity
       BdummyE=BTurbulentViscosity/BViscosity
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-viscosity-ratio', VTK_unit,0)
     endif
     if(LSolveTurbulenceKineticEnergy) then
       if(TurbulenceModel.eq.'kklomega') then         
         dummyE=100.*dsqrt(twothird*(TurbulentKE+TurbulentKL))/dmax1(Uinfinity,tiny)
         BdummyE=100.*dsqrt(twothird*(BTurbulentKE+BTurbulentKL))/dmax1(Uinfinity,tiny)
         call InterpolateToNodesFromElements
       else
         dummyE=100.*dsqrt(twothird*TurbulentKE)/dmax1(Uinfinity,tiny)
         BdummyE=100.*dsqrt(twothird*BTurbulentKE)/dmax1(Uinfinity,tiny)
         call InterpolateToNodesFromElements
       endif
       call printScalar(dummyN,'Turbulent-intensity(%)', VTK_unit,0)
     endif
     if(LTurbulentFlow) then
       if(TurbulenceModel.ne.'spalartallmaras'.and.TurbulenceModel.ne.'wrayagarwal') then
         if(TurbulenceModel.eq.'kklomega') then
           dummyE=Density*ProductionKT   ! Production of k
           do i=1,NumberOfBCSets
             do j=1,NBFaces(i)
               k=NBFaceOwner(i,j)
               BdummyE(i,j)=Density(k)*ProductionKT(k)
             enddo
           enddo
           call InterpolateToNodesFromElements
           call printScalar(dummyN,'Production-of-k(Kg/(m.s3))', VTK_unit,0)
           dummyE=Density*ProductionKL   ! Production of laminar k
           do i=1,NumberOfBCSets
             do j=1,NBFaces(i)
               k=NBFaceOwner(i,j)
               BdummyE(i,j)=Density(k)*ProductionKL(k)
             enddo
           enddo
           call InterpolateToNodesFromElements
           call printScalar(dummyN,'Production-of-laminar-k(Kg/(m.s3))', VTK_unit,0)
         elseif(TurbulenceModel.eq.'sstgamaretheta') then
           dummyE=TurbulenceProduction   ! Production of k
           do i=1,NumberOfBCSets
             do j=1,NBFaces(i)
               k=NBFaceOwner(i,j)
               BdummyE(i,j)=TurbulenceProduction(k)
             enddo
           enddo
           call InterpolateToNodesFromElements
           call printScalar(dummyN,'Production-of-k(Kg/(m.s3))', VTK_unit,0)
         else       
           dummyE=TurbulenceProduction   ! Production of k
           do i=1,NumberOfBCSets
             do j=1,NBFaces(i)
               k=NBFaceOwner(i,j)
               BdummyE(i,j)=TurbulenceProduction(k)
             enddo
           enddo
           call InterpolateToNodesFromElements
           call printScalar(dummyN,'Production-of-k(Kg/(m.s3))', VTK_unit,0)
         endif
       endif
       if(LBuoyancy) then 
         dummyE=TurbulenceProductionB
         BdummyE=BTurbulenceProductionB
         call InterpolateToNodesFromElements
         call printScalar(dummyN,'Production-of-k-by-buoyancy(Kg/(m.s3))', VTK_unit,0)
       endif
     endif
     if(LTurbulentFlow.and.TurbulenceModel.ne.'wrayagarwal'.and.TurbulenceModel.ne.'spalartallmaras') then 
       call CalculateTurbulentReynoldsNumber
       dummyE=ReT
       BdummyE=BReT
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Turbulent-Reynolds-number', VTK_unit,0)
     endif
     if(LTurbulentFlow.and.LtestTurbulenceModel) then
       call CalculateYplusPlotting
       dummyE=yplusPlot
       BdummyE=ByplusPlot
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'yplus', VTK_unit,0)
       dummyE=uplusPlot
       BdummyE=BuplusPlot
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'uplus', VTK_unit,0)
       deallocate(yplusPlot)
       deallocate(uplusPlot)
       deallocate(ByplusPlot)
       deallocate(BuplusPlot)
     endif
     if(LTurbulentFlow) then
       dummyE=WallDistance
       BdummyE=BWallDistance
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Normal-distance-to-wall(m)', VTK_unit,0)
     endif
     if(LSolveEnergy) then
       dummyE=Temperature
       BdummyE=BTemperature
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Temperature(K)', VTK_unit,0)
       dummyE=Temperature+0.5*(uVelocity*uVelocity+vVelocity*vVelocity+wVelocity*wVelocity)/SpecificHeat
       BdummyE=BTemperature+0.5*(BuVelocity*BuVelocity+BvVelocity*BvVelocity+BwVelocity*BwVelocity)/BSpecificHeat
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Total-Temperature(K)', VTK_unit,0)
       dummyE=SpecificHeat*(Temperature-ReferenceTemperature)
       BdummyE=BSpecificHeat*(bTemperature-ReferenceTemperature)
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Static-enthalpy(J/Kg)', VTK_unit,0)
     endif
     if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
       dummyE=Htotal
       BdummyE=BHtotal
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Total-enthalpy(J/Kg)', VTK_unit,0)
     endif
     if(LSolveEnergy.and.LTurbulentFlow) then
       Variable1='temp'
       call CalculateEffectiveDiffusionCoefficient(Variable1)
       dummyE=eDiffCoefficient
       BdummyE=BeDiffCoefficient
       call InterpolateToNodesFromElements
       call printScalar(dummyN,'Effective-thermal-conductivity(W/(m.K))', VTK_unit,0)
     endif
!     
      do irField=1,NumberOfrFieldsToSolve
        if(LSolverField(irField)) then
          dummyE(:)=rField(:,irField)
          BdummyE(:,:)=BrField(:,:,irField)
          call InterpolateToNodesFromElements
          call printScalar(dummyN,rFieldName(irField), VTK_unit,0)
        endif
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LSolveScalar(iScalar)) then
          dummyE(:)=Scalar(:,iScalar)
          BdummyE(:,:)=BScalar(:,:,iScalar)
          call InterpolateToNodesFromElements
          call printScalar(dummyN,ScalarName(iScalar), VTK_unit,0)
        endif
      enddo
      if(LCompressible) then
        dummyE=Density
        BdummyE=BDensity
        call InterpolateToNodesFromElements
        call printScalar(dummyN,'Density(Kg/m3)', VTK_unit,0)
      endif
      if(.not.LCompressible) then
        if(LFreeSurfaceFlow.and.LSolveMomentum) then
          dummyE=Density
          BdummyE=BDensity
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'Density(Kg/m3)', VTK_unit,0)
        endif
      endif
      if(LCompressible) then
        dummyE=MachNumber
        BdummyE=BMachNumber
        call InterpolateToNodesFromElements
        call printScalar(dummyN,'Mach-Number', VTK_unit,0)
      endif
!
      if(LPrintGradients) then
!
        if(LSolveMomentum) then
          dummyE=uVelGradx
          BdummyE=BuVelGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'uVelGradx(1/s)', VTK_unit,0)
          dummyE=uVelGrady
          BdummyE=BuVelGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'uVelGrady(1/s)', VTK_unit,0)
          dummyE=uVelGradz
          BdummyE=BuVelGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'uVelGradz(1/s)', VTK_unit,0)
          dummyE=vVelGradx
          BdummyE=BvVelGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'vVelGradx(1/s)', VTK_unit,0)
          dummyE=vVelGrady
          BdummyE=BvVelGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'vVelGrady(1/s)', VTK_unit,0)
          dummyE=vVelGradz
          BdummyE=BvVelGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'vVelGradz(1/s)', VTK_unit,0)
          dummyE=wVelGradx
          BdummyE=BwVelGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'wVelGradx(1/s)', VTK_unit,0)
          dummyE=wVelGrady
          BdummyE=BwVelGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'wVelGrady(1/s)', VTK_unit,0)
          dummyE=wVelGradz
          BdummyE=BwVelGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'wVelGradz(1/s)', VTK_unit,0)
          if(.not.LTurbulentFlow) then
            call AllocateStrainRateTensor
            call CalculateStrainRateTensor
          endif
          dummyE=StrainRate
          BdummyE=BStrainRate
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'Strain-rate(1/s)', VTK_unit,0)
          if(.not.LTurbulentFlow) call deAllocateStrainRateTensor
          if(.not.LTurbulentFlow) then 
            call AllocateVorticityTensor
            call CalculateVorticityTensor
          endif
          dummyE=Vorticity
          BdummyE=BVorticity
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'Vorticity(1/s)', VTK_unit,0)
          if(.not.LTurbulentFlow) call deAllocateVorticityTensor
        endif
        if(LSolveLambdaELEEquation) then
          dummyE=uVelGradx
          BdummyE=BuVelGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'uVelGradx(1/s)', VTK_unit,0)
          dummyE=uVelGrady
          BdummyE=BuVelGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'uVelGrady(1/s)', VTK_unit,0)
          dummyE=uVelGradz
          BdummyE=BuVelGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'uVelGradz(1/s)', VTK_unit,0)
          dummyE=vVelGradx
          BdummyE=BvVelGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'vVelGradx(1/s)', VTK_unit,0)
          dummyE=vVelGrady
          BdummyE=BvVelGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'vVelGrady(1/s)', VTK_unit,0)
          dummyE=vVelGradz
          BdummyE=BvVelGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'vVelGradz(1/s)', VTK_unit,0)
          dummyE=wVelGradx
          BdummyE=BwVelGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'wVelGradx(1/s)', VTK_unit,0)
          dummyE=wVelGrady
          BdummyE=BwVelGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'wVelGrady(1/s)', VTK_unit,0)
          dummyE=wVelGradz
          BdummyE=BwVelGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'wVelGradz(1/s)', VTK_unit,0)
          dummyE=LambdaELEGradx
          BdummyE=BLambdaELEGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'LambdaGradx', VTK_unit,0)
          dummyE=LambdaELEGrady
          BdummyE=BLambdaELEGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'LambdaGrady', VTK_unit,0)
          dummyE=LambdaELEGradz
          BdummyE=BLambdaELEGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'LambdaGradz', VTK_unit,0)
        endif
        if(LSolveContinuity) then
          dummyE=PressGradx
          BdummyE=BPressGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'PressureGradx(Pa/m)', VTK_unit,0)
          dummyE=PressGrady
          BdummyE=BPressGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'PressureGrady(Pa/m)', VTK_unit,0)
          dummyE=PressGradz
          BdummyE=BPressGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'PressureGradz(Pa/m)', VTK_unit,0)
        endif
        if(LSolveTurbulenceKineticEnergy) then
          dummyE=TKEGradx
          BdummyE=BTKEGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TKEGradx(m/s2)', VTK_unit,0)
          dummyE=TKEGrady
          BdummyE=BTKEGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TKEGrady(m/s2)', VTK_unit,0)
          dummyE=TKEGradz
          BdummyE=BTKEGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TKEGradz(m/s2)', VTK_unit,0)
        endif
        if(LSolveTurbulenceDissipationRate) then
          dummyE=TEDGradx
          BdummyE=BTEDGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TEDGradx(m/s2)', VTK_unit,0)
          dummyE=TEDGrady
          BdummyE=BTEDGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TEDGrady(m/s2)', VTK_unit,0)
          dummyE=TEDGradz
          BdummyE=BTEDGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TEDGradz(m/s2)', VTK_unit,0)
        endif
        if(LSolveTurbulenceV2Equation) then
          dummyE=TurbulentV2Gradx
          BdummyE=BTurbulentV2Gradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'v2Gradx(m/s2)', VTK_unit,0)
          dummyE=TurbulentV2Grady
          BdummyE=BTurbulentV2Grady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'v2Grady(m/s2)', VTK_unit,0)
          dummyE=TurbulentV2Gradz
          BdummyE=BTurbulentV2Gradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'v2Gradz(m/s2)', VTK_unit,0)
        endif
        if(LSolveTurbulenceZetaEquation) then
          dummyE=TurbulentZetaGradx
          BdummyE=BTurbulentZetaGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'ZetaGradx(m/s2)', VTK_unit,0)
          dummyE=TurbulentZetaGrady
          BdummyE=BTurbulentZetaGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'ZetaGrady(m/s2)', VTK_unit,0)
          dummyE=TurbulentZetaGradz
          BdummyE=BTurbulentZetaGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'ZetaGradz(m/s2)', VTK_unit,0)
        endif
        if(LSolveTurbulencefRelaxationEquation) then
          dummyE=TfRelaxationGradx
          BdummyE=BTfRelaxationGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'fGradx(m/s2)', VTK_unit,0)
          dummyE=TfRelaxationGrady
          BdummyE=BTfRelaxationGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'fGrady(m/s2)', VTK_unit,0)
          dummyE=TfRelaxationGradz
          BdummyE=BTfRelaxationGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'fGradz(m/s2)', VTK_unit,0)
        endif
        if(LSolveTurbulenceSpecificDissipationRate) then
          dummyE=TOmegaGradx
          BdummyE=BTOmegaGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TOmegaGradx(1/(m.s))', VTK_unit,0)
          dummyE=TOmegaGrady
          BdummyE=BTOmegaGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TOmegaGrady(1/(m.s))', VTK_unit,0)
          dummyE=TOmegaGradz
          BdummyE=BTOmegaGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TOmegaGradz(1/(m.s))', VTK_unit,0)
        endif
        if(LSolveTurbulenceGammaEquation) then
          dummyE=TGammaGradx
          BdummyE=BTGammaGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TGammaGradx(1/m)', VTK_unit,0)
          dummyE=TGammaGrady
          BdummyE=BTGammaGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TGammaGrady(1/m)', VTK_unit,0)
          dummyE=TGammaGradz
          BdummyE=BTGammaGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TGammaGradz(1/m)', VTK_unit,0)
        endif
        if(LSolveTurbulenceReynoldsThetaEquation) then
          dummyE=TReThetaGradx
          BdummyE=BTReThetaGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TReThetaGradx(1/m)', VTK_unit,0)
          dummyE=TReThetaGrady
          BdummyE=BTReThetaGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TReThetaGrady(1/m)', VTK_unit,0)
          dummyE=TReThetaGradz
          BdummyE=BTReThetaGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TReThetaGradz(1/m)', VTK_unit,0)
        endif
        if(LSolveTurbulentKL) then
          dummyE=TurbulentKLGradx
          BdummyE=BTurbulentKLGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TKLGradx(m/s2)', VTK_unit,0)
          dummyE=TurbulentKLGrady
          BdummyE=BTurbulentKLGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TKLGrady(m/s2)', VTK_unit,0)
          dummyE=TurbulentKLGradz
          BdummyE=BTurbulentKLGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TKLGradz(m/s2)', VTK_unit,0)
        endif
        if(LSolveModifiedED) then
          dummyE=ModifiedEDGradx
          BdummyE=BModifiedEDGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'MEDGradx(m/s)', VTK_unit,0)
          dummyE=ModifiedEDGrady
          BdummyE=BModifiedEDGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'MEDGrady(m/s)', VTK_unit,0)
          dummyE=ModifiedEDGradz
          BdummyE=BModifiedEDGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'MEDGradz(m/s)', VTK_unit,0)
        endif
        if(LSolveEnergy) then
          dummyE=TempGradx
          BdummyE=BTempGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TempGradx(K/m)', VTK_unit,0)
          dummyE=TempGrady
          BdummyE=BTempGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TempGrady(K/m)', VTK_unit,0)
          dummyE=TempGradz
          BdummyE=BTempGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'TempGradz(K/m)', VTK_unit,0)
        endif
        if(LSolveEnergy.and.EnergyEquation.eq.'htotal') then
          dummyE=HtotalGradx
          BdummyE=BHtotalGradx
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'HTotalGradx(J/(Kg.m))', VTK_unit,0)
          dummyE=HtotalGrady
          BdummyE=BHtotalGrady
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'HTotalGrady(J/(Kg.m))', VTK_unit,0)
          dummyE=HtotalGradz
          BdummyE=BHtotalGradz
          call InterpolateToNodesFromElements
          call printScalar(dummyN,'HTotalGradz(J/(Kg.m))', VTK_unit,0)
        endif
        do irField=1,NumberOfrFieldsToSolve
          if(LSolverField(irField)) then
            part1=trim(rFieldName(irField))
            part2="Gradx"
            part3="Grady"
            part4="Gradz"
            namesize=len(part1)+5
            name1=trim(part1)//part2
            name2=trim(part1)//part3
            name3=trim(part1)//part4
            write(12,1,advance='no') name1(1:namesize)
            write(12,1,advance='no') name2(1:namesize)
            write(12,1,advance='no') name3(1:namesize)
            dummyE=rFieldGradx(:,irField) 
            BdummyE=BrFieldGradx(:,:,irField) 
            call InterpolateToNodesFromElements
            call printScalar(dummyN,name1, VTK_unit,0)
            dummyE=rFieldGrady(:,irField) 
            BdummyE=BrFieldGrady(:,:,irField) 
            call InterpolateToNodesFromElements
            call printScalar(dummyN,name2, VTK_unit,0)
            dummyE=rFieldGradz(:,irField) 
            BdummyE=BrFieldGradz(:,:,irField) 
            call InterpolateToNodesFromElements
            call printScalar(dummyN,name3, VTK_unit,0)
          endif
        enddo
        do iScalar=1,NumberOfScalarsToSolve
          if(LSolveScalar(iScalar)) then
            part1=trim(ScalarName(iScalar))
            part2="Gradx"
            part3="Grady"
            part4="Gradz"
            namesize=len(part1)+5
            name1=trim(part1)//part2
            name2=trim(part1)//part3
            name3=trim(part1)//part4
            write(12,1,advance='no') name1(1:namesize)
            write(12,1,advance='no') name2(1:namesize)
            write(12,1,advance='no') name3(1:namesize)
            dummyE=ScalarGradx(:,iScalar) 
            BdummyE=BScalarGradx(:,:,iScalar) 
            call InterpolateToNodesFromElements
            call printScalar(dummyN,name1, VTK_unit,0)
            dummyE=ScalarGrady(:,iScalar) 
            BdummyE=BScalarGrady(:,:,iScalar) 
            call InterpolateToNodesFromElements
            call printScalar(dummyN,name2, VTK_unit,0)
            dummyE=ScalarGradz(:,iScalar) 
            BdummyE=BScalarGradz(:,:,iScalar) 
            call InterpolateToNodesFromElements
            call printScalar(dummyN,name3, VTK_unit,0)
          endif
        enddo
      endif
!       
      deallocate(dummyE)
      deallocate(BdummyE)
      deallocate(dummyN)
      deallocate(SumInvd)
      deallocate(SumPhiInvd)
!    
      return
      end SUBROUTINE PrintVariablesBasedOnNodes
!
!************************************************************************************************
      SUBROUTINE InterpolateToNodesFromElements
!************************************************************************************************
!
      use Geometry1
      use Geometry3, only:NBFaces,NBFacesMax,NBFaceNodes,NIFaces,NBFaceOwner,&
                          NumberOfElementFaceNodes,GlobalFaceNumberOfNodes
      use Geometry4, only:BFaceCentroidx,BFaceCentroidy,BFaceCentroidz,xc,yc,zc
      use tecplot1
!*********************************************************************************************
      implicit none
!********************************************************************************************
      integer i,j,k,k1,NF,j1,j2
      double precision d
!********************************************************************************************
!
!--- First find for interior nodes
!
      SumInvd=0.
      SumPhiInvd=0.
!
      do i=1,NumberOfElements
        do j=1,NumbOfElementNodes(i)
!
          k=ListOfElementNodes(i,j)

          if(NodeFlag(k).eq.1) then
!
            d=dsqrt((xc(i)-x(k))**2+(yc(i)-y(k))**2+(zc(i)-z(k))**2)
            SumInvd(k)=SumInvd(k)+1./d
            SumPhiInvd(k)=SumPhiInvd(k)+dummyE(i)/d
!
          endif
!
        enddo
      enddo
!
      do i=1,NumberOfElements
        do j=1,NumbOfElementNodes(i)
!
          k=ListOfElementNodes(i,j)

          if(NodeFlag(k).eq.1) then
!
            dummyN(k)=SumPhiInvd(k)/SumInvd(k)
!
          endif
!
        enddo
      enddo
!
!--- Second find for boundary nodes
!
      j1=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
!
          j1=j1+1
          do k=1,GlobalFaceNumberOfNodes(j1)
!
            k1=NBFaceNodes(i,j,k)
            d=dsqrt((BFaceCentroidx(i,j)-x(k1))**2+&
                         (BFaceCentroidy(i,j)-y(k1))**2+&
                               (BFaceCentroidz(i,j)-z(k1))**2)
            SumInvd(k1)=SumInvd(k1)+1./d
            SumPhiInvd(k1)=SumPhiInvd(k1)+BdummyE(i,j)/d        
!               
          enddo
!
        enddo
      enddo
!
      j1=NIFaces
      do i=1,NumberOfBCSets
        do j=1,NBCDataRecords(i)
!
          j1=j1+1
          do k=1,GlobalFaceNumberOfNodes(j1)
!
            k1=NBFaceNodes(i,j,k)
            dummyN(k1)=SumPhiInvd(k1)/SumInvd(k1)
!               
          enddo
        enddo
      enddo
!
      return
      end SUBROUTINE InterpolateToNodesFromElements
!************************************************************************************************
     SUBROUTINE ReadFoamMesh
!************************************************************************************************
     use ReadpolyMesh
!************************************************************************************************
     implicit none 
!************************************************************************************************
!
     call Readit
!    
     return
     end SUBROUTINE ReadFoamMesh
!************************************************************************************************
     SUBROUTINE WriteParaviewFile
!************************************************************************************************
     use ReadpolyMesh
!************************************************************************************************
     implicit none
!************************************************************************************************
     integer :: VTK_unit
!************************************************************************************************
!
     VTK_unit=12
     call makeVTU(VTK_unit) 
!     
     return
     end SUBROUTINE WriteParaviewFile

!--------------------------------------------------------------------------------------!
!Added subroutines
!--------------------------------------------------------------------------------------!
subroutine AltitudeToHight(xc,yc,zc,nCells,PolyMeshDirectory)
    !NumberOfElements is nCells
    integer, intent(in) :: nCells
    integer :: nTPoints, IOstatus
    double precision, dimension(nCells) :: xc, yc, zc
    double precision :: zcZero
    double precision, allocatable :: xT(:), yT(:), zT(:)
    double precision :: Xsurrounding(3), Ysurrounding(3), Zsurrounding(3)
    character (len=*) :: PolyMeshDirectory
    !----------------------------------------------------------------------!
    !Reading terrain coordinates
     open(unit = 7, file = trim(PolyMeshDirectory)//"/TerrainXYZ.txt")
     nTPoints=0
     do 
       read (7,*,iostat=IOstatus) 
       if (IOstatus/=0) exit 
       nTPoints = nTPoints+1 
     end do
     
     allocate(xT(nTPoints))
     allocate(yT(nTPoints))
     allocate(zT(nTPoints))

     close(7)
     open(unit = 7, file = trim(PolyMeshDirectory)//"/TerrainXYZ.txt")
     do i=1,nTpoints
       read (7,*,iostat=IOstatus) xT(i), yT(i), zT(i)
     end do
     close(7)
     !----------------------------------------------------------------------!
     do i=1,nCells
        !Find nearest 3 points from the 3 sides
        call SurroundingPoints(xc(i),yc(i),xT,yT,zT,nTPoints,Xsurrounding,Ysurrounding,Zsurrounding)
        !Find the terrain Z coordinate for the corresponing xc, yc coordinates
        call FindTerrainZ(xc(i),yc(i),zcZero,Xsurrounding,Ysurrounding,Zsurrounding)
        !Convert Altitude to Hight
        zc(i) = zc(i)-zcZero
     end do
end subroutine AltitudeToHight

subroutine SurroundingPoints(xc,yc,Xarray,Yarray,Zarray,nPoints,Xsurrounding,Ysurrounding,Zsurrounding)
    !find Id of points that surround xc, yc point
    integer, intent(in) :: nPoints
    double precision ,intent(out):: Xsurrounding(3), Ysurrounding(3), Zsurrounding(3)
    integer :: PointsID(3)
    double precision :: xc, yc
    double precision, dimension(nPoints) :: Xarray, Yarray, Zarray, Difference2D
    logical :: mk(nPoints)

    mk = .true.

    Difference2D = sqrt((Xarray-xc)**2 &
                       +(Yarray-yc)**2)

    !Search the top right point
    mk = .false.                  
    do j=1,nPoints
       if((xc.le.Xarray(j)).and.(yc.le.Yarray(j))) mk(j) = .true.
    end do
    PointsID(1) = minloc(Difference2D, 1, mk)

    !Search the bottom right point
    mk = .false.                  
    do j=1,nPoints
       if((xc.le.Xarray(j)).and.(yc.ge.Yarray(j))) mk(j) = .true.
    end do
    PointsID(2) = minloc(Difference2D, 1, mk)

    !Search the bottom left point
    mk = .false.                  
    do j=1,nPoints
       if((xc.ge.Xarray(j)).and.(yc.ge.Yarray(j))) mk(j) = .true.
    end do
    PointsID(3) = minloc(Difference2D, 1, mk)

    do i=1,3
        Xsurrounding(i) = Xarray(PointsID(i))
        Ysurrounding(i) = Yarray(PointsID(i))
        Zsurrounding(i) = Zarray(PointsID(i))
    end do

end subroutine

subroutine FindTerrainZ(xc,yc,zcZero,Xsurrounding,Ysurrounding,Zsurrounding)
    !Calculate Terrain Z altitude value for given X Y coordinates
    double precision ,intent(In):: Xsurrounding(3), Ysurrounding(3), Zsurrounding(3)
    double precision ,intent(In):: xc,yc
    double precision ,intent(Out)::zcZero
    double precision :: P12(3), P13(3), planeVector(3)
    double precision :: xo, yo, zo, a, b, c

    !Calculating Vectors that represent the plane
    P12(1) = Xsurrounding(2)-Xsurrounding(1)
    P12(2) = Ysurrounding(2)-Ysurrounding(1)
    P12(3) = Zsurrounding(2)-Zsurrounding(1)

    P13(1) = Xsurrounding(3)-Xsurrounding(1)
    P13(2) = Ysurrounding(3)-Ysurrounding(1)
    P13(3) = Zsurrounding(3)-Zsurrounding(1)

    !Cross product of plane vectors to define the plane normal vector
    planeVector(1) = P12(2)*P13(3)-P12(3)*P13(2)
    planeVector(2) = P12(3)*P13(1)-P12(1)*P13(3)
    planeVector(3) = P12(1)*P13(2)-P12(2)*P13(1)

    !Calculate zcZero
    Xo = Xsurrounding(1)
    Yo = Ysurrounding(1)
    Zo = Zsurrounding(1)
    a = planeVector(1)
    b = planeVector(2)
    c = planeVector(3)

    zcZero = Zo-((a*(xc-Xo)+b*(yc-Yo))/c)
end subroutine
!--------------------------------------------------------------------------------------!