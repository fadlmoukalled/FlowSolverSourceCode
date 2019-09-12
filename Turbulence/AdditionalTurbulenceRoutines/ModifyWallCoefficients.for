c
c#############################################################################################
c
      SUBROUTINE ModifyCoefficientsAlongWalls
c
C#############################################################################################
c
      use User0, only: TurbulenceModel,WallTreatment
      use Geometry1, only: NumberOfElements
      use Geometry3, only: NumberofElementNeighbors
      use Variables1, only: TurbulentED,BTurbulentED,
     *                      TurbulentOmega,BTurbulentOmega,
     *                      ModifiedED,BModifiedED
      use Variables2, only: ac,anb,bc
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                                         IWallTurbulenceOwner,
     *                                  IWallTurbulenceNumberOfBCSets,
     *                                  IWallTurbulenceNBFaces
      use Turbulence1, only: ModelNumber
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,j
c********************************************************************************************
c
c-----------------------------------------------------------------------------
      ModifyCoefficientsWall: select case (ModelNumber)
c-----------------------------------------------------------------------------
        case(1,10) ModifyCoefficientsWall           !kepsilon+realizable
c-----------------------------------------------------------------------------
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            ac(i1)=1.
            bc(i1)=0.
c
            do j=1,NumberofElementNeighbors(i1)
c
              anb(i1,j)=0.
c
            enddo
c
            BTurbulentED(i2,i3)=TurbulentED(i1)
c
          enddo
c
c-----------------------------------------------------------------------------
        case(2:9) ModifyCoefficientsWall           !kepsilonchien+kepsilonsharma+kepsilonchc+kepsilonkasagi
                                                   !kepsilontagawa+kepsilonhishida+kelambremhorst
                                                   !kelambremhorstm
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(11) ModifyCoefficientsWall           !komega
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              ac(i1)=1.
              bc(i1)=0.
c
              do j=1,NumberofElementNeighbors(i1)
c
                anb(i1,j)=0.
c
              enddo
c
              BTurbulentOmega(i2,i3)=TurbulentOmega(i1)
c
            enddo
c
          else
c
            return
c
          endif
c
c-----------------------------------------------------------------------------
        case(12) ModifyCoefficientsWall           !komegaepsilon
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              ac(i1)=1.
              bc(i1)=0.
c
              do j=1,NumberofElementNeighbors(i1)
c
                anb(i1,j)=0.
c
              enddo
c
              BTurbulentOmega(i2,i3)=TurbulentOmega(i1)
c
            enddo
c
          else
c
            return
c
          endif
c
c-----------------------------------------------------------------------------
        case(13) ModifyCoefficientsWall           !komegabsl
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              ac(i1)=1.
              bc(i1)=0.
c
              do j=1,NumberofElementNeighbors(i1)
c
                anb(i1,j)=0.
c
              enddo
c
              BTurbulentOmega(i2,i3)=TurbulentOmega(i1)
c
            enddo
c
          else
c
            return
c
          endif
c
c-----------------------------------------------------------------------------
        case(14) ModifyCoefficientsWall           !komegasst
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              ac(i1)=1.
              bc(i1)=0.
c
              do j=1,NumberofElementNeighbors(i1)
c
                anb(i1,j)=0.
c
              enddo
c
              BTurbulentOmega(i2,i3)=TurbulentOmega(i1)
c
            enddo
c
          else
c
            return
c
          endif
c
c-----------------------------------------------------------------------------
        case(15) ModifyCoefficientsWall           !sstgamaretheta
c-----------------------------------------------------------------------------
c
         return
c
c-----------------------------------------------------------------------------
        case(16) ModifyCoefficientsWall           !komega2006
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'wallfunctions') then
c
            do i=1,IwallTurbulence
c
              i1=IWallTurbulenceOwner(i)
              i2=IWallTurbulenceNumberOfBCSets(i)
              i3=IWallTurbulenceNBFaces(i)
c
              ac(i1)=1.
              bc(i1)=0.
c
              do j=1,NumberofElementNeighbors(i1)
c
                anb(i1,j)=0.
c
              enddo
c
              BTurbulentOmega(i2,i3)=TurbulentOmega(i1)
c
            enddo
c
          else
c
            return
c
          endif
c
c-----------------------------------------------------------------------------
        case(17) ModifyCoefficientsWall           !komega2006lrn
c-----------------------------------------------------------------------------
c
         return
c
c-----------------------------------------------------------------------------
        case(18) ModifyCoefficientsWall           !kklmodel
c-----------------------------------------------------------------------------
c
         return
c
c-----------------------------------------------------------------------------
        case(19) ModifyCoefficientsWall           !spalartallmaras
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'wallfunctions') then
c
            return

!            do i=1,IwallTurbulence
!c
!              i1=IWallTurbulenceOwner(i)
!              i2=IWallTurbulenceNumberOfBCSets(i)
!              i3=IWallTurbulenceNBFaces(i)
!c
!              ac(i1)=1.
!              bc(i1)=0.
!c
!              do j=1,NumberofElementNeighbors(i1)
!c
!                anb(i1,j)=0.
!c
!              enddo
!c
!              ModifiedED(i1)=BModifiedED(i2,i3)
!c
!            enddo
c
          else
c
            return
c
          endif
c
c-----------------------------------------------------------------------------
        case(20) ModifyCoefficientsWall           !wrayagarwal
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'wallfunctions') then
c
            return

!            do i=1,IwallTurbulence
!c
!              i1=IWallTurbulenceOwner(i)
!              i2=IWallTurbulenceNumberOfBCSets(i)
!              i3=IWallTurbulenceNBFaces(i)
!c
!              ac(i1)=1.
!              bc(i1)=0.
!c
!              do j=1,NumberofElementNeighbors(i1)
!c
!                anb(i1,j)=0.
!c
!              enddo
!c
!              ModifiedED(i1)=BModifiedED(i2,i3)
!c
!            enddo
c
          else
c
            return
c
          endif
c
c-----------------------------------------------------------------------------
        case(22) ModifyCoefficientsWall           !kepsilonrt
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(23) ModifyCoefficientsWall           !sstgama
c-----------------------------------------------------------------------------
c
          return
c
c-----------------------------------------------------------------------------
        case(24) ModifyCoefficientsWall           !nut92
c-----------------------------------------------------------------------------
c
          return
c-----------------------------------------------------------------------------
        case(25) ModifyCoefficientsWall           !kepsilonrng
c-----------------------------------------------------------------------------
c
          if(WallTreatment.eq.'lowreynoldsnumber') return
c
          do i=1,IwallTurbulence
c
            i1=IWallTurbulenceOwner(i)
            i2=IWallTurbulenceNumberOfBCSets(i)
            i3=IWallTurbulenceNBFaces(i)
c
            ac(i1)=1.
            bc(i1)=0.
c
            do j=1,NumberofElementNeighbors(i1)
c
              anb(i1,j)=0.
c
            enddo
c
            BTurbulentED(i2,i3)=TurbulentED(i1)
c
          enddo
c
c-----------------------------------------------------------------------------------------------      
      end select ModifyCoefficientsWall 
c-----------------------------------------------------------------------------------------------      
c
      return
      end