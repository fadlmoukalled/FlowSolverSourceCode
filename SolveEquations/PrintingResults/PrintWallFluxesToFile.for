c
C#############################################################################################
c
      Subroutine CalculateFluxes
c
C#############################################################################################
c
      use User0, only: LSolveMomentum,LSolveEnergy,ScalarName,
     *                 LanisotropicDiffusion,NumberOfScalarsToSolve,
     *                 LSolveScalar,LTurbulentFlow,name,
     *                 LcalculateReTheta,Filcf,LRough,SolutionDirectory
      use ReferenceValues1
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NBFaces,NBFaceOwner
      use Geometry4, only: BDistanceCFx,BDistanceCFy,BDistanceCFz,
     *                     BFaceAreanx,BFaceAreany,BFaceAreanz,
     *                     BFaceArea,BFaceCentroidx,BFaceCentroidy,
     *                     BFaceCentroidz,BFaceAreax,BFaceAreay,
     *                     BFaceAreaz,BgDiffp,
     *                     BFaceTxp,BFaceTyp,BFaceTzp
      use Variables1, only:uVelocity,vVelocity,wVelocity,BuVelocity,
     *                     BvVelocity,BwVelocity,BPressure,Temperature,
     *                     BTemperature,BTempGradx,BTempGrady,BTempGradz
      use Scalar1, only: Scalar,BScalar,
     *                   BScalarGradx,BScalarGrady,BScalarGradz
      use Scalar2
      use BoundaryConditions1, only: BoundaryType,wallTypeM,wallTypeE
      use PhysicalProperties1, only: BViscosity,BConductivity,
     *                               DiffusionCoefficient,
     *                               BDiffusionCoefficient,
     *                               BeDiffCoefficient,eDiffCoefficient,
     *                               TurbulentViscosity,BDensity,RGas
      use BoundaryConditions2, only: SpecificHeatFarField,
     *                               pressureFarField,MachFarField,
     *                               TemperatureFarField
      use BoundaryFluxes, only: HeatFlux,ScalarFlux
      use Turbulence1, only: yplus,ystar,uplus,ustar,uTau,KsPlus,
     *                       WallViscosity
      use Constants1, only: tiny
      use Extra, only: ReTheta
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      integer :: i,j,k,k1,iScalar,namesize
      double precision :: gDiff,dNorm,nx,ny,nz,dotproduct,
     *                    ShearStressX,ShearStressY,ShearStressZ,
     *                    ShearStress,ForceX,ForceY,ForceZ,Force,
     *                    areaX,areaY,areaZ,area,
     *                    MomentX,MomentY,MomentZ,Moment,
     *                    HeatFluxl,TotalHeatFlux,ratio,cf,cp,
     *                    uWall1,vWall1,wWall1,WallVelocity,Rex
      logical, save, dimension(:), allocatable :: LPrintBoundary
      character(len=8) :: fmt ! format descriptor
      character(len=8) :: x1
      character(len=8) :: x2
c********************************************************************************************
      interface
c********************************************************************************************
        FUNCTION TangentialVelocity(i1)
c********************************************************************************************
          integer :: i1
          double precision :: TangentialVelocity
c********************************************************************************************
        end FUNCTION TangentialVelocity
c********************************************************************************************
        FUNCTION TangentialVelocityLaminar(i1,i2,i3)
c********************************************************************************************
          integer :: i1,i2,i3
          double precision :: TangentialVelocityLaminar
c********************************************************************************************
        end FUNCTION TangentialVelocityLaminar
c********************************************************************************************
      end interface
c********************************************************************************************
      allocate(LPrintBoundary(NumberOfBCSets))
      LPrintBoundary=.true.
c
      namesize=len(trim(name))
c
      fmt = '(I5.5)' ! an integer of width 5 with zeros at the left
c
c--- Print wall shear stress
c
      if(LSolveMomentum) then
c
        Variable='velx'
        call CalculateEffectiveDiffusionCoefficient(Variable)
c
c---  Calculate ReTheta 
c
        if(LTurbulentFlow.and.LcalculateReTheta) then
          call CalculateReTheta
        endif
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'wall'.and.
     *                      wallTypeM(i).eq.'noslip') then
c
            write(11,*) '  '
            write(11,*) '    Shear force along wall number ',i
            write(11,*) '    xc           yc           zc         Tau-x 
     *      Tau-y      Tau-z        Tau       Cf'
c
            if(LPrintBoundary(i)) then
c
              write (x1,fmt) i ! converting integer to string using a 'internal file'
              Filcf=trim(name)//'Cfwall'//trim(x1)//'.dat'
c
              open (unit=17,
     *          file=trim(SolutionDirectory)//'/'//trim(Filcf))
              rewind 17
c
              write(17,*) 'Title = "cf" '
              if(LTurbulentFlow.and.LcalculateReTheta) then
                write(17,*) 
     *              'Variables=x,y,z,ReTheta,Rex,Taux,Tauy,Tauz,Tau,Cf'
              else
                write(17,*) 'Variables=x,y,z,Taux,Tauy,Tauz,Tau,Cf'
              endif
              write(17,*) 'zone T="t1", i=',NBFaces(i)
c              
              LPrintBoundary(i)=.false.
c
            endif
c
            ForceX=0.
            ForceY=0.
            Forcez=0.    
            Force=0.
c
            MomentX=0.
            MomentY=0.
            MomentZ=0.
            Moment=0.
c
            ShearStressX=0.
            ShearStressY=0.
            ShearStressZ=0.
            ShearStress=0.
c
            k1=0
            do j=1,NBFaces(i)
c
              k=NBFaceOwner(i,j)
c
              k1=k1+1
              dNorm=BDistanceCFx(i,j)*BFaceAreanx(i,j)+
     *                BDistanceCFy(i,j)*BFaceAreany(i,j)+
     *                   BDistanceCFz(i,j)*BFaceAreanz(i,j)
c
              area=BFaceArea(i,j)
              areaX=BFaceAreax(i,j)
              areaY=BFaceAreay(i,j)
              areaZ=BFaceAreaz(i,j)
c
              nx=BFaceAreanx(i,j)
              ny=BFaceAreany(i,j)
              nz=BFaceAreanz(i,j)
c
              gDiff=BeDiffCoefficient(i,j)/dNorm
              ShearStressX=gDiff*(
     *               (uVelocity(k)-BuVelocity(i,j))*(1-nx**2)-
     *                        (vVelocity(k)-BvVelocity(i,j))*nx*ny-
     *                           (wVelocity(k)-BwVelocity(i,j))*nx*nz)
              ShearStressY=gDiff*(
     *               (vVelocity(k)-BvVelocity(i,j))*(1-ny**2)-
     *                        (uVelocity(k)-BuVelocity(i,j))*nx*ny-
     *                        (wVelocity(k)-BwVelocity(i,j))*ny*nz)
              ShearStressZ=gDiff*(
     *               (wVelocity(k)-BwVelocity(i,j))*(1-nz**2)-
     *                        (uVelocity(k)-BuVelocity(i,j))*nx*nz-
     *                        (vVelocity(k)-BvVelocity(i,j))*ny*nz)
c
              ForceX=ForceX+ShearStressX*area
              ForceY=ForceY+ShearStressY*area
              Forcez=Forcez+ShearStressZ*area
c
              MomentX=MomentX+BfaceCentroidy(i,j)*ShearStressZ*area-
     *                     BfaceCentroidz(i,j)*ShearStressY*area
              MomentY=MomentY+BfaceCentroidz(i,j)*ShearStressX*area-
     *                     BfaceCentroidx(i,j)*ShearStressZ*area
              MomentZ=MomentZ+BfaceCentroidx(i,j)*ShearStressY*area-
     *                     BfaceCentroidy(i,j)*ShearStressX*area
c
              if(LTurbulentFlow) then
                  WallVelocity=TangentialVelocity(k)
              else
                  WallVelocity=TangentialVelocityLaminar(k,i,j)
              endif
c
              ShearStress=gDiff*WallVelocity
              
              
c              ShearStress=dsqrt(ShearStressX**2+
c     *                             ShearStressY**2+ShearStressZ**2)
              
              cf=ShearStress/
     *              (0.5*dmax1(Rhoinfinity*Uinfinity**2,tiny))
c
              if(LTurbulentFlow.and.LcalculateReTheta) then
c
                Rex=Rhoinfinity*Uinfinity*
     *                    BFaceCentroidx(i,j)/BViscosity(i,j)
                write(11,'(1p10e12.4)') BFaceCentroidx(i,j),
     *                  BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                   ReTheta(k1),Rex,ShearStressX,ShearStressY,
     *                      ShearStressZ,ShearStress,cf
                write(17,*) BFaceCentroidx(i,j),
     *                  BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                   ReTheta(k1),Rex,ShearStressX,ShearStressY,
     *                      ShearStressZ,ShearStress,cf
c
              else
c
                write(11,'(1p8e12.4)') BFaceCentroidx(i,j),
     *                  BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                  ShearStressX,ShearStressY,ShearStressZ,
     *                      ShearStress,cf
                write(17,*) BFaceCentroidx(i,j),
     *                  BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                  ShearStressX,ShearStressY,ShearStressZ,
     *                      ShearStress,cf
c
              endif
c
            enddo
c
            close(17)
c            
            Force=dsqrt(ForceX**2+ForceY**2+ForceZ**2)
            Moment=dsqrt(MomentX**2+MomentY**2+MomentZ**2)
            write(11,*) '  '
            write(11,*)'   Total shear force in x-direction =',ForceX
            write(11,*)'   Total shear force in y-direction =',ForceY
            write(11,*)'   Total shear force in z-direction =',ForceZ
            write(11,*)'   Total shear force =',Force
     *                           
            write(11,*) 
     *        ' Moment due to viscous forces around x-axis = ',MomentX
            write(11,*) 
     *        ' Moment due to viscous forces around y-axis = ',MomentY
            write(11,*) 
     *        ' Moment due to viscous forces around z-axis = ',MomentZ
            write(11,*)'   Moment due to viscous forces     =',Moment
            write(11,*) '  '
c
          endif
c        
        enddo
c
        LPrintBoundary=.true.
c
c--- Print normal forces due to pressure along walls
c
        do i=1,NumberOfBCSets
c
          if(BoundaryType(i).eq.'wall') then
c
            write(11,*) '  '
            write(11,*) 'Pressure force along wall number ',i
            write(11,*) '    xc         yc         zc          Pressure 
     *FPx          FPy          FPy          FP         cp'
c
            if(LPrintBoundary(i)) then
c
              write (x1,fmt) i ! converting integer to string using a 'internal file'
              Filcf=trim(name)//'Cpwall'//trim(x1)//'.dat'
c
              open (unit=17,
     *                file=trim(SolutionDirectory)//'/'//trim(Filcf))
              rewind 17
c
              write(17,*) 'Title = "cp" '
              write(17,*) 'Variables=x,y,z,Pressure,FPx,FPy,FPz,FP,Cp'
              write(17,*) 'zone T="t1", i=',NBFaces(i)
c              
              LPrintBoundary(i)=.false.
c
            endif
c
            ForceX=0.
            ForceY=0.
            ForceZ=0.
            Force=0.
c
            MomentX=0.
            MomentY=0.
            MomentZ=0.
            Moment=0.
c
            do j=1,NBFaces(i)
c
              area=BFaceArea(i,j)
              areaX=BFaceAreax(i,j)
              areaY=BFaceAreay(i,j)
              areaZ=BFaceAreaz(i,j)
c
              ForceX=ForceX+BPressure(i,j)*areaX
              ForceY=ForceY+BPressure(i,j)*areaY
              Forcez=Forcez+BPressure(i,j)*areaZ
              Force=Force+BPressure(i,j)*area
c
              MomentX=MomentX+BfaceCentroidy(i,j)*BPressure(i,j)*areaZ-
     *                     BfaceCentroidz(i,j)*BPressure(i,j)*areaY
              MomentY=MomentY+BfaceCentroidz(i,j)*BPressure(i,j)*areaX-
     *                     BfaceCentroidx(i,j)*BPressure(i,j)*areaZ
              MomentZ=MomentZ+BfaceCentroidx(i,j)*BPressure(i,j)*areaY-
     *                     BfaceCentroidy(i,j)*BPressure(i,j)*areaX
c
              cp=(BPressure(i,j)-pressureFarField(i))/
     *                  (0.5*dmax1(Rhoinfinity*Uinfinity**2,tiny))
c
              write(11,'(1p9e12.4)') BFaceCentroidx(i,j),
     *                  BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                  BPressure(i,j),BPressure(i,j)*areaX,
     *                  BPressure(i,j)*areaY,BPressure(i,j)*areaZ,
     *                  BPressure(i,j)*area,cp
              write(17,*) BFaceCentroidx(i,j),
     *                  BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                  BPressure(i,j),BPressure(i,j)*areaX,
     *                  BPressure(i,j)*areaY,BPressure(i,j)*areaZ,
     *                  BPressure(i,j)*area,cp
c
            enddo
c
            close(17)
c            
            write(11,*) '  '
            write(11,*)'   Total normal force in x-direction =',ForceX
            write(11,*)'   Total normal force in y-direction =',ForceY
            write(11,*)'   Total normal force in z-direction =',ForceZ
            write(11,*)'   Total normal force ',Force
            write(11,*)'   Total normalized normal force in x-direction
     * =',ForceX/(0.5*dmax1(Rhoinfinity*Uinfinity**2,tiny))
            write(11,*)'   Total normalized normal force in y-direction 
     * =',ForceY/(0.5*dmax1(Rhoinfinity*Uinfinity**2,tiny))
            write(11,*)'   Total normalized normal force in z-direction 
     * =',ForceZ/(0.5*dmax1(Rhoinfinity*Uinfinity**2,tiny))
            write(11,*)'   Total normalized normal force  
     * =',Force/(0.5*dmax1(Rhoinfinity*Uinfinity**2,tiny))
            write(11,*)
     *       '   Moment due to pressure forces around x-axis =',MomentX
            write(11,*)
     *       '   Moment due to pressure forces around y-axis =',MomentY
            write(11,*)
     *       '   Moment due to pressure forces around z-axis =',MomentZ
            write(11,*)
     *       '   Total Moment due to pressure forces =',Moment
            write(11,*) '  '
c
          endif
c        
        enddo
c
        LPrintBoundary=.true.
c
        if(LTurbulentFlow) then
c
          k=0
          do i=1,NumberOfBCSets
c
            if(BoundaryType(i).eq.'wall'.and.
     *                    wallTypeM(i).eq.'noslip') then
c
              write(11,*) '  '
              write(11,*) 'Yplus and viscosity along wall number ',i
c
            if(.not.LRough) then
            write(11,*) '    xc         yc         zc          Yplus    
     *Ystar        uplus        ustar        uTau       MuWall'
            else
            write(11,*) '    xc         yc         zc          Yplus    
     *Ystar        uplus        ustar        uTau       MuWall   KsPlus'
            endif
c
            if(LPrintBoundary(i)) then
c
              write (x1,fmt) i ! converting integer to string using a 'internal file'
              Filcf=trim(name)//'YPluswall'//trim(x1)//'.dat'
c
              open (unit=17,
     *            file=trim(SolutionDirectory)//'/'//trim(Filcf))
              rewind 17
c
              write(17,*) 'Title = "wallTurbulence" '
              if(LRough) then
                write(17,*) 
     *'Variables=x,y,z,Yplus,Ystar,uplus,ustar,uTau,MuWall,KsPlus,TVisc'
              else
                write(17,*) 
     *'Variables=x,y,z,Yplus,Ystar,uplus,ustar,uTau,MuWall,TVisc'
             endif
              write(17,*) 'zone T="t1", i=',NBFaces(i)
c              
              LPrintBoundary(i)=.false.
c
            endif
c
            do j=1,NBFaces(i)
c
              k=k+1
              k1=NBFaceOwner(i,j)
c              
              if(LRough) then
                write(11,'(1p11e12.4)') BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                yplus(k),ystar(k),uplus(k),ustar(k),uTau(k),
     *                KsPlus(k),WallViscosity(k),TurbulentViscosity(k1)
                write(17,*) BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                yplus(k),ystar(k),uplus(k),ustar(k),uTau(k),
     *                KsPlus(k),WallViscosity(k),TurbulentViscosity(k1)
              else
c
                write(11,'(1p10e12.4)') BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                yplus(k),ystar(k),uplus(k),ustar(k),uTau(k),
     *                WallViscosity(k),TurbulentViscosity(k1)
                write(17,*) BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),
     *                yplus(k),ystar(k),uplus(k),ustar(k),uTau(k),
     *                WallViscosity(k),TurbulentViscosity(k1)
              endif
            enddo
c
            close(17)
c            
            write(11,*) '  '
c
          endif
c        
        enddo
c
        endif
c
      endif
c
c--- Print wall heat flux
c
      LPrintBoundary=.true.
c
      if(LSolveEnergy) then
c
        if(LanisotropicDiffusion) then
 
          Variable='temp'
          call calculateSprime(Variable)
c
          do i=1,NumberOfBCSets
c
            if(BoundaryType(i).eq.'wall') then
c
              write(11,*) '  '
              write(11,*) '    Heat flux along wall number ',i
              write(11,*) '     xc          yc          zc
     *      Heat Flux  '
c
              if(LPrintBoundary(i)) then
c
                write (x1,fmt) i ! converting integer to string using a 'internal file'
                Filcf=trim(name)//'HeatFluxwall'//trim(x1)//'.dat'
c
                open (unit=17,
     *              file=trim(SolutionDirectory)//'/'//trim(Filcf))
                rewind 17
c
                write(17,*) 'Title = "HeatFlux" '
                write(17,*) 'Variables=x,y,z,HeatFlux'
                write(17,*) 'zone T="t1", i=',NBFaces(i)
c              
                LPrintBoundary(i)=.false.
c
              endif
c
              TotalHeatFlux=0.
c
              do j=1,NBFaces(i)
c
                k=NBFaceOwner(i,j)
c
                if(wallTypeE(i).eq.'vonneumann') then
c              
                  HeatFluxl=HeatFlux(i,j)              
                  TotalHeatFlux=TotalHeatFlux+HeatFluxl*BFaceArea(i,j)
c              
                else
c
                  gDiff=BgDiffp(i,j)
                  HeatFluxl=gDiff*(BTemperature(i,j)-Temperature(k))+
     *                             BTempGradx(i,j)*BFaceTxp(i,j)+
     *                                BTempGrady(i,j)*BFaceTyp(i,j)+
     *                                  BTempGradz(i,j)*BFaceTzp(i,j)
                  TotalHeatFlux=TotalHeatFlux+HeatFluxl
                  HeatFluxl=HeatFluxl/BFaceArea(i,j)
c
                endif
c
                write(11,'(1p4e12.4)') BFaceCentroidx(i,j),
     *              BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
                write(17,*) BFaceCentroidx(i,j),
     *              BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
              enddo
c
              close(17)
c            
              write(11,*) '  '
              write(11,*)'   Total Heat Flux =',TotalHeatFlux
              write(11,*) '  '
c
            endif
c
          enddo
c
        else
c
          Variable='temp'
          call CalculateEffectiveDiffusionCoefficient(Variable)
c
          do i=1,NumberOfBCSets
c
            if(BoundaryType(i).eq.'wall') then
c
              write(11,*) '  '
              write(11,*) '    Heat flux along wall number ',i
              write(11,*) '     xc          yc          zc
     *      Heat Flux  '
c
              if(LPrintBoundary(i)) then
c
                write (x1,fmt) i ! converting integer to string using a 'internal file'
                Filcf=trim(name)//'HeatFluxwall'//trim(x1)//'.dat'
c
                open (unit=17,
     *              file=trim(SolutionDirectory)//'/'//trim(Filcf))
                rewind 17
c
                write(17,*) 'Title = "HeatFlux" '
                write(17,*) 'Variables=x,y,z,HeatFlux'
                write(17,*) 'zone T="t1", i=',NBFaces(i)
c              
                LPrintBoundary(i)=.false.
c
              endif
c
              TotalHeatFlux=0.
c
              do j=1,NBFaces(i)
c
                k=NBFaceOwner(i,j)
c
                if(wallTypeE(i).eq.'vonneumann') then
c              
                  HeatFluxl=HeatFlux(i,j)              
c              
                else
c
                  dNorm=BDistanceCFx(i,j)*BFaceAreanx(i,j)+
     *                       BDistanceCFy(i,j)*BFaceAreany(i,j)+
     *                          BDistanceCFz(i,j)*BFaceAreanz(i,j)
                  gDiff=BeDiffCoefficient(i,j)/dNorm
                  HeatFluxl=gDiff*(BTemperature(i,j)-Temperature(k))
c
                endif
c 
                TotalHeatFlux=TotalHeatFlux+HeatFluxl*BFaceArea(i,j)
c
                write(11,'(1p4e12.4)') BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
                write(17,*) BFaceCentroidx(i,j),
     *              BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
              enddo
c
              close(17)
c            
              write(11,*) '  '
              write(11,*)'   Total Heat Flux =',TotalHeatFlux
              write(11,*) '  '
c
            endif
c
          enddo
c
        endif
c
      endif
c
      if(NumberOfScalarsToSolve.gt.0) then
c
        do iScalar=1,NumberOfScalarsToSolve
c
          if(LSolveScalar(iScalar)) then
c
            LPrintBoundary=.true.
c
            if(LanisotropicDiffusion) then
 
              Variable=ScalarName(iScalar)
              iScalarVariable=iScalar
              call calculateSprime(Variable)
c
              do i=1,NumberOfBCSets
c
                if(BoundaryType(i).eq.'wall') then
c
                  write(11,*) '  '
                  write(11,*) '    Scalar ',iScalar,
     *                               ' flux along wall number ',i
                  write(11,*) '     xc          yc          zc
     *     Scalar Flux  '
c
                  if(LPrintBoundary(i)) then
c
                    write (x1,fmt) iscalar ! converting integer to string using a 'internal file'
                    write (x2,fmt) i ! converting integer to string using a 'internal file'
                    Filcf=trim(name)//'Scalar'//trim(x1)//
     *                                  'Fluxwall'//trim(x2)//'.dat'
c
                    open (unit=17,
     *                 file=trim(SolutionDirectory)//'/'//trim(Filcf))
                    rewind 17
c
                    write(17,*) 'Title = "HeatFlux" '
                    write(17,*) 'Variables=x,y,z,HeatFlux'
                    write(17,*) 'zone T="t1", i=',NBFaces(i)
c              
                    LPrintBoundary(i)=.false.
c
                  endif
c
                  TotalHeatFlux=0.
c
                  do j=1,NBFaces(i)
c
                    k=NBFaceOwner(i,j)
c
                    if(wallTypeE(i).eq.'vonneumannscalar') then
c              
                      HeatFluxl=ScalarFlux(i,j,iScalar)              
                      TotalHeatFlux=TotalHeatFlux+
     *                               HeatFluxl*BFaceArea(i,j)
c              
                    else
c
                      gDiff=BgDiffp(i,j)
                      HeatFluxl=
     *                 gDiff*(BScalar(i,j,iScalar)-Scalar(k,iScalar))+
     *                      BScalarGradx(i,j,iScalar)*BFaceTxp(i,j)+
     *                        BScalarGrady(i,j,iScalar)*BFaceTyp(i,j)+
     *                          BScalarGradz(i,j,iScalar)*BFaceTzp(i,j)
                      TotalHeatFlux=TotalHeatFlux+HeatFluxl
                      HeatFluxl=HeatFluxl/BFaceArea(i,j)
c
                    endif
c
                    write(11,'(1p4e12.4)') BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
                    write(17,*) BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
                  enddo
c
                  close(17)
c            
                  write(11,*) '  '
                  write(11,*)'   Total Scalar Flux =',TotalHeatFlux
                  write(11,*) '  '
c
                endif
c
              enddo
c
            else
c
              Variable=ScalarName(iScalar)
              call CalculateEffectiveDiffusionCoefficient(Variable)
              do i=1,NumberOfBCSets
c
                if(BoundaryType(i).eq.'wall') then
c
                  write(11,*) '  '
                  write(11,*) '    Heat flux along wall number ',i
                  write(11,*) '     xc          yc          zc
     *      Scalar Flux  '
c
                  if(LPrintBoundary(i)) then
c
                    write (x1,fmt) iscalar ! converting integer to string using a 'internal file'
                    write (x2,fmt) i ! converting integer to string using a 'internal file'
                    Filcf=trim(name)//'Scalar'//trim(x1)//
     *                             'Fluxwall'//trim(x2)//'.dat'
c
                    open (unit=17,
     *                file=trim(SolutionDirectory)//'/'//trim(Filcf))
                    rewind 17
c
                    write(17,*) 'Title = "HeatFlux" '
                    write(17,*) 'Variables=x,y,z,HeatFlux'
                    write(17,*) 'zone T="t1", i=',NBFaces(i)
c              
                    LPrintBoundary(i)=.false.
c
                  endif
c
                  TotalHeatFlux=0.
c
                  do j=1,NBFaces(i)
c
                    k=NBFaceOwner(i,j)
c
                    if(wallTypeE(i).eq.'vonneumannscalar') then
c              
                      HeatFluxl=ScalarFlux(i,j,iScalar)              
c              
                    else
c
                      dNorm=BDistanceCFx(i,j)*BFaceAreanx(i,j)+
     *                       BDistanceCFy(i,j)*BFaceAreany(i,j)+
     *                          BDistanceCFz(i,j)*BFaceAreanz(i,j)
                      gDiff=BeDiffCoefficient(i,j)/dNorm
                      HeatFluxl=gDiff*
     *                       (BScalar(i,j,iScalar)-Scalar(k,iScalar))
c
                    endif
c 
                    TotalHeatFlux=TotalHeatFlux+HeatFluxl*BFaceArea(i,j)
c
                    write(11,'(1p3e12.4)') BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
                    write(17,*) BFaceCentroidx(i,j),
     *                BFaceCentroidy(i,j),BFaceCentroidz(i,j),HeatFluxl
                  enddo
c
                  close(17)
c            
                  write(11,*) '  '
                  write(11,*)'   Total Scalar Flux =',TotalHeatFlux
                  write(11,*) '  '
c
                endif
c
              enddo
c
            endif
c
          endif
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
      SUBROUTINE CalculateReTheta
c
C#############################################################################################
c
      use User0, only: LcalculateReTheta,ConstantViscosity
      use ReferenceValues1
      use WallDistance1, only: WallDistance,iTau,jTau
      use Geometry1, only: NumberOfElements
      use BoundaryConditionsTurbulence2, only: IwallTurbulence,
     *                       IWallTurbulenceOwner,
     *                       IWallTurbulenceNumberOfBCSets,
     *                       IWallTurbulenceNBFaces
      use PhysicalProperties1, only: Density
      use Extra, only:ReTheta
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      integer :: i,i1,i2,i3,i4,j,j1,j2,k,b
      double precision :: term1,term2,term3,term4,a,WallVelocity
      double precision, save, dimension(:), allocatable :: dy
      integer, save, dimension(:), allocatable :: yindex
c********************************************************************************************
      interface
c********************************************************************************************
        FUNCTION TangentialVelocity(i1)
c********************************************************************************************
          integer :: i1
          double precision :: TangentialVelocity
c********************************************************************************************
        end FUNCTION TangentialVelocity
c********************************************************************************************
      end interface
c********************************************************************************************
c
      allocate(ReTheta(IwallTurbulence))
      allocate(dy(NumberOfElements))
      allocate(yindex(NumberOfElements))
c
      do i=1,IwallTurbulence
c
        i1=IWallTurbulenceOwner(i)
        i2=IWallTurbulenceNumberOfBCSets(i)
        i3=IWallTurbulenceNBFaces(i)
c
        k=0
        do j=1,NumberOfElements
c
          WallVelocity=TangentialVelocity(j)
          term1=WallVelocity/Uinfinity           
          if(iTau(j).eq.i2.and.jTau(j).eq.i3.and.term1.le.0.99) then
c                      
            k=k+1       
            dy(k)=WallDistance(j)
            yindex(k)=j
c
          endif
        enddo
c
c---- Sort by distance form wall
c
        do j=2,k                    !Pick out each element in turn.
          a=dy(j)
          b=yindex(j)
          do i4=j-1,1,-1             !Look for the place to insert it.
            if (dy(i4) <= a) exit
            dy(i4+1)=dy(i4)
            yindex(i4+1)=yindex(i4)
          enddo
          dy(i4+1)=a                !Insert it.
          yindex(i4+1)=b
        enddo          
c
c---- Calculate ReTheta
c
        ReTheta(i)=0.
        do j=1,k         
c
          WallVelocity=TangentialVelocity(yindex(j))
          term1=Density(yindex(j))/Rhoinfinity          
          term2=WallVelocity/Uinfinity           
          term3=1.-term2
          if(j.eq.1) then
            term4=dy(j)
          else
            term4=dy(j)-dy(j-1)
          endif
c
          ReTheta(i)=ReTheta(i)+term1*term2*term3*term4           
c
        enddo
c
        ReTheta(i)=Rhoinfinity*Uinfinity*
     *                      ReTheta(i)/ConstantViscosity
c
      enddo
c
      deallocate(dy)
      deallocate(yindex)
c
      return
      end