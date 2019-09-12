c
C#############################################################################################
c
      SUBROUTINE PrintConvergenceHistory(nIter,nIterTotal,source)
c
C#############################################################################################
      use User0, only: IMonitor,NumberOfScalarsToSolve,LsolveContinuity,
     *                 LsolveMomentum,LsolveEnergy,LsolveScalar,
     *                 ScalarName,NstopType,LsolveModifiedED,
     *                 LSolveTurbulenceKineticEnergy,
     *                 LSolveTurbulenceDissipationRate,
     *                 LSolveTurbulenceSpecificDissipationRate,
     *                 LSolveTurbulentKL,EnergyEquation,
     *                 LSolveTurbulenceGammaEquation,
     *                 LSolveTurbulenceReynoldsThetaEquation,
     *                 LSolveTurbulencefRelaxationEquation,
     *                 LSolveTurbulenceV2Equation,
     *                 LSolveTurbulenceZetaEquation,
     *                 NumberOfrFieldsToSolve,LsolverField,rFieldName,
     *                 LUnsteady,time,LSolveLambdaELEEquation
      use Transient1, only: ndt
      use Residuals1
      use Variables1, only: uVelocity,vVelocity,wVelocity,Pressure,
     *                      Temperature,TurbulentKE,TurbulentED,
     *                      TurbulentOmega,ModifiedED,TurbulentKL,
     *                      Htotal,TGamma,TReTheta,TfRelaxation,
     *                      TurbulentV2,TurbulentZeta,LambdaELE
      use Scalar1, only: Scalar
      use VolumeOfFluid1, only: rField
      use Geometry4, only: xc,yc,zc
      use PhysicalProperties1, only: SpecificHeat,ReferenceTemperature
c
c********************************************************************************************
c
      implicit none
c********************************************************************************************
      double precision source,Tref
      Integer nIter,iScalar,irField,nIterTotal
c
c********************************************************************************************
c
c--- Start with headings
c
      Tref=ReferenceTemperature
c
      write(*,10) IMonitor  
  10  format(2x,'*--Residuals and Field Values at Monitoring Location
     *(',1X,I6')--*')
      write(*,11) xc(IMonitor),yc(IMonitor),zc(IMonitor)
  11  format(2x,'xMonitor=',1P1E10.3,2x,'yMonitor=',1P1E10.3,2x,
     *           'zMonitor=',1P1E10.3)
      if(LUnsteady) then
        write(*,9) time,ndt,nIterTotal
      endif
   9  format(2x,'  time= ',1P1E10.3,2x,' timestep=  ',I7,
     *2x,' Total # of iterations= ',I7)
      write(*,20,advance='no')
  20  format(15X,'  Iter')   
      if(LsolveContinuity) write(*,30,advance='no')
  30  format(2X,'   Mass   ')
      if(LsolveMomentum) write(*,40,advance='no')
  40  format(2X,'   xMom   ',2X,'   yMom   ',2X,'   zMom   ')
      if(LSolveTurbulenceKineticEnergy) write(*,42,advance='no') 
  42  format(2X,'   TKE    ')
      if(LSolveTurbulenceDissipationRate) write(*,43,advance='no') 
  43  format(2X,'   TED    ')
      if(LSolveTurbulenceV2Equation) write(*,61,advance='no') 
  61  format(2X,'   v2     ')
      if(LSolveTurbulenceZetaEquation) write(*,62,advance='no') 
  62  format(2X,'   zeta   ')
      if(LSolveTurbulencefRelaxationEquation) write(*,63,advance='no') 
  63  format(2X,'   f      ')
      if(LSolveTurbulenceSpecificDissipationRate)
     *                                    write(*,44,advance='no') 
  44  format(2X,'   Omega  ')
      if(LSolveTurbulenceGammaEquation) write(*,52,advance='no') 
  52  format(2X,'   Gamma   ')
      if(LSolveTurbulenceReynoldsThetaEquation) write(*,53,advance='no')
  53  format(2X,'  ReTheta  ')
      if(LSolveModifiedED) write(*,45,advance='no') 
  45  format(2X,'ModifiedED')
      if(LSolveTurbulentKL) write(*,46,advance='no') 
  46  format(2X,'   TKL    ')
      if(LsolveEnergy) write(*,50,advance='no')
  50  format(2X,'  Energy  ')
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) then
          write(*,51,advance='no') rFieldName(irField)
        endif
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) then
          write(*,51,advance='no') ScalarName(iScalar)
        endif
      enddo
  51  format(2x,A10)
      if(LSolveLambdaELEEquation) then
        write(*,71,advance='no')
      endif
  71  format(2x,'  Lambda  ')
      write(*,60)
  60  format(1x)
c
c---  Write values at the Monitoring location
c
      write(*,70,advance='no') 
  70  format(1x,'Values at MP')
      write(*,80,advance='no') nIter
  80  format(2x,I6) 
      if(LsolveContinuity) write(*,90,advance='no') Pressure(IMonitor)
  90  format(2x,1P1E10.3) 
      if(LsolveMomentum) write(*,100,advance='no') uVelocity(IMonitor),
     *       vVelocity(IMonitor),wVelocity(IMonitor)
 100  format(2x,1P1E10.3,2x,1P1E10.3,2x,1P1E10.3) 
      if(LSolveTurbulenceKineticEnergy)  
     *               write(*,110,advance='no') TurbulentKE(IMonitor)
      if(LSolveTurbulenceDissipationRate)  
     *               write(*,110,advance='no') TurbulentED(IMonitor)
      if(LSolveTurbulenceV2Equation)  
     *               write(*,110,advance='no') TurbulentV2(IMonitor)
      if(LSolveTurbulenceZetaEquation)  
     *               write(*,110,advance='no') TurbulentZeta(IMonitor)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(*,110,advance='no') TfRelaxation(IMonitor)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(*,110,advance='no') TurbulentOmega(IMonitor)
      if(LSolveTurbulenceGammaEquation)  
     *               write(*,110,advance='no') TGamma(IMonitor)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(*,110,advance='no') TReTheta(IMonitor)
      if(LsolveModifiedED) 
     *        write(*,110,advance='no') ModifiedED(IMonitor)
      if(LsolveTurbulentKL) 
     *        write(*,110,advance='no') TurbulentKL(IMonitor)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *        write(*,110,advance='no') Temperature(IMonitor)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *          write(*,110,advance='no') Htotal(IMonitor)
 110  format(2x,1P1E10.3) 
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(*,111,advance='no') 
     *                                        rField(IMonitor,irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(*,111,advance='no') 
     *                                        Scalar(IMonitor,iScalar)
      enddo
 111  format(2x,1P1E10.3) 
      if(LSolveLambdaELEEquation) then
        write(*,110,advance='no') LambdaELE(IMonitor)
      endif
c      
      write(*,60)
c
c---  Write absolute residuals
c
      write(*,120,advance='no') 
 120  format(1x,'   absRes   ')
      write(*,130,advance='no') nIter
 130  format(2x,I6) 
      if(LsolveContinuity) write(*,140,advance='no') ResorAbs(5)
 140  format(2x,1P1E10.3) 
      if(LsolveMomentum) write(*,150,advance='no') ResorAbs(1),
     *       ResorAbs(2),ResorAbs(3)
 150  format(2x,1P1E10.3,2x,1P1E10.3,2x,1P1E10.3) 
      if(LSolveTurbulenceKineticEnergy)  
     *               write(*,160,advance='no') ResorAbs(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(*,160,advance='no') ResorAbs(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(*,160,advance='no') ResorAbs(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(*,160,advance='no') ResorAbs(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(*,160,advance='no') ResorAbs(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(*,160,advance='no') ResorAbs(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(*,160,advance='no') ResorAbs(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(*,160,advance='no') ResorAbs(14)
      if(LSolveModifiedED) write(*,160,advance='no') ResorAbs(11)
      if(LSolveTurbulentKL) write(*,160,advance='no') ResorAbs(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                  write(*,160,advance='no') ResorAbs(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                  write(*,160,advance='no') ResorAbs(7)
 160  format(2x,1P1E10.3) 
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(*,161,advance='no') 
     *                                            ResorAbs(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(*,161,advance='no') 
     *                 ResorAbs(18+NumberOfrFieldsToSolve+iScalar)
      enddo
 161  format(2x,1P1E10.3) 
      if(LSolveLambdaELEEquation)
     *                  write(*,160,advance='no') ResorAbs(18)
      write(*,60)
c
c---  Write maximum residuals
c
      write(*,170,advance='no') 
 170  format(1x,'   maxRes   ')
      write(*,180,advance='no') nIter
 180  format(2x,I6) 
      if(LsolveContinuity) write(*,190,advance='no') ResorMax(5)
 190  format(2x,1P1E10.3) 
      if(LsolveMomentum) write(*,200,advance='no') ResorMax(1),
     *       ResorMax(2),ResorMax(3)
 200  format(2x,1P1E10.3,2x,1P1E10.3,2x,1P1E10.3) 
      if(LSolveTurbulenceKineticEnergy)  
     *               write(*,210,advance='no') ResorMax(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(*,210,advance='no') ResorMax(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(*,210,advance='no') ResorMax(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(*,210,advance='no') ResorMax(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(*,210,advance='no') ResorMax(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(*,210,advance='no') ResorMax(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(*,210,advance='no') ResorMax(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(*,210,advance='no') ResorMax(14)
      if(LSolveModifiedED) write(*,210,advance='no') ResorMax(11)
      if(LSolveTurbulentKL) write(*,210,advance='no') ResorMax(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                       write(*,210,advance='no') ResorMax(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                       write(*,210,advance='no') ResorMax(7)
 210  format(2x,1P1E10.3) 
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(*,211,advance='no') 
     *                                            ResorMax(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(*,211,advance='no') 
     *                     ResorMax(18+NumberOfrFieldsToSolve+iScalar)
      enddo
 211  format(2x,1P1E10.3) 
      if(LSolveLambdaELEEquation)
     *                       write(*,210,advance='no') ResorMax(18)
      write(*,60)
c
c---  Write rms residuals
c
      write(*,220,advance='no') 
 220  format(1x,'   rmsRes   ')
      write(*,230,advance='no') nIter
 230  format(2x,I6) 
      if(LsolveContinuity) write(*,240,advance='no') ResorRMS(5)
 240  format(2x,1P1E10.3) 
      if(LsolveMomentum) write(*,250,advance='no') ResorRMS(1),
     *       ResorRMS(2),ResorRMS(3)
 250  format(2x,1P1E10.3,2x,1P1E10.3,2x,1P1E10.3) 
      if(LSolveTurbulenceKineticEnergy)  
     *               write(*,260,advance='no') ResorRMS(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(*,260,advance='no') ResorRMS(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(*,260,advance='no') ResorRMS(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(*,260,advance='no') ResorRMS(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(*,260,advance='no') ResorRMS(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(*,260,advance='no') ResorRMS(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(*,260,advance='no') ResorRMS(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(*,260,advance='no') ResorRMS(14)
      if(LSolveModifiedED)  
     *               write(*,260,advance='no') ResorRMS(11)
      if(LSolveTurbulentKL)  
     *               write(*,260,advance='no') ResorRMS(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                       write(*,260,advance='no') ResorRMS(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                       write(*,260,advance='no') ResorRMS(7)
 260  format(2x,1P1E10.3) 
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(*,261,advance='no') 
     *                                          ResorRMS(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(*,261,advance='no') 
     *                   ResorRMS(18+NumberOfrFieldsToSolve+iScalar)
      enddo
 261  format(2x,1P1E10.3) 
      if(LSolveLambdaELEEquation)
     *                       write(*,260,advance='no') ResorRMS(18)
      write(*,60)
c
c---  Write scaled residuals
c
      write(*,270,advance='no') 
 270  format(1x,'  scaledRes ')
      write(*,280,advance='no') nIter
 280  format(2x,I6) 
      if(LsolveContinuity) write(*,290,advance='no') ResorScaled(5)
 290  format(2x,1P1E10.3) 
      if(LsolveMomentum) write(*,300,advance='no') ResorScaled(1),
     *       ResorScaled(2),ResorScaled(3)
 300  format(2x,1P1E10.3,2x,1P1E10.3,2x,1P1E10.3) 
      if(LSolveTurbulenceKineticEnergy)  
     *               write(*,310,advance='no') ResorScaled(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(*,310,advance='no') ResorScaled(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(*,310,advance='no') ResorScaled(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(*,310,advance='no') ResorScaled(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(*,310,advance='no') ResorScaled(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(*,310,advance='no') ResorScaled(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(*,310,advance='no') ResorScaled(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(*,310,advance='no') ResorScaled(14)
      if(LSolveModifiedED)  
     *               write(*,310,advance='no') ResorScaled(11)
      if(LSolveTurbulentKL)  
     *               write(*,310,advance='no') ResorScaled(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                write(*,310,advance='no') ResorScaled(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                write(*,310,advance='no') ResorScaled(7)
 310  format(2x,1P1E10.3) 
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(*,311,advance='no') 
     *                                          ResorScaled(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(*,311,advance='no') 
     *                   ResorScaled(18+NumberOfrFieldsToSolve+iScalar)
      enddo
 311  format(2x,1P1E10.3) 
      if(LSolveLambdaELEEquation)
     *                write(*,310,advance='no') ResorScaled(18)
      write(*,320)
 320  format(1x,/)
c
c---  Write convergence history to file; Start with headings 
c
      if(nIter.eq.1) then
c
        write(13,10) IMonitor    
        write(13,11) xc(IMonitor),yc(IMonitor),zc(IMonitor)
        if(LUnsteady) then
          write(13,9) time,ndt,nIterTotal
        endif
        write(13,20,advance='no')
        if(LsolveContinuity) write(13,30,advance='no')
        if(LsolveMomentum) write(13,40,advance='no')
        if(LSolveTurbulenceKineticEnergy) write(13,42,advance='no') 
        if(LSolveTurbulenceDissipationRate) write(13,43,advance='no') 
        if(LSolveTurbulenceV2Equation) write(13,61,advance='no') 
        if(LSolveTurbulenceZetaEquation) write(13,62,advance='no') 
        if(LSolveTurbulencefRelaxationEquation)write(13,63,advance='no')
        if(LSolveTurbulenceSpecificDissipationRate)
     *                                    write(13,44,advance='no') 
        if(LSolveTurbulenceGammaEquation) write(13,52,advance='no') 
        if(LSolveTurbulenceReynoldsThetaEquation) 
     *                                    write(13,53,advance='no')
        if(LSolveModifiedED) write(13,45,advance='no') 
        if(LSolveTurbulentKL) write(13,46,advance='no') 
        if(LsolveEnergy) write(13,50,advance='no')
        do irField=1,NumberOfrFieldsToSolve
          if(LsolverField(irField)) then
            write(13,51,advance='no') rFieldName(irField)
          endif
        enddo
        do iScalar=1,NumberOfScalarsToSolve
          if(LsolveScalar(iScalar)) then
            write(13,51,advance='no') ScalarName(iScalar)
          endif
        enddo
        if(LSolveLambdaELEEquation) then
          write(13,71,advance='no')
          write(13,40,advance='no')
        endif
c       
        write(13,60)
c
      endif
c
c---  Write values at the Monitoring location
c
      write(13,70,advance='no') 
      write(13,80,advance='no') nIter
      if(LsolveContinuity) write(13,90,advance='no') Pressure(IMonitor)
      if(LsolveMomentum) write(13,100,advance='no') uVelocity(IMonitor),
     *       vVelocity(IMonitor),wVelocity(IMonitor)
      if(LSolveTurbulenceKineticEnergy)  
     *               write(13,110,advance='no') TurbulentKE(IMonitor)
      if(LSolveTurbulenceDissipationRate)  
     *               write(13,110,advance='no') TurbulentED(IMonitor)
      if(LSolveTurbulenceV2Equation)  
     *               write(13,110,advance='no') TurbulentV2(IMonitor)
      if(LSolveTurbulenceZetaEquation)  
     *               write(13,110,advance='no') TurbulentZeta(IMonitor)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(13,110,advance='no') TfRelaxation(IMonitor)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(13,110,advance='no') TurbulentOmega(IMonitor)
      if(LSolveTurbulenceGammaEquation)  
     *               write(13,110,advance='no') TGamma(IMonitor)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(13,110,advance='no') TReTheta(IMonitor)
      if(LSolveModifiedED)  
     *               write(13,110,advance='no') ModifiedED(IMonitor)
      if(LSolveTurbulentKL)  
     *               write(13,110,advance='no') TurbulentKL(IMonitor)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                write(13,110,advance='no') Temperature(IMonitor)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *          write(13,110,advance='no') Htotal(IMonitor)
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(13,111,advance='no') 
     *                                        rField(IMonitor,irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(13,111,advance='no') 
     *                                        Scalar(IMonitor,iScalar)
      enddo
      if(LSolveLambdaELEEquation) then
        write(13,110,advance='no') LambdaELE(IMonitor)
        write(13,100,advance='no') uVelocity(IMonitor),
     *             vVelocity(IMonitor),wVelocity(IMonitor)
      endif
c      
      write(13,60)
c
c---  Write absolute residuals
c
      write(13,120,advance='no') 
      write(13,130,advance='no') nIter
      if(LsolveContinuity) write(13,140,advance='no') ResorAbs(5)
      if(LsolveMomentum) write(13,150,advance='no') ResorAbs(1),
     *       ResorAbs(2),ResorAbs(3)
      if(LSolveTurbulenceKineticEnergy)  
     *               write(13,160,advance='no') ResorAbs(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(13,160,advance='no') ResorAbs(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(13,160,advance='no') ResorAbs(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(13,160,advance='no') ResorAbs(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(13,160,advance='no') ResorAbs(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(13,160,advance='no') ResorAbs(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(13,160,advance='no') ResorAbs(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(13,160,advance='no') ResorAbs(14)
      if(LSolveModifiedED)  
     *               write(13,160,advance='no') ResorAbs(11)
      if(LSolveTurbulentKL)  
     *               write(13,160,advance='no') ResorAbs(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                write(13,160,advance='no') ResorAbs(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                write(13,160,advance='no') ResorAbs(7)
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(13,161,advance='no') 
     *                                            ResorAbs(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(13,161,advance='no') 
     *                    ResorAbs(18+NumberOfrFieldsToSolve+iScalar)
      enddo
      if(LSolveLambdaELEEquation) 
     *               write(13,160,advance='no') ResorAbs(18)
      write(13,60)
c
c---  Write maximum residuals
c
      write(13,170,advance='no') 
      write(13,180,advance='no') nIter
      if(LsolveContinuity) write(13,190,advance='no') ResorMax(5)
      if(LsolveMomentum) write(13,200,advance='no') ResorMax(1),
     *       ResorMax(2),ResorMax(3)
      if(LSolveTurbulenceKineticEnergy)  
     *               write(13,210,advance='no') ResorMax(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(13,210,advance='no') ResorMax(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(13,210,advance='no') ResorMax(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(13,210,advance='no') ResorMax(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(13,210,advance='no') ResorMax(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(13,210,advance='no') ResorMax(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(13,210,advance='no') ResorMax(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(13,210,advance='no') ResorMax(14)
      if(LSolveModifiedED)  
     *               write(13,210,advance='no') ResorMax(11)
      if(LSolveTurbulentKL)  
     *               write(13,210,advance='no') ResorMax(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                 write(13,210,advance='no') ResorMax(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                 write(13,210,advance='no') ResorMax(7)
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(13,211,advance='no') 
     *                                          ResorMax(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(13,211,advance='no') 
     *                   ResorMax(18+NumberOfrFieldsToSolve+iScalar)
      enddo
      if(LSolveLambdaELEEquation) 
     *               write(13,210,advance='no') ResorMax(18)
      write(13,60)
c
c---  Write rms residuals
c
      write(13,220,advance='no') 
      write(13,230,advance='no') nIter
      if(LsolveContinuity) write(13,240,advance='no') ResorRMS(5)
      if(LsolveMomentum) write(13,250,advance='no') ResorRMS(1),
     *       ResorRMS(2),ResorRMS(3)
      if(LSolveTurbulenceKineticEnergy)  
     *               write(13,260,advance='no') ResorRMS(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(13,260,advance='no') ResorRMS(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(13,260,advance='no') ResorRMS(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(13,260,advance='no') ResorRMS(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(13,260,advance='no') ResorRMS(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(13,260,advance='no') ResorRMS(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(13,260,advance='no') ResorRMS(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(13,260,advance='no') ResorRMS(14)
      if(LSolveModifiedED)  
     *               write(13,260,advance='no') ResorRMS(11)
      if(LSolveTurbulentKL)  
     *               write(13,260,advance='no') ResorRMS(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                write(13,260,advance='no') ResorRMS(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                write(13,260,advance='no') ResorRMS(7)
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(13,261,advance='no') 
     *                                          ResorRMS(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(13,261,advance='no') 
     *                   ResorRMS(18+NumberOfrFieldsToSolve+iScalar)
      enddo
      if(LSolveLambdaELEEquation) 
     *                write(13,260,advance='no') ResorRMS(18)
      write(13,60)
c
c---  Write scaled residuals
c
      write(13,270,advance='no') 
      write(13,280,advance='no') nIter
      if(LsolveContinuity) write(13,290,advance='no') ResorScaled(5)
      if(LsolveMomentum) write(13,300,advance='no') ResorScaled(1),
     *       ResorScaled(2),ResorScaled(3)
      if(LSolveTurbulenceKineticEnergy)  
     *               write(13,310,advance='no') ResorScaled(8)
      if(LSolveTurbulenceDissipationRate)  
     *               write(13,310,advance='no') ResorScaled(9)
      if(LSolveTurbulenceV2Equation)  
     *               write(13,310,advance='no') ResorScaled(15)
      if(LSolveTurbulenceZetaEquation)  
     *               write(13,310,advance='no') ResorScaled(16)
      if(LSolveTurbulencefRelaxationEquation)  
     *               write(13,310,advance='no') ResorScaled(17)
      if(LSolveTurbulenceSpecificDissipationRate)  
     *               write(13,310,advance='no') ResorScaled(10)
      if(LSolveTurbulenceGammaEquation)  
     *               write(13,310,advance='no') ResorScaled(13)
      if(LSolveTurbulenceReynoldsThetaEquation)  
     *               write(13,310,advance='no') ResorScaled(14)
      if(LSolveModifiedED)  
     *               write(13,310,advance='no') ResorScaled(11)
      if(LSolveTurbulentKL)  
     *               write(13,310,advance='no') ResorScaled(12)
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                 write(13,310,advance='no') ResorScaled(6)
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                 write(13,310,advance='no') ResorScaled(7)
      do irField=1,NumberOfrFieldsToSolve
        if(LsolverField(irField)) write(13,311,advance='no') 
     *                                         ResorScaled(18+irField)
      enddo
      do iScalar=1,NumberOfScalarsToSolve
        if(LsolveScalar(iScalar)) write(13,311,advance='no') 
     *                  ResorScaled(18+NumberOfrFieldsToSolve+iScalar)
      enddo
      if(LSolveLambdaELEEquation) 
     *                 write(13,310,advance='no') ResorScaled(18)
      write(13,320)
c
c--- Check for convergence
c
      source=0.
c
      if(NstopType.eq.1) then
c
        if(LsolveMomentum) source=
     *        DMAX1(source,dabs(ResorAbs(1)),
     *                      dabs(ResorAbs(2)),dabs(ResorAbs(3)))
        if(LsolveContinuity) source=DMAX1(source,dabs(ResorAbs(5)))
        if(LSolveTurbulenceKineticEnergy) 
     *              source=DMAX1(source,dabs(ResorAbs(8)))
        if(LSolveTurbulenceDissipationRate) 
     *              source=DMAX1(source,dabs(ResorAbs(9))) 
        if(LSolveTurbulenceV2Equation)  
     *              source=DMAX1(source,dabs(ResorAbs(15))) 
        if(LSolveTurbulenceZetaEquation)  
     *              source=DMAX1(source,dabs(ResorAbs(16))) 
        if(LSolveTurbulencefRelaxationEquation)  
     *              source=DMAX1(source,dabs(ResorAbs(17))) 
        if(LSolveTurbulenceSpecificDissipationRate) 
     *              source=DMAX1(source,dabs(ResorAbs(10)))
        if(LSolveTurbulenceGammaEquation)  
     *              source=DMAX1(source,dabs(ResorAbs(13)))
        if(LSolveTurbulenceReynoldsThetaEquation)  
     *              source=DMAX1(source,dabs(ResorAbs(14)))
        if(LSolveModifiedED) 
     *              source=DMAX1(source,dabs(ResorAbs(11)))
        if(LSolveTurbulentKL) 
     *              source=DMAX1(source,dabs(ResorAbs(12)))
        if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *               source=DMAX1(source,dabs(ResorAbs(6)))
        if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *               source=DMAX1(source,dabs(ResorAbs(7)))
      if(LSolveLambdaELEEquation) 
     *               source=DMAX1(source,dabs(ResorAbs(18)))
c
        do irField=1,NumberOfrFieldsToSolve
          if(LsolverField(irField)) 
     *                 source=DMAX1(source,dabs(ResorAbs(18+irField)))
        enddo
c
        do iScalar=1,NumberOfScalarsToSolve
          if(LsolveScalar(iScalar)) 
     *       source=DMAX1(source,
     *            dabs(ResorAbs(18+NumberOfrFieldsToSolve+iScalar)))
        enddo
c
      elseif(NstopType.eq.2) then
c
        if(LsolveMomentum) 
     *         source=DMAX1(source,dabs(ResorMax(1)),
     *                        dabs(ResorMax(2)),dabs(ResorMax(3)))
        if(LsolveContinuity) 
     *         source=DMAX1(source,dabs(ResorMax(5)))
        if(LSolveTurbulenceKineticEnergy) 
     *              source=DMAX1(source,dabs(ResorMax(8)))
        if(LSolveTurbulenceDissipationRate) 
     *              source=DMAX1(source,dabs(ResorMax(9))) 
        if(LSolveTurbulenceV2Equation)  
     *              source=DMAX1(source,dabs(ResorMax(15))) 
        if(LSolveTurbulenceZetaEquation)  
     *              source=DMAX1(source,dabs(ResorMax(16))) 
        if(LSolveTurbulencefRelaxationEquation)  
     *              source=DMAX1(source,dabs(ResorMax(17))) 
        if(LSolveTurbulenceSpecificDissipationRate) 
     *              source=DMAX1(source,dabs(ResorMax(10)))
        if(LSolveTurbulenceGammaEquation)  
     *              source=DMAX1(source,dabs(ResorMax(13)))
        if(LSolveTurbulenceReynoldsThetaEquation)  
     *              source=DMAX1(source,dabs(ResorMax(14)))
        if(LSolveModifiedED) 
     *              source=DMAX1(source,dabs(ResorMax(11)))
        if(LSolveTurbulentKL) 
     *              source=DMAX1(source,dabs(ResorMax(12)))
      if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                source=DMAX1(source,dabs(ResorMax(6)))
      if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                source=DMAX1(source,dabs(ResorMax(7)))
      if(LSolveLambdaELEEquation) 
     *              source=DMAX1(source,dabs(ResorMax(18)))
c
        do irField=1,NumberOfrFieldsToSolve
          if(LsolverField(irField)) 
     *                 source=DMAX1(source,dabs(ResorMax(18+irField)))
        enddo
c
        do iScalar=1,NumberOfScalarsToSolve
          if(LsolveScalar(iScalar)) 
     *         source=DMAX1(source,
     *            dabs(ResorMax(18+NumberOfrFieldsToSolve+iScalar)))
        enddo
c
      elseif(NstopType.eq.3) then
c
        if(LsolveMomentum) 
     *         source=DMAX1(source,dabs(ResorRMS(1)),
     *                         dabs(ResorRMS(2)),dabs(ResorRMS(3)))
        if(LsolveContinuity) 
     *         source=DMAX1(source,dabs(ResorRMS(5)))
        if(LSolveTurbulenceKineticEnergy) 
     *              source=DMAX1(source,dabs(ResorRMS(8)))
        if(LSolveTurbulenceDissipationRate) 
     *              source=DMAX1(source,dabs(ResorRMS(9))) 
        if(LSolveTurbulenceV2Equation)  
     *              source=DMAX1(source,dabs(ResorRMS(15))) 
        if(LSolveTurbulenceZetaEquation)  
     *              source=DMAX1(source,dabs(ResorRMS(16))) 
        if(LSolveTurbulencefRelaxationEquation)  
     *              source=DMAX1(source,dabs(ResorRMS(17))) 
        if(LSolveTurbulenceSpecificDissipationRate) 
     *              source=DMAX1(source,dabs(ResorRMS(10)))
        if(LSolveTurbulenceGammaEquation)  
     *              source=DMAX1(source,dabs(ResorRMS(13)))
        if(LSolveTurbulenceReynoldsThetaEquation)  
     *              source=DMAX1(source,dabs(ResorRMS(14)))
        if(LSolveModifiedED) 
     *              source=DMAX1(source,dabs(ResorRMS(11)))
        if(LSolveTurbulentKL) 
     *              source=DMAX1(source,dabs(ResorRMS(12)))
        if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *                source=DMAX1(source,dabs(ResorRMS(6)))
        if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *                source=DMAX1(source,dabs(ResorRMS(7)))
      if(LSolveLambdaELEEquation) 
     *              source=DMAX1(source,dabs(ResorRMS(18)))
c
        do irField=1,NumberOfrFieldsToSolve
          if(LsolverField(irField)) 
     *                 source=DMAX1(source,dabs(ResorRMS(18+irField)))
        enddo
c
        do iScalar=1,NumberOfScalarsToSolve
          if(LsolveScalar(iScalar)) 
     *            source=DMAX1(source,
     *               dabs(ResorRMS(18+NumberOfrFieldsToSolve+iScalar)))
        enddo
c
      elseif(NstopType.eq.4) then
c
        if(LsolveMomentum) 
     *             source=DMAX1(source,dabs(ResorScaled(1)),
     *                    dabs(ResorScaled(2)),dabs(ResorScaled(3)))
        if(LsolveContinuity) 
     *         source=DMAX1(source,dabs(ResorScaled(5)))
        if(LSolveTurbulenceKineticEnergy) 
     *              source=DMAX1(source,dabs(ResorScaled(8)))
        if(LSolveTurbulenceDissipationRate) 
     *              source=DMAX1(source,dabs(ResorScaled(9))) 
        if(LSolveTurbulenceV2Equation)  
     *              source=DMAX1(source,dabs(ResorScaled(15))) 
        if(LSolveTurbulenceZetaEquation)  
     *              source=DMAX1(source,dabs(ResorScaled(16))) 
        if(LSolveTurbulencefRelaxationEquation)  
     *              source=DMAX1(source,dabs(ResorScaled(17))) 
        if(LSolveTurbulenceSpecificDissipationRate) 
     *              source=DMAX1(source,dabs(ResorScaled(10)))
        if(LSolveTurbulenceGammaEquation)  
     *              source=DMAX1(source,dabs(ResorScaled(13)))
        if(LSolveTurbulenceReynoldsThetaEquation)  
     *              source=DMAX1(source,dabs(ResorScaled(14)))
        if(LSolveModifiedED) 
     *              source=DMAX1(source,dabs(ResorScaled(11)))
        if(LSolveTurbulentKL) 
     *              source=DMAX1(source,dabs(ResorScaled(12)))
        if(LsolveEnergy.and.EnergyEquation.eq.'temperature') 
     *               source=DMAX1(source,dabs(ResorScaled(6)))
        if(LsolveEnergy.and.EnergyEquation.eq.'htotal') 
     *               source=DMAX1(source,dabs(ResorScaled(7)))
      if(LSolveLambdaELEEquation) 
     *              source=DMAX1(source,dabs(ResorScaled(18)))
c
        do irField=1,NumberOfrFieldsToSolve
          if(LsolverField(irField)) 
     *              source=DMAX1(source,dabs(ResorScaled(18+irField)))
        enddo
c
        do iScalar=1,NumberOfScalarsToSolve
          if(LsolveScalar(iScalar)) 
     *        source=DMAX1(source,
     *           dabs(ResorScaled(18+NumberOfrFieldsToSolve+iScalar)))
        enddo
c
      endif
c
      return
      end