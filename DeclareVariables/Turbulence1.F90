MODULE Turbulence1
implicit none

  integer, save :: ModelNumber
  double precision, save :: sigT,sigTKE,sigTED,ce1,ce2,ce3
  double precision, save :: cmu,cmu25,cmu75,cmu50,ctrans,cappa,cc
  double precision, save :: alpha,betta,sigB,sigC,c1limiter,a1sst
  double precision, save :: sigTKE1,sigTED1,alpha1,betta1
  double precision, save :: sigTKE2,sigTED2,alpha2,betta2
  double precision, save :: betta0,sigdo,Clim,CrC,alphaStar0,alpha0
  double precision, save :: ReB,ReK,ReW,c2,c3,ct3,ct4
  double precision, save :: cb1,cb2,cw1,cw2,cw3,cv1,sigMED,c5,cn1,crot
  double precision, save :: cr1,cr2,cr3,Cr1SA,Cr1WA
  double precision, save :: c1KW,c1KE,sigKW,sigKE,cw,cm,c2KW,c2KE
  double precision, save :: sigKL,xi1,xi2,xi3,c11,c12,cd1
  double precision, save :: A0,ca1,ca2,cot,s1,sigTGamma,sigTReTheta
  double precision, save :: sigRT,AmuRT,alphaRT,BettaRT,C1,phiRT
!
  double precision, save :: As,Av,Abp,Anat,Ats,Cbpcrit
  double precision, save :: Cnc,Cnatcrit,Cint,Ctscrit,Crnat
  double precision, save :: Cr,Ca0,Css,Ctl,Cwr,Clambda
!
  double precision, save :: FlengthGama,Ctu1,Ctu2,Ctu3
  double precision, save :: Cpg1,Cpg2,Cpg3,Cpg1lim,Cpg2lim
  double precision, save :: Ck,Csep,Reoclim
!
  double precision, save :: Anut1,Anut2,Cnut0,Cnut1,Cnut2,Cnut3
  double precision, save :: Cnut4,Cnut5,Cnut6,Cnut7,Cnut8
!
  double precision, save :: C1V2f,C2V2f,CLV2f,CetaV2f
  double precision, save :: C1Zeta,C2Zeta,CLZeta,CetaZeta
!
  double precision, save :: LengthScale
!
  double precision, save :: alpha0RNG,eta0RNG,betaRNG
  double precision, save, dimension(:), allocatable :: C2eRNG
  double precision, save, dimension(:), allocatable :: sigTKERNG
  double precision, save, dimension(:), allocatable :: sigTEDRNG
  double precision, save, dimension(:,:), allocatable :: BsigTKERNG
  double precision, save, dimension(:,:), allocatable :: BsigTEDRNG
!
  double precision, save, dimension(:), allocatable :: SigScalar
  double precision, save, dimension(:), allocatable :: rhok
  double precision, save, dimension(:,:), allocatable :: Brhok
  double precision, save, dimension(:), allocatable :: drhokdx
  double precision, save, dimension(:), allocatable :: drhokdy
  double precision, save, dimension(:), allocatable :: drhokdz
  double precision, save, dimension(:,:), allocatable :: Bdrhokdx
  double precision, save, dimension(:,:), allocatable :: Bdrhokdy
  double precision, save, dimension(:,:), allocatable :: Bdrhokdz
  double precision, save, dimension(:), allocatable :: rhoTED
  double precision, save, dimension(:,:), allocatable :: BrhoTED
  double precision, save, dimension(:), allocatable :: yplus
  double precision, save, dimension(:), allocatable :: yplusT
  double precision, save, dimension(:,:), allocatable :: yplusS
  double precision, save, dimension(:,:), allocatable :: yplusR
  
  double precision, save, dimension(:), allocatable :: ustar
  double precision, save, dimension(:), allocatable :: uplus
  double precision, save, dimension(:), allocatable :: duplusdyplus
  double precision, save, dimension(:), allocatable :: KsPlus
  double precision, save, dimension(:), allocatable :: uTau
  double precision, save, dimension(:), allocatable :: ystar

  double precision, save, dimension(:), allocatable :: RoughccM
  double precision, save, dimension(:), allocatable :: RoughccE
  double precision, save, dimension(:,:), allocatable :: RoughccS
  
  double precision, save, dimension(:), allocatable :: WallViscosity
  double precision, save, dimension(:), allocatable :: WallViscosityT
  double precision, save, dimension(:,:), allocatable :: WallViscosityS
  double precision, save, dimension(:,:), allocatable :: WallViscosityR
  
  double precision, save, dimension(:), allocatable :: F1factor
  double precision, save, dimension(:,:), allocatable :: BF1factor
  double precision, save, dimension(:), allocatable :: F2factor
  double precision, save, dimension(:,:), allocatable :: BF2factor
  double precision, save, dimension(:), allocatable :: F3factor
  double precision, save, dimension(:,:), allocatable :: BF3factor
  double precision, save, dimension(:), allocatable :: F4factor
  double precision, save, dimension(:,:), allocatable :: BF4factor

  double precision, save, dimension(:), allocatable :: S11
  double precision, save, dimension(:), allocatable :: S12
  double precision, save, dimension(:), allocatable :: S13
  double precision, save, dimension(:), allocatable :: S22
  double precision, save, dimension(:), allocatable :: S23
  double precision, save, dimension(:), allocatable :: S33
  double precision, save, dimension(:,:), allocatable :: BS11
  double precision, save, dimension(:,:), allocatable :: BS12
  double precision, save, dimension(:,:), allocatable :: BS13
  double precision, save, dimension(:,:), allocatable :: BS22
  double precision, save, dimension(:,:), allocatable :: BS23
  double precision, save, dimension(:,:), allocatable :: BS33
!  
  double precision, save, dimension(:), allocatable :: S11Old
  double precision, save, dimension(:), allocatable :: S12Old
  double precision, save, dimension(:), allocatable :: S13Old
  double precision, save, dimension(:), allocatable :: S22Old
  double precision, save, dimension(:), allocatable :: S23Old
  double precision, save, dimension(:), allocatable :: S33Old
!  
  double precision, save, dimension(:), allocatable :: S11OldOld
  double precision, save, dimension(:), allocatable :: S12OldOld
  double precision, save, dimension(:), allocatable :: S13OldOld
  double precision, save, dimension(:), allocatable :: S22OldOld
  double precision, save, dimension(:), allocatable :: S23OldOld
  double precision, save, dimension(:), allocatable :: S33OldOld
!  
  double precision, save, dimension(:), allocatable :: DS11Dt
  double precision, save, dimension(:), allocatable :: DS12Dt
  double precision, save, dimension(:), allocatable :: DS13Dt
  double precision, save, dimension(:), allocatable :: DS22Dt
  double precision, save, dimension(:), allocatable :: DS23Dt
  double precision, save, dimension(:), allocatable :: DS33Dt
!
  double precision, save, dimension(:), allocatable :: StrainRate
  double precision, save, dimension(:,:), allocatable :: BStrainRate
!  
!
  double precision, save, dimension(:), allocatable :: LScale
  double precision, save, dimension(:,:), allocatable :: BLScale
  double precision, save, dimension(:), allocatable :: TScale
  double precision, save, dimension(:,:), allocatable :: BTScale
  double precision, save, dimension(:), allocatable :: Ce1Coefficient
!
  double precision, save, dimension(:), allocatable :: W11
  double precision, save, dimension(:), allocatable :: W12
  double precision, save, dimension(:), allocatable :: W13
  double precision, save, dimension(:), allocatable :: W22
  double precision, save, dimension(:), allocatable :: W23
  double precision, save, dimension(:), allocatable :: W33
  double precision, save, dimension(:,:), allocatable :: BW11
  double precision, save, dimension(:,:), allocatable :: BW12
  double precision, save, dimension(:,:), allocatable :: BW13
  double precision, save, dimension(:,:), allocatable :: BW22
  double precision, save, dimension(:,:), allocatable :: BW23
  double precision, save, dimension(:,:), allocatable :: BW33
  double precision, save, dimension(:), allocatable :: Vorticity
  double precision, save, dimension(:,:), allocatable :: BVorticity

  double precision, save, dimension(:), allocatable :: Tau11
  double precision, save, dimension(:), allocatable :: Tau12
  double precision, save, dimension(:), allocatable :: Tau13
  double precision, save, dimension(:), allocatable :: Tau22
  double precision, save, dimension(:), allocatable :: Tau23
  double precision, save, dimension(:), allocatable :: Tau33

  double precision, save, dimension(:), allocatable :: TauWall
  
  double precision, save, dimension(:), allocatable :: alfa
  double precision, save, dimension(:), allocatable :: alfaStar
  double precision, save, dimension(:), allocatable :: bettaStar
  double precision, save, dimension(:,:), allocatable :: Balfa
  double precision, save, dimension(:,:), allocatable :: BalfaStar
  double precision, save, dimension(:,:), allocatable :: BbettaStar

  double precision, save, dimension(:), allocatable :: fmuCoefficient
  double precision, save, dimension(:,:), allocatable :: BfmuCoefficient
  double precision, save, dimension(:), allocatable :: f1Coefficient
  double precision, save, dimension(:), allocatable :: f2Coefficient
  double precision, save, dimension(:), allocatable :: LTKE
  double precision, save, dimension(:), allocatable :: LTED
  double precision, save, dimension(:), allocatable :: ReT
  double precision, save, dimension(:,:), allocatable :: BReT

  double precision, save, dimension(:), allocatable :: sqrtTurbulentKE
  double precision, save, dimension(:,:), allocatable :: BsqrtTurbulentKE
  double precision, save, dimension(:), allocatable :: sqrtTKEGradx
  double precision, save, dimension(:), allocatable :: sqrtTKEGrady
  double precision, save, dimension(:), allocatable :: sqrtTKEGradz
  double precision, save, dimension(:,:), allocatable :: BsqrtTKEGradx
  double precision, save, dimension(:,:), allocatable :: BsqrtTKEGrady
  double precision, save, dimension(:,:), allocatable :: BsqrtTKEGradz
  double precision, save, dimension(:), allocatable :: SRateGradx
  double precision, save, dimension(:), allocatable :: SRateGrady
  double precision, save, dimension(:), allocatable :: SRateGradz
  double precision, save, dimension(:,:), allocatable :: BSRateGradx
  double precision, save, dimension(:,:), allocatable :: BSRateGrady
  double precision, save, dimension(:,:), allocatable :: BSRateGradz

  double precision, save, dimension(:), allocatable :: Stelda
  double precision, save, dimension(:), allocatable :: fwCoefficient
  double precision, save, dimension(:), allocatable :: ft2Coefficient
  double precision, save, dimension(:), allocatable :: fr1Coefficient
  double precision, save, dimension(:), allocatable :: fv1Coefficient
  double precision, save, dimension(:,:), allocatable :: Bfv1Coefficient
  double precision, save, dimension(:), allocatable :: fnCoefficient
  double precision, save, dimension(:,:), allocatable :: BfnCoefficient
 
  double precision, save, dimension(:), allocatable :: f1WA
  double precision, save, dimension(:,:), allocatable :: Bf1WA
  double precision, save, dimension(:), allocatable :: fmuWA
  double precision, save, dimension(:,:), allocatable :: BfmuWA
  
  double precision, save, dimension(:), allocatable :: cmuR
  double precision, save, dimension(:,:), allocatable :: BcmuR
  double precision, save, dimension(:), allocatable :: c1R
  double precision, save, dimension(:,:), allocatable :: Bc1R


  double precision, save, dimension(:), allocatable :: LambdaT
  double precision, save, dimension(:,:), allocatable :: BLambdaT
  double precision, save, dimension(:), allocatable :: LambdaEff
  double precision, save, dimension(:,:), allocatable :: BLambdaEff
  double precision, save, dimension(:), allocatable :: coefficientFW
  double precision, save, dimension(:,:), allocatable :: BcoefficientFW
  double precision, save, dimension(:), allocatable :: TurbulentKTs
  double precision, save, dimension(:,:), allocatable :: BTurbulentKTs
  double precision, save, dimension(:), allocatable :: TurbulentViscosityTs
  double precision, save, dimension(:,:), allocatable :: BTurbulentViscosityTs
  double precision, save, dimension(:), allocatable :: TurbulentViscosityTl
  double precision, save, dimension(:,:), allocatable :: BTurbulentViscosityTl
  double precision, save, dimension(:), allocatable :: AlfaT
  double precision, save, dimension(:,:), allocatable :: BAlfaT
  double precision, save, dimension(:), allocatable :: AlfaTheta
  double precision, save, dimension(:,:), allocatable :: BAlfaTheta
!
  double precision, save, dimension(:), allocatable :: ProductionKT
  double precision, save, dimension(:), allocatable :: ProductionKL
  double precision, save, dimension(:,:), allocatable :: BProductionKT
  double precision, save, dimension(:,:), allocatable :: BProductionKL
  double precision, save, dimension(:), allocatable :: SourceRbp
  double precision, save, dimension(:), allocatable :: SourceRnat
!
  double precision, save, dimension(:), allocatable :: f2RT
  double precision, save, dimension(:,:), allocatable :: Bf2RT
  double precision, save, dimension(:), allocatable :: velRT
  double precision, save, dimension(:,:), allocatable :: BvelRT
  double precision, save, dimension(:), allocatable :: velRTGradx
  double precision, save, dimension(:,:), allocatable :: BvelRTGradx
  double precision, save, dimension(:), allocatable :: velRTGrady
  double precision, save, dimension(:,:), allocatable :: BvelRTGrady
  double precision, save, dimension(:), allocatable :: velRTGradz
  double precision, save, dimension(:,:), allocatable :: BvelRTGradz
  double precision, save, dimension(:), allocatable :: fmuKE
  double precision, save, dimension(:,:), allocatable :: BfmuKE
  double precision, save, dimension(:), allocatable :: fmuRT
  double precision, save, dimension(:,:), allocatable :: BfmuRT
  double precision, save, dimension(:), allocatable :: TurbulentViscosityKE
  double precision, save, dimension(:,:), allocatable :: BTurbulentViscosityKE
  double precision, save, dimension(:), allocatable :: TurbulentViscosityRT
  double precision, save, dimension(:,:), allocatable :: BTurbulentViscosityRT


  double precision, dimension(:), allocatable :: ModifiedEDGrad2x
  double precision, dimension(:), allocatable :: ModifiedEDGrad2y
  double precision, dimension(:), allocatable :: ModifiedEDGrad2z
  double precision, dimension(:,:), allocatable :: BModifiedEDGrad2x
  double precision, dimension(:,:), allocatable :: BModifiedEDGrad2y
  double precision, dimension(:,:), allocatable :: BModifiedEDGrad2z
  double precision, dimension(:), allocatable :: ModifiedMut
  double precision, dimension(:,:), allocatable :: BModifiedMut
  double precision, dimension(:), allocatable :: ModifiedMutGradx
  double precision, dimension(:), allocatable :: ModifiedMutGrady
  double precision, dimension(:), allocatable :: ModifiedMutGradz
  double precision, dimension(:,:), allocatable :: BModifiedMutGradx
  double precision, dimension(:,:), allocatable :: BModifiedMutGrady
  double precision, dimension(:,:), allocatable :: BModifiedMutGradz
  double precision, dimension(:), allocatable :: G1Nut92
  double precision, dimension(:), allocatable :: G2Nut92
  double precision, dimension(:), allocatable :: F1Nut92
  double precision, dimension(:), allocatable :: F2Nut92
  double precision, dimension(:), allocatable :: N1Nut92
  double precision, dimension(:), allocatable :: N2Nut92
  double precision, dimension(:,:), allocatable :: BN1Nut92
  double precision, dimension(:), allocatable :: N1Nut92Gradx
  double precision, dimension(:), allocatable :: N1Nut92Grady
  double precision, dimension(:), allocatable :: N1Nut92Gradz
  double precision, dimension(:,:), allocatable :: BN1Nut92Gradx
  double precision, dimension(:,:), allocatable :: BN1Nut92Grady
  double precision, dimension(:,:), allocatable :: BN1Nut92Gradz
  double precision, dimension(:), allocatable :: ModifiedDensGradx
  double precision, dimension(:), allocatable :: ModifiedDensGrady
  double precision, dimension(:), allocatable :: ModifiedDensGradz
  double precision, dimension(:,:), allocatable :: BModifiedDensGradx
  double precision, dimension(:,:), allocatable :: BModifiedDensGrady
  double precision, dimension(:,:), allocatable :: BModifiedDensGradz
  double precision, dimension(:), allocatable :: ModifiedDensGrad2x
  double precision, dimension(:), allocatable :: ModifiedDensGrad2y
  double precision, dimension(:), allocatable :: ModifiedDensGrad2z
  double precision, dimension(:,:), allocatable :: BModifiedDensGrad2x
  double precision, dimension(:,:), allocatable :: BModifiedDensGrad2y
  double precision, dimension(:,:), allocatable :: BModifiedDensGrad2z

  double precision, dimension(:), allocatable :: XiStarFmt

  double precision, dimension(:), allocatable :: TurbulentZetaGrad2x
  double precision, dimension(:), allocatable :: TurbulentZetaGrad2y
  double precision, dimension(:), allocatable :: TurbulentZetaGrad2z
  double precision, dimension(:,:), allocatable :: BTurbulentZetaGrad2x
  double precision, dimension(:,:), allocatable :: BTurbulentZetaGrad2y
  double precision, dimension(:,:), allocatable :: BTurbulentZetaGrad2z
  double precision, dimension(:), allocatable :: TurbulentZetaGradxy
  double precision, dimension(:), allocatable :: TurbulentZetaGradxz
  double precision, dimension(:,:), allocatable :: BTurbulentZetaGradxy
  double precision, dimension(:,:), allocatable :: BTurbulentZetaGradxz

end MODULE Turbulence1