MODULE Variables4

implicit none
double precision, save, dimension(:), allocatable :: Sc
double precision, save, dimension(:), allocatable :: Sb
double precision, save, dimension(:), allocatable :: ScMomentumx
double precision, save, dimension(:), allocatable :: SbMomentumx
double precision, save, dimension(:), allocatable :: ScMomentumy
double precision, save, dimension(:), allocatable :: SbMomentumy
double precision, save, dimension(:), allocatable :: ScMomentumz
double precision, save, dimension(:), allocatable :: SbMomentumz
double precision, save, dimension(:), allocatable :: ScEnergy
double precision, save, dimension(:), allocatable :: SbEnergy
double precision, save, dimension(:), allocatable :: ScTKE
double precision, save, dimension(:), allocatable :: SbTKE
double precision, save, dimension(:), allocatable :: ScTED
double precision, save, dimension(:), allocatable :: SbTED
double precision, save, dimension(:), allocatable :: ScTOmega
double precision, save, dimension(:), allocatable :: SbTOmega
double precision, save, dimension(:), allocatable :: ScTurbulentKL
double precision, save, dimension(:), allocatable :: SbTurbulentKL
double precision, save, dimension(:), allocatable :: ScModifiedED
double precision, save, dimension(:), allocatable :: SbModifiedED
double precision, save, dimension(:), allocatable :: ScTGamma
double precision, save, dimension(:), allocatable :: SbTGamma
double precision, save, dimension(:), allocatable :: ScTReTheta
double precision, save, dimension(:), allocatable :: SbTReTheta
double precision, save, dimension(:), allocatable :: ScTfRelaxation
double precision, save, dimension(:), allocatable :: SbTfRelaxation
double precision, save, dimension(:), allocatable :: ScTurbulentV2
double precision, save, dimension(:), allocatable :: SbTurbulentV2
double precision, save, dimension(:), allocatable :: ScTurbulentZeta
double precision, save, dimension(:), allocatable :: SbTurbulentZeta
double precision, save, dimension(:), allocatable :: ScLambdaELE
double precision, save, dimension(:), allocatable :: SbLambdaELE

end MODULE Variables4
