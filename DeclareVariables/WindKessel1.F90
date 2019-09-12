MODULE WindKessel1
!
  implicit none
!
  double precision, save, dimension(:), allocatable :: ResistanceToBloodFlow
  double precision, save, dimension(:), allocatable :: TotalPeripheralResistance
  double precision, save, dimension(:), allocatable :: ComplianceC
  double precision, save, dimension(:), allocatable :: InertiaL
  double precision, save, dimension(:), allocatable :: OutletPressure
!
  double precision, save, dimension(:), allocatable :: MassFlowRateOld
  double precision, save, dimension(:), allocatable :: MassFlowRateOldOld
  double precision, save, dimension(:), allocatable :: OutletPressureOld
  double precision, save, dimension(:), allocatable :: OutletPressureOldOld
!  
end MODULE WindKessel1
