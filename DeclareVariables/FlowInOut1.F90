MODULE FlowInOut1
!
  implicit none
  double precision, save :: mdotIn
  double precision, save, dimension(:), allocatable :: MassFlowFraction
  double precision, save, dimension(:), allocatable :: mdotOut
  double precision, save, dimension(:), allocatable :: mdot1
  double precision, save, dimension(:), allocatable :: MassFlowRate
  Logical, save, dimension(:), allocatable ::LPrintMassFlowRate
!
end MODULE FlowInOut1
