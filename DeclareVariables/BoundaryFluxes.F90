MODULE BoundaryFluxes
      implicit none
!
      double precision, save, dimension(:,:), allocatable :: HeatFlux
      double precision, save, dimension (:,:), allocatable :: Hinfinity
      double precision, save, dimension (:,:), allocatable :: Tinfinity
      
      double precision, save, dimension(:,:,:), allocatable :: ScalarFlux
      double precision, save, dimension (:,:,:), allocatable :: ScalarConvectionCoefficient
      double precision, save, dimension (:,:,:), allocatable :: ScalarPhiInfinity

      double precision, save, dimension(:,:,:), allocatable :: rFieldFlux
      double precision, save, dimension (:,:,:), allocatable :: rFieldConvectionCoefficient
      double precision, save, dimension (:,:,:), allocatable :: rFieldPhiInfinity

!
end MODULE BoundaryFluxes
