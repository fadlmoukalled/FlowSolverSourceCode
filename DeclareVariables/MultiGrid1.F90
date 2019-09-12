MODULE MultiGrid1


    implicit none

        integer, save, dimension(:), allocatable :: NumberOfElementsMG
        integer, save, dimension(:), allocatable :: NElementFacesMG
        integer, save, dimension(:), allocatable :: NBChildrenMaxMG
        integer, save, dimension(:), allocatable :: ijbeginE
        integer, save, dimension(:), allocatable :: ijbeginN
        integer, save, dimension(:), allocatable :: ijbeginNOld
        integer, save, dimension(:), allocatable :: ElementNeighborMG
        integer, save, dimension(:), allocatable :: NumberofElementNeighborsMG
        integer, save, dimension(:), allocatable :: Parents
        integer, save, dimension(:), allocatable :: Children
        integer, save, dimension(:), allocatable :: NumberOfChildren
        integer, save, dimension(:), allocatable :: NeighborOfParent
        integer, save, dimension(:), allocatable :: NumberOfNeighborOfParent
        integer, save, dimension(:), allocatable :: ElementParent

        double precision, save, dimension(:), allocatable :: maxanb
        double precision, save, dimension(:), allocatable :: anbMG
        double precision, save, dimension(:), allocatable :: acMG
        double precision, save, dimension(:), allocatable :: bcMG
        double precision, save, dimension(:), allocatable :: Residuals
        double precision, save, dimension(:), allocatable :: dphiMG




end MODULE MultiGrid1