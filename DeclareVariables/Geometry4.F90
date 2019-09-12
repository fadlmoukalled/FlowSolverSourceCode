      MODULE Geometry4
      Implicit none
      double precision, save, dimension(:), allocatable :: xc,yc,zc
      double precision, save, dimension(:), allocatable :: Volume
!      double precision, save, dimension(:), allocatable :: EdgeCentroidx
!      double precision, save, dimension(:), allocatable :: EdgeCentroidy
!      double precision, save, dimension(:), allocatable :: EdgeCentroidz

      double precision, save, dimension(:), allocatable :: FaceCentroidx
      double precision, save, dimension(:), allocatable :: FaceCentroidy
      double precision, save, dimension(:), allocatable :: FaceCentroidz
      double precision, save, dimension(:), allocatable :: FaceArea
      double precision, save, dimension(:), allocatable :: FaceAreax
      double precision, save, dimension(:), allocatable :: FaceAreay
      double precision, save, dimension(:), allocatable :: FaceAreaz
      double precision, save, dimension(:), allocatable :: FaceAreanx
      double precision, save, dimension(:), allocatable :: FaceAreany
      double precision, save, dimension(:), allocatable :: FaceAreanz
      double precision, save, dimension(:), allocatable :: DistanceCF
      double precision, save, dimension(:), allocatable :: DistanceCFx
      double precision, save, dimension(:), allocatable :: DistanceCFy
      double precision, save, dimension(:), allocatable :: DistanceCFz
      double precision, save, dimension(:), allocatable :: DistanceCFux
      double precision, save, dimension(:), allocatable :: DistanceCFuy
      double precision, save, dimension(:), allocatable :: DistanceCFuz
      double precision, save, dimension(:), allocatable :: GFactorCF
      double precision, save, dimension(:,:), allocatable :: BFaceCentroidx
      double precision, save, dimension(:,:), allocatable :: BFaceCentroidy
      double precision, save, dimension(:,:), allocatable :: BFaceCentroidz
      double precision, save, dimension(:,:), allocatable :: BFaceArea
      double precision, save, dimension(:,:), allocatable :: BFaceAreax
      double precision, save, dimension(:,:), allocatable :: BFaceAreay
      double precision, save, dimension(:,:), allocatable :: BFaceAreaz
      double precision, save, dimension(:,:), allocatable :: BFaceAreanx
      double precision, save, dimension(:,:), allocatable :: BFaceAreany
      double precision, save, dimension(:,:), allocatable :: BFaceAreanz
      double precision, save, dimension(:,:), allocatable :: BDistanceCF
      double precision, save, dimension(:,:), allocatable :: BDistanceCFx
      double precision, save, dimension(:,:), allocatable :: BDistanceCFy
      double precision, save, dimension(:,:), allocatable :: BDistanceCFz
      double precision, save, dimension(:,:), allocatable :: BDistanceCFux
      double precision, save, dimension(:,:), allocatable :: BDistanceCFuy
      double precision, save, dimension(:,:), allocatable :: BDistanceCFuz

      double precision, save, dimension(:), allocatable :: FaceE
      double precision, save, dimension(:), allocatable :: FaceEx
      double precision, save, dimension(:), allocatable :: FaceEy
      double precision, save, dimension(:), allocatable :: FaceEz
      double precision, save, dimension(:), allocatable :: gDiff
      double precision, save, dimension(:), allocatable :: FaceT
      double precision, save, dimension(:), allocatable :: FaceTx
      double precision, save, dimension(:), allocatable :: FaceTy
      double precision, save, dimension(:), allocatable :: FaceTz
      double precision, save, dimension(:,:), allocatable :: BFaceE
      double precision, save, dimension(:,:), allocatable :: BFaceEx
      double precision, save, dimension(:,:), allocatable :: BFaceEy
      double precision, save, dimension(:,:), allocatable :: BFaceEz
      double precision, save, dimension(:,:), allocatable :: BgDiff
      double precision, save, dimension(:,:), allocatable :: BFaceT
      double precision, save, dimension(:,:), allocatable :: BFaceTx
      double precision, save, dimension(:,:), allocatable :: BFaceTy
      double precision, save, dimension(:,:), allocatable :: BFaceTz



      double precision, save, dimension(:), allocatable :: FaceAreap
      double precision, save, dimension(:), allocatable :: FaceAreaxp
      double precision, save, dimension(:), allocatable :: FaceAreayp
      double precision, save, dimension(:), allocatable :: FaceAreazp
      double precision, save, dimension(:,:), allocatable :: BFaceAreap
      double precision, save, dimension(:,:), allocatable :: BFaceAreaxp
      double precision, save, dimension(:,:), allocatable :: BFaceAreayp
      double precision, save, dimension(:,:), allocatable :: BFaceAreazp
      double precision, save, dimension(:), allocatable :: FaceEp
      double precision, save, dimension(:), allocatable :: FaceExp
      double precision, save, dimension(:), allocatable :: FaceEyp
      double precision, save, dimension(:), allocatable :: FaceEzp
      double precision, save, dimension(:), allocatable :: gDiffp
      double precision, save, dimension(:), allocatable :: FaceTp
      double precision, save, dimension(:), allocatable :: FaceTxp
      double precision, save, dimension(:), allocatable :: FaceTyp
      double precision, save, dimension(:), allocatable :: FaceTzp
      double precision, save, dimension(:,:), allocatable :: BFaceEp
      double precision, save, dimension(:,:), allocatable :: BFaceExp
      double precision, save, dimension(:,:), allocatable :: BFaceEyp
      double precision, save, dimension(:,:), allocatable :: BFaceEzp
      double precision, save, dimension(:,:), allocatable :: BgDiffp
      double precision, save, dimension(:,:), allocatable :: BFaceTp
      double precision, save, dimension(:,:), allocatable :: BFaceTxp
      double precision, save, dimension(:,:), allocatable :: BFaceTyp
      double precision, save, dimension(:,:), allocatable :: BFaceTzp

      end MODULE Geometry4