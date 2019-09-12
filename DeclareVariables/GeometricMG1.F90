MODULE GeometricMG1

implicit none
integer, save :: AllocateMGGeometricStorage=50
integer, save, dimension(:),allocatable :: ijbeginNode
integer, save, dimension(:),allocatable :: ijbeginENode
integer, save, dimension(:),allocatable :: ijbeginNElement
integer, save, dimension(:),allocatable :: FlagNodesMG
integer, save, dimension(:),allocatable :: FlagNodesMGTemp
integer, save, dimension(:),allocatable :: NumberOfNodesMG
integer, save, dimension(:),allocatable :: NumberOfElementsConnectedToNode
integer, save, dimension(:),allocatable :: ListOfElementNodesMG
integer, save, dimension(:),allocatable :: ListOfElementsConnectedToNodeMG
integer, save, dimension(:),allocatable :: NodesMG
integer, save, dimension(:),allocatable :: NodesMGTemp
integer, save, dimension(:),allocatable :: NumbOfElementNodesMG
integer, save, dimension(:),allocatable :: NumberOfElementsConnectedToNodeMG
double precision, save, dimension(:),allocatable :: VolumeMG
double precision, save, dimension(:),allocatable :: VolumeMGTemp
integer, save, dimension(:),allocatable :: ListOfElementNodesMGTemp
integer, save, dimension(:),allocatable :: NumbOfElementNodesMGTemp

end MODULE GeometricMG1