c
c#############################################################################################
c
      SUBROUTINE ComputeEffectiveDivergence(Variable)
c
c#############################################################################################
c
      use Variables1, only: mdot,Bmdot,effdiv
      use Geometry1, only: NumberOfBCSets
      use Geometry3, only: NIFaces,NIFaceOwner,NIFaceNeighbor,NBFaces,
     *                     NBFaceOwner
      use PhysicalProperties1, only: SpecificHeat,BspecificHeat,
     *                        SpecificHeatScalar,BSpecificHeatScalar
      use Geometry4, only: GFactorCF
      use Scalar2
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,k
      character*10 Variable
      double precision cpf,gf
c********************************************************************************************
c
      effdiv=0.
c
      if(Variable.eq.'velx'.or.Variable.eq.'vely'.or.
     *          Variable.eq.'velz'.or.Variable.eq.'rfield') then
c
c--- Internal faces
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          effdiv(i)=effdiv(i)+mdot(k)
          effdiv(j)=effdiv(j)-mdot(k)
c
        enddo
c
c--- Boundary faces
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            k=NBFaceOwner(i,j)           
c  
            effdiv(k)=effdiv(k)+Bmdot(i,j)
c
          enddo
        enddo
c
      elseif(Variable.eq.'tke'.or.Variable.eq.'ted'.or.
     *         Variable.eq.'tomega'.or.Variable.eq.'med'.or.
     *           Variable.eq.'tkl'.or.Variable.eq.'htotal'.or.
     *             Variable.eq.'tv2'.or.Variable.eq.'tzeta'.or.
     *              Variable.eq.'tgamma'.or.Variable.eq.'tretheta') then
c
c--- Internal faces
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          effdiv(i)=effdiv(i)+mdot(k)
          effdiv(j)=effdiv(j)-mdot(k)
c
        enddo
c
c--- Boundary faces
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            k=NBFaceOwner(i,j)           
c  
            effdiv(k)=effdiv(k)+Bmdot(i,j)
c
          enddo
        enddo
c
      elseif(Variable.eq.'temp') then
c
c--- Internal faces
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gf=GFactorCF(k)
          cpf=gf*SpecificHeat(i)+(1.-gf)*SpecificHeat(j)
c
          effdiv(i)=effdiv(i)+mdot(k)*cpf
          effdiv(j)=effdiv(j)-mdot(k)*cpf
c
        enddo
c
c--- Boundary faces
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            k=NBFaceOwner(i,j)           
c
            cpf=BSpecificHeat(i,j)
c  
            effdiv(k)=effdiv(k)+Bmdot(i,j)*cpf
c
          enddo
        enddo
c
      else
c
c--- Internal faces
c
        do k=1,NIFaces
c        
          i=NIFaceOwner(k)
          j=NIFaceNeighbor(k)
c
          gf=GFactorCF(k)
          cpf=gf*SpecificHeatScalar(i,iScalarVariable)+
     *            (1.-gf)*SpecificHeatScalar(j,iScalarVariable)
c
          effdiv(i)=effdiv(i)+mdot(k)*cpf
          effdiv(j)=effdiv(j)-mdot(k)*cpf
c
        enddo
c
c--- Boundary faces
c
        do i=1,NumberOfBCSets
          do j=1,NBFaces(i)
c   
            k=NBFaceOwner(i,j)           
c
            cpf=BSpecificHeatScalar(i,j,iScalarVariable)
c  
            effdiv(k)=effdiv(k)+Bmdot(i,j)*cpf
c
          enddo
        enddo
c
      endif
c
      return
      end