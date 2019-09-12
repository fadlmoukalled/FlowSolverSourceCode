c
c#############################################################################################
c
      SUBROUTINE calculatesigRNG
c
C#############################################################################################
c
      use Geometry1, only: NumberOfElements,NumberOfBCSets
      use Geometry3, only: NBFaces
      use Turbulence1, only: sigTKERNG,sigTEDRNG,BsigTKERNG,BsigTEDRNG
      use PhysicalProperties1, only: Viscosity,TurbulentViscosity,
     *                               BViscosity,BTurbulentViscosity
c
c********************************************************************************************
      implicit none
c********************************************************************************************
      integer :: i,j,iter
      integer :: itermax=100
      double precision :: epsilon=1.e-6
      double precision :: a1,a2,am,rhs,diff,term1,term2
c********************************************************************************************
c
      do i=1,NumberOfElements
c
        a1=1.392899999999
        a2=2.5
        rhs=Viscosity(i)/(TurbulentViscosity(i)+Viscosity(i))
        diff=1.
        iter=1 
c
        do while (dabs(diff).gt.epsilon.and.iter.lt.itermax)
c
          am=0.5*(a1+a2)
          term1=dabs((am-1.3929)/(-0.3929))**0.6321
          term2=dabs((am+2.3929)/(3.3929))**0.3679
c
          diff=term1*term2-rhs
c        
          if(diff.lt.0.) then
            a1=am
          else
            a2=am
          endif       
c
          iter=iter+1        
c
        enddo
c
        if(iter.ge.100) then
c
          print*,'maximum iterations reached'
          print*,'residuals= ',dabs(diff)
c
        endif
c        
        sigTKERNG(i)=am
        sigTEDRNG(i)=am
c
      enddo
c
      do i=1,NumberOfBCSets
        do j=1,NBFaces(i)
c
          a1=1.392899999999
          a2=2.5
          rhs=BViscosity(i,j)/(BTurbulentViscosity(i,j)+BViscosity(i,j))
          diff=1.
          iter=1 
c
          do while (dabs(diff).gt.epsilon.and.iter.lt.itermax)
c
            am=0.5*(a1+a2)
            term1=dabs((am-1.3929)/(-0.3929))**0.6321
            term2=dabs((am+2.3929)/(3.3929))**0.3679
c
            diff=term1*term2-rhs
c        
            if(diff.lt.0.) then
              a1=am
            else
              a2=am
            endif       
c
            iter=iter+1        
c
          enddo
c
          if(iter.ge.100) then
c
            print*,'maximum iterations reached'
            print*,'residuals= ',dabs(diff)
c
          endif
c        
          BsigTKERNG(i,j)=am
          BsigTEDRNG(i,j)=am
c
        enddo
      enddo
c
      return
      end