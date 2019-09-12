!-------------------------------------------------------------------------------
!    GnuPlot Interface
!-------------------------------------------------------------------------------
!    Purpose:   Object Based Interface to GnuPlot from Fortran (ogpf)
!    Platform:  Windows XP/Vista/7/10
!               (It should work on other platforms, see the Write2GnuPlot subroutine below)
!    Language:  Fortran 2003 and 2008
!    Requires:  1. Fortran 2003 compiler (e.g gfortran 4.7, IVF 12.1, ...) or 2008
!               2. gnuplot 4.5 and higher (other previous version can be used
!    Author:    Mohammad Rahmani
!               Chem Eng Dep., Amirkabir Uni of Tech
!               Tehran, Ir
!               url: aut.ac.ir/m.rahmani
!               email: m[dot]rahmani[at]aut[dot]ac[dot]ir
!   License:    MIT

! This file demonstrate the capability of ogpf module
! An object based Fortran interface to gnuplot
!
! Acknowledgement:
! Special thanks to Hagen Wierstorf (http://www.gnuplotting.org)
! For vluable codes and examples on using gnuplot
! Some examples and color palletes are provided by gnuplotting.
!
!
! Revision 0.22
! Date: Mar 9th, 2018
! - see ogpf.f90 for details
! - more examples to reflect the new features

! Revision:  0.20
! Date:     Feb 20th, 2018
!  - more examples
!  - animation of 2D and 3D plots

! Revision:  0.19
! Date:     Jan 15th, 2018
!  - new contour plot procedure


! Revision:  0.18
! Date:     Dec 22th, 2017
! More example based on ogpf 0.18

! Version:  0.17
! Date:     Dec 18th, 2017
!   Minor corrections
! - Multi window plots using script



! Version:  0.16
! Date:     Feb 11th, 2016
!   Minor corrections
!   Correct the lspec processing in plot2D_matrix_vs_vector
!   Some examples revised!


! Version:  0.15
! Date:     Apr 20th, 2012
!   Minor corrections
!   Use of select_precision module and working precision: wp


! Version:  0.14
! Date:     Mar 28th, 2012
!   Minor corrections
!   use of import keyboard and  removing ogpf precision module


! Version:  0.13
! Date:     Feb 12th, 2012
!   Minor corrections
!   Added more samples
!   Use ogpf precision module



! Version:  0.12
! Date:     Feb 9th, 2012
! New examples for semilogx, semilogy, loglog
! New set options method

! Version:  0.11
! Date:     Feb 9th, 2012

subroutine PlotConvergenceHistory

    use ogpf
    use MultiGrid2, only: nIter
    use Residuals1
    implicit none
    ! parameters
    ! local variables
    integer:: i=0
    type(gpf):: gp
    double precision,save, dimension(:), allocatable :: x,y
    double precision,save, dimension(:), allocatable :: xT,yT
        
    if(.not.allocated(xT)) then
      allocate(x(1))
      allocate(y(1))
      x(1)=1
      y(1)=ResorMax(19)
      allocate(xT(1))
      allocate(yT(1))
      xT=x
      yT=y
    else
      allocate(x(size(xT)+1))
      allocate(y(size(xT)+1))
      x(1:size(xT))=xT(1:size(xT))
      y(1:size(xT))=yT(1:size(xT))
      x(size(xT)+1)=nIter
      y(size(xT)+1)=ResorMax(19)
      deallocate(xT,yT)
      allocate(xT(size(x)))
      allocate(yT(size(y)))
      xT=x
      yT=Y
    endif        
        ! Annotation: set title, xlabel, ylabel
        call gp%title('Example 1. A simple xy plot','#990011')
        !call gp%xlabel('my x axis ...','#99aa33',font_name="Tahoma")
        !call gp%ylabel('my y axis ...')
        call gp%options('set border lc "#99aa33"; set ylabel "my label..." tc "#99aa33"')
  !
        call gp%ylabel("y,logarithmic scale")
        call gp%xlabel("x, normal scale")
 !     call gp%options('set logscale y2')
        call gp%semilogy(x,y)
!        write ( command, * ) 'gnuplot -persist ' // trim ( command_filename ) // ' &'
     deallocate(x,y)
    return
    end subroutine PlotConvergenceHistory
    
    
    
    
    
    
    

