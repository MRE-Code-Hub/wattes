! WATTES - Wind And Tidal Turbine Embedded Simulator
! An actuator disc/line turbine model for coupling with CFD software
!
! Copyright (C) 2017 Heriot Watt University and the University of Edinburgh.
!
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
!
!	  Dr. A Creech
!	  Institute of Energy Systems
!	  School of Engineering
!	  University of Edinburgh
!	  
!	  angus_creech@hotmail.com
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
! -----------------------------------------------------------------------------


#include <turbdefs.h>
  
subroutine turbineFarmFluidityInterface(numXNodes, numUNodes, &
     xNode, yNode, zNode, &
     uNode, vNode, wNode, &
     referenceDensity, &
     xSource, ySource, zSource, &
     t, dt, &
     currentNonLinearIteration, numIterations, &
     remesh, bladeDistrib )

  use main
  implicit none

  logical remesh

  integer :: numXNodes, numUNodes, &
       currentNonLinearIteration, numIterations
  real, dimension(numUNodes), intent(in) :: xNode, yNode, zNode
  real, dimension(numUNodes), intent(in) :: uNode, vNode, wNode

  ! Currently we only deal with prescribed, constant density fluids.
  ! For Boussinesq, we'll need the density field instead.
  real :: referenceDensity
  real, dimension(numUNodes) :: xSource, ySource, zSource
  real, dimension(numUNodes), optional :: bladeDistrib
  real t, dt

  integer i
  logical updateTurbine, outputTurbineData


  
  print*, "turbineFarmFluidityInterface()"
  print*, "xnonod:", numXNodes, "  nonods:", numUNodes
  print*, "iteration, nonlin:", currentNonLinearIteration, numIterations

  do i=1, numUNodes
!     write(*,"(1I5, 3e5.3)"),i, xNode(i), yNode(i), zNode(i)
  end do


  ! WARNING: DIRTY HACK AHEAD (THIS IS OUTDATED)
  !
  ! The following conditional statements get around some odd
  ! behaviour on the part of Fluidity:
  !
  ! 1) Fluidity calls the momentum source function for each nonlinear
  !     iteration -- which could be several times each timestep.
  !
  ! 2) Fluidity does not conserve mesh values after an adapt, but uses
  !     interpolation to 'hope things work out'. This causes stability
  !     problems with parallel (MPI) simulations, since velocity values
  !     at nodes within the halo may be suspect.
  !
  ! THE WORK-AROUND:
  !
  ! 1) If Fluidity is in last iteration in a timestep and not right after an adapt sweep,
  !     choose to recalculate the turbine values. Also output the turbine values to
  !      the screen/file.
  !
  ! 2) Fluidity does not conserve mesh values after an adapt, but uses
  !     interpolation to 'hope things work out'. This causes stability
  !     problems with parallel (MPI) simulations, since velocity values
  !     at nodes within the halo may be suspect.
  !
  !
  ! FOOTNOTE
  !
  ! It is doubtful this level of complexity would be needed for other 
  ! CFD programs, which is why it is included in the Fluidity interface 
  ! rather than the solver itself.
  !


  if ( currentNonLinearIteration == 1 ) then
     updateTurbine = .true.
     outputTurbineData = .true.
  else
     updateTurbine = .false.
     outputTurbineData = .false.
  end if


  ! print*,"@ currentNonLinearIteration:", currentNonLinearIteration
  ! print*,"@ numIterations:", numIterations
  ! print*,"@ remesh:", remesh
  ! print*,"@ updateTurbineFlag:", updateTurbine
  ! print*,"@ outputTurbineData:", outputTurbineData


  ! Notice dirty hack for incompressible fluids. Only
  ! passing first node's density.

  if(.not. present(bladeDistrib)) then
      call turbineFarmSolver( numUNodes, &
           xNode, yNode, zNode, &
           uNode, vNode, wNode, &
           referenceDensity,  &
           xSource, ySource, zSource, &
           t, dt, &
           updateTurbine, outputTurbineData )
  else
      call turbineFarmSolver( numUNodes, &
           xNode, yNode, zNode, &
           uNode, vNode, wNode, &
           referenceDensity,  &
           xSource, ySource, zSource, &
           t, dt, &
           updateTurbine, outputTurbineData, bladeDistrib )
  end if

  ! print*,"End turbineFarmFluidityInterface"

end subroutine turbineFarmFluidityInterface


!------------------------------------------------------------------------------


! Turbine Farm OpenFOAM Interface
!
! Unlike Fluidity, we can't pass bladeDistrib currently. C++ confuses the
! Fortran interface


subroutine turbineFarmOpenFOAMInterface(numXNodes, numUNodes, &
     xNode, yNode, zNode, &
     uNode, vNode, wNode, &
     referenceDensity, &
     xSource, ySource, zSource, &
     t, dt, &
     currentNonLinearIteration, numIterations, &
     remesh )

  use main
  implicit none

  logical remesh

  integer :: numXNodes, numUNodes, &
       currentNonLinearIteration, numIterations
  real, dimension(numUNodes), intent(in) :: xNode, yNode, zNode
  real, dimension(numUNodes), intent(in) :: uNode, vNode, wNode

  ! Currently we only deal with prescribed, constant density fluids.
  ! For Boussinesq, we'll need the density field instead.
  real :: referenceDensity
  real, dimension(numUNodes) :: xSource, ySource, zSource
  real t, dt

  integer i
  logical updateTurbine, outputTurbineData

  !currentNonLinearIteration = 1
  !numIterations = 2
  print*, ""
  print*, ""
  print*, "turbineFarmOpenFOAMInterface()"
  print*, "DEBUG EDITION"
  print*, "xnonod:", numXNodes, "  nonods:", numUNodes
  print*, "iteration, nonlin:", currentNonLinearIteration, numIterations

  if ( currentNonLinearIteration == 1 ) then
     updateTurbine = .true.
     outputTurbineData = .true.
  else
     updateTurbine = .false.
     outputTurbineData = .false.
  end if

  

  print*, "X"
  do i=1, numXNodes
     print 755, i, xNode(i), yNode(i), zNode(i)
755  format(I5, 3(F8.3, " "))
  end do

  print*, ""
  print*, "U"
  do i=1, numXNodes
     print 756, i, uNode(i), vNode(i), wNode(i)
756  format(I5, 3(F8.3, " "))
  end do

  call turbineFarmSolver( numUNodes, &
       xNode, yNode, zNode, &
       uNode, vNode, wNode, &
       referenceDensity,  &
       xSource, ySource, zSource, &
       t, dt, &
       updateTurbine, outputTurbineData )

  print*,"End turbineFarmOpenFOAMInterface"
  print*, ""

end subroutine turbineFarmOpenFOAMInterface
