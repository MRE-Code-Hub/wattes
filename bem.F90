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

module bem
  use turbtypes
  implicit none

contains

  ! --------------------------------------------------------------------------
  ! Routine to calculate momenta of inertia of turbine

  function calculateBemTurbineInertia( turb )
    implicit none

    type(ModelTurbine) :: turb
    real calculateBemTurbineInertia

    type(BladeElement) :: el1, el2

    integer i
    real inertia, innerSum

    ! Currently broken

    if( .false. ) then
       ! Calculation of surf. area of each blade segment. Assumes dm(r) linear
       innerSum = 0
       do i=1, size(turb%bladeElements)-1
          el1 = turb%bladeElements(i)
          el2 = turb%bladeElements(i+1)

          innerSum = innerSum + &
               ( el1%chord + el2%chord ) &
               * abs( el2%radius - el1%radius )
       end do

       inertia = 0.5 * turb%massPerUnitArea * turb%numBlades * innerSum

       ! default case just now. need to fix 
    else
       inertia = turb%numBlades * turb%radius*3. * turb%massPerUnitArea
    endif

    calculateBemTurbineInertia = inertia

  end function calculateBemTurbineInertia

end module bem
