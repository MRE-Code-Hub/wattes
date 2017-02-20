! WATTES - Wind And Tidal Turbine Embedded Simulator
! An actuator disc/line turbine model for coupling with CFD software
!
! Copyright (C) 2017 Heriot Watt University and the University of Edinburgh.
!
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
!
!         Dr. A Creech
!         Institute of Energy Systems
!         School of Engineering
!         University of Edinburgh
!
!         angus_creech@hotmail.com
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
! --------------------------------------------------------------------------

#include <turbdefs.h>

module mathsmisc
  use turbtypes
  implicit none

contains

  ! --------------------------------------------------------------------------
  ! Adams-Bashforth two-step technique

  real function adamsBashforthIntegrate( y0, ydot1, ydot0, h )
    implicit none
    real y0, ydot1, ydot0, h
    real y1

    y1 = y0 + 1.5*h*ydot1 - 0.5*h*ydot0

    adamsBashforthIntegrate = y1

  end function adamsBashforthIntegrate



  ! --------------------------------------------------------------------------
  ! Quick and dirty routine used in several places

  subroutine horizRotate( componentx, componenty, componentz, rad)
    implicit none
    real componentx,  componenty,  componentz
    real rad

    real tempx, tempy

    tempx = componentx*cos(rad) - componenty*sin(rad)
    tempy = componentx*sin(rad) + componenty*cos(rad)

    componentx=tempx
    componenty=tempy

  end subroutine horizRotate


  ! --------------------------------------------------------------------------
  ! Another quick and dirty routine used in several places

  subroutine tiltRotate( componentx, componenty, componentz, tiltRad)
    implicit none
    real componentx,  componenty,  componentz
    real tiltRad

    real tempx, tempz

    tempx = componentx*cos(tiltRad) - componentz*sin(tiltRad)
    tempz = componentx*sin(tiltRad) + componentz*cos(tiltRad)

    componentx=tempx
    componentz=tempz

  end subroutine tiltRotate

  ! --------------------------------------------------------------------------
  ! Rotate vector about specified axis
  !
  ! CURRENTLY UNUSED !!!

  subroutine rotateVector( vec, rad, axis )
    implicit none

    type(Vector) vec
    real rad

    real tempx, tempy, tempz
    integer axis

    select case (axis)
       ! clockwise
    case (X_AXIS)
       tempx=vec%x
       tempy=vec%y*cos(rad) - vec%z*sin(rad)
       tempz=vec%y*sin(rad) + vec%z*cos(rad)

       ! clockwise
    case (Y_AXIS)
       tempx=vec%x*cos(rad) + vec%z*sin(rad)
       tempy=vec%y
       tempz=( -vec%x*sin(rad) ) + vec%z*cos(rad)

       ! anti-clockwise
    case (Z_AXIS)
       tempx=vec%x*cos(rad) - vec%y*sin(rad)
       tempy=vec%x*sin(rad) + vec%y*cos(rad)
       tempz=vec%z

    end select

    vec%x=tempx
    vec%y=tempy
    vec%z=tempz

  end subroutine rotateVector

  ! --------------------------------------------------------------------------

  function calculateU0( turb )
    implicit none

    type(ModelTurbine) :: turb

    real calculateU0

    select case( turb%status )

	    case( BEM_DISC_MODEL, BEM_LINE_MODEL )
	        calculateU0 = (2./3.) * turb%uMax

        case default
            print*, "*** WARNING: can't identify turbine type"
            calculateU0 = turb%uMax

    end select

  end function calculateU0

end module mathsmisc
