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

module chordlift
  use turbtypes
  implicit none

contains

  ! --------------------------------------------------------------------------
  ! Calculates twist and chord length at a particular node for BEM model

  subroutine calculateBemChordTwist(chord, twist, turb, r)
    real chord, twist, r
    type (ModelTurbine) turb

    real dr, rStep, rq
    integer i
    type(BladeElement) :: be1, be2

    chord = 0
    twist = 0

    ! Linear interpolation
    do i=1, size(turb%bladeElements)-1

       be1 = turb%bladeElements(i)
       be2 = turb%bladeElements(i+1)

       if( be1%radius .le. r .and. r .lt. be2%radius ) then
          dr = r - be1%radius
          rStep = be2%radius - be1%radius
          rq = (dr/rStep)

          chord = (1.-rq)*be1%chord + rq*be2%chord
          twist = (1.-rq)*be1%twist + rq*be2%twist
       end if

    end do

  end subroutine calculateBemChordTwist



  ! --------------------------------------------------------------------------
  ! Calculates lift and drag coefficients for a given turbine and attack angle
  ! (Reynolds number current unused)

  subroutine calculateBemLiftDragCoeffs(lift, drag, turb, attack, Re)
    real :: lift, drag, attack, Re
    type (ModelTurbine) :: turb

    real :: dTheta, thetaStep, thetaQ
    integer :: i

    ! if( turbineDebug ) print*,"@ ... calcBemLiftDragCoeffs()"
    ! if( turbineDebug ) print*,"@ bladeCoeffs_sz:", size(turb%bladeCoeffs)

    lift = 0.
    drag = 0.

    ! Linear interpolation
    do i=1, size(turb%bladeCoeffs)-1

       if( turb%bladeCoeffs(i)%angle .le. attack &
            .and. attack .lt. turb%bladeCoeffs(i+1)%angle ) then

          dTheta = attack - turb%bladeCoeffs(i)%angle
          thetaStep = turb%bladeCoeffs(i+1)%angle - turb%bladeCoeffs(i)%angle
          thetaQ = (dTheta/thetaStep)

          lift = (1-thetaQ)*turb%bladeCoeffs(i)%lift + thetaQ*turb%bladeCoeffs(i+1)%lift
          drag = (1-thetaQ)*turb%bladeCoeffs(i)%drag + thetaQ*turb%bladeCoeffs(i+1)%drag
       end if

    end do


  end subroutine calculateBemLiftDragCoeffs



  ! --------------------------------------------------------------------------


  subroutine calculateAttackMaxLiftPerf( turb )
    type(ModelTurbine) :: turb

    integer :: i
    real :: attackMax, optLift, diff, diffMax


    attackMax=0
    optLift=0.
    diffMax=0

    do i=1, size(turb%bladeCoeffs)
       if(turb%traditionalAttackOpt .eqv. .false.) then
            diff = turb%bladeCoeffs(i)%lift - turb%bladeCoeffs(i)%drag
       else
            if( abs(turb%bladeCoeffs(i)%drag) >  verySmall) then
                diff = turb%bladeCoeffs(i)%lift/turb%bladeCoeffs(i)%drag
            else
                diff = veryBig
            end if
       end if

       if(diff > diffMax) then
          attackMax = turb%bladeCoeffs(i)%angle
          optLift =  turb%bladeCoeffs(i)%lift
          diffMax = diff
       end if
       
    end do

    turb%attackMaxPerf = attackMax
    turb%optLift = optLift

    print*, "@ attackMaxPerf:", attackMax

  end subroutine calculateAttackMaxLiftPerf


  ! --------------------------------------------------------------------------

end module chordlift
