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

module bladepitch
  use turbtypes
  implicit none

contains

  ! --------------------------------------------------------------------------
  ! correctTurbineBladePitch() :
  ! For thesis model -
  !     This corrects the blade pitch to make sure that omegaTurbOpt
  !     is not exceeded.
  ! For BEM - 
  !     Currently slowly turns from initial angle to target angle
  !     (needed for high performance aerofoils)

  subroutine correctTurbineBladePitch( turb, t, dt )
    implicit none

    type(ModelTurbine) :: turb
    real :: t, dt
    real :: powerCurrentMax
    real :: u0_Current, deltaPower, deltaU0, powerGrad
    real :: k, omegaTurbMax

    real :: newAlpha, deltaAlpha
    real :: l, m
    real :: deltaOmegaT
    real :: omegaTDotTarget, omegaTDotMax
    real :: alphaDotTarget, alphaDotMax
    real :: attackTarg, attackDotMax, attackDotTarg
    real :: normOmegaDiff


    ! For each active turbine...
    select case ( turb%status )
       
    ! Blade element momentum model.
    case( BEM_DISC_MODEL, BEM_LINE_MODEL )

        if( abs(turb%wAttack) < veryBig ) then
		    alphaDotMax = 2*pi / turb%alphaRelaxTime
		    attackDotMax = -alphaDotMax

            k = -turb%maxPower / (pi/2. - turb%attackMaxPerf)

            if( turb%power < turb%maxPower ) then
                attackTarg = turb%attackMaxPerf
            else
                attackTarg = 0.
            end if

		    attackDotTarg = (attackTarg - turb%wAttack) / turb%alphaRelaxTime
		    m = (attackDotTarg - turb%wAttackDot) / attackDotMax

            deltaAlpha = dt * ( ( turb%alphaDot + m * alphaDotMax )/2. )

            ! Limit deltaAlpha
            if( abs(deltaAlpha) > abs(alphaDotMax*dt) ) then
               deltaAlpha=alphaDotMax*dt * (deltaAlpha/abs(deltaAlpha))
            end if

            newAlpha = turb%alpha + deltaAlpha
            
            ! Catch-all for bad pitch values
            if( newAlpha < -pi/2. ) then
                newAlpha = -pi/2.
            elseif( newAlpha > pi/2. ) then
                newAlpha = pi/2.
            end if

		    turb%alpha = newAlpha
		    turb%alphaDot = (turb%alpha - turb%oldAlpha) / dt
		    turb%oldAlpha = turb%alpha

        end if
       
       ! Nothing sensible to do here
    case default
       print*, "<!> Awooga awooga! Unidentified model type!"
       
    end select


  end subroutine correctTurbineBladePitch


  ! --------------------------------------------------------------------------


end module bladepitch
