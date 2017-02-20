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

! --------------------------------------------------------------------------

module calcglobals
  use turbtypes
  use parallel
  use mathsmisc
  use interpnodes
  use chordlift

  implicit none

  private
  public :: calculateTurbineGlobals, calculateBemNodeAverages

contains


  ! --------------------------------------------------------------------------


  subroutine calculateTurbineGlobals( turb, t, dt )
    implicit none

    integer :: notFirstFlag

    type (ModelTurbine) :: turb
    type (turbineNode) :: thisNode

    real :: t, dt
    real :: turbRelax

    real :: omegaFLMean, omegaFLSqMean, omegaFree, omegaTurb, omegaSqMean
    real :: BEffective
    real :: halfInertiaMeanU, umax, uMean, vMean, wMean
    real :: uHubMean, uSqMean, vSqMean, wSqMean, deltaOmegaSqMean
    real :: relFlowAngle, absFlowDir, turnAngle, maxTurnAngle
    real :: kappaOmegaT, newPower, densitySum
    real :: calculateMeanDeltaOmegaSq

    ! BEM variables
    real :: inertia, dotOmegaT, newOmegaTurb, newBladeFirstAngle


    ! Can't use nice node structures since MPI won't support it
    integer :: i

    logical, parameter :: adamsBashforthMethod = .true.

    !    real calculateEffectiveNetSolidity

    turbRelax = turb%turbineRelax

    if( turbineDebug ) print*, "--- calculateTurbineGlobals()"


    ! Fluidity doesn't update it's halo values properly after an adapt.
    ! So, we have to make a guess.

!    if( fluidityStabilityHack .and. isParallelRun() ) then
!       call correctSuspiciousFluidityValues( turb, dt )
!    end if


    ! Calculate various means of turbine
    if( turb%numberNodes.gt.0 ) then

       ! mean element dx (assuming elements all same size)
       turb%meanElementDx = ( &
            pi * (turb%radius**2.0) * turb%length &
            / turb%numberNodes )**(1/3)

       call calculateTurbineVelocityMeans( &
            uMean, vMean, wMean, uHubMean, &
            uSqMean, vSqMean, wSqMean, &
            omegaFLMean, &
            omegaFLSqMean, &
            omegaFree, &
            deltaOmegaSqMean, &
            turb )


       call calculateInterpolatedNodes( turb )

       if( turbineDebug ) print*,"@uMean:", uMean
       if( turbineDebug ) print*,"@deltaOmegaSqMean:", deltaOmegaSqMean


       turb%instantUMean = uMean
       turb%instantVMean = vMean
       turb%instantWMean = wMean

       turb%instantUSqMean = uSqMean
       turb%instantVSqMean = vSqMean
       turb%instantWSqMean = wSqMean

       turb%uHubMean = uHubMean * (1-turbRelax) + turb%oldUHubMean * turbRelax
       turb%oldUHubMean = turb%uHubMean


       turb%uMean = uMean * (1-turbRelax) + turb%oldUMean * turbRelax
       turb%oldUMean = turb%uMean

       turb%vMean = vMean * (1-turbRelax) + turb%oldVMean * turbRelax
       turb%oldVMean = turb%vMean

       turb%wMean = wMean * (1-turbRelax) + turb%oldWMean * turbRelax
       turb%oldWMean = turb%wMean


       turb%uSqMean = uSqMean * (1-turbRelax) + turbRelax*turb%oldUSqMean
       turb%oldUSqMean = turb%uSqMean

       turb%omegaFLSqMean = (1-turbRelax)*omegaFLSqMean &
            + turbRelax * turb%oldOmegaFLSqMean
       turb%oldOmegaFLSqMean = turb%omegaFLSqMean

       if( turbineDebug ) then
          print*, "@ uMean, vMean:", turb%uMean, turb%vMean
          print*, "@ relFlowAngle:", radToDeg * relFlowAngle
       end if


       ! Now to the turbine orientation. The turbine can rotate towards 
       ! the flow. uMean and vMean are relaxed values; no further
       ! relaxation necessary.

       relFlowAngle = atan2(turb%vMean, turb%uMean)

       ! The ideal turn angle is the orientation of the flow relative
       ! to the turbine
       turnAngle = relFlowAngle

       ! Diagnostics: give the absolute flow direction 
       ! (not relative to the turbine)
       absFlowDir = relFlowAngle + turb%orientation
       turb%flowDir = absFlowDir
       turb%oldFlowDir = turb%flowDir ! Redundant - get rid of!

       ! We're limited by the speed at which the turbine can rotate.
       ! if orientRevTime = 0, this fixes the turbine at the initial
       ! angle.

       if(turb%orientationFixed .eqv. .false. ) then
           if( abs(turb%orientRevTime) .lt. verySmall ) then
              maxTurnAngle = 0.
           else
              maxTurnAngle = (2*pi/turb%orientRevTime) * dt
           end if

           if( abs(turnAngle) > maxTurnAngle ) then
              turnAngle = sign( maxTurnAngle, turnAngle )
           end if

           turb%orientation = &
                mod(turb%orientation+turnAngle, 2.*pi)
       end if


       ! maximum velocity in turbine nodes
       uMax=0
       do i=1, turb%numberNodes
          if( abs(turb%nodes(i)%u) .gt. uMax) then
             uMax = abs(turb%nodes(i)%u)
          end if
       end do

       turb%instantUMax = uMax

       turb%uMax = uMax * (1-turbRelax) &
            + turb%oldUMax * turbRelax

       turb%oldUMax = turb%uMax

       ! ----------------------------------------
       ! Calculations based on model type

       select case (turb%status)

       case (BEM_DISC_MODEL, BEM_LINE_MODEL)

          turb%fluidTorque = (1.-0.25*turbRelax) &
               * turb%instFluidTorque &
               + 0.25*turbRelax * turb%oldFluidTorque
          turb%oldFluidTorque = turb%fluidTorque

          ! power torque is in same direction as blade torque
          turb%powerTorque = turb%k_Tau &
               *  sign(turb%omegaTurb**2., turb%omegaTurb)

          turb%bladesTorque = -turb%fluidTorque &
               - turb%powerTorque


          dotOmegaT = turb%bladesTorque / turb%inertia

          ! Two-step Adams Bashforth integration
          if( adamsBashforthMethod ) then

             ! First order Euler needed to kick start solution
             if( (t-turb%windRampTime) .le. dt) then

                newOmegaTurb = turb%omegaTurb + dotOmegaT * dt

                if(turb%status == BEM_LINE_MODEL) then
                    newBladeFirstAngle=turb%bladeFirstAngle+turb%omegaTurb*dt
                end if

             else
                ! Otherwise...
                newOmegaTurb = adamsBashforthIntegrate( &
                     turb%omegaTurb, &
                     turb%omegaTurbDot, &
                     turb%oldOmegaTurbDot, &
                     dt)

                if(turb%status == BEM_LINE_MODEL) then
                    newBladeFirstAngle = adamsBashforthIntegrate( &
                         turb%bladeFirstAngle, &
                         turb%omegaTurb, &
                         turb%oldOmegaTurb, &
                         dt)
                end if
             end if

          else
             ! If Adams-Bashforth not selected, 1st order Euler integration
             newOmegaTurb = turb%omegaTurb + dotOmegaT * dt
             if(turb%status == BEM_LINE_MODEL) then
                 newBladeFirstAngle=turb%bladeFirstAngle+turb%omegaTurb*dt
             end if

          end if

          if(turb%status == BEM_LINE_MODEL) then
            if (newBladeFirstAngle >= 0 ) then
              newBladeFirstAngle = modulo(newBladeFirstAngle, 2*pi)
            else
              newBladeFirstAngle=2*pi &
                    + modulo(newBladeFirstAngle, -(2*pi))
             end if
          end if

          ! mean angular velocity of fluid in turbine
          turb%omegaFLMean = (1-turbRelax)*omegaFLMean &
               + turbRelax * turb%oldOmegaFLMean
          turb%oldOmegaFLMean = turb%omegaFLMean

          ! BEM relaxed - for testing purposes
          !            turb%omegaTurb = (1-turbRelax) * newOmegaTurb &
          !                  + turbRelax * turb%oldOmegaTurb


          turb%oldOmegaTurb = turb%omegaTurb
          turb%omegaTurb = newOmegaTurb

          if(turb%status == BEM_LINE_MODEL) then
             turb%bladeFirstAngle = newBladeFirstAngle
          end if

          ! Milarky for Runge-Kutta
          turb%oldOmegaTurbDot = turb%omegaTurbDot
          turb%omegaTurbDot = dotOmegaT

          turb%power = turb%powerEfficiency &
               * overPredictFactor*abs(turb%k_Tau  * turb%omegaTurb**3.)

       end select


    end if

  end subroutine calculateTurbineGlobals



  ! --------------------------------------------------------------------------


  subroutine calculateTurbineVelocityMeans( &
       uMean, vMean, wMean, uHubMean, &
       uSqMean, vSqMean, wSqMean, &
       omegaFLMean, &
       omegaFLSqMean, &
       omegaFreeMean, &
       deltaOmegaSqMean, &
       turb )

    implicit none

    real uMean, vMean, wMean, uHubMean, uSqMean, vSqMean, wSqMean
    real omegaFL, omegaFreeMean, deltaOmegaSqMean
    type (ModelTurbine) :: turb
    integer i, nNodes

    real r, localBeta
    real uSum, vSum, wSum, uHubSum, uSqSum, vSqSum, wSqSum
    integer nBladeNodes, nHubNodes
    real omegaFreeSum, omegaFLSum, omegaFLSqSum, deltaOmegaSqSum
    real omegaFree, omegaFLMean, omegaFLSqMean, deltaOmegaSq

    type (turbineNode) node
    !    real calculateEffectiveLocalSolidity

    print*,"--- calculateTurbineVelocityMeans()"

    nNodes = turb%numberNodes
    uSum = 0
    vSum = 0
    wSum = 0
    uHubSum = 0
    uSqSum = 0
    vSqSum = 0
    wSqSum = 0

    nHubNodes = 0
    nBladeNodes = 0

    omegaFLSum = 0
    omegaFLSqSum = 0
    omegaFreeSum = 0
    deltaOmegaSqSum = 0

    do i=1, nNodes
       node = turb%nodes(i)

       r = sqrt ( node%y**2.0 + node%z**2.0 )

       !        print*, "uValues:", uValues(i)
       uSum = uSum + node%u
       vSum = vSum + node%v
       wSum = wSum + node%w

       uSqSum = uSqSum + (node%u**2.)
       vSqSum = vSqSum + (node%v**2.)
       wSqSum = wSqSum + (node%w**2.)

       if( r .lt. turb%hubRadius ) then
          uHubSum = uHubSum + node%u

          nHubNodes = nHubNodes + 1
       else

          if( turb%status==BEM_DISC_MODEL .or. &
                turb%status==BEM_LINE_MODEL) then
             omegaFLSum = omegaFLSum &
                  + (node%y * node%w - node%z * node%v) / r**2.
             omegaFLSqSum = omegaFLSqSum &
                  + ((node%y * node%w - node%z * node%v) / r**2.)**2.
          end if

          nBladeNodes = nBladeNodes + 1
       end if


    end do

    uMean = uSum / nNodes
    vMean = vSum / nNodes
    wMean = wSum / nNodes

    uSqMean = uSqSum / nNodes
    vSqMean = vSqSum / nNodes
    wSqMean = wSqSum / nNodes

    if (nHubNodes .gt. 0 ) then
       uHubMean = uHubSum / nHubNodes
    else
       uHubMean = 0.
    end if

    if (nBladeNodes .gt. 0 .and. &
        ( turb%status==BEM_DISC_MODEL .or. turb%status==BEM_LINE_MODEL) ) then
       omegaFLMean = omegaFLSum / nBladeNodes
       omegaFLSqMean = omegaFLSqSum / nBladeNodes
    end if


  end subroutine calculateTurbineVelocityMeans



  ! --------------------------------------------------------------------------
  ! What's the local average contribution to the torque,
  ! and average angle of attack

  subroutine calculateBemNodeAverages( turb, dt )
    implicit none

    type(ModelTurbine) turb
    real :: dt

    integer :: numNodes, numBladeNodes, i
    real :: r, fluidTorqueSum, forceSum
    real :: volume
    real :: chord, twist
    real :: wAttack, attackSum, weight, weightsSum, wAttackSum, eta
    real :: relax

    forceSum = 0.
    fluidTorqueSum = 0.
    attackSum = 0.
    weightsSum = 0.
    wAttackSum = 0.

    relax = turb%turbineRelax

    numBladeNodes = 0
    do i=1, turb%numberNodes
       r = sqrt(turb%nodes(i)%y**2. + turb%nodes(i)%z**2.)

       ! Angle of attack, thrust etc. only valid for blades section
       if(r > turb%hubRadius ) then
           attackSum = attackSum + turb%nodes(i)%attack
		       call calculateBemChordTwist(chord, twist, turb, r)

		       ! Weight = chord * uRel^2 (why? measure of contribution to lift)
		       select case(turb%status)
                    case(BEM_LINE_MODEL)
                        weight = chord * turb%nodes(i)%uRel**2. &
                            * turb%nodes(i)%bladeDistrib

                    case(BEM_DISC_MODEL)
                        weight = chord * turb%nodes(i)%uRel**2.

		       end select

		       weightsSum = weightsSum + weight
		       wAttackSum = wAttackSum + turb%nodes(i)%attack * weight

		       forceSum = forceSum + turb%nodes(i)%bodyThrust
		       fluidTorqueSum = fluidTorqueSum + turb%nodes(i)%fluidBodyTorque

		       numBladeNodes = numBladeNodes+1
		   end if

    end do

    volume = pi* (turb%radius**2.0-turb%hubRadius**2.0) * turb%length
    if( numBladeNodes > 0 ) then

       turb%thrust = - overPredictFactor * (forceSum / numBladeNodes) * volume

       turb%bladeLoading = turb%thrust / turb%bladesArea
       turb%instFluidTorque = (fluidTorqueSum / numBladeNodes) * volume
       turb%attackAv = attackSum / numBladeNodes

       if ( weightsSum > verySmall .and. weightsSum < veryBig ) then
            ! print*, "attURelSqSum:", attackURelSqSum
            ! print*, "uRelSqlSum:", uRelSqSum
            wAttack = wAttackSum / weightsSum
            turb%wAttack = wAttack * (1.-relax) + relax * turb%oldWAttack
            turb%wAttackDot = (turb%wAttack - turb%oldWAttack)/dt
            turb%oldWAttack = turb%wAttack
       end if
    end if


  end subroutine calculateBemNodeAverages


  ! --------------------------------------------------------------------------


end module calcglobals
