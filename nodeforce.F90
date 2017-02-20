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


module nodeforce
  use turbtypes
  use turbulence
  use chordlift
  use interpnodes
  use parallel
  use mathsmisc

  implicit none

contains



  ! --------------------------------------------------------------------------
  ! Calculates force terms for BEM actuator line
  !

  subroutine calculateBemLineNodeForce( node, turbine, updateTurbine )

    type (TurbineNode) :: node
    type (ModelTurbine) :: turbine
    logical :: updateTurbine

    real :: a, r, R_T, mu, hubSolidity
    real :: du, dv, dw
    real :: tDu, tDv, tDw
    real :: bladeOrbVel, fluidOrbVel, relOrbVel, uRel
    real :: mag

    ! orbital aerodynamic force on fluid
    real :: orbComp, axialComp

    real :: N, rho, chord, twist, nodeOrbAngle
    real :: bladeRotAngle, bladePitchAngle, attack
    real :: velAngle, theta
    real :: fakeReynolds, c_L, c_D, liftBody, dragBody, ncf, lambda
    real :: x, y, z,  u, v, w
    real :: rotfactor, rotdr

    real :: tipLoss, expN, expo


    real, parameter :: defaultForceFactor = 1.0
    real, parameter :: hubFactor = 0.5

    real :: inv_dt, forceFactor

    ! Gaussian convolution filter for forces: see Sorensen and Shen
    ! (actuator line theory) - sort of ! Been tweaked
    ! real :: eta
    real :: epsilon

    integer :: i, numInterps, ix, procID

    ! Calc eta up here.
    real :: mux, hubEta, etaSum
    real, allocatable ::  eta(:)
    ! This standard dev gives ~95% of the force within -2*sigma < x < 2*sigma
    real :: sigma

    real :: lift, drag, sense, dTheta, dc, funcDistSq

    allocate( eta(turbine%numBlades) )

    ! Variables to initialise
    if(turbine%sigma <= 0) then
        sigma = turbine%length/4
    else
        sigma = turbine%sigma
    end if

    if(turbine%forceFactor <= 0) then
        forceFactor = defaultForceFactor
    else
        forceFactor = turbine%forceFactor
    end if


    ! -1 for anticlockwise, +1 for clockwise
    sense = 1.0 * turbine%bladeSense

    if( turbineDebug ) print*,"% ------ calcBemNodeAbsSrc()"

    node%source%x = 0.0
    node%source%y = 0.0
    node%source%z = 0.0

    node%absorption%x = 0.0
    node%absorption%y = 0.0
    node%absorption%z = 0.0

    node%fluidBodyTorque = 0.
    node%bodyThrust = 0.0
    axialComp = 0.0

    du = 0.
    dv = 0.
    dw = 0.

    tDu = 0.0
    tDv = 0.0
    tDw = 0.0

    x = node%x
    y = node%y
    z = node%z

    u = node%u
    v = node%v
    w = node%w

    N = 1. * turbine%numBlades
    rho = turbine%densityFluid

    ! Work out perpendicular distance from node to axis of rotation
    r = sqrt( node%y**2.0 + node%z**2.0 )

    ! Calculate weights for forces from each blade:
    ! we have multiple blades in the model now (as opposed to actuator
    ! disc rings)

    hubEta = (1/(sqrt(2*pi)*sigma)) * exp(-0.5* ((node%x/sigma)**2.) )

    etaSum=0

    if( r > turbine%hubRadius) then
        do i=1, turbine%numBlades
            bladeRotAngle = modulo(turbine%bladeFirstAngle &
                            + turbine%bladeRelativeAngle(i), 2*pi)

            ! This bit finds out the orbital difference between the actuator
            ! line and the mesh node
            nodeOrbAngle = modulo(atan2(z, y)+2*pi, 2*pi)
            dTheta = abs( modulo( (nodeOrbAngle - bladeRotAngle) &
            	+pi, 2*pi)-pi )
            dc = dTheta * r

            ! Functional distance for weighting function eta
            funcDistSq = dc**2. + node%x**2.

            ! Weight from each blade for this node
            eta(i) = (1/(2*pi*sigma**2.)) * exp(-0.5*(funcDistSq/(sigma**2.)) )
            etaSum = etaSum + eta(i)

        end do

        node%bladeDistrib = etaSum
    else
        node%bladeDistrib = 0
    end if

    node%bladeDistrib = etaSum

    ! Work out perpendicular distance from node to axis of rotation
    r = sqrt( node%y**2.0 + node%z**2.0 )
    rho = turbine%densityFluid

    ! Check if point inside the general volume (you never know)
    if( r <= turbine%radius ) then


       ! If in general blade volume, not hub, and if mean u +ve...
       ! (do nothing in reverse flow!)
       if( r > turbine%hubRadius .and. turbine%uMean > 0) then

          ! Angle of the blade (not attack angle)
          call calculateBemChordTwist(chord, twist, turbine, r)

          bladePitchAngle = turbine%alpha + twist

          ! Some simplifications (make calculations stable)

          bladeOrbVel = r * turbine%omegaTurb

          if( r > verySmall ) then
             fluidOrbVel = (node%y * node%w - node%z * node%v ) / r
          else
            fluidOrbVel = 0.
          end if

          relOrbVel = bladeOrbVel - fluidOrbVel

          uRel = sqrt(u**2.+relOrbVel**2.)

          if( abs(u) < verySmall &
               .and. abs(relOrbVel) < verySmall ) then
             velAngle = 0.
          else
             velAngle = atan2( sense*relOrbVel, u )
          end if

          attack =   0.5*pi - velAngle - bladePitchAngle

          R_T = turbine%radius
          mu = r / R_T


          ! Only if there is some flow through the turbine, do something

          if( abs(uRel) > verySmall ) then

                 call calculateBemLiftDragCoeffs( c_L, c_D, &
                      turbine, attack, fakeReynolds )

                 ! Tip-speed ratio for wind turbines
                 lambda = (turbine%omegaTurb/turbine%omegaTurbMax) &
                      * turbine%tipSpeedRatioMax

                 ! flow induction factor. Picked for 'optimised operation',
                 ! cf. Wind Energy Handbook.

                 ! Mean axial induction factor
                 a = 0.333333333333333

                 ! NB. Localised axial induction factor can be ~=0.5
                 ! see fig. 3.36, pp89.
                 ! a=0.5

                 ! Calculate tip loss factor

                 expo = - 0.5* N * ((1-mu)/mu) &
                      * sqrt( 1+ (lambda*mu)**2.0 / ((1.0-a)**2.0) )
                 expN = exp( expo )

                 if ( abs(expN) .gt. 1 ) then
                    expN = expN / abs(expN)
                 end if

                 ! tip-loss
                 if (turbine%tipLoss .eqv. .true.) then
                    tiploss = (2/pi) * acos( expN )
                 else
                    tiploss = 1.0
                 end if

                 ! orbital force component
                 orbComp = - sense * (tiploss * 0.5 * rho* uRel**2. &
                      * (c_L * cos(velAngle) - c_D * sin(velAngle) ) &
                      * chord )

                 ! axial component
                 axialComp =  -  tiploss * 0.5 * rho * uRel**2. &
                      * (c_L * sin(velAngle)  + c_D * cos(velAngle) ) &
                      * chord

                 ! Aerodynamic forces on fluid without power extraction
                 du = forceFactor * etaSum * axialComp

                 dv = - forceFactor * etaSum * (node%z/r) * orbComp
                 dw =   forceFactor * etaSum * (node%y/r) * orbComp

                 ! Torque per-unit-volume at node and angle of attack (used later)

                 node%fluidBodyTorque = forceFactor * etaSum * r * orbComp
                 node%uRel = uRel
                 node%attack = attack

              else

                 ! The blade has no relative motion to the fluid. Do nothing.
                 du = 0.
                 dv = 0.
                 dw = 0.

                 node%fluidBodyTorque = 0.
                 node%uRel = 0.
                 node%attack = 0.

              end if
       else

          if ( abs(turbine%uHubMean)  .lt. verySmall) then
             inv_dt = 0.
          else
             ! Removing fluid density for now - does Fluidity require it?
             ! inv_dt = turbine%densityFluid * abs(turbine%uHubMean)
             inv_dt = abs(turbine%uHubMean)
          end if

          ! Hub resistive force
          du = - hubEta * forceFactor * hubFactor * turbine%uHubMean * inv_dt
          ! du = 0.
          dv = 0.
          dw = 0.

          node%fluidBodyTorque = 0.
          node%uRel = 0.
          node%attack = 0.

       end if

       call calculateNodeTurbulence( turbine, node, &
                turbine%uMean, r, tDu,tDv,tDw )

       ! Spit out source values (notice the turbulence terms)
       ! Divide out fluid density, as Fluidity wants velocity source
       node%source%x = (du + etaSum*tDu) / rho
       node%source%y = (dv + etaSum*tDv) / rho
       node%source%z = (dw + etaSum*tDw) / rho

       ! thrust / unit vol for this node
       ! node%bodyThrust = eta * axialComp

       node%bodyThrust = du

       ! Rotate source and absorption terms back into original reference frame

       call tiltRotate(node%source%x, node%source%y, &
            node%source%z, turbine%upwardTilt)

       call horizRotate(node%source%x, node%source%y, &
            node%source%z, turbine%orientation)

       call tiltRotate(node%absorption%x, node%absorption%y, &
            node%absorption%z, turbine%upwardTilt)

       call horizRotate(node%absorption%x, node%absorption%y, &
            node%absorption%z, turbine%orientation)

    end if

    if( turbineDebug ) print*,"% ------ end calcBemNodeAbs()"

    deallocate(eta)

  end subroutine calculateBemLineNodeForce


  ! --------------------------------------------------------------------------
  ! Calculates force terms for BEM disc
  !

  subroutine calculateBemDiscNodeForce( node, turbine, updateTurbine )

    type (TurbineNode) :: node
    type (ModelTurbine) :: turbine
    logical :: updateTurbine

    real :: a, r, R_T, mu, hubSolidity
    real :: du, dv, dw
    real :: tDu, tDv, tDw
    real :: bladeOrbVel, fluidOrbVel, relOrbVel, magURel
    real :: mag

    ! orbital aerodynamic force on fluid
    real :: aF_theta, axialComp

    real :: N, dr, rho, chord, twist, bladeAngle, velAngle, attack
    real :: theta
    real :: fakeReynolds, c_L, c_D, lift, drag, ncf, ls, lambda
    real :: expo, expN, tipLoss

    real :: uMean, uSqMean
    real :: uInt, vInt, wInt
    real :: usqInt, vsqInt, wsqInt

    real :: rotfactor, rotdr

    ! Accounts for loss due to use of global turbine properties being used,
    ! not time-averaged spatially varying velocities
    ! real, parameter :: forceFactor = 0.3333333333333
    ! real, parameter :: forceFactor = 0.1
    ! real, parameter :: forceFactor = 1.0
    ! real, parameter :: forceFactor = 0.9
    real, parameter :: defaultForceFactor = 0.75
    ! real, parameter :: forceFactor = 0.6666666666666
    ! real, parameter :: forceFactor = 0.5
    ! real, parameter :: forceFactor = 0.4
    ! real, parameter :: forceFactor = 0.3333333333333
    ! real, parameter :: forceFactor = 0.05

    ! real, parameter :: rotateScale = 0.5

    ! real, parameter :: hubFactor = 0.25
    real, parameter :: hubFactor = 0.5
    ! real, parameter :: hubFactor = 0.75
    ! real, parameter :: hubFactor = 1.0

    ! real, parameter :: axialInductionFactor=1.0
    real, parameter :: axialInductionFactor=4./3.
    ! real, parameter :: axialInductionFactor=5./3.
    ! real, parameter :: torqueFactor = 1.0
    real, parameter :: torqueFactor = 1.00

    real :: inv_dt, forceFactor

    ! Gaussian convolution filter for forces: see Sorensen and Shen 
    ! (actuator line theory) - sort of ! Been tweaked
    ! real :: eta
    real :: epsilon

    integer :: numInterps, ix
    integer procID


    ! Calc eta up here.
    real :: mux, eta
    ! This standard dev gives ~95% of the force within -2*sigma < x < 2*sigma
    real :: sigma

    if(turbine%sigma <= 0) then
        sigma = turbine%length/2
    else
        sigma = turbine%sigma
    end if

    if(turbine%forceFactor <= 0) then
        forceFactor = defaultForceFactor
    else
        forceFactor = turbine%forceFactor
    end if

    eta = (1/(sqrt(2*pi)*sigma)) * exp(-0.5* ((node%x/sigma)**2.) )

    if( turbineDebug ) print*,"% ------ calcBemNodeAbsSrc()"

    node%source%x = 0.0
    node%source%y = 0.0
    node%source%z = 0.0

    node%absorption%x = 0.0
    node%absorption%y = 0.0
    node%absorption%z = 0.0

    node%bodyThrust = 0.0
    axialComp = 0.0

    tDu = 0.0
    tDv = 0.0
    tDw = 0.0


    ! Work out perpendicular distance from node to axis of rotation
    r = sqrt( node%y**2.0 + node%z**2.0 )
    rho = turbine%densityFluid

    ! Check if point inside the general volume (you never know)
    if( r <= turbine%radius ) then


       ! If in general blade volume, not hub, and if mean u +ve...
       ! (do nothing in reverse flow!)
       if( r > turbine%hubRadius .and. turbine%uMean > 0) then

#ifdef OLD_MAG_U_REL
          uMean = turbine%uMean
          uSqMean = turbine%uSqMean
#else
          ! Calculate relative velocity in the new way
          call calculateNodalVelFromInterp(turbine, node, &
               uInt, vInt, wInt, &
               usqInt, vsqInt, wsqInt )

          uMean = abs(uInt)
          uSqMean = abs(usqInt)
#endif


          ! Rotate to span reference frame (around x-axis)

          ! Angle of the blade (not attack angle)
          call calculateBemChordTwist(chord, twist, turbine, r)

          bladeAngle = turbine%alpha + twist

          ! Some simplifications (make calculations stable)

          bladeOrbVel = r * turbine%omegaTurb

          if( r > verySmall ) then
             ! fluidOrbVel = turbine%omegaFLMean
             fluidOrbVel = (node%y * node%w - node%z * node%v ) / r

             !fluidOrbVel = 0.
          else
            fluidOrbVel = 0.
          end if

          ! Not sure if this assumption is strictly correct...
          !            relOrbVel = fluidOrbVel + bladeOrbVel
          ! print*, "bladeOrbVel:", bladeOrbVel
          ! print*, "fluidOrbVel:", fluidOrbVel
          relOrbVel = bladeOrbVel - fluidOrbVel

          magURel = sqrt(uSqMean*(axialInductionFactor**2.)+relOrbVel**2.)

          if( abs(uMean) .lt. verySmall &
               .and. abs(relOrbVel) .lt. verySmall ) then
             velAngle = 0.
          else
             velAngle = atan2( relOrbVel, uMean )
          end if

          attack =   0.5*pi - velAngle - bladeAngle

          R_T = turbine%radius
          mu = r / R_T


          ! Print out diagnostics if in debug mode
          if( turbineDebug .and. r<turbine%radius*0.25) then
             print*,"% uMean:", uMean, "r:", r
             print*,"% uSqMean:", uSqMean, "r:", r
             print*,"% bladeOrbVel:", bladeOrbVel, "r:", r
             print*,"% fluidOrbVel:", fluidOrbVel, "r:", r
             print*,"% relOrbVel:", relOrbVel, "r:", r

             print*,"% bladeAngle:", bladeAngle*radToDeg, "r:", r
             print*,"% velAngle:", velAngle*radToDeg, "r:", r
             print*,"% attack:", attack*radToDeg, "r:", r

             if (abs(attack) > pi) then
                print*, "*** Error: |attack| > 180"
                call parallelHaltSimulation()
                stop
             end if

          end if



          ! Only if there is some flow through the turbine, do something

          if( abs(magURel) > verySmall ) then

             call calculateBemLiftDragCoeffs( c_L, c_D, &
                  turbine, attack, fakeReynolds )

             dr = turbine%meanElementDx
             N = turbine%numBlades
             ls = N * chord / (2. * pi * r)


             ! Tip-speed ratio for wind turbines
             lambda = (turbine%omegaTurb/turbine%omegaTurbMax) &
                  * turbine%tipSpeedRatioMax

             ! flow induction factor. Picked for 'optimised operation',
             ! cf. Wind Energy Handbook.

             ! Mean axial induction factor
             a = 0.333333333333333

             ! NB. Localised axial induction factor can be ~=0.5
             ! see fig. 3.36, pp89.
             ! a=0.5

             ! Calculate tip loss factor

             expo = - 0.5* N * ((1-mu)/mu) &
                  * sqrt( 1+ (lambda*mu)**2.0 / ((1.0-a)**2.0) )
             expN = exp( expo )

             if ( abs(expN) .gt. 1 ) then
                expN = expN / abs(expN)
             end if

             ! tip-loss
             if (turbine%tipLoss .eqv. .true.) then
                tiploss = (2/pi) * acos( expN )
             else
                tiploss = 1.0
             end if

             ! orbital force component
             aF_theta = - (tiploss * 0.5 * rho* magURel**2. &
                  * (c_L * cos(velAngle) - c_D * sin(velAngle) ) &
                  * ls )

             axialComp =  -  tiploss * 0.5 * rho * magURel**2. &
                  * (c_L * sin(velAngle)  + c_D * cos(velAngle) ) &
                  * ls

             ! Aerodynamic forces on fluid without power extraction
             du = forceFactor * eta * axialComp

             ! FIX TO DAMP ROTATION
             ! KEEP UNTIL ROTATIONAL AUGMENTATION / ACTUATOR LINE IMPLEMENTED

             !if(r >= turbine%radius) then
             !   rotfactor = 1
             !else
             !   rotdr = r-turbine%hubRadius
             !   rotfactor = (rotdr/(turbine%radius-turbine%hubRadius))**3.0
             !end if

             dv = - forceFactor * eta * (node%z/r) * aF_theta
             dw =   forceFactor * eta * (node%y/r) * aF_theta

             ! Set to 0 for now
             ! dv = 0.
             ! dw = 0.

             ! Torque per-unit-volume at node and
             ! angle of attack (used later)

             node%fluidBodyTorque = torqueFactor * forceFactor * eta * r * aF_theta
             node%uRel = magURel
             node%attack = attack

          else

             ! The blade has no relative motion to the fluid. Do nothing.
             du = 0.
             dv = 0.
             dw = 0.

             node%fluidBodyTorque = 0.
             node%uRel = 0.
             node%attack = 0.

          end if
       else

          if ( abs(turbine%uHubMean)  .lt. verySmall) then
             inv_dt = 0.
          else
            ! Removing fluid density for now - does Fluidity require it?
             ! inv_dt = turbine%densityFluid * abs(turbine%uHubMean)
             inv_dt = abs(turbine%uHubMean)
          end if

          ! Hub resistive force
          du = - eta * forceFactor * hubFactor * turbine%uHubMean * inv_dt
          ! du = 0.
          dv = 0.
          dw = 0.

          node%fluidBodyTorque = 0.
          node%uRel = 0.
          node%attack = 0.

       end if

       call calculateNodeTurbulence( turbine, node, uMean, r, tDu,tDv,tDw )

       ! Spit out source values (notice the turbulence terms)
       node%source%x = (du / turbine%densityFluid + tDu) / rho
       node%source%y = (dv / turbine%densityFluid+ tDv) / rho
       node%source%z = (dw / turbine%densityFluid+ tDw) / rho

       ! thrust / unit vol for this node
       node%bodyThrust = eta * axialComp

       ! Rotate source and absorption terms back into original reference frame

       call tiltRotate(node%source%x, node%source%y, &
            node%source%z, turbine%upwardTilt)

       call horizRotate(node%source%x, node%source%y, &
            node%source%z, turbine%orientation)

       call tiltRotate(node%absorption%x, node%absorption%y, &
            node%absorption%z, turbine%upwardTilt)

       call horizRotate(node%absorption%x, node%absorption%y, &
            node%absorption%z, turbine%orientation)

    end if

    if( turbineDebug ) print*,"% ------ end calcBemNodeAbs()"

  end subroutine calculateBemDiscNodeForce

end module nodeforce
