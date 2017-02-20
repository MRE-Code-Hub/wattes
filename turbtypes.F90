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

module turbtypes
  implicit none


  ! MPI support is assumed.

#ifndef nompi
include 'mpif.h'
#endif

  type Vector
     real :: x, y, z
  end type Vector


  ! Element mesh nodes
  type TurbineNode
     integer :: indexRef
     real :: x, y, z
     real :: u, v, w
     real :: bladeDistrib
     real :: fluidBodyTorque, bodyThrust, attack, uRel

     type(Vector) :: source, absorption
  end type TurbineNode


  type InterpolatedNodes
     integer numberNodes

     real, pointer :: r(:)

     real, pointer :: u(:), rOmega(:)
     real, pointer :: uOld(:), rOmegaOld(:)

     real, pointer :: usq(:), rOmegaSq(:)
     real, pointer :: usqOld(:), rOmegaSqOld(:)
  end type InterpolatedNodes



type BladeElement
     real :: radius, chord, twist
end type BladeElement


type BladeCoefficients
     real :: angle, lift, drag
end type BladeCoefficients



  ! Turbine type
  type ModelTurbine
     ! file variables
     integer :: status
     character(len=256) :: id
     real :: x, y, z
     real :: orientation, oldOrientation ! in radians
     real :: orientRevTime     ! How long for one complete revolution?
     real :: orientOmegaDotMax ! How long to slow down from full orient rotation?
     real :: orientOmega
     real :: upwardTilt	    ! Angle that turbine points towards the sky
     real :: flowDir, oldFlowDir ! Whither the wind blows
     real :: radius
     real :: length

     ! ========== BEM model-specific variables

     type (BladeElement), pointer :: bladeElements(:)
     type (BladeCoefficients), pointer :: bladeCoeffs(:)

     logical :: autoTwist, autoChord, orientationFixed
     logical :: traditionalAttackOpt, tipLoss

     real :: outputPeriod, timeUntilOutput

     integer :: numBlades
     real :: maxPower, powerEfficiency, tipSpeedRatioMax, u0_MaxPower, k_tau
     real :: massPerUnitArea
     real :: omegaFLMean, oldOmegaFLMean
     real :: omegaFLSqMean, oldOmegaFLSqMean
     real :: inertia, fluidTorque, oldFluidTorque
     real :: powerTorque, bladesTorque
     real :: attackAv, wAttack, oldWAttack, wAttackDot
     real :: instFluidTorque
     real :: alpha, oldAlpha, alphaDot, oldAlphaDot
     real :: attackMaxPerf, optLift

     real :: thrust, bladeLoading, bladesArea

     real :: windRampTime, forceFactor, sigma

     integer :: bladeSense

    ! Actuator line / surface specific

    ! Angular position of first blade
    real :: bladeFirstAngle
    real, pointer :: bladeRelativeAngle(:)

    logical :: haveBladeDistrib
     
     ! ========== Thesis model and general variables

     real :: initialAlpha, targetAlpha, alphaRelaxTime
     ! Net solidity
     real :: B, BEffectiveOpt
     ! F = flow induction factor. Makes the blades spin.
     real :: F
     ! Q = power efficiency of turbine: proportion of input power extracted,
     ! sort of.
     real :: Q
     ! What is kappa? A scaling factor for net effective solidity.
     real :: kappa
     ! Scaling constant for local effective solidity
     real :: k
     ! The optimum and maximum operating angular velocity
     ! of this particular turbine
     real :: omegaTurbOpt, omegaTurbCutout
     real :: powerOpt, powerCoeffOpt
     real :: powerCutout, powerCoeffCutout
     real :: u0_Opt, u0_Cutout, uCutIn, uCutOut

     real :: turbulenceOpt
     real :: uMaxOpt
     real :: tipWidthFraction, tipRadiusFraction
     real :: tipStartRadius

     ! very roughly, how many seconds we want to relax over.
     real :: turbineRelaxSecs, bladeRelaxSecs
     ! and the resultant relaxation parameters
     real :: turbineRelax, bladeRelax

     ! ================= calculated variables ======================

     real :: densityFluid, fluidDensityScaling

     type (TurbineNode), pointer :: nodes(:)

     ! These are for the interpolated grid
     type(InterpolatedNodes) :: intNodes

     real :: hubRadius, hubLocalSolidity

     real :: omega, oldOmega, omegaFree, oldOmegaFree
     ! omegaTurb not relaxed. Just here for convienience
     real :: omegaTurb, oldOmegaTurb, omegaTurbDot, oldOmegaTurbDot
     real :: omegaTurbMax, oldOmegaTurbMax

     real :: deltaOmegaSqMean, oldDeltaOmegaSqMean

     real :: uHubMean, oldUHubMean, instantUHubMeany
     real :: uSqMean, oldUSqMean
     real :: uMax, oldUMax
     real :: instantUMean, instantVMean, instantWMean
     real :: instantUSqMean, instantVSqMean, instantWSqMean
     real :: instantUMax
     real :: uMean, oldUMean, vMean, oldVMean, wMean, oldWMean
     real :: meanElementDx

     ! Effective net solidity, and old value for it
     real :: BEffective
     real :: power

     ! turbulence
     real :: turbulenceTheta
     real :: turbulencePhi
     real :: turbulenceTI

     real :: volume

     ! Need this, since size incorrectly returns a value for 1 for
     ! an non-initialised pointer array
     integer :: numberNodes, numberNativeNodes

  end type ModelTurbine

end module turbtypes
