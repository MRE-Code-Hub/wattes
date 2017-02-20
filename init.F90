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

module init
  use turbtypes
  use bem
  use chordlift
  use interpnodes
  use parallel

  implicit none

  private
  public :: initialiseTurbines

contains


  ! --------------------------------------------------------------------


  subroutine initialiseTurbines( turbines, haveBladeDistrib, densityFluid, dt)

    type (ModelTurbine), pointer :: turbines(:)
    type (ModelTurbine) :: turb

    integer :: i, j, numBEMs
    real :: densityFluid, dt
    real :: a, b, m, S, NcMean, Nc1, R_T, R_H
    real :: bladeLength, B_blades, B_hub
    real :: turbRelax, bladeRelax

    real :: areaSum, chordAv, dr
    integer :: numBlades

    logical :: haveBladeDistrib

    print*, "--* initialiseTurbines()"

    ! Set up the random number generator
    call initGaussGen( dt, 0 )

    do i=1, size(turbines)

       turb = turbines(i)

       R_H = turb%hubRadius
       R_T = turb%radius

       turb%densityFluid = densityFluid * turb%fluidDensityScaling

       turb%uMean = 0.
       turb%vMean = 0.
       turb%wMean = 0.
       
       turb%oldUMean = 0.
       turb%oldVMean = 0.
       turb%oldWMean = 0.
       
       turb%uSqMean = 0.
       turb%oldUSqMean = 0.
       
       turb%omegaTurb = 0.
       turb%oldOmegaTurb = 0.
       turb%omegaTurbDot = 0.
       
       turb%fluidTorque = 0.
       turb%oldFluidTorque = 0.
       
       turb%omega = 0.
       turb%omegaFree = 0.
       turb%oldOmegaFree = 0.

       turb%oldOmegaTurbMax=0.
       turb%bEffective=0.
       turb%uMax=0.
       turb%oldUMax=0.
       turb%uMean=0.
       turb%oldUMean=0.
       turb%deltaOmegaSqMean=0.
       turb%oldDeltaOmegaSqMean=0.
       turb%power=0.
       turb%instantUMean=0.
       turb%instantUMax=0.
       turb%uHubMean=0.
       turb%oldUHubMean=0.
       turb%instFluidTorque=0.
       turb%omegaFLMean=0.
       turb%oldOmegaFLMean=0.
       turb%flowDir=0.
       turb%oldFlowDir=0.
       turb%uSqMean=0.
       turb%oldUSqMean=0.
       turb%omegaFLSqMean=0.
       turb%oldOmegaFLSqMean=0.
       turb%timeUntilOutput=0.


       select case( turb%status )

       case( BEM_DISC_MODEL, BEM_LINE_MODEL )
          turb%alpha = turb%initialAlpha
          turb%oldAlpha = turb%alpha
          turb%alphaDot = 0.
          turb%wAttack = 0.
          turb%oldWAttack = 0.
          turb%wAttackDot = 0.
          turb%thrust = 0.
          turb%bladeLoading = 0.

          turb%omegaTurbMax =  turb%tipSpeedRatioMax &
               * turb%u0_MaxPower / turb%radius

          turb%k_Tau = turb%maxPower &
               / (turb%omegaTurbMax**3.)

         ! This sets the relative angles of each blade. Currently fixed
         ! delta, but this could be expanded upon at a later date for
         ! irregular blade arrangements.
         if (turb%status==BEM_LINE_MODEL ) then
            numBlades = turb%numBlades
            allocate( turb%bladeRelativeAngle(numBlades) )
            do j=1, numBlades
                turb%bladeRelativeAngle(j) = (j-1) * 2. * pi / (1.*numBlades )
                print*, "calc rel angle:", (j-1) * 2. * pi / (1.*numBlades )
                print*, "@@ calc rel angle:", turb%bladeRelativeAngle(j)
            end do

            turb%haveBladeDistrib = haveBladeDistrib

            print*, "@ bladeSense: ", turb%bladeSense

         end if


          ! Calculate blade area (1/2 total surface area)
          numBEMs = size(turb%bladeElements)
          areaSum = 0.

          do j=1, numBEMs-1
            chordAv = ( turb%bladeElements(j)%chord &
                + turb%bladeElements(j+1)%chord ) / 2.
            dr = abs( turb%bladeElements(j+1)%radius &
                - turb%bladeElements(j)%radius )
            areaSum = areaSum  + chordAv * dr
          end do

          turb%bladesArea = areaSum * turb%numBlades

          ! turb%inertia = calculateBemTurbineInertia( turbines(i) )
          turb%inertia = turb%numBlades &
               * turb%radius**3.0 * turb%massPerUnitArea

          ! If relaxation time is 0, don't relax.
          if( abs(turb%turbineRelaxSecs) .lt. verySmall) then
             turbRelax = 0.

          else

             ! Otherwise proceed as normal
             turbRelax = 1 - dt / turb%turbineRelaxSecs

             ! If less than relaxDefault, set value to that
             if(turbRelax .lt. relaxDefault ) then
                print*,"@ WARNING: turbRelax falling back to relaxDefault=", &
                     relaxDefault
                turbRelax = relaxDefault
             end if
          end if

          turb%turbineRelax = turbRelax
          ! Reasonable estimate for now
          turb%tipRadiusFraction = 0.25
          turb%tipStartRadius = (1-turb%tipRadiusFraction) &
               * turb%radius

          call calculateAttackMaxLiftPerf( turbines(i) )

       end select



       turb%numberNodes=0

       ! Incompressible for now.
       turb%densityFluid=densityFluid
       turb%volume = pi * turb%radius**2. * turb%length

       write(*,*) "-----------------------------"
       write(*,*) "Turbine(",trim(turb%id),")"
       select case(turb%status)
       case (THESIS_MODEL)
          ! Should never get here, but just in case.
          print*,"turbine status: thesis model - no longer supported: exiting"
          stop
       case (BEM_DISC_MODEL)
          print*,"turbine status: BEM disc model"
       case (BEM_LINE_MODEL)
          print*,"turbine status: BEM line model"
       case (TURBINE_OFF)
          print*,"turbine status: off"
       case default
          print*,"turbine status unrecognised"
          stop
       end select

       write(*,*) "Position: (", &
            turb%x, turb%y, turb%z, ")"
       write(*,*) "Orientation: ", radToWind(turb%orientation)
       write(*,*) "Orient rev time:", turb%orientRevTime
       write(*,*) "Radius: ", turb%radius
       write(*,*) "Hub radius: ", turb%hubRadius
       write(*,*) "Tip radius fraction: ", turb%tipRadiusFraction
       write(*,*) "Length: ", turb%length
       write(*,*) "Cut-in speed: ", turb%uCutIn
       write(*,*) "Cut-out speed: ", turb%uCutOut
       write(*,*) "Turbulence opt: ", turb%turbulenceOpt
       write(*,*) "Turbine relax time(s): ", turb%turbineRelaxSecs
       
       if( turb%status==BEM_DISC_MODEL &
            .or. turb%status==BEM_LINE_MODEL ) then
        write(*,*) "Omega_turb_max: ", turb%omegaTurbMax
        write(*,*) "Blades area: ", turb%bladesArea
        write(*,*) "Tip loss: ", turb%tipLoss
        write(*,*) "Traditional attack opt: ", turb%traditionalAttackOpt
        write(*,*) "Output period: ", turb%outputPeriod
       end if

       turbines(i) = turb

       ! Lastly, create interpolated nodes arrays

    end do

    call initInterpolatedNodes( turbines )

    write(*,*) "Number of turbines: ", size(turbines)
    write(*,*) "--* end initialiseTurbines()"

  end subroutine initialiseTurbines


  ! --------------------------------------------------------------------


end module init
