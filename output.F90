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

module output
  use turbtypes
  use parallel
  use volume

  implicit none

contains

  subroutine outputTurbineDiagnostics(turbines, t, tEffective, dt)
    implicit none

    type (ModelTurbine), pointer :: turbines(:), turb
    real :: t, tEffective, dt
    integer :: i, j
    logical :: allBemTurbines
    integer :: processID

    character (len=1000) :: formatString


    print*, "---- outputTurbineDiagnostics()"

    if(turbineStdOut) then
		    print*, "@ --------------------------------------------------------------------------"
		    print*, "@ Turbine diagnostics at time=", t
    end if


    if( isParallelRun() ) call parallelGetProcID(processID)

    ! If this is a serial run, or this is a parallel run and this is the
    ! master process

    if( isParallelRun() .eqv. .false. .or. processID.eq.1) then

       ! Is this the first time-step?
       if( tEffective .lt. dt ) then

          ! Are all the turbines BEM? If so, special format string
          allBemTurbines=.true.
          do i=1, size(turbines)
            if( (turbines(i)%status /= BEM_DISC_MODEL .and. &
                turbines(i)%status /= BEM_LINE_MODEL ) &
              .and. turbines(i)%status /= 0 ) allBemTurbines = .false.
          end do

          ! If we have just BEM turbines (likely), output a reduced set of data.
          ! Skips any thesis model diagnostics.
          if(allBemTurbines) then
            formatString='t, id, turbine_nodes, orientation, ' // &
               'mean_u, hub_mean_u, max_u, flow_dir, alpha, weighted_attack_average, '// &
               'bem_meanOmegaFL, turbine_omega, turbine_omega_max, ' // &
               'blade_first_angle, thrust, blade_loading, fluid_density, power, ' // &
               'opt_power_ratio, fluid_torque, inst_mean_u, inst_max_u'
          else
            formatString='t, turbine_no, turbine_nodes, orientation, ' // &
               'mean_u, hub_mean_u, max_u, flow_dir, alpha, weighted_attack_average, '// &
               'thesis_omegaFL, bem_meanOmegaFL, turbine_omega, ' // &
               'turbine_omega_max, thrust, blade_loading, B_effective, fluid_density, power, ' // &
               'opt_power_ratio, fluid_torque, inst_mean_u, inst_max_u'
          end if
          ! Then open the file for writing (overwrite or create)
          open(1, file=trim(turbineOutputFile), status='replace')

          ! The CSV headers
          ! for the record, I despise format statements.
1010      write(1,*) trim(formatString)

       else
          ! Open up previously-existing file for appending entries
          open(1, file=trim(turbineOutputFile), status='old', position='append')
       end if

       print*, "size(turbines):", size(turbines)


       ! If we've found the marine turbine file, load it in.
       do i=1, size(turbines)
          turb => turbines(i)

          ! Print out turbine diagnostics
          if(turbineStdOut) then
              if(i.gt.1) print*, "@"

              print*, "@ Turbine id: ", trim(turb%id)
              print*, "@ native nodes=", turb%numberNativeNodes
              print*, "@ total nodes=", getGlobalTurbineNodes(turb,i)
              print*, "@ mean u=", turb%uMean
              print*, "@ mean u_sq=", turb%uSqMean
              print*, "@ hub mean u=", turb%uHubMean
              print*, "@ max u =", turb%uMax
              print*, "@ orientation=", radToWind( turb%orientation )
              print*, "@ flow dir=", radToWind( turb%flowDir )
              
              print*, "@ omegaFLMean=", turb%omegaFLMean
              print*, "@ omegaFLSqMean=", turb%omegaFLSqMean
              
              
              select case (turb%status)
                 
              case (BEM_DISC_MODEL, BEM_LINE_MODEL)
                 write(*, "(A, F16.12)") " @ blade pitch=", turb%alpha*radToDeg
                 write(*, "(A, F16.12)") " @ attack av.=", turb%attackAv*radToDeg
                 write(*, "(A, F16.12)") " @ weighted attack av.=", turb%wAttack*radToDeg
                 write(*, "(A, F16.12)") " @ turbineOmega=", turb%omegaTurb
                 write(*, "(A, F16.12)") " @ blade_first_angle=", turb%bladeFirstAngle*radToDeg
                 
                 write(*, "(A, F16.12)") " @ mean fluid omega=", turb%omegaFLMean
                 print*, "@ thrust=", turb%thrust
                 print*, "@ blade loading=", turb%bladeLoading
                 print*, "@ fluid density=", turb%densityFluid
                 print*, "@ fluidTorque=", turb%fluidTorque
                 print*, "@ power=", turb%power
                 write(*, "(A, F16.12)") " @ opt power ratio =", &
                      turb%power / (turb%maxPower)
                 print*, "@ instantaneous mean u=", turb%instantUMean
                 print*, "@ instantaneous max u =", turb%instantUMax
                 
                 ! Under BEM, this is the same thing
                 turb%powerOpt = turb%maxPower
                 
              end select
          end if
          ! Write to turbines file

          ! If all the turbines are BEM, output only BEM data

          if( turb%timeUntilOutput <= verySmall .or. &
               turb%outputPeriod < dt ) then
             if(allBemTurbines) then
3142         format ( f10.4, ",", a, ",", 1(i10.1, ","), &
                   17(f20.5,  ","), &
                    e14.8, ",", e14.8 )


              write(1, 3142 ) t, trim(turb%id), getGlobalTurbineNodes(turb,i), &
                   radToWind( turb%orientation ), &
                   turb%uMean, &
                   turb%uHubMean, turb%uMax, &
                   radToWind( turb%flowDir ), &
                   turb%alpha*radToDeg, &
                   turb%wAttack*radToDeg, &
                   turb%omegaFLMean, &
                   turb%omegaTurb, turb%omegaTurbMax, &
                   turb%bladeFirstAngle*radToDeg, &
                   turb%thrust, &
                   turb%bladeLoading, &
                   turb%densityFluid, &
                   turb%power, turb%power / turb%maxPower, &
                   turb%fluidTorque, &
                   turb%instantUMean, turb%instantUMax
              else
                  ! Otherwise print out both BEM and thesis diagnostics
3143         format ( f10.4, ",",  2(i10.1, ","), &
                   14(f20.5,  ","), &
                   f20.5, ",", f20.5, ",", &

                   e14.8, ",", e14.8 )


                 write(1, 3143 ) t, i, turb%numberNodes, &
                       radToWind( turb%orientation ), &
                       turb%uMean, &
                       turb%uHubMean, turb%uMax, &
                       radToWind( turb%flowDir ), &
                       turb%alpha*radToDeg, &
                       turb%wAttack*radToDeg, &
                       turb%omega, turb%omegaFLMean, &
                       turb%omegaTurb, turb%omegaTurbMax, &
                       turb%thrust, turb%bladeLoading, &
                       turb%BEffective, turb%densityFluid, &
                       turb%power, turb%power / (turb%powerOpt), &
                       turb%fluidTorque, &
                       turb%instantUMean, turb%instantUMax
                 end if
            turb%timeUntilOutput=turb%outputPeriod
          end if

          if(turb%outputPeriod >= dt) then
             turb%timeUntilOutput = turb%timeUntilOutput - dt
          end if

          if( turbineDebug ) then
             print *, "==INTERP=="
             do j=1, size(turb%intNodes%r)
                print*, "r:", turb%intNodes%r(j), "u:", turb%intNodes%u(j)
             end do
          end if
       end do

       ! Now that we've finished writing out entries for turbines, close file.
       close(1)

       print*, "@ --------------------------------------------------------------------------"
    end if

  end subroutine outputTurbineDiagnostics

end module output
