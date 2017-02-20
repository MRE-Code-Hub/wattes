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

module main
  use turbtypes
  use load
  use init
  use volume
  use bladepitch
  use absorbsrc
  use memory
  use parallel
  use output
  use transform
  use checkpoint

  implicit none

  private
  public :: turbineFarmSolver

  ! Module-wide persistent variables
  integer :: turbinesDisabled, turbinesLoaded
  type (ModelTurbine), pointer :: turbines(:)

contains


  ! -----------------------------------------------------------------------------


  subroutine turbineFarmSolver( numMeshNodes, &
       xNode, yNode, zNode, &
       uNode, vNode, wNode, &
       densityFluid, &
       xSource, ySource, zSource, &
       tActual, dt, &
       updateTurbine, outputTurbineData, bladeDistrib )

    real, parameter :: errorTol = verySmall

    integer numMeshNodes
    real, dimension(numMeshNodes), intent(in) ::  xNode, yNode, zNode
    real, dimension(numMeshNodes), intent(in) ::  uNode, vNode, wNode
    real, intent(in) ::  densityFluid
    real, dimension(numMeshNodes) ::  xSource, ySource, zSource
    real, dimension(numMeshNodes), optional :: bladeDistrib

    real :: tActual, tEffective, dt, dummyVar

    real, save :: tStart

    logical :: updateTurbine, outputTurbineData, stateFileFound
    logical :: haveBladeDistrib

    integer :: i, procID, stateFileErr

    print*, "Entered turbineFarmSolver()"

    ! Set module flag
    if( present(bladeDistrib) ) then
        print*, "bladeDistributionField present"
        haveBladeDistrib = .true.
    else
        haveBladeDistrib = .false.
    end if


    stateFileFound = .false.

    ! Is this the first time we've come here? Then try and load in the 
    ! turbines

    tEffective=tActual

    if(turbinesLoaded.eq.0 .and. turbinesDisabled.eq.0) then

      open(4721, file="turbines.state", status="old", &
            action="read", iostat=stateFileErr)

      if(stateFileErr ==  0) then
        stateFileFound=.true.
        close(4721)

      else
        print*, "Not found turbines.state"
        tStart = tActual
        print*, "tStart:", tStart

        tEffective = tActual - tStart
      end if

       call loadTurbines(turbines, turbinesDisabled, turbinesLoaded)
       if(turbinesDisabled.eq.1) then
          print*, "---* turbines disabled!"
       else
          ! Notice the dirty hack for incompresible fluid density
          call initialiseTurbines(turbines, haveBladeDistrib, densityFluid, dt)

          if( isParallelRun() ) call parallelInitProcNodeArrays(turbines)

          ! If we've got a non-zero time-step, assume we're reading
          ! in a checkpoint -- for the prior time-step.
          if( abs(tEffective) > verySmall ) then
            call loadCheckpoint(tActual-dt, dt, turbines)
          end if

       end if
    end if



    ! If the turbine routines have have been enabled (typically if
    ! turbines.dat exists and has been processed correcty),
    ! then run them.
    if(turbinesDisabled==0) then

       call collectTurbineNodes( numMeshNodes, &
            xNode, yNode, zNode, &
            uNode, vNode, wNode, &
            turbines, updateTurbine )

       ! For this section, only do something if we have nodes
       do i=1, size(turbines)

          ! This should never happen, but if it does, best to halt the
          ! simulation as continuing can give unpredictable results.
          if( isTurbineNodeless(turbines(i), i) ) then

            call parallelGetProcID(procID)
            print*,""
            write(*,"(2(a, i3))") "*** STOP: Processor ", procID, " reports no nodes in any processor for turbine ", i
            write(*,"(a)") "*** Undefined behaviour. Please check initial mesh and TurbinePresenceFn field"

            if ( isParallelRun() ) call parallelHaltSimulation()
            stop
          end if


          ! Only do the following if we have nodes - saves on CPU/memory usage
          if( turbines(i)%numberNativeNodes > 0 ) then
             call transformToTurbineCoords(turbines(i))
             
             if( updateTurbine ) then
                call correctTurbineBladePitch(turbines(i), tEffective, dt)
             end if
             
             call calculateAbsorptionSourceTerms(turbines(i), tEffective, dt, &
                  updateTurbine)
             
           end if

       end do

       ! Now update the absorption/source terms for the CFD solver

      if( .not. present(bladeDistrib) ) then
           print*, "bladeDistrib not present"
           call updateGlobalAbsorptionSources(turbines, &
                      xSource, ySource, zSource, &
                      numMeshNodes)
      else
           call updateGlobalAbsorptionSources(turbines, &
                      xSource, ySource, zSource, &
                      numMeshNodes, &
                      bladeDistrib )
      end if

       ! If running in parallel mode, turbine data needs to be gathered.
       if( isParallelRun() ) then
          call parallelGatherTurbineData(turbines)
       end if
       
       if( outputTurbineData ) then
          call outputTurbineDiagnostics(turbines, tActual, tEffective, dt )

           ! Checkpointing turbine states.
          call saveCheckpoint(tActual, tEffective, turbines)
       end if


       call deallocateTurbineNodes(turbines)

    end if


  end subroutine turbineFarmSolver

end module main
