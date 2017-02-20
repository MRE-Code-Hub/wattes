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

! ---------------------------------------------------------------------
! Module to save / load checkpoints
! (necessary for Cray / HECToR)

module checkpoint
    use turbtypes
    use parallel

    implicit none

    character(len=256), parameter :: stateFilename = "turbines.state"
    integer, parameter :: fileUnit = 8104

contains



    ! ---------------------------------------------------------------------

    subroutine saveCheckpoint(time, tEffective, turbines)
        real :: time, tEffective
        type(ModelTurbine), pointer :: turbines(:)
        integer :: i, numTurbines, procID
        real, allocatable :: turbinesGlobalState(:,:), rWork(:)

        call parallelGetProcID(procID)

        ! If this is the master process, checkpoint.
        if(procID==1) then

            numTurbines = size(turbines)

            ! If first time step, open and wipe old turbine state file
            if(abs(tEffective) < verySmall) then
                open( fileUnit, file=trim(stateFilename), &
                    action="write", form="formatted", &
                    status="replace", err=101 )
            else
             ! Otherwise append to old one
                open( fileUnit, file=trim(stateFilename), &
                    action="write", form="formatted", &
                    position="append", access="stream", err=101 )
            end if

            allocate(turbinesGlobalState(numTurbines, parWorkArraySize))
            allocate(rWork(parWorkArraySize))

            ! Populate turbines state array
            print*, "writing turbines state for time:", time, ", no. turbines:", numTurbines
            do i=1, numTurbines
                call turbinetoWorkArray(turbines(i), rWork)
                turbinesGlobalState(i, :) = rWork(:)
            end do
            write(fileUnit, *) "time ", time
            write(fileUnit, *) numTurbines, turbinesGlobalState

            ! Save interpnodes state (necessary!)
            do i=1, numTurbines
               write(fileUnit, *) turbines(i)%intNodes%u, turbines(i)%intNodes%rOmega,&
                    turbines(i)%intNodes%uOld, turbines(i)%intNodes%rOmegaOld, &
                    turbines(i)%intNodes%usq, turbines(i)%intNodes%rOmegaSq, &
                    turbines(i)%intNodes%usqOld, turbines(i)%intNodes%rOmegaSqOld
            end do

            print*, "Saving turbines checkpoint data at t=", time

            deallocate(turbinesGlobalState)
            deallocate(rWork)

            close(fileUnit)

        end if

        return
101 write(0,*) "Error opening turbine state file '", trim(stateFilename), "' for writing. Exiting."
        call exit(1)

    end subroutine saveCheckpoint



    ! ---------------------------------------------------------------------



    subroutine loadCheckpoint(checkpointTime, dt, turbines)
        real :: checkpointTime, dt, thisTime
        type(ModelTurbine), pointer :: turbines(:)
        integer :: i, j, numTurbines, thisNumTurbines, numIntNodes
        integer :: filePosition, procID
        real, allocatable :: turbinesGlobalState(:,:), rWork(:)
        logical :: notFoundCheckpoint
        character(len=256) :: label

        real :: savedVal

        print*, "loadCheckpoint()"

        numTurbines = size(turbines)

        allocate(turbinesGlobalState(numTurbines, parWorkArraySize))
        allocate(rWork(parWorkArraySize))

        open( fileUnit, file=trim(stateFilename), &
            form="formatted", status="old", access="stream", &
            action="read", err=102 )

        ! Can't conceive of more than 10^10 time steps!
        notFoundCheckpoint=.true.
        i=1
        do while(i<veryBig .and. notFoundCheckpoint)
            inquire(fileUnit, POS=filePosition)
            read(fileUnit, *, err=103) label, thisTime
            read(fileUnit, *) thisNumTurbines, turbinesGlobalState

            ! Load in turbine interpolated nodes
            do j=1, numTurbines
                read(fileUnit, *, err=103) &
                      turbines(j)%intNodes%u, turbines(j)%intNodes%rOmega,&
                      turbines(j)%intNodes%uOld, turbines(j)%intNodes%rOmegaOld, &
                      turbines(j)%intNodes%usq, turbines(j)%intNodes%rOmegaSq,&
                      turbines(j)%intNodes%usqOld, turbines(j)%intNodes%rOmegaSqOld
            end do

            ! Sanity check
            if(thisNumTurbines .ne. numTurbines) then
                write(0,*) "Error. Specifed", numTurbines, &
                    "turbines, but checkpoint entry has", thisNumTurbines
                call exit(1)
            end if

            ! If we're at the right time-step (or close enough), load turbines
            if( abs(checkpointTime-thisTime) < dt .and. &
                thisNumTurbines==numTurbines) then
                notFoundCheckpoint=.false.

                print*, "Loading turbines checkpoint data for t=", checkpointTime
                do j=1, numTurbines
                    ! If we've fixed the orientation in turbines.dat, don't use the checkpointed value
                    if( turbines(j)%orientationFixed .eqv. .true. ) then
                        savedVal = turbines(j)%orientation
                    end if

                    rWork(:) = turbinesGlobalState(j,:)
                    call workArraytoTurbine(rWork, turbines(j))

                    ! For fixed orientation, restore  initial value
                    if( turbines(j)%orientationFixed .eqv. .true. ) then
                        turbines(j)%orientation = savedVal
                        turbines(j)%oldOrientation = savedVal
                    end if
                end do

                call parallelGetProcID(procID)
                ! If master, truncate files: wipe later turbine states
                if(procID==1) then
                    ! Deleted turbines.state truncation code. Its
                    ! functionality is not really needed.

                    call truncateCSVFileFromNow(thisTime)
                end if
            end if
            i=i+1

        end do
        close(fileUnit)

        deallocate(turbinesGlobalState)
        deallocate(rWork)

        return

        ! I/O errors
102 write(0,*)  "Error opening old turbine state file '", trim(stateFilename), "'. Exiting. "
        call exit(1)

103 write(0,*) "Error: reached end of '", trim(stateFilename), "' but not found checkpoint time", &
            checkpointTime
        call exit(1)

    end subroutine loadCheckpoint



    ! ---------------------------------------------------------------------
    ! Truncates all turbine CSV entries after nowTime.

    subroutine truncateCSVFileFromNow(nowTime)
        real :: nowTime

        integer, parameter :: csvUnit=9851
        integer :: i, thisPosition, lastPosition, fileErr
        integer :: comPosition
        real :: thisTime
        logical :: foundComma, justPastNow
        character(len=10000) :: thisLine, lastLine, colWord

        open( csvUnit, file=trim(turbineOutputFile), &
            form="formatted", status="old", access="stream", &
            action="readwrite", err=202 )

        justPastNow=.false.

        ! Sift through CSV file, looking for time just /after/ checkpoint time.
        ! Not exactly fast, but only done once per continuous simulation
        ! run. If nothing is found, do nothing (as there is nothign to truncate).
        i=1
        thisPosition=1
        lastPosition=1
        thisLine=""
        lastLine=""
        do while (i<veryBig .and. justPastNow .eqv. .false.)

            ! Save position, then read in a line
            inquire(csvUnit, POS=thisPosition)
            read(csvUnit, iostat=fileErr, fmt="(A)") thisLine

            ! If we've reached the end of the file, exit loop.
            if(fileErr<0) exit

            ! If there's a disasterous IO error, abort.
            if(fileErr>0) then
                write(0,*) "Error: fatal I/O error reading '", trim(turbineOutputFile), "'."
                call exit(1)
            end if

            ! Get first column for this row
            comPosition = scan(thisLine, ",", foundComma)
            if(comPosition /= 0 ) then
                colWord=""
                colWord=adjustl(thisLine(1:comPosition-1))
            end if

            ! If it's not empty, and we're not on the first line, it must be a
            ! floating point number.
            if(trim(colWord) /= "" .and. i>1 ) then
                read(colWord, "(F10.4)") thisTime

		            if(thisTime-nowTime > verySmall) then
		                justPastNow = .true.
		            else
		                lastLine(:)=thisLine(:)
		                lastPosition=thisPosition
		            end if
		        end if

            i=i+1
        end do

        ! If we're just past the time we were going to truncate after, rewind
        ! to there and write it out. This has the effect of truncating the file.
        if(justPastNow .and. fileErr==0) then
            print*, "Truncating turbine CSV file '", trim(turbineOutputFile), "' after time=", thisTime
            write(csvUnit, POS=lastPosition, fmt="(A)") trim(lastLine)
        end if

        close(csvUnit)
        return

        ! I/O error

202 write(0,*) "truncateCSVFileFromNow(): failed to open '", &
            trim(turbineOutputFile), "'."
        call exit(1)


    end subroutine truncateCSVFileFromNow

end module checkpoint
