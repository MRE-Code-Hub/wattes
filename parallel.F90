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

module parallel
  use turbtypes
  use memory

  implicit none

  integer, pointer :: numberProcTurbNodes(:,:)

  integer, parameter :: parWorkArraySize=54


contains


  subroutine parallelInitProcNodeArrays(turbs)
    type(ModelTurbine), pointer :: turbs(:)
    integer :: numProcs, i

    call getNumberOfProcessors( numProcs )

    allocate( numberProcTurbNodes(numProcs, size(turbs)) )

    numberProcTurbNodes(:,:) = 0

  end subroutine parallelInitProcNodeArrays


! --------------------------------------------------------------------------




  subroutine parallelSyncProcNodeArrays(turbs)
    type(ModelTurbine), pointer :: turbs(:)
    integer :: numProcs, procID, t, p

    real, pointer :: rProcTurbNodes(:), rWork(:)

    call parallelGetProcID( procID )
    call getNumberOfProcessors( numProcs )

    allocate( rProcTurbNodes(numProcs) )
    allocate( rWork(numProcs) )

    do t=1, size(turbs)
        rProcTurbNodes(:) = 0.
        rProcTurbNodes(procID) = turbs(t)%numberNativeNodes

        call parallelSumArrays( rProcTurbNodes, rWork, numProcs )

        do p=1, numProcs
            numberProcTurbNodes(p, t) = int(rProcTurbNodes(p))
        end do

    end do

    deallocate(rProcTurbNodes)
    deallocate(rWork)


  end subroutine parallelSyncProcNodeArrays



! --------------------------------------------------------------------------



  subroutine parallelPassDummyMPI()

    integer numProcs
    real, pointer :: numTurbProcs(:), alaWork(:)

    ! if(turbineDebug) print*,"    --- parallelPassDummyMPI()"

    call getNumberOfProcessors( numProcs )

    allocate( numTurbProcs(numProcs) )
    allocate( alaWork(numProcs) )

    numTurbProcs = 0
    alaWork = 0

    call parallelSumArrays(numTurbProcs, alaWork, numProcs)

    deallocate( numTurbProcs )
    deallocate ( alaWork )

    ! if(turbineDebug) print*,"    --- end parallelPassDummyMPI()"

  end subroutine parallelPassDummyMPI


  ! --------------------------------------------------------------------------

  subroutine parallelGetProcID( procID )
    integer procID, iError

    if( isParallelRun() ) then
        call mpi_comm_rank( MPI_COMM_WORLD, procID, iError )
    else
        procID = 0
    end if

    procID = procID+1

  end subroutine parallelGetProcID



  ! --------------------------------------------------------------------------


  subroutine parallelHaltSimulation()
    integer iError

    call mpi_abort( MPI_COMM_WORLD, iError )

  end subroutine parallelHaltSimulation



  ! --------------------------------------------------------------------------

  subroutine getNumberOfProcessors( returnNum )
    integer iError, returnNum

    call mpi_comm_size( MPI_COMM_WORLD, returnNum, iError )
  end subroutine getNumberOfProcessors



  ! --------------------------------------------------------------------------

  subroutine parallelSumArrays(mainArray, workArray, arraySize)

    integer arraySize
    real, intent(inout) :: mainArray(arraySize)
    real, intent(inout) :: workArray(arraySize)
    integer n

    integer iError

    ! if( turbineDebug ) print*,"       ::::: parallelSumArrays()"

    if( isParallelRun() ) then

       ! if( turbineDebug ) print*,"       ::::: arraySize:", arraySize

       call mpi_allreduce( mainArray, workArray, arraySize, &
            MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, iError )

       do n=1, arraySize
          mainArray(n) = workArray(n)
       end do
    end if

    ! if( turbineDebug ) print*,"       ::::: end parallelSumArrays()"

  end subroutine parallelSumArrays


  ! --------------------------------------------------------------------------

  function isParallelRun()
    implicit none

    logical isParallelRun
    integer numProcs

    call getNumberOfProcessors( numProcs )
    if( numProcs .le. 1 ) then
       isParallelRun = .false.
    else
       isParallelRun = .true.
    end if

  end function isParallelRun



  ! --------------------------------------------------------------------------


  subroutine parallelGatherTurbineNodes( turb, turbID )
    type (ModelTurbine) :: turb
    integer :: turbID

    real, pointer :: xOrig(:), yOrig(:), zOrig(:), &
         uOrig(:), vOrig(:), wOrig(:)
    real, pointer :: xWork(:), yWork(:), zWork(:), &
         uWork(:), vWork(:), wWork(:)
    real, pointer :: numProcNodes(:)

    integer :: numOrigNodes, numNewNodes
    integer :: iError, MPI_status(MPI_STATUS_SIZE)

    integer :: procID, numProcs, numTotalNodes, numAddedNodes
    integer :: newNumNodes
    integer :: i, p, sp, newix, numWorkNodes

    integer MPI_realType

    if( turbineDebug ) print*, "    == parallelGatherTurbineNodes()"


    ! The number of turbine nodes will change over time, since nodes from
    ! other processes will get merged in
    numOrigNodes = turb%numberNativeNodes

    ! Get information about number of nodes each process has in its space.
    call parallelGetProcID(procID)

    ! if( turbineDebug ) print*,"procID:", procID


    call getNumberOfProcessors(numProcs)
    MPI_realType = MPI_DOUBLE_PRECISION

    allocate( numProcNodes(numProcs) )
    numProcNodes(:) = numberProcTurbNodes(:, turbID)
    numTotalNodes = sum(numProcNodes)

    ! allocate memory for node arrays

    allocate( xWork(numTotalNodes) )
    allocate( yWork(numTotalNodes) )
    allocate( zWork(numTotalNodes) )
    allocate( uWork(numTotalNodes) )
    allocate( vWork(numTotalNodes) )
    allocate( wWork(numTotalNodes) )

    allocate( xOrig(numOrigNodes) )
    allocate( yOrig(numOrigNodes) )
    allocate( zOrig(numOrigNodes) )
    allocate( uOrig(numOrigNodes) )
    allocate( vOrig(numOrigNodes) )
    allocate( wOrig(numOrigNodes) )


    ! Copy the original to send via MPI
    xOrig( 1:numOrigNodes ) = turb%nodes( 1:numOrigNodes )%x
    yOrig( 1:numOrigNodes ) = turb%nodes( 1:numOrigNodes )%y
    zOrig( 1:numOrigNodes ) = turb%nodes( 1:numOrigNodes )%z
    uOrig( 1:numOrigNodes ) = turb%nodes( 1:numOrigNodes )%u
    vOrig( 1:numOrigNodes ) = turb%nodes( 1:numOrigNodes )%v
    wOrig( 1:numOrigNodes ) = turb%nodes( 1:numOrigNodes )%w


    ! Index in xValues, etc for adding more nodes
    newix = numOrigNodes+1

    ! Go through list of processors
    do p=1, numProcs

       ! if( turbineDebug ) print*,"p:", p, " npn(p):", int(numProcNodes(p))

       ! If the list process has turbine nodes
       if( int(numProcNodes(p)).gt.0 ) then

          ! If the list process is /this/ process
          if( p.eq.procID ) then

             ! Send to others, that have turbine nodes themselves
             do sp=1, numProcs
                if (sp .ne. p .and. numProcNodes(sp) > 0) then

                   call MPI_Send( xOrig, numOrigNodes, MPI_realType, &
	                sp-1,1000, MPI_COMM_WORLD, iError )

                   call MPI_Send( yOrig, numOrigNodes, MPI_realType, &
	                sp-1,2000, MPI_COMM_WORLD, iError )

                   call MPI_Send(zOrig, numOrigNodes, MPI_realType, &
	                sp-1,3000, MPI_COMM_WORLD, iError )

                   call MPI_Send( uOrig, numOrigNodes, MPI_realType, &
	                sp-1,4000, MPI_COMM_WORLD, iError )

                   call MPI_Send( vOrig, numOrigNodes, MPI_realType, &
	                sp-1,5000, MPI_COMM_WORLD, iError )

                   call MPI_Send( wOrig, numOrigNodes, MPI_realType, &
	                sp-1,6000, MPI_COMM_WORLD, iError )

                end if

             end do ! end send to others

          else
             ! If this list process is another process, receive nodes
             ! and merge
             numWorkNodes = numProcNodes(p)

             call MPI_Recv( xWork, numWorkNodes, MPI_realType, &
                  p-1, 1000, MPI_COMM_WORLD, MPI_status, iError )

             call MPI_Recv( yWork, numWorkNodes, MPI_realType, &
                  p-1, 2000, MPI_COMM_WORLD, MPI_status, iError )

             call MPI_Recv( zWork, numWorkNodes, MPI_realType, &
                  p-1, 3000, MPI_COMM_WORLD, MPI_status, iError )

             call MPI_Recv( uWork, numWorkNodes, MPI_realType, &
                  p-1, 4000, MPI_COMM_WORLD, MPI_status, iError )

             call MPI_Recv( vWork, numWorkNodes, MPI_realType, &
                  p-1, 5000, MPI_COMM_WORLD, MPI_status, iError )

             call MPI_Recv( wWork, numWorkNodes, MPI_realType, &
                  p-1, 6000, MPI_COMM_WORLD, MPI_status, iError )

             call parallelMergeNodes( turb, &
                  xWork, yWork, zWork, &
                  uWork, vWork, wWork, &
                  numWorkNodes )

          end if
       end if
    end do


    deallocate( xWork )
    deallocate( yWork )
    deallocate( zWork )
    deallocate( uWork )
    deallocate( vWork )
    deallocate( wWork )

    deallocate( xOrig )
    deallocate( yOrig )
    deallocate( zOrig )
    deallocate( uOrig )
    deallocate( vOrig )
    deallocate( wOrig )

    deallocate( numProcNodes )

    if( turbineDebug ) print*, "    == end parallelGatherTurbineNodes()"

  end subroutine parallelGatherTurbineNodes


  ! --------------------------------------------------------------------------
  ! This merges the turbine nodes and external node arrays together. Any point that is
  ! within the error tolerance is assumed to be identical, and so cast aside.
  !
  ! This could be in parallelGatherTurbineNodes(), but I think it's a
  ! little tidier here. Not much performance overhead either, since it's
  ! only called once for each process, per process.
  !

  subroutine parallelMergeNodes( turb, &
       xSrc, ySrc, zSrc, uSrc, vSrc, wSrc, &
       numSrcNodes )

    type(ModelTurbine) :: turb
    type(TurbineNode), pointer :: dnodes(:)

    real, pointer :: xSrc(:), ySrc(:), zSrc(:), uSrc(:), vSrc(:), wSrc(:)
    integer :: numSrcNodes

    real :: dr, speedSrc, speedDest
    real, parameter :: errorTol=10e-7
    integer :: numOrigDestNodes, numAddedNodes, oversized, finalsize
    integer :: s, d, newix
    logical :: nodeIsUnique

    numOrigDestNodes = turb%numberNodes
    oversized = turb%numberNodes+numSrcNodes
    finalsize = turb%numberNodes

    allocate( dnodes(oversized) )
    dnodes(1:numOrigDestNodes) = turb%nodes(1:numOrigDestNodes)

    newix = turb%numberNodes+1


    call resizeNodeArray(dnodes, numOrigDestNodes, oversized)

    do s = 1, numSrcNodes
       nodeIsUnique = .true.

       do d =1, numOrigDestNodes

          dr = sqrt( ( xSrc(s) - dnodes(d)%x )**2.0 &
               + ( ySrc(s) - dnodes(d)%y )**2.0 &
               + ( zSrc(s) - dnodes(d)%z )**2.0 )

          speedSrc = sqrt( uSrc(s)**2.0 + vSrc(s)**2.0 + wSrc(s)**2.0 )
          speedDest = sqrt( dnodes(d)%u**2.0 + dnodes(d)%v**2.0 + dnodes(d)%w**2.0 )

          ! Houston, we have a duplicate point
          if(dr.lt.errorTol) then
             nodeIsUnique = .false.

             ! Always choose the bigger value
             if ( speedSrc.gt.speedDest ) then
                dnodes(d)%u = uSrc(s)
                dnodes(d)%v = vSrc(s)
                dnodes(d)%w = wSrc(s)
             end if

          end if
       end do

       ! If point (i) is unique, add it to the list
       if( nodeIsUnique ) then

          dnodes(newix)%x = xSrc(s)
          dnodes(newix)%y = ySrc(s)
          dnodes(newix)%z = zSrc(s)

          dnodes(newix)%u = uSrc(s)
          dnodes(newix)%v = vSrc(s)
          dnodes(newix)%w = wSrc(s)
          dnodes(newix)%indexRef = -1

          newix = newix + 1
          finalsize = finalsize +1
       end if

    end do

    call resizeNodeArray( turb%nodes, numOrigDestNodes, finalsize )
    turb%nodes(1:finalsize) = dnodes(1:finalsize)

    turb%numberNodes = finalsize
    deallocate(dnodes)

  end subroutine parallelMergeNodes



  ! --------------------------------------------------------------------------
  ! parallelGatherTurbineStats()
  ! Master process gets turbine stats


  subroutine parallelGatherTurbineData(turbines)
    implicit none

    type (ModelTurbine), pointer :: turbines(:)
    real, pointer :: rWork(:), rintWork(:)

    integer t, p, processID, numProcs, nNodes
    integer, pointer :: nNodesArray(:)
    integer iError, MPI_Status(MPI_STATUS_SIZE)
    integer numTotalTurbNodes

    ! Can't really stick this anywhere else.. a bit dirty
    real newPowerTorque

    allocate( rWork(parWorkArraySize) )
    allocate (rintWork(size(turbines(1)%intNodes%u)))

    if( turbineDebug ) print*,"-------- parallelGatherTurbineStats()"


    call parallelGetProcID(processID)
    call getNumberOfProcessors(numProcs)

    allocate( nNodesArray(numProcs) )

    do t=1, size(turbines)

       nNodesArray = 0
       nNodes = turbines(t)%numberNativeNodes
       nNodesArray(processID) = nNodes
       rWork = 0.

       call turbineToWorkArray(turbines(t), rWork)

        if( turbineDebug ) then
	       print*,"pregather==="
	       print*,""
	       print*,"wAttack:", turbines(t)%wAttack
	       print*,"oldWAttack:", turbines(t)%oldWAttack
	       print*,"wAttackDot:", turbines(t)%wAttackDot
	       print*,""
	       print*,"turbOmega:", turbines(t)%omegaTurb
	       print*,"oldTurbOmega:", turbines(t)%oldOmegaTurb
       end if


       call parallelGatherAveragedArray( rWork, turbines(t)%numberNativeNodes )
       call workArrayToTurbine(rWork, turbines(t))

       ! The interpolated nodes data
       rintWork(:) = turbines(t)%intNodes%u
       call parallelGatherAveragedArray( rintWork, turbines(t)%numberNativeNodes )
       turbines(t)%intNodes%u = rintWork(:)

       rintWork(:) = turbines(t)%intNodes%uOld
       call parallelGatherAveragedArray( rintWork, turbines(t)%numberNativeNodes )
       turbines(t)%intNodes%uOld = rintWork(:)

       rintWork(:) = turbines(t)%intNodes%usq
       call parallelGatherAveragedArray( rintWork, turbines(t)%numberNativeNodes )
       turbines(t)%intNodes%usq = rintWork(:)

       rintWork(:) = turbines(t)%intNodes%usqOld
       call parallelGatherAveragedArray( rintWork, turbines(t)%numberNativeNodes )
       turbines(t)%intNodes%usqOld = rintWork(:)


         if( turbineDebug ) then
	       print*,"post-gather==="
	       print*,""
	       print*,"wAttack:", turbines(t)%wAttack
	       print*,"oldWAttack:", turbines(t)%oldWAttack
	       print*,"wAttackDot:", turbines(t)%wAttackDot
	       print*,""
	       print*,"turbOmega:", turbines(t)%omegaTurb
	       print*,"oldTurbOmega:", turbines(t)%oldOmegaTurb
        end if

    end do

    deallocate( rWork )
    deallocate( rintWork )


    if( turbineDebug ) print*,"---- end parallelGatherTurbineStats()"


  end subroutine parallelGatherTurbineData




  ! --------------------------------------------------------------------------
  ! parallelGatherInterpNodes()
  ! Master process gets turbine interpolated nodes

  subroutine parallelGatherInterpNodes(turbs)
    type (ModelTurbine), pointer :: turbs(:)
    integer :: i, numLocal

    do i=1, size(turbs)
       
       numLocal = turbs(i)%numberNativeNodes

       call parallelGatherAveragedArray(turbs(i)%intNodes%u, numLocal)
       call parallelGatherAveragedArray(turbs(i)%intNodes%rOmega, numLocal)

       call parallelGatherAveragedArray(turbs(i)%intNodes%uOld, numLocal)
       call parallelGatherAveragedArray(turbs(i)%intNodes%rOmegaOld, numLocal)

       call parallelGatherAveragedArray(turbs(i)%intNodes%usq, numLocal)
       call parallelGatherAveragedArray(turbs(i)%intNodes%rOmegaSq, numLocal)

       call parallelGatherAveragedArray(turbs(i)%intNodes%usqOld, numLocal)
       call parallelGatherAveragedArray(turbs(i)%intNodes%rOmegaSqOld, numLocal)

    end do


  end subroutine parallelGatherInterpNodes



  subroutine parallelGatherAveragedArray(sentArray, nNodes)
    real, pointer :: sentArray(:)
    integer :: arrsz, i
    integer :: nNodes

    real, allocatable :: rWork(:), sumWork(:), avWork(:), mWork(:)

    integer :: t, p, processID, numProcs
    integer, pointer :: nNodesArray(:)
    integer :: iError, MPI_Status(MPI_STATUS_SIZE)
    integer :: numTotalTurbNodes

    arrsz = size(sentArray)

    allocate( rWork(arrsz) )
    allocate( sumWork(arrsz) )
    allocate( avWork(arrsz) )
    allocate( mWork(arrsz) )

  
    call parallelGetProcID(processID)
    call getNumberOfProcessors(numProcs)

    allocate( nNodesArray(numProcs) )
    
    if(processID.gt.1) then
       ! If we are a slave process, send out information to master
       rWork(:) = sentArray(:)

       call MPI_Send(  nNodes, 1, &
            MPI_INTEGER, 0, 10000, &
            MPI_COMM_WORLD, iError )

       
       if(nNodes.gt.0) then
          ! If we have nodes, send our turbine information to the master
          if( turbineDebug ) then
             do i=1, arrsz
                if(rWork(i) > veryBig .or. IsNaN(rWork(i)) ) then
                   print*, "SLAVE SEND: weighted array: rWork element", i, "is infinite or NaN"
                   call parallelHaltSimulation()
                   stop
                end if
             end do
          end if
       
          call MPI_Send( rWork, arrsz, MPI_DOUBLE_PRECISION, 0, 10001, &
               MPI_COMM_WORLD, iError )
       end if

       ! The clever bit: the master knows about the the turbine, so send all
       ! the processes that don't have turbine nodes information, so that if 
       ! they do at some point have turbine nodes (possibly through adapted
       ! mesh), they'll have up-to-date relaxed values etc.

       call MPI_Recv( rWork, arrsz, MPI_DOUBLE_PRECISION, 0, 10002, &
            MPI_COMM_WORLD, MPI_status, iError )

       if( turbineDebug ) then
          do i=1, arrsz
             if(rWork(i) > veryBig .or. IsNaN(rWork(i)) ) then
                print*, "SLAVE RECV: weighted array: rWork element", i, "is infinite or NaN"
                call parallelHaltSimulation()
                stop
             end if
          end do
       end if
       

       sentArray(:) = rWork(:)

    else
       ! This process is the master, initialise values for summing
       sumWork = 0
       numTotalTurbNodes = 0

       if( nNodes .gt. 0 ) then
          mWork(:) = sentArray(:)

          numTotalTurbNodes = nNodesArray(1)
          sumWork = mWork * nNodesArray(1)
       end if

       ! The master process gathers turbine information from slaves
       do p=2, numProcs
          ! For each process

          call MPI_Recv( nNodesArray(p), 1, MPI_INTEGER, p-1, 10000, &
               MPI_COMM_WORLD, MPI_status, iError )

          ! If process p contains turbine nodes
          if( nNodesArray(p) .gt.0 ) then
             call MPI_Recv( rWork, arrsz, MPI_DOUBLE_PRECISION, p-1, 10001, &
                  MPI_COMM_WORLD, MPI_status, iError )
             
             if( turbineDebug ) then
                do i=1, arrsz
                   if(rWork(i) > veryBig .or. IsNaN(rWork(i)) ) then
                      print*, "MASTER RECV: weighted array: rWork element", i, "recv'd from proc", p, "is infinite or NaN"
                      call parallelHaltSimulation()
                      stop
                   end if
                end do
             end if

             ! Each contribution to these pan-process values are weighted
             ! according to the number of native nodes each process holds
             
             sumWork = sumWork + rWork * nNodesArray(p)
             numTotalTurbNodes = numTotalTurbNodes + nNodesArray(p)
          end if
       end do
       
       
       ! Averaged weighted array contributions
       
       if( numTotalTurbNodes .eq. 0 ) then
          sumWork = mWork
          numTotalTurbNodes = 1
       end if
       avWork = sumWork / numTotalTurbNodes

       sentArray(:) = avWork(:)
       
       ! The master now sends to each slave processor
       
       do p=2, numProcs
          call MPI_Send( avWork, arrsz, MPI_DOUBLE_PRECISION, p-1, 10002, &
               MPI_COMM_WORLD, iError )

          if( turbineDebug ) then
             do i=1, arrsz
                if(rWork(i) > veryBig .or. IsNaN(rWork(i)) ) then
                   print*, "MASTER SEND: weighted array: rWork element", i, "sent to proc", p, "is infinite or NaN"
                   call parallelHaltSimulation()
                   stop
                end if
             end do
          end if

       end do
       
    end if
    
    
    deallocate( nNodesArray )
    
    deallocate( rWork )
    deallocate( sumWork )
    deallocate( avWork )
    deallocate( mWork )
    
    
  end subroutine parallelGatherAveragedArray


  ! --------------------------------------------------------------------------
  ! turbineToWorkArray()

  subroutine turbineToWorkArray(turb, wArray)
   real :: wArray(:)
   type(ModelTurbine) :: turb

   wArray(1) = turb%omega
   wArray(2) = turb%omegaFree
   wArray(3) = turb%oldOmegaFree
   wArray(4) = turb%omegaTurb
   wArray(5) = turb%omegaTurbMax
   wArray(6) = turb%oldOmegaTurbMax
   wArray(7) = turb%bEffective
   wArray(8) = turb%alpha
   wArray(9) = turb%oldAlpha
   wArray(10) = turb%uMax
   wArray(11) = turb%oldUMax
   wArray(12) = turb%uMean
   wArray(13) = turb%oldUMean
   wArray(14) = turb%deltaOmegaSqMean
   wArray(15) = turb%oldDeltaOmegaSqMean
   wArray(16) = turb%power
   ! wArray(17) = turb%numberNativeNodes
   wArray(18) = turb%instantUMean
   wArray(19) = turb%instantUMax
   wArray(20) = turb%omegaTurb
   wArray(21) = turb%oldOmegaTurb
   wArray(22) = turb%alphaDot
   wArray(23) = turb%omegaTurbDot
   wArray(24) = turb%uHubMean
   wArray(25) = turb%oldUHubMean
   wArray(26) = turb%instFluidTorque
   wArray(27) = turb%attackAv
   wArray(28) = turb%omegaFLMean
   wArray(29) = turb%oldOmegaFLMean
   wArray(30) = turb%fluidTorque
   wArray(31) = turb%oldFluidTorque
   wArray(32) = turb%flowDir
   wArray(33) = turb%oldFlowDir
   wArray(34) = turb%orientation
   wArray(35) = turb%uSqMean
   wArray(36) = turb%oldUSqMean
   wArray(37) = turb%omegaFLSqMean
   wArray(38) = turb%oldOmegaFLSqMean
   wArray(39) = turb%vMean
   wArray(40) = turb%oldVMean
   wArray(41) = turb%orientation
   wArray(42) = turb%oldOrientation
   wArray(43) = turb%orientOmega

   wArray(44) = turb%oldOmegaTurb
   wArray(45) = turb%omegaTurbDot
   wArray(46) = turb%oldOmegaTurbDot

   wArray(47) = turb%wAttack
   wArray(48) = turb%oldWAttack
   wArray(49) = turb%wAttackDot

   wArray(50) = turb%alphaDot
   wArray(51) = turb%oldAlphaDot

   wArray(52) = turb%thrust
   wArray(53) = turb%bladeLoading
   wArray(54) = turb%bladeFirstAngle

end subroutine turbineToWorkArray



  ! --------------------------------------------------------------------------
  ! workArrayToTurbine()


  subroutine workArrayToTurbine(wArray, turb)
    real :: wArray(:)
    type(ModelTurbine) :: turb

   turb%omega = wArray(1)
   turb%omegaFree = wArray(2)
   turb%oldOmegaFree = wArray(3)
   turb%omegaTurb = wArray(4)
   turb%omegaTurbMax = wArray(5)
   turb%oldOmegaTurbMax = wArray(6)
   turb%bEffective = wArray(7)
   turb%alpha = wArray(8)
   turb%oldAlpha = wArray(9)
   turb%uMax = wArray(10)
   turb%oldUMax = wArray(11)
   turb%uMean = wArray(12)
   turb%oldUMean = wArray(13)
   turb%deltaOmegaSqMean = wArray(14)
   turb%oldDeltaOmegaSqMean = wArray(15)
   turb%power = wArray(16)
   ! turb%globalNumberNodes = int( wArray(17) )
   turb%instantUMean = wArray(18)
   turb%instantUMax = wArray(19)
   turb%omegaTurb = wArray(20)
   turb%oldOmegaTurb = wArray(21)
   turb%alphaDot = wArray(22)
   turb%omegaTurbDot = wArray(23)
   turb%uHubMean = wArray(24)
   turb%oldUHubMean =wArray(25)
   turb%instFluidTorque = wArray(26)
   turb%attackAv = wArray(27)
   turb%omegaFLMean = wArray(28)
   turb%oldOmegaFLMean = wArray(29)
   turb%fluidTorque = wArray(30)
   turb%oldFluidTorque = wArray(31)
   turb%flowDir = wArray(32)
   turb%oldFlowDir = wArray(33)
   turb%orientation = wArray(34)
   turb%uSqMean = wArray(35)
   turb%oldUSqMean = wArray(36)
   turb%omegaFLSqMean = wArray(37)
   turb%oldOmegaFLSqMean = wArray(38)
   turb%vMean = wArray(39)
   turb%oldVMean = wArray(40)
   turb%orientation = wArray(41)
   turb%oldOrientation = wArray(42)
   turb%orientOmega = wArray(43)

   turb%oldOmegaTurb = wArray(44)
   turb%omegaTurbDot = wArray(45)
   turb%oldOmegaTurbDot = wArray(46)

   turb%wAttack = wArray(47)
   turb%oldWAttack = wArray(48)
   turb%wAttackDot = wArray(49)

   turb%alphaDot = wArray(50)
   turb%oldAlphaDot = wArray(51)

   turb%thrust = wArray(52)
   turb%bladeLoading = wArray(53)
   turb%bladeFirstAngle = wArray(54)

end subroutine workArrayToTurbine



end module parallel
