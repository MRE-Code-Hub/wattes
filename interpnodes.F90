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

module interpnodes
  use turbtypes
  implicit none

contains


  ! --------------------------------------------------------------------------

  subroutine initInterpolatedNodes( turbs )

    type(ModelTurbine), pointer :: turbs(:)
    type(InterpolatedNodes) :: inodes

    real :: r, dr
    integer :: t, n, nint

    do t=1, size(turbs)

       if( turbs(t)%status > 0 ) then

          inodes = turbs(t)%intNodes 
          nint = NUMBER_INTERP
          
          inodes%numberNodes = nint
          
          allocate( inodes%r(nint) )
          
          allocate( inodes%u(nint) )
          allocate( inodes%rOmega(nint) )
          
          allocate( inodes%uOld(nint) )
          allocate( inodes%rOmegaOld(nint) )
          
          allocate( inodes%usq(nint) )
          allocate( inodes%rOmegaSq(nint) )
          
          allocate( inodes%usqOld(nint) )
          allocate( inodes%rOmegaSqOld(nint) )
          

          ! Initialising the interpolated node values
          
          dr = turbs(t)%radius / (1.*(nint-1))

          do n=1, nint
             inodes%r(n) = (n-1) * dr
             
             inodes%u(n) = 0.
             inodes%uOld(n) = 0.
             
             inodes%usq(n) = 0.
             inodes%usqOld(n) = 0.

             inodes%rOmega(n) = 0.
             inodes%rOmegaOld(n) = 0.

             inodes%rOmegaSqOld(n) = 0.
             inodes%rOmegaSq(n) = 0.
          end do
          
          turbs(t)%intNodes = inodes
       end if

    end do
       
  end subroutine initInterpolatedNodes



  ! --------------------------------------------------------------------------
  ! Use weighted distance to calculate interpolated values
  ! This really REALLY needs to be replaced with proper interpolation

  subroutine calculateInterpolatedNodes ( turb )
    type(ModelTurbine) :: turb

    real, allocatable :: rworking(:), runiq(:)
    real, allocatable :: uuniq(:), rOmegaUniq(:)
    real, allocatable :: usquniq(:), rOmegaSqUniq(:)
    integer :: numNodes, numInt, numUniq, i
    
    integer :: addHub, addTip

    numNodes = turb%numberNodes
    allocate( rworking(turb%numberNodes) )

    call calculateOrderedUniqueR( turb, rworking, numUniq )

    addHub = 0
    addTip = 0
    if( rworking(1) > verySmall ) addHub = 1
    if( rworking(numUniq) < turb%radius ) addTip = 1

    numUniq = numUniq + addHub + addTip

    allocate( runiq(numUniq) )

    allocate( uuniq(numUniq) )
    allocate( rOmegaUniq(numUniq) )

    allocate( usquniq(numUniq) )
    allocate( rOmegaSqUniq(numUniq) )

    runiq(1+addHub:numUniq-addTip) = rworking(1:numUniq-addHub-addTip)

    call calculateOrderedVels(turb, runiq, &
         uuniq, rOmegaUniq, &
         usquniq, rOmegaSqUniq, &
         addHub, addTip)

    call calculateSmoothRelaxInterpVels( turb, runiq, &
         uuniq, rOmegaUniq, &
         usquniq, rOmegaSqUniq )


    deallocate(runiq)

    deallocate(uuniq)
    deallocate(rOmegaUniq)

    deallocate(usquniq)
    deallocate(rOmegaSqUniq)

    deallocate(rworking)

  end subroutine calculateInterpolatedNodes


  ! -----------------------------------------------------------------------------


  subroutine calculateOrderedUniqueR( turb, runiq, numUniq )
    type (ModelTurbine) :: turb
    real, allocatable :: runiq(:)
    integer :: numUniq

    real :: r, rtemp
    integer :: i, j, ixnew
    logical :: ixexists

    runiq = -1
    ixnew=0


    ! Find unique values of r
    do i=1, turb%numberNodes
       r = sqrt( turb%nodes(i)%y**2. + turb%nodes(i)%z**2. )
       
       ixexists = .false.

       do j=1, ixnew
          if( abs(runiq(j)-r) < verySmall ) then
             ixexists = .true.
          end if
       end do 

       if( ixexists .eqv. .false. ) then
          ixnew = ixnew+1
          runiq(ixnew) = r
       end if

    end do

    numUniq = ixnew

    ! Order runiq(:)
    do i=1, numUniq
       do j=i, numUniq
          if( runiq(i) > runiq(j) ) then
             rtemp = runiq(i)
             runiq(i) = runiq(j)
             runiq(j) = rtemp
          end if
       end do
    end do


  end subroutine calculateOrderedUniqueR


  ! -----------------------------------------------------------------------------


  subroutine calculateOrderedVels( turb, runiq, uuniq, rOmegaUniq, &
       usquniq, rOmegaSqUniq, &
       addHub, addTip )
    type (ModelTurbine) :: turb
    real, allocatable :: runiq(:), uuniq(:), rOmegaUniq(:)
    real, allocatable :: usquniq(:), rOmegaSqUniq(:)
    integer :: addHub, addTip

    type (TurbineNode) :: node
    real, allocatable :: usum(:), rOmegaSum(:)
    real, allocatable :: usqsum(:), rOmegaSqSum(:)
    real :: r, rnode, rOmega, rOmegaSq

    integer, allocatable :: ct(:)
    integer :: i, j, numUniq

    numUniq = size(runiq)

    allocate( usum(numUniq) )
    allocate( usqsum(numUniq) )

    allocate( rOmegaSum(numUniq) )
    allocate( rOmegaSqSum(numUniq) )

    allocate( ct(numUniq) )

    usum=0
    rOmegaSum=0

    usqsum=0
    rOmegaSqSum = 0

    ct=0

    ! Look for nodes with same radii as in runiq(:)
    do i=1, numUniq
       r = runiq(i)

       do j=1, turb%numberNodes
          node = turb%nodes(j)
          rnode = sqrt( node%y**2. + node%z**2. )
          rOmega = rnode * 0.5 * (-node%v/node%z + node%w/node%y)
          rOmegaSq = rOmega**2.

          if( abs(r-rnode) < verySmall ) then
             usum(i) = usum(i) + node%u
             usqsum(i) = usqsum(i) + node%u**2.

             rOmegaSum(i) = rOmegaSum(i) + rOmega
             rOmegaSqSum(i) = rOmegaSqSum(i) + rOmegaSq

             ct(i) = ct(i) + 1
          end if
       end do
    end do


    do i=1, numUniq
       uuniq(i) = usum(i) / (1.*ct(i))
       usquniq(i) = usqsum(i) / (1.*ct(i))

       rOmegaUniq(i) = rOmegaUniq(i) / (1.*ct(i))
       rOmegaSqUniq(i) = rOmegaSqUniq(i) / (1.*ct(i))
    end do


    ! If needed, fill in empty values at hub and tip
    ! No interpolation here: simply copy from nearest neighbour
    if( addHub==1 ) then
       runiq(1) = 0.

       uuniq(1) = uuniq(2)
       usquniq(1) = usquniq(2)

       rOmegaUniq(1) = rOmegaUniq(2)
       rOmegaSqUniq(1) = rOmegaSqUniq(2)
    end if

    if( addTip==1 ) then
       runiq(numUniq) = turb%radius

       uuniq(numUniq) = uuniq(numUniq-1)
       usquniq(numUniq) = usquniq(numUniq-1)

       rOmegaUniq(numUniq) = rOmegaUniq(numUniq-1)
       rOmegaSqUniq(numUniq) = rOmegaSqUniq(numUniq-1)
    end if


    !    print*,""
    !    print*,"=== UNIQUE ==="

    ! do i=1, numUniq
    !    write(*, "(a, f6.3, a, f6.3, a, f6.3, a, f6.3)") &
    !         "r: ", runiq(i), "  u: ", uuniq(i), &
    !         "  v: ", vuniq(i), "  w: ", wuniq(i)
    ! end do
    ! print*,""

    deallocate( usum )
    deallocate( usqsum )

    deallocate( rOmegaSum )
    deallocate( rOmegaSqSum )

    deallocate( ct )

  end subroutine calculateOrderedVels


  
  subroutine calculateSmoothRelaxInterpVels( turb, runiq, &
       uuniq, romuniq, &
       usquniq, romsquniq )

    type(ModelTurbine) :: turb
    type(InterpolatedNodes) :: inodes
    real, allocatable :: runiq(:), uuniq(:), romuniq(:)
    real, allocatable :: usquniq(:), romsquniq(:)

    integer :: numInt, numUniq, i, j
    real :: rInner, rOuter, rLeft, rRight
    real :: leftdr, rightdr, rdelta, rdeltasum
    real :: r, length
    real :: relax

    real :: uGrad, romGrad
    real :: uArea, romArea
    real :: uLeft, romLeft
    real :: uRight, romRight

    real :: usqGrad, romsqGrad
    real :: usqArea, romsqArea
    real :: usqLeft, romsqLeft
    real :: usqRight, romsqRight

    real :: u, rom
    real :: usq, romsq

    logical :: foundArea
    
    inodes = turb%intNodes
    numInt = size(inodes%r)
    numUniq = size(runiq)

 !   print*,""

    ! Step through each interpolated node
    do i=1, numInt
       r = inodes%r(i)

       ! Find inner and outer radii to average values over
       if(i<=1) then
          rInner = inodes%r(1)
       else 
          rInner = inodes%r(i-1)
       end if

       if(i>=numInt) then
          rOuter = inodes%r(numInt)
       else
          rOuter = inodes%r(i+1)
       end if

       length = rOuter-rInner

!       write(*, "(A, F6.3, A, F6.3, A, F6.3)") "r:", r, &
!            "   rInner: ", rInner, "   rOuter: ", rOuter

       ! Find total area
       uArea = 0
       romArea = 0

       usqArea = 0
       romsqArea = 0

       rdeltasum=0

       ! Integration range for u, v, w in terms of r
       do j=2, numUniq
          foundArea=.false.


          ! leftmost block
          if( rInner >= runiq(j-1) .and. rInner < runiq(j) ) then
             foundArea = .true.
             rLeft = rInner

             if( rOuter > runiq(j-1) ) then
                if ( rOuter <= runiq(j) ) then
                   rRight = rOuter
                else
                   rRight = runiq(j)
                end if
             end if

          ! middle blocks
          elseif( rInner < runiq(j-1) .and. rOuter > runiq(j) ) then
             foundArea = .true.

             rLeft = runiq(j-1)
             rRight = runiq(j)

          ! rightmost block
          elseif( rOuter >= runiq(j-1) .and. rOuter <= runiq(j) ) then
             foundArea = .true.
             rRight = rOuter

             if( rInner < runiq(j) ) then
                if ( rInner >= runiq(j-1) ) then
                   rLeft = rInner
                else
                   rLeft = runiq(j-1)
                end if
             end if
          end if


          if( foundArea ) then
             rdelta = rRight-rLeft
             leftdr = rLeft-runiq(j-1)
             rightdr = rRight-runiq(j-1)

             uGrad = ( uuniq(j)-uuniq(j-1) )/(runiq(j) - runiq(j-1))
             romGrad = ( romuniq(j)-romuniq(j-1) )/(runiq(j) - runiq(j-1))
             
             uLeft = uuniq(j-1) + uGrad*leftdr
             romLeft = romuniq(j-1) + romGrad*leftdr
             
             uRight = uuniq(j-1) + uGrad*rightdr
             romRight = romuniq(j-1) + romGrad*rightdr

             usqGrad = ( usquniq(j)-usquniq(j-1) )/(runiq(j) - runiq(j-1))
             romsqGrad = ( romsquniq(j)-romsquniq(j-1) )/(runiq(j) - runiq(j-1))
             
             usqLeft = usquniq(j-1) + usqGrad*leftdr
             romsqLeft = romsquniq(j-1) + romsqGrad*leftdr
             
             usqRight = usquniq(j-1) + usqGrad*rightdr
             romsqRight = romsquniq(j-1) + romsqGrad*rightdr

             rdeltasum = rdeltasum + rdelta
!             write(*, "(A, f6.3, A, f6.3, A, f6.3)") "rdelta: ", rdelta, &
!                  "   uLeft: ",uLeft, &
!                  "   uRight: ",uRight

             uArea = uArea + (uLeft+uRight)*0.5*rdelta
             romArea = romArea + (romLeft+romRight)*0.5*rdelta

             usqArea = usqArea + (usqLeft+usqRight)*0.5*rdelta
             romsqArea = romsqArea + (romsqLeft+romsqRight)*0.5*rdelta
          end if
       end do

!       write(*,"(A, f12.8)"),"rdeltasum: ", rdeltasum
       
       u = uArea / length
       rom = romArea / length

       usq = usqArea / length
       romsq = romsqArea / length

       ! Now put into inodes in relaxed form
       relax = turb%turbineRelax

       inodes%u(i) = (1-relax)*u + relax*inodes%uOld(i)
       inodes%rOmega(i) = (1-relax)*rom + relax*inodes%rOmegaOld(i)

       inodes%uOld(i) = inodes%u(i)
       inodes%rOmegaOld(i) = inodes%rOmega(i)

       inodes%usq(i) = (1-relax)*usq + relax*inodes%usqOld(i)
       inodes%rOmegaSq(i) = (1-relax)*romsq + relax*inodes%rOmegaSqOld(i)

       inodes%usqOld(i) = inodes%usq(i)
       inodes%rOmegaSqOld(i) = inodes%rOmegaSqOld(i)

!       print*,""

    end do

    turb%intNodes = inodes

!    print*,""
!    print*,""

  end subroutine calculateSmoothRelaxInterpVels


  ! -----------------------------------------------------------------------------


  subroutine calculateNodalVelFromInterp(turb, node, &
       uInt, vInt, wInt, &
       usqInt, vsqInt, wsqInt )

    type(ModelTurbine) :: turb
    type(TurbineNode) :: node
    type(InterpolatedNodes) :: inodes
    real :: r, dr, intdr, uInt, vInt, wInt, usqInt, vsqInt, wsqInt
    real :: rwInt, rwSqInt, rnode
    real :: w1, w2
    integer :: ix1, ix2

    integer :: i
    logical :: foundR

    ! print*, "calculateNodalVelFromInterp():"

    foundR = .false.

    inodes = turb%intNodes

    do i=1, size(inodes%r)-1

       intdr = inodes%r(i+1) - inodes%r(i)

       if( r>=inodes%r(i) .and. r<=inodes%r(i+1) &
            .and. foundR .eqv. .false.) then
          dr = r-inodes%r(i)

          w1 = 1.-dr/intdr
          w2 = dr/intdr

          ix1 = i
          ix2 = i+1
       end if
    end do


    uInt = w1 * inodes%u(ix1) + w2 * inodes%u(ix2)
    usqInt = w1 * inodes%usq(ix1) + w2 * inodes%usq(ix2)

    ! v, w components extrapolated from r omega
    rwInt = w1 * inodes%rOmega(ix1) + w2 * inodes%rOmega(ix2)
    rwsqInt = w1 * inodes%rOmegaSq(ix1) + w2 * inodes%rOmegaSq(ix2)

    rnode = sqrt(node%y**2. + node%z**2.)
    if(rnode .gt. verySmall ) then

	   		vInt = node%z * rwInt / rnode
	   		wInt = -node%y * rwsqInt / rnode

	   		vsqInt = -node%z * rwInt / rnode
	   		wsqInt = node%y * rwsqInt / rnode
	  else
	      vInt = 0.
	      wInt = 0.
	      vsqInt = 0.
	      wsqInt = 0.
	  end if


    ! print*, "end calculateNodalVelFromInterp():"

  end subroutine calculateNodalVelFromInterp


  ! -----------------------------------------------------------------------------


end module interpnodes
