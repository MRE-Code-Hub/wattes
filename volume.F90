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

module volume
  use turbtypes
  use memory
  use parallel
  use mathsmisc

  implicit none

  private
  public :: collectTurbineNodes, isTurbineNodeless, getGlobalTurbineNodes

contains

  ! --------------------------------------------------------------------------
  ! This routine collects nodes, and associates them with turbines if they
  ! are inside the turbine volume.

  subroutine collectTurbineNodes(numMeshNodes, &
       xNode, yNode, zNode, &
       uNode, vNode, wNode, &
       turbines, updateTurbine )

    ! Passed variables
    integer, intent(in) :: numMeshNodes
    real, dimension(numMeshNodes), intent(in) :: xNode, yNode, zNode
    real, dimension(numMeshNodes), intent(in) :: uNode, vNode, wNode

    type (ModelTurbine), pointer :: turbines(:)
    logical :: updateTurbine

    ! internal vars
    integer i, j, nodesSize, newNodeIndex
    logical insideFlag
    type (TurbineNode) :: newNode
    type (TurbineNode), pointer :: nodesCopy(:)

    turbines(:)%numberNodes = 0
    turbines(:)%numberNativeNodes = 0


    ! print*, "in collectTurbineNodes()"

    ! For each point...
    ! (Makes sense to do it this way: the theory being, no point can be inside
    ! several turbines at once)
    !
    ! Footnote:
    ! We have to search through each node in each element, to insure we map
    ! position nodes to velocity nodes (see ndglno stuff in Fluidity code)
    ! This may duplication of effort (shared nodes between elements), but
    ! there is not really much choice in the matter...
    !
    ! Bootnote to the footnote:
    ! Ach, I'm dumping this -- I see no reason I can't just use nonods;
    ! simulation runs seem to show indentical indicies.


    ! Look for nodes inside each turbine
    do i=1, numMeshNodes
       ! print*,"i=", i

       ! For each turbine
       do j=1, size(turbines)

          ! print*,"turb%omegaTurb:", turbines(j)%omegaTurb
          ! print*,"j=", j

          ! Note the title, 'isPoint...BladeVolume(): this is important.

          insideFlag= &
               isPointInsideBladeVolume( xNode(i), yNode(i), zNode(i), &
               turbines(j) )

          ! If inside turbine, add to turbine node list
          ! (Probably better implemented as a linked list)
          if( insideFlag ) then

             ! print*,"insideFlag=true"

             newNode%indexRef = i
             newNode%x = xNode(i)
             newNode%y = yNode(i)
             newNode%z = zNode(i)

             newNode%u = uNode(i)
             newNode%v = vNode(i)
             newNode%w = wNode(i)

             newNode%source%x = 0.0
             newNode%source%y = 0.0
             newNode%source%z = 0.0

             newNode%absorption%x = 0.0
             newNode%absorption%y = 0.0
             newNode%absorption%z = 0.0

             nodesSize=turbines(j)%numberNodes

             ! If we have not actually allocated memory for nodes before,
             ! do clumsy array resizing

             call resizeNodeArray( turbines(j)%nodes, &
                  turbines(j)%numberNodes, turbines(j)%numberNodes+1 )

             turbines(j)%numberNodes = turbines(j)%numberNodes+1
             newNodeIndex=nodesSize+1

             turbines(j)%nodes(newNodeIndex)=newNode

          end if

       end do
    end do


    ! At this point, this is valid for both serial and parallel 
    turbines(:)%numberNativeNodes = turbines(:)%numberNodes

    ! Parallel merge nodes
    if( isParallelRun() .and. updateTurbine ) then

       call parallelSyncProcNodeArrays(turbines)

       do i=1, size(turbines)
          if( turbines(i)%numberNativeNodes > 0 ) then
                call parallelGatherTurbineNodes( turbines(i), i )
          end if
       end do

    end if

    ! print*, "end collectTurbineNodes()"

  end subroutine collectTurbineNodes



  ! --------------------------------------------------------------------------
  ! Is a point inside a particular turbine?

  function isPointInsideBladeVolume(xtest, ytest, ztest, turbine)
    use turbtypes

    logical isPointInsideBladeVolume, outTheBox
    real, intent(in) :: xtest, ytest, ztest
    type (ModelTurbine) :: turbine

    ! Calculation variables
    real rp, radius
    real transx, transy, transz

    ! Assume point is not inside turbine initially
    isPointInsideBladeVolume=.false.
    outTheBox=.false.

	radius = turbine%radius

    ! Transform and rotate coords to turbine coordinates

    transx = xtest - turbine%x
    transy = ytest - turbine%y
    transz = ztest - turbine%z

	! A quick way of ruling out a lot of points
	if(abs(transx) >radius) then
		outTheBox = .true.
	elseif(abs(transy) >radius) then
		outTheBox = .true.
	elseif(abs(transz) >radius) then
		outTheBox = .true.
	end if

	! If not out the box, then ..
	if(outTheBox .eqv. .false.) then
		call horizRotate( transx, transy, transz, -turbine%orientation )
		call tiltRotate( transx, transy, transz, -turbine%upwardTilt )

		! Radial distance from hub centre
		rp = sqrt(transy**2 + transz**2)

		! Check y is within circle r^2 = y^2 + z^2
		! ie. is point inside blade/axle volume?
		if( abs(transx) <= 0.5*turbine%length .and. rp < radius ) then
		       isPointInsideBladeVolume=.true.
		end if
	end if

  end function IsPointInsideBladeVolume



  ! --------------------------------------------------------------------------


  logical function isTurbineNodeless(turb, n)
    type (ModelTurbine) :: turb
    integer :: n

    if ( getGlobalTurbineNodes(turb,n)==0 ) then
        isTurbineNodeless = .true.
    else
        isTurbineNodeless = .false.
    end if

  end function isTurbineNodeless


  ! --------------------------------------------------------------------------


  integer function getGlobalTurbineNodes(turb, n)
    type (ModelTurbine) :: turb
    integer :: n

    if( isParallelRun() ) then
        getGlobalTurbineNodes = sum(numberProcTurbNodes(:, n))
    else
        getGlobalTurbineNodes = turb%numberNodes
    end if

  end function getGlobalTurbineNodes


  ! --------------------------------------------------------------------------



end module volume
