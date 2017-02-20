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

module memory
  use turbtypes
  implicit none

contains

  ! --------------------------------------------------------------------------
  ! Subroutine to resize the turbine node arrays dynamically

  subroutine resizeNodeArray(nodesArray, oldSize, newSize)

    integer :: oldSize, newSize, copySize
    integer :: i

    type(TurbineNode), pointer :: nodesArray(:)
    type(TurbineNode), pointer :: nodesCopy(:)

    ! print*,"in resizeNodeArray()"

    if(oldSize.ne.newSize) then
       if(oldSize > newSize) then
          copySize=newSize
       else
          copySize=oldSize
       end if

       if(oldSize == 0 .and. newSize >= 1) then
          allocate(nodesArray(newSize))
       else
          allocate (nodesCopy(newSize))

          do i=1, copySize
             nodesCopy(i)=nodesArray(i)
          end do

          deallocate(nodesArray)
          allocate(nodesArray(newSize))

          nodesArray(:)=nodesCopy(:)
          deallocate(nodesCopy)
       end if
    end if

    ! print*,"end resizeNodeArray()"

  end subroutine resizeNodeArray



  ! --------------------------------------------------------------------------
  ! Needs to be done every iteration, since the nodes inside the turbines
  ! can easily change with hr-adaptivity

  subroutine deallocateTurbineNodes(turbines)

    type(ModelTurbine), pointer :: turbines(:)

    integer i

    if( turbineDebug ) print*, "deallocateTurbineNodes()"

    do i=1, size(turbines)

       if( turbines(i)%numberNodes.gt.0 ) then
          deallocate( turbines(i)%nodes )
          turbines(i)%numberNodes=0
          turbines(i)%numberNativeNodes=0
       end if

    end do


  end subroutine deallocateTurbineNodes

end module memory
