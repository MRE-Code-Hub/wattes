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

  ! --------------------------------------------------------------------------

module transform
    use turbtypes
    use mathsmisc

    implicit none

contains

  subroutine transformToTurbineCoords(turb)

    type (ModelTurbine) :: turb
    type (TurbineNode) :: node

    integer i
    real transx, transy, transz
    real rotatx, rotaty, rotatz
    real tu, tv, tw
    real rad

    ! print*,"transformToTurbineCoords()"

    rad=turb%orientation
    
    do i=1, turb%numberNodes
       node=turb%nodes(i)
       
       ! translate
       transx=node%x-turb%x
       transy=node%y-turb%y
       transz=node%z-turb%z
       
       ! rotate
       call horizRotate(transx, transy, transz, -rad)
       call tiltRotate(transx, transy, transz, -turb%upwardTilt)

       node%x=transx
       node%y=transy
       node%z=transz
       
       tu = node%u
       tv = node%v
       tw = node%w
       
       call horizRotate(tu, tv, tw, -rad)
       call tiltRotate(tu, tv, tw, -turb%upwardTilt)
       
       node%u = tu
       node%v = tv
       node%w = tw
       
       turb%nodes(i)=node
       
    end do

    ! print*,"end transformToTurbineCoords()"

  end subroutine transformToTurbineCoords

end module transform
