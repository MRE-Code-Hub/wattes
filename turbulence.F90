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

module turbulence
  use turbtypes
  use mathsmisc

  implicit none

contains

  ! --------------------------------------------------------------------------
  ! Calculates (fractional) turbulence to add for orbital components
  ! These are passed back, and balanced by backthrust


  subroutine calculateNodeTurbulence( turb, node, uMean, r, tDu, tDv, tDw )

    type(ModelTurbine) :: turb
    type(TurbineNode) :: node
    real :: uMean, r, tDu, tOrb, tRad, tDv, tDw

    real :: pl, dxTol, factor, tiTip, tiInner, uMag, uRadMag, u0,flowDir, tipCr
    real :: a, b, f, g, rDelta, meanDx

    real :: rand, theta, phi
    real :: gaussGen

    if (turbineDebug) then
       print*,""
       print*, "@ calculateNodeTurbulence:"
    end if

    ! Calculate turbulence intensity

    meanDx = turb%meanElementDx
    dxTol = 2 * meanDx

    factor = turb%length / meanDx

     ! distance from leading edge?
     if ( turb%uMean .gt. 0 ) then
        pl= node%x + turb%length/2
     else
        pl = turb%length/2 - node%x
     end if

    tDu = 0
    tDv = 0
    tDw = 0
    tRad = 0
    tOrb = 0
    tiTip = 0
    tiInner=0

    flowDir = turb%uMean / abs(turb%uMean)

    if (turbineDebug) then
       print*, "pl:", pl
       print*,"turb%length-dxTol:", turb%length-dxTol
    end if

!    if ( pl .gt. turb%length-dxTol) then
    if ( pl <= dxTol) then

        select case( turb%status )

            case( BEM_DISC_MODEL, BEM_LINE_MODEL )
                tiTip = turb%turbulenceOpt * &
                    ( (turb%power / turb%maxPower)**0.3333 )

                u0 = turb%uMean

        end select


        tiInner = tiTip / 2.
        tDu = 0

        ! Turbulence for blades
        if(r > turb%hubRadius) then

            if (r > turb%tipStartRadius .or. r <=turb%hubRadius ) then
                uMag = abs(u0) * tiTip * factor
            else
                uMag = abs(u0) * tiInner * factor
            end if

            tDu = uMag * gaussGen()
            tRad = 0.5 * uMag * gaussGen()
            tOrb = 0.0 

            tDv = tOrb * (node%z/r) + tRad * (node%y/r)
            tDw = (- tOrb * (node%y/r)) + tRad * (node%z/r)

        ! Turbulence for hub
        else
            uMag = turb%uHubMean * tiInner * factor

            tDu = uMag * gaussGen()
            tDv = uMag * gaussGen()
            tDw = uMag * gaussGen()
        end if

    end if

    if (turbineDebug) then
       print*, "@ r: ", r
       print*, "@ status: ", turb%status
       print*, "@ uMean:", uMean
       print*, "@ uHubMean:", turb%uHubMean
       print*, "@ u0:", u0
       print*, "@ uMag:", uMag
       print*, "@ tDu/tDv/tDw:", tDu, tDv, tDw
    end if

  end subroutine calculateNodeTurbulence

end module turbulence
