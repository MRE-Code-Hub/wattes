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

module absorbsrc
  use turbtypes
  use calcglobals
  use parallel
  use nodeforce

  implicit none

  private
  public :: calculateAbsorptionSourceTerms, updateGlobalAbsorptionSources

contains


  subroutine calculateAbsorptionSourceTerms(turb, t, dt, updateTurbine)
    type (ModelTurbine) turb

    real t, dt
    real r, uRatio
    integer j, numberNodes

    logical updateTurbine

    if( turbineDebug ) print*,"calculateAbsorptionTerms()"

    ! Calculate mean turbine inflow and thus angular velocity

    if ( updateTurbine ) then
       call calculateTurbineGlobals( turb, t, dt )
    end if

    ! calculate forces, and thus absorption/source terms 
    ! for Fluidity
    do j=1, turb%numberNodes
       select case( turb%status )
          
       case (BEM_DISC_MODEL)
          call calculateBemDiscNodeForce(&
               turb%nodes(j), &
               turb, updateTurbine )
          
       case (BEM_LINE_MODEL)
          call calculateBemLineNodeForce(&
               turb%nodes(j), &
               turb, updateTurbine )

       end select
    end do

    if ( (turb%status==BEM_DISC_MODEL .or. &
            turb%status==BEM_LINE_MODEL ) &
         .and. updateTurbine .eqv. .true. ) then
       call calculateBemNodeAverages( turb, dt )
    end if


    !    print*,"----- end calculateAbsorptionTerms()"

  end subroutine calculateAbsorptionSourceTerms


  ! --------------------------------------------------------------------------

  subroutine updateGlobalAbsorptionSources( turbines, &
       xSource, ySource, zSource, &
       numGlobalNodes, &
       bladeDistrib )

    type (ModelTurbine), pointer :: turbines(:)
    real :: xSource(:), ySource(:), zSource(:)
    real :: sourceScale
    real, optional :: bladeDistrib(:)

    integer :: numGlobalNodes

    integer t, n, ix
    type (TurbineNode) :: node

    ! Only called once (finally quashed the disappearing turbines problem)

    if(.not. present(bladeDistrib) ) then
        call resetAllAbsorptionSource (xSource, ySource, zSource, &
                      numGlobalNodes)
    else
        call resetAllAbsorptionSource (xSource, ySource, zSource, &
                      numGlobalNodes, &
                      bladeDistrib )
    end if

    do t=1, size(turbines)
       do n=1, turbines(t)%numberNodes
          
          node = turbines(t)%nodes(n)
          ix = node%indexRef
          
          ! Occasionally, the CFD fluid density and the real fluid
          ! density will differ. For power and thrust values to be
          ! correct, the CFD density needs to be overridden with
          ! a real density value, and source terms for the CFD
          ! solver need to be scaled back
    
          sourceScale = turbines(t)%fluidDensityScaling
          
          if(ix >= 1 ) then
             xSource(ix) = sourceScale * node%source%x
             ySource(ix) = sourceScale * node%source%y
             zSource(ix) = sourceScale * node%source%z
             
             if(turbines(t)%haveBladeDistrib) then
                bladeDistrib(ix) = node%bladeDistrib
                ! print*,"bladeDistrib(",ix,"):", bladeDistrib(ix)
             end if
             
          end if
          
       end do
    end do
    
  end subroutine updateGlobalAbsorptionSources



  ! --------------------------------------------------------------------------
  ! Any nodes that drift out of turbine volume must have zero absorption

  subroutine resetAllAbsorptionSource(xSource, ySource, zSource, &
                        arrSize, &
                        bladeDistrib )

    real :: xSource(:),  ySource(:),  zSource(:)
    real, optional :: bladeDistrib(:)
    integer :: arrSize

    xSource(1:arrSize)=0
    ySource(1:arrSize)=0
    zSource(1:arrSize)=0


    if( present(bladeDistrib) ) then
        bladeDistrib(1:arrSize) = 0
    end if

  end subroutine resetAllAbsorptionSource

end module absorbsrc
