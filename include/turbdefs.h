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


! Defaults for various variables - should never need to change this 

! For scaling power and thrust
#define overPredictFactor 1.0

! This outputs all the diagnostics to standard output
! (turn off for large parallel simulations...)
#define turbineStdOut .false.


! Turns out lots of debugging out if .true.
#define turbineDebug .false.

#define tprint write(678,*)

! Number of radial interp points
# define NUMBER_INTERP 11

#define pi 3.1415926535897932384626433832795

! Defaults to this if calculated relax is less 
#define relaxDefault  0.5

! turns on hacks to stability under Fluidity - DEFUNCT 
#define fluidityStabilityHack .false.

! Hub solidity - DEFUNCT 
#define hubSolidityDensity 0.05

! turbine data log file 
#define turbineOutputFile 'turbines-output.csv'


! -----------------------------------------------------------------
!  Below - shouldnt need to touch of any of this, but well, 
!  delve in at your own peril ... 


#define waterViscosity 1.5e-3
#define vonKarman 0.41

#define verySmall 10e-5
#define veryBig 10e10

#define TURBINE_OFF 0
! THESIS_MODEL no longer used, just here for check at start.
#define THESIS_MODEL 1
#define BEM_DISC_MODEL 2
#define BEM_LINE_MODEL 3

#define X_AXIS 1
#define Y_AXIS 2
#define Z_AXIS 3


#define ANTICLOCKWISE 1
#define CLOCKWISE -1


#define MAX_BLADE_REFS 100000

#define orientFractionTime 0.1


! ------------- FUNCTION HACKS ---------------------------------------


#define radToDeg (360/(2*pi))
#define degToRad ((2*pi)/360)

#define radToWind(x) mod((270.-(x*radToDeg)),360.)
#define windToRad(x) (270.-x)*degToRad



! -------- Thesis stuff. Largely redundant ---------------------------

! Function for calculating B effective (see turbine model section in thesis)

#define omegaBFunc(x,y)  (y* ((abs(x/y))**0.333333) )
