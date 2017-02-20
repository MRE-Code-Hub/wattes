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

module load
  use turbtypes
  use parallel
  use chordlift

  implicit none

  private
  public :: loadTurbines

  type ConfigKey
     character(len=256) :: name
     integer :: count
     character(len=256) :: type, flags, helpline
  end type ConfigKey

  type KeyReturn
     character(len=256) :: name
     character(len=256) :: type
     integer :: ival
     real :: rval
     real :: rvals(3)
     character(len=256) :: strval
  end type KeyReturn

  ! Array for the configuration keys used to configure the BEM model,
  ! and looked for in the turbine configuration file.

  type(ConfigKey) ::  bemKeyList(34)


contains

  subroutine loadTurbines(turbines, turbinesDisabled, turbinesLoaded)
    type(ModelTurbine), pointer :: turbines(:)

    integer :: turbinesDisabled, turbinesLoaded

    integer :: numberOfTurbines, fileError
    integer :: i, j, k, iTemp, turbineStatusFlag

    real :: rTemp, hubRadiusFraction, bladePitchDegrees
    real :: windDegrees, sigma

    character*256 :: bemfilename, coeffsFilename, strTemp

    type(KeyReturn) :: retkey
    logical, allocatable :: foundKeyFlags(:)
    logical :: allFound, sectionEnd, clonedTurbine


    ! Open turbine spec file - 'Turbines.dat' and read specs
    ! in, up to maxTurbines turbines.
    !
    ! Format is:
    !     statusFlag (0=off, 1=thesis model, 2=BEM model),

    ! ================= BEM MODEL ==============
    ! THIS HAS SINCE BEEN ADDED TO - SEE bemKeyList() below
    !
    !     turbine status 
    !     x, y, z
    !     initial angle of orientation (degrees: converted to radians)
    !     orientRevTime (secs),
    !     radius
    !     hub radius fraction
    !     turbine length
    !     initial blade pitch (degrees)
    !     blade relax time
    !     number of blades
    !     Tip speed ratio
    !     mass per unit area
    !     blade element filename (for radius/chord/twist angle)
    !     lift/drag coeff. filename
    !     cut in u_0
    !     cut out u_0
    !     turbulence opt
    !     maxPower
    !     velocity relax time

    ! These are the configuration keys used to configure the BEM model,
    ! and looked for in the turbine configuration file.

    bemKeyList = (/ConfigKey("clone", 0, "none", "optional", &
         "Clone previous turbine (defined properties will overwrite cloned values)"), &
       ConfigKey("type", 1, "string", "", &
         "Type of model: thesis=thesis model, disc=BEM actuator disc, line=BEM actuator line, off=model off"), &
       ConfigKey("id", 1, "string", "optional", &
         "String identifier for turbine"), &
       ConfigKey("coords", 3, "real", "", &
         "Physical co-ordinations. Format: X,Y,Z"), &
       ConfigKey("orientation", 1, "real", "", &
         "Orientation of turbine front in compass degrees"), &
       ConfigKey("orientation time", 1, "real", "", &
         "Time for turbine to rotate 360 about z-axis. Set to -ve to remain fixed."), &
       ConfigKey("upward tilt", 1, "real", "optional", &
         "Degrees upward tilt of front side of turbine"), &
       ConfigKey("radius", 1, "real", "", &
         "Radius of turbine blades in simulation units (eg. metres)"), &
       ConfigKey("hub fraction", 1, "real", "", &
         "Hub radius as fraction of turbine radius"), &
       ConfigKey("length", 1, "real", "", &
          "Length of turbine model cylinder volume"), &
       ConfigKey("traditional attack opt", 1, "string", "optional", &
         "yes = alpha @ max(cL/cD), no/omitted = alpha @ max(cL/cD)"), &
       ConfigKey("tip loss", 1, "string", "", &
         "yes or no (default)"), &
       ConfigKey("blade pitch", 1, "real", "", &
         "Initial pitch of blades, degrees clockwise from blade plane"), &
       ConfigKey("blade sense", 1, "string", "", &
         "Aerofoil sense of rotatation: clockwise or anticlockwise"), &
       ConfigKey("blade relax time", 1, "real", "", &
         "Relaxation time for blade rotation. Higher means slow rotation"), &
       ConfigKey("blades", 1, "integer", "", &
         "Number of blades"), &
       ConfigKey("first blade angle", 1, "real", "optional", &
         "Starting angle of first blade (degrees)"), &
       ConfigKey("tip speed ratio", 1, "real", "", &
         "Ratio of tip speed to freestream wind speed"), &
       ConfigKey("mass density", 1, "real", "", &
         "Mass per unit area of blades"), &
       ConfigKey("auto twist", 1, "string", "optional", &
         "If set to 'yes', will automatically calculate blade twist angles"), &
       ConfigKey("auto chord", 1, "string", "optional", &
         "If set to 'yes', will automatically calculate chord lengths"), &
       ConfigKey("bem filename", 1, "string", "", &
         "Filename for BEM geometry. if='auto', calculated automagically"), &
       ConfigKey("liftdrag filename", 1, "string", "", &
         "Filename for lift-drag data, eg. CSV w/ angle of attack versus C_L, C_D"), &
       ConfigKey("cut in u0", 1, "real", "", &
         "Cut-in wind speed for model turbine"), &
       ConfigKey("cut out u0", 1, "real", "", &
         "Cut-out wind speed for model turbine"), &
       ConfigKey("max power u0", 1, "real", "", &
         "Wind speed at which maximum power first occurs"), &
       ConfigKey("power efficiency", 1, "real", "optional", &
         "Mechanical -> electrical power efficiency (1.0=100%)"), &
       ConfigKey("turbulence opt", 1, "real", "", &
         "Turbulence at maximum angular velocity"), &
       ConfigKey("max power", 1, "real", "", &
         "Maximum rated power of turbine"), &
       ConfigKey("velocity relax time", 1, "real", "", &
         "Relaxation time for velocity values"), &
       ConfigKey("sigma", 1, "real", "optional", &
         "standard deviation of blade distribution"), &
       ConfigKey("force factor", 1, "real", "optional", &
         "body force scaling factor - if omitted, resorts to default values"), &
       ConfigKey("fluid density scaling", 1, "real", "optional", &
         "Scales the CFD-supplied fluid density (NB. forces are scaled accordingly)"), &
       ConfigKey("output frequency", 1, "real", "optional", &
         "Frequency of outputted data (in 1 / time units)")/)


    print*, "--- loadTurbines()"

    ! Open turbine file
    open(1, file='turbines.dat', status='old', iostat=fileError)

    ! If we've found the turbine file, load it in.
    if(fileError.ne.0) then
       print*,"**** Cannot find turbines.dat ! Disabling turbine model."
       turbinesDisabled=1
       return
    end if


    ! Count the number of turbines by looking for 'turbine'. Yes it's crude,
    ! but if the format's hosed it'll be picked up later. Just for now, just 
    ! for memory allocation purposes.

    ! Initial values

    numberOfTurbines=0
    do i=1, 10000000
       ! Read in a line. If it's the start of a turbine block, count it
       read(1, "(A)", end=101 ) strTemp
       if(trim(strTemp) == "turbine") numberOfTurbines = numberOfTurbines+1
    end do

101 if(numberOfTurbines==0) then
        print*, "error: found no 'turbine' sections in turbines.dat"
        call printFileFormat()
        stop
    end if

    ! Allocate memory, rewind file, and load values in.
    allocate( turbines(numberOfTurbines) )

    ! All the optional keys that have a default value when omitted,
    ! set here.

    turbines(:)%outputPeriod = -1
    turbines(:)%powerEfficiency = 1.0
    turbines(:)%upwardTilt=0.0
    turbines(:)%autoTwist=.false.
    turbines(:)%autoChord=.false.
    turbines(:)%traditionalAttackOpt=.true.
    turbines(:)%tipLoss=.false.
    turbines(:)%bladeFirstAngle=0.0
    turbines(:)%forceFactor=-1.0
    turbines(:)%sigma = -1.0
    turbines(:)%fluidDensityScaling = 1.0

    do i = 1, numberOfTurbines
      write ( strTemp, * ) i
      turbines(i)%id = trim(adjustl(strTemp))
    end do

    rewind(1)
    do i=1, numberOfTurbines
       ! End marker flag
       sectionEnd=.false.

       clonedTurbine=.false.

       ! Read in a line ("turbine")
       read(1, "(A)", end=102) strTemp

       if(trim(strTemp) .ne. "turbine") then
          print*, "**** can't find 'turbine' label for entry", i
          call printFileFormat()
          stop
       end if

       ! Find key values for each key in list.
       allocate(foundKeyFlags(size(bemKeyList)))
       foundKeyFlags(:)=.false.
       do j=1, size(bemKeyList)
          ! read in one line.
          read(1, "(A)") strTemp

          ! If we've reached the end of this turbine section, break the loop
          if( trim(strTemp)=="end turbine" ) then
             sectionEnd=.true.
             exit
          end if

          ! Else we continue and scan for a key
          call scanLineForKey(strTemp, retkey, bemKeyList, foundKeyFlags)

          select case(trim(retkey%name))
          case("clone")
             ! We can't clone a turbine if we have no previous turbine
             if(i==1) then
                print*, "error: 'clone' keyword cannot be used for first turbine entry"
                print*, "(clone copies previous turbine)"
                stop
             else
                ! Clone values from previous turbine, but save 'id' key.
                clonedTurbine = .true.
                strTemp = turbines(i)%id
                turbines(i) = turbines(i-1)
                turbines(i)%id = strTemp
             end if

          case("type")
             select case(trim(retkey%strval))
             case ("thesis")
                turbines(i)%status = THESIS_MODEL
                print*, "error: thesis model is no longer supported"
                stop
             case("disc")
                turbines(i)%status = BEM_DISC_MODEL
             case("line")
                turbines(i)%status = BEM_LINE_MODEL
             case("off")
                turbines(i)%status = 0
             case default
                print*, "error: unrecognised model type ", trim(retkey%strval)
                stop
             end select
             
          case("id")
             ! Is this ID tag already in use?
             if(i>1) then
                do k=1, i-1
                  if( trim(turbines(k)%id) .eq. trim(retkey%strval) ) then
                    print*, "error: turbine", i, "has same id as turbine", k, &
                        " - '", trim(turbines(j)%id), "'"
                    stop
                  end if
                end do
             end if
             turbines(i)%id = trim(retkey%strval)

          case("coords")
             turbines(i)%x = retkey%rvals(1)
             turbines(i)%y = retkey%rvals(2)
             turbines(i)%z = retkey%rvals(3)

          case("orientation")
            print*, "setting orientation"
             turbines(i)%orientation = windToRad(retkey%rval)
             turbines(i)%oldOrientation = windToRad(retkey%rval)
             print*, "orientation:", radToWind(turbines(i)%orientation)
             
          case("orientation time")
             if(retkey%rval < 0) then
                turbines(i)%orientationFixed = .true.
                turbines(i)%orientRevTime = 0
             else
                turbines(i)%orientationFixed = .false.
                turbines(i)%orientRevTime = retkey%rval
             end if

          case("upward tilt")
             turbines(i)%upwardTilt = -degToRad * retKey%rval
             
          case("radius")
             turbines(i)%radius = retkey%rval
             
          case("hub fraction")
             hubRadiusFraction = retkey%rval
             
          case("length")
             turbines(i)%length = retkey%rval
             
          case("traditional attack opt")
            select case(trim(retkey%strval))
                case("yes")
                    turbines(i)%traditionalAttackOpt = .true.
            end select

          case("tip loss")
            select case(trim(retkey%strval))
                case("yes")
                    turbines(i)%tipLoss = .true.
            end select

          case("blade pitch")
             turbines(i)%initialAlpha = retkey%rval * degToRad

          case("blade sense")
            select case(trim(retkey%strval))
                case("clockwise")
                    turbines(i)%bladeSense = CLOCKWISE
                case("anticlockwise")
                    turbines(i)%bladeSense = ANTICLOCKWISE
            end select
             
          case("blade relax time")
             turbines(i)%alphaRelaxTime = retkey%rval
             
          case("blades")
             if(retkey%ival < 1) then
                print*, "error: turbine must have number of blades >0"
                stop
             end if

             turbines(i)%numBlades = retkey%ival
             
          case("first blade angle")
            turbines(i)%bladeFirstAngle=retkey%rval * degToRad

          case("tip speed ratio")
             turbines(i)%tipSpeedRatioMax = retkey%rval
             
          case("mass density")
             turbines(i)%massPerUnitArea = retkey%rval
             
          case("auto twist")
             if(trim(retkey%strval)=="yes") turbines(i)%autoTwist = .true.
             
          case("auto chord")
             if(trim(retkey%strval)=="yes") turbines(i)%autoChord = .true.
             
          case("bem filename")
             bemFilename = trim(retkey%strval)
             
          case("liftdrag filename")
             coeffsFilename = trim(retkey%strval)
             
          case("cut in u0")
             turbines(i)%uCutin = retkey%rval
             
          case("cut out u0")
             turbines(i)%uCutout = retkey%rval

          case("max power u0")
             turbines(i)%u0_MaxPower = retkey%rval
             
          case("turbulence opt")
             turbines(i)%turbulenceOpt = retkey%rval
             
          case("max power")
             turbines(i)%maxPower = retkey%rval

          case("power efficiency")
             turbines(i)%powerEfficiency = retkey%rval
             
          case("velocity relax time")
             turbines(i)%turbineRelaxSecs = retkey%rval
             
          case("force factor")
             turbines(i)%forceFactor = retkey%rval

          case("sigma")
             turbines(i)%sigma = retkey%rval

          case("fluid density scaling")
             turbines(i)%fluidDensityScaling = retkey%rval

          case("output frequency")
             if(retkey%rval > verySmall) then
                turbines(i)%outputPeriod = 1.0/retkey%rval
             end if

          end select
          
       end do

       ! Abort if not clone and we didn't get all the key values we needed.

       if(.not. clonedTurbine) then
           allFound = .true.
		       do j=1, size(bemKeyList)

		          if(foundKeyFlags(j) .eqv. .false. ) then
		             if(scan(trim(bemKeyList(j)%flags), "optional")==0) then
		                allFound = .false.
		             end if
		          end if
		       end do

		       if( allFound .eqv. .false. ) then
		          print*, "Error: missing keys for config file -"

		          do j=1, size(bemKeyList)
		             if(foundKeyFlags(j) .eqv. .false. &
		                  .and. scan(trim(bemKeyList(j)%flags), "optional")<1) then
		                if(trim(bemKeylist(j)%flags)=="") then
		                    print*,"   ", trim(bemKeyList(j)%name)
		                else
				                print*,"   ", trim(bemKeyList(j)%name), &
				                    " <",trim(bemKeyList(j)%flags), ">"
				            end if
		             end if
		          end do

		          call printFileFormat()
		          stop
		       end if
		   end if
       deallocate(foundKeyFlags)

       ! Additional calculations for turbine structure, load geometry and
       ! lift/drag coefficient files

       turbines(i)%hubRadius = hubRadiusFraction * turbines(i)%radius
       turbines(i)%windRampTime=0.

       turbines(i)%orientOmegaDotMax &
            = 2.*pi /( orientFractionTime*turbines(i)%orientRevTime**2. )

       ! Geometry calculations. If we set the filename to 'auto', this will
       ! automagically calculate the values for twist and chord width in
       ! initialiseTurbines(). See init.F90 for more details

       ! Load in blade coefficients
       call loadBladeCoeffs( turbines(i), coeffsFilename )
       call calculateAttackMaxLiftPerf( turbines(i) )

       call loadBladeElements( turbines(i), bemFilename)

       ! Read in the section end marker if we haven't already
       if(.not. sectionEnd) then
          ! Read in a line ("end turbine")
          read(1, "(A)") strTemp
          
          if( trim(strTemp) .ne. "end turbine" ) then
             print*, "**** can't find 'end turbine' label for entry", i
             call printFileFormat()
             stop
          end if
       end if

    end do
102 close(1)

    turbinesDisabled=1
    do i=1, size(turbines)
       if(turbines(i)%status .ne. TURBINE_OFF) &
            turbinesDisabled=0
    end do

    if(turbinesDisabled==1) &
         print*, "All turbines switched off: disabling calculations"

    turbinesLoaded=1

    print*, "--------------------------------------------------------"


  end subroutine loadTurbines




  ! --------------------------------------------------------------------------
  ! Read in blade element radius, chord length, and twist angle from file.
  !
  ! If filename is specified as 'auto' (or autotwist=yes AND autochord=yes are
  ! defined in config file), then calculate values based upon equations for
  ! constant TSR turbines.  If only one of autotwist=yes OR autochord is true,
  ! the program will read in two columns and re-create either the chord or 
  ! twist values.
  !
  ! See Wind Energy Handbook 
  ! (pp.72-73, Section 3.7.2 Optimal design for variable-speed operation &
  ! Section 3.7.3 A practical blade design)

 
  subroutine loadBladeElements(turb, filename)
    implicit none

    type(ModelTurbine) :: turb
    character*256 filename
    real, pointer, dimension(:) :: a, b, c
    integer nRows, i, fileError

    real :: mu, r, dr, lambda, alpha, beta, beta0, chord, maxChord
    real :: lambdamu, lowerfrac, Cl, betaRT, rootTerm
    real, parameter :: twothirds=2./3.


    lambda = turb%tipSpeedRatioMax
    alpha  = turb%attackMaxPerf
    Cl     = turb%optLift

    ! If filename is set to 'auto', set both twist and chord to 'auto'
    if( trim(filename) == "auto" ) then
       turb%autoTwist=.true.
       turb%autoChord=.true.
    end if


    if( turb%autoTwist .and. turb%autoChord ) then

       print*, "automagically generating both blade twist and chord length..."

       nRows = 21
       allocate( a(nRows) )
       allocate( b(nRows) )
       allocate( c(nRows) )

       dr = (turb%radius-turb%hubRadius) / (nRows-1.)
       do i=1, nRows
          mu = ((i-1.)*dr + turb%hubRadius)/turb%radius
          a(i) = mu
       end do
    else
       if((turb%autoTwist.eqv..false.) .and. (turb%autoChord.eqv..false.)) then
          call loadAndOrderNColumns(a, b, c, filename, 3, nRows)

       elseif( turb%autoTwist ) then
          print*, "automagically generating blade twist..."
          call loadAndOrderNColumns(a, b, c, filename, 2, nRows)

       elseif( turb%autoChord ) then
          print*, "automagically generating chord length..."
          call loadAndOrderNColumns(a, b, c, filename, 2, nRows)
       end if

    end if
    
    allocate( turb%bladeElements(nRows) )

    ! Blade twist at R_T (so that twist = 0 at tip)
    betaRT = atan(twothirds / (lambda*(1+twothirds/(lambda**2.))))

    write(*,*) "@ r, twist, chord"
    write(*,*) "----------------------------------------------------------------------"

    do i=1, nRows
       mu = a(i)
       lambdamu = lambda * mu
       maxChord = (2.*pi*mu*turb%radius) / turb%numBlades

       turb%bladeElements(i)%radius = mu * turb%radius

       ! Set blade twist
       if( turb%autoTwist ) then
!          if( a(i) > verySmall .and. a(i) >= turb%hubRadius/turb%radius ) then
          if( a(i) > verySmall ) then
             lowerfrac = lambdamu * ( 1. + twothirds/(lambdamu**2.) )
             beta = atan( twothirds / lowerfrac ) - betaRT
          else
             beta = 0.
          end if
       else
          beta = (b(i) / 360.0) * 2* pi
       end if
       turb%bladeElements(i)%twist = beta

       ! Set chord length
       if( turb%autoChord ) then
!          if( a(i) > verySmall .and. a(i) >= turb%hubRadius/turb%radius ) then
          if( a(i) > verySmall ) then

#ifdef LINEAR_TAPER
             ! Linearly tapering blade
             chord = ( 10./(9.*lambda) ) &
                  * (2 - 1.25 * mu) &
                  * (2*pi/(Cl * lambda *turb%numBlades)) &
                  * turb%radius
#else
             ! Curved taper
             rootTerm = 4./9. + lambdamu**2. * &
                  (1 + 2./(9.* lambdamu**2.))**2.

             chord = 16.*pi/( 9.*turb%numBlades*lambda*Cl*sqrt(rootTerm) ) &
                            * turb%radius
#endif
          else
             chord = 0.
          end if
       else
          if( turb%autoTwist ) then
             chord = b(i) * turb%radius
          else
             chord = c(i) * turb%radius
          end if
       end if

       ! It's doubtful that blades would overlap each other
       if(chord > maxChord) chord=maxChord

       turb%bladeElements(i)%chord = chord

       write(*,"(F6.3, A, F6.3, A, F6.3)") &
            turb%bladeElements(i)%radius, ", ", &
            turb%bladeElements(i)%twist*radToDeg, ", ", &
            turb%bladeElements(i)%chord

    end do

    deallocate(a)
    deallocate(b)
    deallocate(c)

  end subroutine loadBladeElements


  ! --------------------------------------------------------------------------
  ! read in angle, blade lift and drag coefficients
  !

  subroutine loadBladeCoeffs(turb, filename)
    implicit none

    type(ModelTurbine) :: turb
    character*256 filename
    real, pointer :: a(:), b(:), c(:)
    integer nRows, i

    call loadAndOrderNColumns(a, b, c, &
         filename, 3, nRows)

    allocate( turb%bladeCoeffs(nRows) )

    do i=1, nRows
       turb%bladeCoeffs(i)%angle = (a(i) / 360.0) * 2 * pi
       turb%bladeCoeffs(i)%lift = b(i)
       turb%bladeCoeffs(i)%drag = c(i)
    end do

    deallocate(a)
    deallocate(b)
    deallocate(c)

  end subroutine loadBladeCoeffs



  ! --------------------------------------------------------------------------
  ! load, read in to 2-3 arrays and order
  !

  subroutine loadAndOrderNColumns(a, b, c, filename, nCols, nRows)
    implicit none

    real, dimension(:), pointer :: a, b, c
    real :: aTemp, bTemp, cTemp
    character*256 :: filename
    integer :: fileError
    integer :: nCols, nRows, orderColumn, i, j


    if(nCols<2 .and. nCols>3) then
       print*, "ERROR: loadAndOrderNColumns() only supports 2 or 3 columns"
       stop
    end if

    open(5, file=trim(filename), status='old', iostat=fileError)

    if( fileError .ne. 0 ) then
       print*,"**** Error: can't open file '", trim(filename), "'!"
       stop
    else


       do i=1, MAX_BLADE_REFS
          if(nCols==2) then
             read(5, end=201, fmt=*) aTemp, bTemp
          elseif(nCols==3) then
             read(5, end=201, fmt=*) aTemp, bTemp, cTemp
          end if
       end do

201    nRows=i-1        
       rewind(5)

       ! Allocate memory for column arrays and reset values
       allocate( a(nRows) )
       allocate( b(nRows) )
       allocate( c(nRows) )
       a=0.
       b=0.
       c=0.

       ! Read in and order rows based on first column (in array a)
       do i=1, nRows
          if(nCols==2) then
             read(5, end=201, fmt=*) aTemp, bTemp
          elseif(nCols==3) then
             read(5, end=201, fmt=*) aTemp, bTemp, cTemp
          end if

          ! Sort rule: a(i) < a(i+1)
          if( aTemp .lt. a(i) .and. i .ne. 1) then
             do j=nRows-1, i, -1
                a(j+1) = a(j)
                b(j+1) = b(j)
                c(j+1) = c(j)
             end do
          end if

          a(i) = aTemp
          b(i) = bTemp
          c(i) = cTemp
       end do

       close(5)
    end if

  end subroutine loadAndOrderNColumns


  ! --------------------------------------------------------------------------

  subroutine printFileFormat()
    integer :: i, j
    character(len=256) :: line, subline

    type(ConfigKey) :: key

    print*, ""
    print*, "-----------------------------------------------------------------"
    print*, "**** FOR BEM models:"
    print*, ""

    do i=1, size(bemKeyList)
       key = bemKeyList(i)

       subline=""
       write(subline, "(2A)") "<", trim(key%type)

       if(key%count > 1) then
          do j=2, key%count
             write(subline, "(3A)") trim(subline), ", ", trim(key%type)
          end do

       end if
       write(subline, "(2A)") trim(subline), ">"

       if( trim(key%type) == "none" .or. key%count==0 ) then
          write(line, "(2A)") trim(key%name), " <no arguments>"
       else
          write(line, "(3A)") trim(key%name), " = ", trim(subline)
       end if
       print*, trim(line)

       write(line, "(2A)") "   ", trim(key%helpline)
       print*, trim(line)

       if(trim(key%flags) .ne. "") then
          write(line, "(3A)") "   [ ", trim(key%flags)," ]"
          print*, trim(line)
       end if
       
    end do
  end subroutine printFileFormat


  ! --------------------------------------------------------------------------
  ! Scan a single line character string for a config parameter key

  subroutine scanLineForKey( line, keyRet, keyList, foundKeyFlags )
    character (len=256) :: line
    type(KeyReturn) :: keyRet
    type(ConfigKey) :: keyList(:), foundKey
    logical :: foundKeyFlags(:)
    logical :: endSection=.false.

    integer :: i, j, position, foundIx
    integer :: fileError
    character(len=256) :: type, rSplit="", formatString=""
    logical :: foundSign, validKey, singleKeyNoValue

    ! Reset key values
    keyRet  = KeyReturn("", "", 0, 0., (/0.,0.,0./), "")

    ! Hunt for '=' sign
    foundSign = .false.
    singleKeyNoValue = .false.

    ! Scan for an '=', and see if we require it
    position = scan(line, "=", foundSign)

    if(position==0) then
       ! Check we don't have a key that doesn't require a value
       do i=1, size(keyList)
          if( trim(keyList(i)%type)=="none" &
            .and. trim(adjustl(line)) .eq. trim(keyList(i)%name)) then
              singleKeyNoValue = .true.
          end if
       end do

       if(singleKeyNoValue .eqv. .false.) then
          print*, "**** format error: must have 'key=value' on each line"
          call printFileFormat()
          stop
       end if
    end if

    ! If we require no value, just grab line and trim for key.
    if(singleKeyNoValue) then
      keyRet%name = trim(adjustl(line))
    else
      ! Otherwise take from before split
      keyRet%name = trim(adjustl(line(1:position-1)))
    end if

    ! Check against valid key names in key list
    foundIx=0
    do i=1, size(keyList)
       if(trim(keyRet%name) == trim(keyList(i)%name)) foundIx=i
    end do

    if(foundIx==0) then
       print*, "**** Error: invalid key name '", trim(keyRet%name), "'"
       call printFileFormat()
       stop
    end if

    if( foundKeyFlags(foundIx) ) then
       print*, "**** Error: already defined key '", &
            trim(keyList(foundIx)%name), "'"
       stop
    end if

    ! Check against type definition for matching entry in list
    rSplit = adjustl(line(position+1 : len_trim(line)))


    foundKey=keyList(foundIx)
    select case (trim(foundKey%type))
    case("none")
       keyRet%strval = ""
       keyRet%rval = 0.0
       keyRet%ival = 0
       keyRet%rvals(:) = 0.0

    case ("string")
       if(foundKey%count > 1) then
          print*,"Error: unsupported multiple string values for config key"
          stop
       end if
       keyRet%strval = trim(rSplit)

    case("integer")
       if(foundKey%count > 1) then
          print*, "Error: unsupported multiple integer values for config key"
          stop
       end if
       
       read(rSplit, *, iostat=fileError) keyRet%ival
       if(fileError .ne. 0) then
          print*, "Error: can't convert key value '",trim(rSplit), "' to int"
          stop
       end if

    case("real")
       if(foundKey%count==1) then
          read(rSplit, *, iostat=fileError) keyRet%rval

          if(fileError .ne. 0) then
             print*, "Error: can't convert key value '",trim(rSplit), "' to real"
             stop
          end if

       elseif(foundKey%count==3) then
          read(rSplit, *, iostat=fileError) &
               keyRet%rvals(1), &
               keyRet%rvals(2), &
               keyRet%rvals(3)

          if(fileError .ne. 0) then
             print*, "Error: can't convert key value '",trim(rSplit), "' to real(3)"
             stop
          end if
       else
          print*, "Error: unsupported", foundKey%count, "values for key"
          stop
       end if
          
    end select

    ! Make sure to set type of returned key
    foundKeyFlags(foundIx) = .true.
    keyRet%type = trim(foundKey%type)

    
  end subroutine scanLineForKey
    

end module load
