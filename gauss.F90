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
!
! Gaussian noise generation routines
!
! --------------------------------------------------------------------------


! ==========================================================================
! initialize the RNG's
! and set the color of gaussian noise
! DT is the time step used in whatever process the colored Gaussian noise
!   is used.
! CORTIM is correlation time in the same units as time step DT.
! WHITE=.true. means generate white gaussian noise which happens when
!   CORTIM=0. This flag is used in CGAUSS.
! Here we use the flat distribution RAN1 also taken from Numerical Recipe
! but any other good flat distribution random number generator will do.

subroutine initGaussGen(dt,cortim)
  real cape,dt,l1me2,cgauss
  real cortim,x
  logical white
  common /color/l1me2,cape,white

  if (cortim.eq.0.) then
     white=.true.
     l1me2=-2.000                        !white noise
     cape=0.0
  else
     white=.false.
     cape=dexp(-dt/dble(cortim))
     !parameter needed in CGAUSS
     l1me2=-(dble(1.)-cape*cape)*dble(2./cortim)
  endif
  !        idum=-1
  !        x=ran1(idum)            !initialize flat rng
  x=gaussGen()            !initialize CGAUSS value
  return
end subroutine initGaussGen


! ==========================================================================
! Program to produce exponentially correlated colored (Gaussian) noise.
! based on Fox et al Physical Review A vol.38(1988)5938 and
! modification of GASDEV from Numerical Recipes for Fortran(2nd ed.pg279)

! CAPE is capital E in the article by Fox et. al.
! PREV is the previous value of CGAUSS used in the next iteration
! L1ME2 is the main parameters causing colored noise in Fox et al
!       and represents (lamda*(1-E)^2). Ditto for H in that article.

! routine is illustrated in Double Precision in case it is needed in this
! mode, otherwise all Double Precision variables maybe changed to REAL
! but the corresponding changes must be made to CGAUS0 and the calling
! programs.


real function gaussGen()
  integer iset
  logical white
  real fac,gset,rsq,v1,v2,l1me2,h,cape
  common /color/l1me2,cape,white

  save iset,gset,prev

  data iset/0/
  data prev/0.0d0/

  if (iset.eq.0) then
     !1       v1=2.*ran1(idum)-1.
1    call random_number(v1)
     v1=2.*v1-1
     !        v2=2.*ran1(idum)-1.
     call random_number(v2)
     v2=2.*v2-1
     rsq=v1**2.+v2**2.

     if(rsq.ge.1..or.rsq.eq.0.) goto 1

     !took out sqrt(2) vs eq(28) Fox etal
     fac=dsqrt(l1me2*dlog( dble(rsq) )/rsq)
     gset=v1*fac
     h=v2*fac
     iset=1
  else
     h=gset
     iset=0
  endif

  if(white)then  !please note that the time step vs its sqrt
     gaussGen=h      !in integration is previously set in PARAM
  else
     gaussGen=prev*cape+h
     prev=gaussGen
  endif

  return
end function gaussGen


! --------------------------------------------------------------------------


subroutine initGaussGenInlet(dt,cortim)
  real cape,dt,l1me2,cgauss
  real cortim,x
  logical white
  common /colorInlet/l1me2,cape,white

  if (cortim.eq.0.) then
     white=.true.
     l1me2=-2.000                        !white noise
     cape=0.0
  else
     white=.false.
     cape=dexp(-dt/dble(cortim))
     !parameter needed in CGAUSS
     l1me2=-(dble(1.)-cape*cape)*dble(2./cortim)
  endif
  !        idum=-1
  !        x=ran1(idum)            !initialize flat rng
  x=gaussGen()            !initialize CGAUSS value
  return
end subroutine initGaussGenInlet



! --------------------------------------------------------------------------

real function gaussGenInlet()
  integer iset
  logical white
  real fac,gset,rsq,v1,v2,l1me2,h,cape
  common /colorInlet/l1me2,cape,white

  !  save iset,gset,prev
  data iset/0/
  data prev/0.0d0/

  if (iset.eq.0) then
     !1       v1=2.*ran1(idum)-1.
1    call random_number(v1)
     v1=2.*v1-1
     !        v2=2.*ran1(idum)-1.
     call random_number(v2)
     v2=2.*v2-1
     rsq=v1**2.+v2**2.

     if(rsq.ge.1..or.rsq.eq.0.) goto 1

     !took out sqrt(2) vs eq(28) Fox etal
     fac=dsqrt(l1me2*dlog( dble(rsq) )/rsq)
     gset=v1*fac
     h=v2*fac
     iset=1
  else
     h=gset
     iset=0
  endif

  if(white)then  !please note that the time step vs its sqrt
     gaussGenInlet=h      !in integration is previously set in PARAM
  else
     gaussGenInlet=prev*cape+h
     prev=gaussGenInlet
  endif

  return
end function gaussGenInlet
