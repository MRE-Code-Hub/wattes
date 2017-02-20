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

! program to test unresolved dependencies

program testprog
  implicit none

  integer, parameter :: nx=250, ny=250, nz=30
  integer, parameter :: ndim=(nx*ny*nz)
  real, parameter :: xsize=1000, ysize=1000, zsize=150
  integer numXNodes, numUNodes, numElements, numLocalNodes

  real :: density

  real, pointer :: xNode(:), yNode(:), zNode(:)
  real, pointer :: uNode(:), vNode(:), wNode(:)

  real, pointer :: xAbsorb(:), yAbsorb(:), zAbsorb(:)
  real, pointer :: xSource(:), ySource(:), zSource(:)

  real t, dt, ltime, x,y,z, dx, dy, dz
  integer i, j, k, n

  integer currentNonLinearIteration, numIterations
  logical remesh

  integer ierr, mpiID
  character*256 logfile, mpich

  integer ncpus
  real offsetx, largedx

  include 'mpif.h'

  allocate(xNode(ndim))
  allocate(yNode(ndim))
  allocate(zNode(ndim))

  allocate(uNode(ndim))
  allocate(vNode(ndim))
  allocate(wNode(ndim))

  allocate(xAbsorb(ndim))
  allocate(yAbsorb(ndim))
  allocate(zAbsorb(ndim))

  allocate(xSource(ndim))
  allocate(ySource(ndim))
  allocate(zSource(ndim))
  

  print*,"begin testprog"

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, ncpus, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpiID, ierr)


  largedx = 1000./ncpus
  offsetx = largedx * mpiID


  write( mpich, '(I0)') mpiID

  logfile='output-'//trim(mpich)//'.log'

  open(678, file=trim(logfile), status='replace')


  numXNodes=ndim
  numUNodes=ndim

  dx=xsize/(nx-1)
  dy=ysize/(ny-1)
  dz=zsize/(nz-1)


  ! Initialise the arrays
  n=1
  do k=0, nz-1
     do j=0, ny-1
        do i=0, nx-1
           x= i*dx + offsetx
           y= j*dy
           z= k*dz

           !           tprint  x,y,z

           xNode(n)=x
           yNode(n)=y
           zNode(n)=z
           uNode(n)=5.
           vNode(n)=0.
           wNode(n)=0.

           n=n+1
        end do
     end do
  end do

  density=1.227


  t=0
  dt=1
  remesh=.false.
  currentNonLinearIteration=1

  ! Last time
  ltime = 100000
  ! do some pretend loops



  do i=1, int(ltime/dt)+1
     do j=1, 2
        call turbineFarmFluidityInterface( &
             numXNodes, numUNodes, &
             xNode, yNode, zNode, &
             uNode, vNode, wNode, &
             density, &
             xAbsorb, yAbsorb, zAbsorb, &
             xSource, ySource, zSource, &
             t, dt, &
             1, 2, &
             remesh )
     end do
     t=t+dt
  end do


  call MPI_Finalize()

  deallocate(xNode)
  deallocate(yNode)
  deallocate(zNode)

  deallocate(uNode)
  deallocate(vNode)
  deallocate(wNode)

  deallocate(xAbsorb)
  deallocate(yAbsorb)
  deallocate(zAbsorb)

  deallocate(xSource)
  deallocate(ySource)
  deallocate(zSource)
  
  print*, "end testprog"

end program testprog
