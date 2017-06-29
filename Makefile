# WATTES - Wind And Tidal Turbine Embedded Simulator
# An actuator disc/line turbine model for coupling with CFD software
#
# Copyright (C) 2017 Heriot Watt University and the University of Edinburgh.
#
# Please see the AUTHORS file in the main source directory for a full list
# of copyright holders.
#
#	  Dr. A Creech
#	  Institute of Energy Systems
#	  School of Engineering
#	  University of Edinburgh
#	  
#	  angus_creech@hotmail.com
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
# -----------------------------------------------------------------------------

# This builds the turbine library and associated test program

SHELL=/bin/sh

libdir=lib
sharedlib=$(libdir)/libwattes.so
staticlib=$(libdir)/libwattes.a
linktest=-L$(libdir) -lwattes

libobjs=\
	turbtypes.o memory.o \
	parallel.o checkpoint.o \
	mathsmisc.o \
	bem.o \
	gauss.o turbulence.o \
	interpnodes.o \
	chordlift.o \
	load.o init.o nodeforce.o \
	bladepitch.o \
	volume.o \
	calcglobals.o \
	absorbsrc.o \
	transform.o \
	output.o \
	main.o

testprog=testprog
testobjs=testprog.o interfaceturb.o 

moddir=modules

inc=-I./include

debugflags=

linkflags=$(debugflags) -shared 

# Install library to $(prefix)/lib
prefix=$(HOME)/local


# define your Fortran compiler

# Need to have Fortran-enabled MPI installed
FC=mpif90
# FC=ftn

# Gfortran
FFLAGS=$(inc) $(debugflags) -fPIC -O2 -fbacktrace -J$(moddir) -fdefault-real-8

# FFLAGS=$(inc) $(debugflags) -fPIC -O2 -J$(moddir) -fdefault-real-8

#FFLAGS=$(inc) $(debugflags) -fPIC -Ofast -J$(moddir) -fdefault-real-8 \
-march=bdver1 -Ofast -mprefer-avx128 -funroll-all-loops \
-fprefetch-loop-arrays --param prefetch-latency=300 \
-minline-all-stringops \
-fno-tree-pre -ftree-vectorize -ffast-math


# Ifort
# FFLAGS=-O3 $(inc) -g -module $(moddir) -vec-report0 -traceback -fPIC

all: directories $(sharedlib)
# all: directories $(staticlib)
$(sharedlib): $(libobjs)
	$(CC) $(linkflags) -o $(sharedlib) $(libobjs)

$(staticlib): $(libobjs)
	ar cr $(staticlib) $(libobjs)

flinstall: staticinstall
	@echo "Installing changes to Fluidity source"
	@echo "... copying interfaceturb.o to main/"
	@cp -f interfaceturb.o ~/fluidity/main/

	@echo "--* Done"

interfaceturb.o: turbtypes.o

install: $(sharedlib) interfaceturb.o 
	@echo "Copying $(sharedlib) to $(prefix)/lib"
	@cp -f $(sharedlib) $(prefix)/lib

staticinstall: $(staticlib) interfaceturb.o
	@echo "Copying $(staticlib) to $(prefix)/lib"
	@cp -f $(staticlib) $(prefix)/lib

directories: $(moddir) $(libdir)
$(moddir):
	@mkdir -p $(moddir)
$(libdir):
	@mkdir -p $(libdir)


test: $(sharedlib) $(testobjs) 
	$(FC) $(FFLAGS) -o $(testprog) $(testobjs) -Llib -lturb

modtest:  $(libobjs) $(testobjs)
	$(FC) $(FFLAGS) -o $(testprog) $(testobjs) $(libobjs)

gitclean: svnclean
svnclean: clean
	rm -rf .svn */.svn .git */.git

clean:
	rm -f *.o *.mod $(moddir)/*.mod *~ \#*\# */*~
	rm -f $(testprog) $(sharedlib) $(staticlib)
	rm -f test/output test/output*.log test/turbines-output.csv


.SUFFIXES: .f90 .F90 .o .so .a

.f90.o:
	@echo "  FC $<"
	@$(FC) -c (FFLAGS) $<
.F90.o:
	@echo "  FC $<"
	@$(FC) -c  $(FFLAGS) $<


.PHONY: directories

depends:
	@echo "Making dependencies (requires makedepf90) ..."
	@makedepf90 $(inc) *.F90 > .make-depends


# Include dependencies
include .make-depends
