
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

testobjs=Ctestprog.o interfaceturb.o 
testcxx=testcxx

moddir=modules

inc=-I./include

debugflags=

#linkflags=$(debugflags) -shared 
# linkflags=-g $(debugflags)

# Install library to $(prefix)/lib
prefix=$(HOME)/local
linkflags=-lgfortran

CXX=mpicxx
FC=mpif90
#CC=...
#CFLAGS=...
# FC=ftn

# Gfortran
CXXFLAGS=$(inc) $(debugflags) -O2 
FFLAGS=$(inc) $(debugflags) -fPIC -O2 -fbacktrace -J$(moddir) -fdefault-real-8

# FFLAGS=$(inc) $(debugflags) -fPIC -O2 -J$(moddir) -fdefault-real-8

#FFLAGS=$(inc) $(debugflags) -fPIC -Ofast -J$(moddir) -fdefault-real-8 \
-march=bdver1 -Ofast -mprefer-avx128 -funroll-all-loops \
-fprefetch-loop-arrays --param prefetch-latency=300 \
-minline-all-stringops \
-fno-tree-pre -ftree-vectorize -ffast-math


# Ifort
# FFLAGS=-O3 $(inc) -g -module $(moddir) -vec-report0 -traceback -fPIC

#all: directories $(sharedlib) # $(staticlib)
all: directories $(staticlib)
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

#$(EXEC): $(testobjs)
#	$(CXX) $(linkflags) -o $(XEC) $(testobjs)
#testprog: $(staticlib) $(testobjs) 
#	@echo "Building test program"
#	$(FC) $(FFLAGS) -o $(Ctestprog) $(testobjs) -Llib -lwattes


$(testcxx): $(staticlib) $(testobjs) 
	@echo "Building c++ test program"
	$(CXX) $(CXXFLAGS) -o $(testcxx) $(testobjs) \
	-Llib -lmpi_mpifh -lwattes -lgfortran \
	-lstdc++ -lquadmath -lm


gitclean: svnclean
svnclean: clean
	rm -rf .svn */.svn .git */.git

clean:
	rm -f *.o *.mod $(moddir)/*.mod *~ \#*\# */*~
	rm -f $(testcxx) $(sharedlib) $(staticlib) $(testobjs)
	rm -f test/output test/output*.log test/turbines-output.csv


.SUFFIXES: .f90 .F90 .o .so .a

.f90.o:
	@echo "  FC $<"
	@$(FC) (FFLAGS) -c $<
.F90.o:
	@echo "  FC $<"
	@$(FC) $(FFLAGS) -c $<
.cpp.o:
	@echo "  CXX $<"
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: directories

depends:
	@echo "Making dependencies (requires makedepf90) ..."
	@makedepf90 $(inc) *.F90 > .make-depends


# Include dependencies
include .make-depends
