# Top level Makefile (included by those in subdirectories)

# These seem to have no effect -- perhaps make sets them (to cc and g++) itself?
CC       ?= g++ -g
CXX      ?= g++ -g
# PICFLAG -- set to -fPIC to also have the option to build shared libraries (position independent code)
# This is only needed on linux.
OPTFLAG  ?=  $(PICFLAG) -g -O2 -DNEW_OP_ORDER -DUSE_PARI_FACTORING

# default settings
PARI_PREFIX ?= /usr/local
NTL_PREFIX ?= /usr/local
PARIINCDIR ?= ${PARI_PREFIX}/include
PARILIBDIR ?= ${PARI_PREFIX}/lib
NTLINCDIR ?= ${NTL_PREFIX}/include
NTLLIBDIR ?= ${NTL_PREFIX}/lib

RANLIB = ranlib
CP = cp -p

# these paths are relative to the current path of the Makefile it is included from
INCDIR = ../include
LIBDIR = ../lib
BINDIR = .



# possible values for ARITH are
#  
#  NTL
#  NTL_INTS
#
# override default with "make ARITH=NTL" or "make ARITH=NTL_INTS"
# or set ARITH as a shell environment variable
#
#
ifndef ARITH
   ARITH = NTL
endif


# targets:
#
# all:  builds libraries, puts them into ./lib and includes in ./include
#

all: dirs
	cd procs && ${MAKE} lib install
	cd qcurves && ${MAKE} lib install
	cd qrank && ${MAKE} lib install
	cd g0n && ${MAKE} lib install

dirs: 
	mkdir -p include
	mkdir -p lib

dylib:
	cd procs && ${MAKE} install_dylib && ${MAKE} tests progs
	cd qcurves && ${MAKE} install_dylib && ${MAKE} tests progs
	cd qrank && ${MAKE} install_dylib && ${MAKE} tests progs
	cd g0n && ${MAKE} install_dylib && ${MAKE} tests progs

so:
	cd procs && ${MAKE} install_so && ${MAKE} tests progs
	cd qcurves && ${MAKE} install_so && ${MAKE} tests progs
	cd qrank && ${MAKE} install_so && ${MAKE} tests progs
	cd g0n && ${MAKE} install_so && ${MAKE} tests progs

dll: all
	cd procs && ${MAKE} install_dll && ${MAKE} tests progs
	cd qcurves && ${MAKE} install_dll && ${MAKE} tests progs
	cd qrank && ${MAKE} install_dll && ${MAKE} tests progs
	cd g0n && ${MAKE} install_dll && ${MAKE} tests progs

clean:
	cd procs && ${MAKE} clean
	cd qcurves && ${MAKE} clean
	cd qrank && ${MAKE} clean
	cd g0n && ${MAKE} clean
	cd lib && /bin/rm -f *.a *.so *.dylib *.def
	cd include && /bin/rm -f *.h

veryclean:
	cd procs && ${MAKE} veryclean
	cd qcurves && ${MAKE} veryclean
	cd qrank && ${MAKE} veryclean
	cd g0n && ${MAKE} veryclean
	cd lib && /bin/rm -f *.a *.so *.dylib *.def
	cd include && /bin/rm -f *.h

show:
	echo $(OPTFLAG)

check: 
	cd procs && ${MAKE} check
	cd qcurves && ${MAKE} check
	cd qrank && ${MAKE} check
	cd g0n && ${MAKE} check

#################################################################
# Shared object libraries
# Used for creating a shared library on OS X.
DYN_OPTS = -dynamiclib
# Used for creating a shared library on Linux/Unix
SO_OPTS = -fPIC --shared
# Used for creating a shared library on Cygwin:
DLL_OPTS = -shared

shared_dylib: all
	cd procs && ${MAKE} install_dylib
	cd qcurves && ${MAKE} install_dylib
	cd qrank && ${MAKE} install_dylib
	cd g0n && ${MAKE} install_dylib

shared_so: all
	cd procs && ${MAKE} install_so
	cd qcurves && ${MAKE} install_so
	cd qrank && ${MAKE} install_so
	cd g0n && ${MAKE} install_so


