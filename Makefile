# Top level Makefile (included by those in subdirectories)

# These seem to have no effect -- perhaps make sets them (to cc and g++) itself?
CC       ?= g++ -g
CXX      ?= g++ -g
# PICFLAG -- set to -fPIC to also have the option to build shared libraries (position independent code)
# This is only needed on linux.
OPTFLAG  ?=  $(PICFLAG) -g -O0 -DNEW_OP_ORDER -DUSE_PARI_FACTORING

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

all:
	cd procs && make tests progs install
	cd qcurves && make tests progs install
	cd qrank && make tests progs install
	cd g0n && make tests progs install

dylib:
	cd procs && make install_dylib
	cd qcurves && make install_dylib
	cd qrank && make install_dylib
	cd g0n && make install_dylib

so:
	cd procs && make install_so
	cd qcurves && make install_so
	cd qrank && make install_so
	cd g0n && make install_so

clean:
	cd procs && make clean
	cd qcurves && make clean
	cd qrank && make clean
	cd g0n && make clean
	cd lib && /bin/rm -f *.a *.so *.dylib 
	cd include && /bin/rm -f *.h

veryclean:
	cd procs && make veryclean
	cd qcurves && make veryclean
	cd qrank && make veryclean
	cd g0n && make veryclean
	cd lib && /bin/rm -f *.a *.so *.dylib 
	cd include && /bin/rm -f *.h

show:
	echo $(OPTFLAG)

#################################################################
# Shared object libraries
# Used for creating a shared library on OS X.
DYN_OPTS = -dynamiclib
# Used for creating a shared library on Linux/Unix
SO_OPTS = -fPIC --shared

shared_dylib: all
	cd procs && make install_dylib
	cd qcurves && make install_dylib
	cd qrank && make install_dylib
	cd g0n && make install_dylib

shared_so: all
	cd procs && make install_so
	cd qcurves && make install_so
	cd qrank && make install_so
	cd g0n && make install_so

