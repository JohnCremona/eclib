CC ?= g++ -g
CXX ?= g++ -g
OPTFLAG = -O2 -DNEW_OP_ORDER -DLONG_IS_64_BIT
OPTFLAG2 = -g -O2

# default settings
PARI_PREFIX ?= /usr/local
NTL_PREFIX ?= /usr/local
PARIINCDIR ?= ${PARI_PREFIX}/include
PARILIBDIR ?= ${PARI_PREFIX}/lib
NTLINCDIR ?= ${NTL_PREFIX}/include
NTLLIBDIR ?= ${NTL_PREFIX}/lib

CLEAN = rcsclean
RANLIB = ranlib
CP = cp -p

# these paths are relative to the current path of the Makefile it is included from
INCDIR = ../include
LIBDIR = ../lib
BINDIR = .

#
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
	ARITH=NTL
endif

all:
	cd procs && make install
	cd qcurves && make install
	cd qrank && make
	cd g0n && make

show:
	echo $(OPTFLAG)

