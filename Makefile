# Top level Makefile (included by those in subdirectories)

CC ?= g++ -g
CXX ?= g++ -g
OPTFLAG = -O2 -DNEW_OP_ORDER -DLONG_IS_64_BIT -DUSE_PARI_FACTORING
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

# targets:
#
# all:  builds libraries, puts them into ./lib and includes in ./include
#

all:
	cd procs && make install
	cd qcurves && make install
	cd qrank && make install
	cd g0n && make install

clean:
	cd procs && make clean
	cd qcurves && make clean
	cd qrank && make clean
	cd g0n && make clean

show:
	echo $(OPTFLAG)

