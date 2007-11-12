#Makefile for PROCS directory with test programs

ifndef CC
CC = g++
endif
ifndef CXX
CXX = g++
endif

	# -DCONIC_METHOD=2 # -D DEBUG_LEGENDRE_PARAM  -DTEST_PARAM #-DCHECK_CLAIMS  -DDEBUG_LLL -DTRACE_LLL #-DVERBOSITY=3 \
	  #-DCHECK_LEM3\
          #-DDEBUG_REDUCE \
          #-DWATCH_REDUCTION \
          #-DDEBUG_LEGENDRE \
          #-DCHECK_ALL \
          #-DTEST_PARAM \
          #-DCHECK_RESTRICT
          #-DTRACE_LLL \
	  #-DDEBUG_CUBIC\
          #-DREDUCE_INTERMEDIATES \
          #-DMORDELL_REDUCE \
          #-DTRACE_HOLZER=1 \
          #-DDEBUG_LEM1 \

ALT_OPTFLAG = -g -O2
ALPHA_OPTFLAG = -g -O2 -mcpu=ev6 -Wa,-arch,ev6
OPTFLAG = -g -O3 -Wall -DUSE_PARI_FACTORING

show:
	echo $(OPTFLAG)

INCDIR = ../include
LIBDIR = ../lib
BINDIR = .

IOLIB =  -lstdc++ -lgcc
EFENCE = ../ElectricFence-2.0.5/libefence.a 


LIDIALIBDIR = .
LIDIAINCDIR = .

NTLINCDIR = .
NTLLIBDIR = .

CLEAN = rcsclean
RANLIB = ranlib
CP = cp -p


#
# possible values for ARITH are
#  
#  LIDIA
#  LIDIA_INTS
#  NTL
#  NTL_INTS
#
# override default with "make ARITH=LIDIA" or "make ARITH=LIDIA_INTS"
# or set ARITH as a shell environment variable
#
#
ifndef ARITH
	ARITH=NTL_INTS
endif
#
# Now include the architecture-independent stuff:
# (in fact we could now append the targets file here)

include targets
