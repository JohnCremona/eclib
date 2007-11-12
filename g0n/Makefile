#Makefile for g0n directory

# choose the linear algebra method (default=2):
METHOD = 2

GCC=g++ -g
CC = $(GCC)
OPTFLAG = -g -O2 -DNEW_OP_ORDER -DLONG_IS_64_BIT

INCDIR = ../include
LIBDIR = ../lib

IOLIB = -lstdc++  # for gcc-2.7.2 and above

LIDIAINCDIR = $(HOME)/herainclude
LIDIALIBDIR = $(HOME)/heralib

PARILIBDIR = .
PARIINCDIR = .

NTLINCDIR = /usr/local/sage-2.5.3-64bit-linux-x86_64-Linux/local/include
NTLLIBDIR = /usr/local/sage-2.5.3-64bit-linux-x86_64-Linux/local/lib

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
#

include targets
