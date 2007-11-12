#Makefile for QCURVES directory

ifndef CC
CC = g++ -static
endif
ifndef CXX
CXX = g++ -static
endif

ALPHA_OPTFLAG = -g -O2 -mcpu=ev6 -Wa,-arch,ev6
OPTFLAG = -g -O3 -Wall# -DDEBUG

INCDIR = ../include
LIBDIR = ../lib
BINDIR = .
PROGSDIR = $(HOME)/bin
FTP_DIR = $(HOME)/.html/ftp/progs
FTP_SUBDIR = $(FTP_DIR)/linux-intel
ifeq ($(HOSTNAME),selmer)
	FTP_SUBDIR = $(FTP_DIR)/alpha
endif
ifeq ($(HOSTNAME),vibri)
	FTP_SUBDIR = $(FTP_DIR)/linux-intel
	PROGSDIR = $(HOME)/linbin
endif


LIDIALIBDIR = /usr/local/LiDIA/lib
LIDIAINCDIR = /usr/local/LiDIA/include

NTLINCDIR = /usr/local/NTL/include
NTLLIBDIR = /usr/local/NTL/lib

CLEAN = rcsclean
RANLIB = ranlib
CP = cp -p

%:: RCS/%,v
	$(CO) $(COFLAGS) $< $@

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
