#Makefile for QRANK directory

ifndef CC
CC = g++
endif
ifndef CXX
CXX = g++
endif

#-DSELMER_ONLY#-DNO_PADIC_FILTERING #-DDEBUG_REDUCE -DDEBUG_REDUCE_QUAD

OPTFLAG = -g -Wall -O3
OPTFLAG2 = -g -O2

INCDIR = ../include
LIBDIR = ../lib
BINDIR = .
PROGSDIR = $(HOME)/bin

LIDIALIBDIR = /usr/local/LiDIA/lib
LIDIAINCDIR = /usr/local/LiDIA/include

NTLINCDIR = /usr/local/NTL/include
NTLLIBDIR = /usr/local/NTL/lib

CLEAN = rcsclean
RANLIB = ranlib
CP = cp -p

%:: RCS/%,v
#  commands to checkout RCS files.  This is the default:
#       test -f $@ || $(CO) $(COFLAGS) $< $@
# and here is my version
	$(CO) $(COFLAGS) $< $@


FTP_DIR = $(HOME)/.html/ftp/progs
FTP_SUBDIR = $(FTP_DIR)/linux-intel
ifeq ($(HOSTNAME),selmer)
	FTP_SUBDIR = $(FTP_DIR)/alpha
endif
ifeq ($(HOSTNAME),vibri)
	FTP_SUBDIR = $(FTP_DIR)/linux-intel
	PROGSDIR = $(HOME)/linbin
endif

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
