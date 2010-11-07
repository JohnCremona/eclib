# File Makefile: architecture-independent elements of qrank Makefile

include ../Makefile

# Useful targets:
# 
# progs (makes mwrank, tmrank, ratpoint)
# tests (makes various tests)
# install (makes lib & copies lib & includes to ../lib, ../include)
#

.PHONY = clean

ifeq ($(ARITH),NTL_INTS)
LIBNAME = libjcntlint.a
NTLCFLAGS = -DNTL_INTS  -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -L$(LIBDIR) -lgmp  -lpari
LIBFLAG = -L$(LIBDIR) -lcurvesntlint -ljcntlint
CURVELIBNAME = libcurvesntlint.a
CURVELIBFLAG = -L../libs -lcurvesntlint
RANKLIBBASE = librankntlint
RANKLIBFLAG = -L../libs -lrankntlint
OBJ_SUF = _ntl
EXE_SUFFIX = _ntl
else
ifeq ($(ARITH),NTL)
LIBNAME = libjcntl.a
NTLCFLAGS = -DNTL_ALL -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -L$(LIBDIR) -lgmp  -lpari
LIBFLAG = -L$(LIBDIR) -lcurvesntl -ljcntl
CURVELIBNAME = libcurvesntl.a
CURVELIBFLAG = -L../libs -lcurvesntl
RANKLIBBASE = librankntl
RANKLIBFLAG = -L../libs -lrankntl
OBJ_SUF = _n
EXE_SUFFIX = _n
else
NTLCFLAGS = 
NTLLFLAGS = 
LIBNAME = libjc.a
LIBFLAG = -L$(LIBDIR) -lcurves -ljc
CURVELIBNAME = libcurves.a
CURVELIBFLAG = -lcurves
RANKLIBBASE = librank
RANKLIBFLAG = -L../libs -lrank
OBJ_SUF = _
EXE_SUFFIX = _
endif
endif

RANKLIBNAME=$(RANKLIBBASE).a

CFLAGS = -c $(OPTFLAG) -I$(INCDIR) $(NTLCFLAGS) -I$(PARIINCDIR) #-pg
LFLAGS = $(DYN_FLAGS) $(LIBFLAG) -lpari -L$(PARILIBDIR) $(NTLLFLAGS) $(IOLIB) -lpari -lm #-pg

LIBFLAGS = $(LFLAGS)

%$(OBJ_SUF).o:   %.cc
	$(CXX) $(CFLAGS) $< -o $@

# TARGETS

HEADERS = mrank2.h mequiv.h mlocsol.h mglobsol.h msoluble.h mquartic.h mrank1.h qc.h sqfdiv.h version.h minim.h reduce.h transform.h desc2.h bitspace.h options.h GetOpt.h twoadic.h descent.h
headers: $(HEADERS)

DOTCCS = mequiv.cc  mrank1.cc mlocsol.cc mglobsol.cc mquartic.cc mrank2.cc tmrank.cc ratpoint.cc qc.cc mwrank.cc sqfdiv.cc version.cc minim.cc reduce.cc transform.cc desc2.cc tred.cc bitspace.cc GetOpt.cc twoadic.cc descent.cc # d2.cc tequiv.cc tred.cc

sources: $(HEADERS) $(DOTCCS)

DOCS = mwrank.info mwrank.options mwrank.changes mwrank.readme
docs: $(DOCS)

DISTFILES = $(HEADERS) mequiv.cc mrank1.cc mlocsol.cc mglobsol.cc mquartic.cc mrank2.cc tmrank.cc ratpoint.cc d2.cc qc.cc mwrank.cc sqfdiv.cc version.cc minim.cc reduce.cc transform.cc desc2.cc bitspace.cc twoadic.h twoadic.cc descent.cc $(DOCS) tmrank.in

distfiles: $(DISTFILES)

clean:
	/bin/rm -f *$(OBJ_SUF).o *~ $(RANKLIBNAME) *.a *.so *.dylib

OBJS = mrank1$(OBJ_SUF).o mrank2$(OBJ_SUF).o mequiv$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  qc$(OBJ_SUF).o sqfdiv$(OBJ_SUF).o desc2$(OBJ_SUF).o transform$(OBJ_SUF).o minim$(OBJ_SUF).o reduce$(OBJ_SUF).o bitspace$(OBJ_SUF).o GetOpt$(OBJ_SUF).o twoadic$(OBJ_SUF).o descent$(OBJ_SUF).o

objs: $(OBJS)

lib: $(OBJS)
	rm -f $(RANKLIBNAME)
	ar r $(RANKLIBNAME) $(OBJS)
	$(RANLIB) $(RANKLIBNAME)

install_lib: $(HEADERS) lib
	(cd $(INCDIR); rm -f $(HEADERS))
	for f in $(HEADERS); \
	 do cp $${f} $(INCDIR); \
	done
	chmod a+r $(INCDIR)/*
	rm -f $(LIBDIR)/$(RANKLIBNAME)
	cp $(RANKLIBNAME) $(LIBDIR)
	chmod a+rx $(LIBDIR)/$(RANKLIBNAME)

install: lib install_lib

install_docs: $(DOCS) 
	for f in $(DOCS); \
	do cp $${f} $(FTP_DIR)/$${f}; chmod 644 $(FTP_DIR)/$${f}; done

PROGS = mwrank ratpoint
progs: $(PROGS)
STATIC_PROGS = mwrank_sta ratpoint_sta
static_progs: $(STATIC_PROGS)

TESTS = $(PROGS) tmrank tequiv tred
tests: $(TESTS) 

test1: tmrank
	time ./tmrank < tmrank.in

test2: mwrank
	time awk 'NF==6{print $$1,$$2,$$3,$$4,$$5;}' tmrank.in | ./mwrank -v 0 -q -o -S10

alltests: tests test1 test2

veryclean: clean
	(cd $(BINDIR); /bin/rm -f $(TESTS) $(PROGS))

install_progs: progs install_docs
	for f in $(PROGS); \
	do strip $${f}; \
	   cp $${f} $(FTP_SUBDIR)/$${f}$(EXE_SUFFIX); \
	   /bin/rm -f $(FTP_SUBDIR)/$${f}$(EXE_SUFFIX).gz; \
	   gzip $(FTP_SUBDIR)/$${f}$(EXE_SUFFIX); \
	   chmod 755 $(FTP_SUBDIR)/$${f}$(EXE_SUFFIX).gz; \
	   cp $${f} $(PROGSDIR)/$${f}$(EXE_SUFFIX); \
	   chmod 755 $(PROGSDIR)/$${f}$(EXE_SUFFIX); \
	done

dist: distfiles
	(cd ../rankdist; ${MAKE} clean; /bin/rm -f *.cc *.h $(DOCS) tmrank.in)
	(cd ../procs; ${MAKE} dist;)
	(cd ../qcurves; ${MAKE} dist)
	for f in $(HEADERS) $(DOTCCS) tmrank.in $(DOCS); \
          do $(CP) $${f} $(HOME)/curves/rankdist; done
	chmod 644 $(HOME)/curves/rankdist/*
	chmod 755 $(HOME)/curves/rankdist/configure
	chmod 755 $(HOME)/curves/rankdist/configure.ac

check: $(PROGS) $(TESTS)
	rm -f PRIMES t
	./tequiv < tequiv.in > t && diff t tequiv.out
	./tred < tred.in > t 2>/dev/null && diff t tred.out
	./ratpoint < ratpoint.in > t 2>/dev/null && diff t ratpoint$(OBJ_SUF).out
	./mwrank < mwrank.in > t 2>/dev/null && diff t mwrank.out
	rm -f PRIMES t

mwrank: version$(OBJ_SUF).o mwrank$(OBJ_SUF).o descent$(OBJ_SUF).o mrank1$(OBJ_SUF).o mrank2$(OBJ_SUF).o mequiv$(OBJ_SUF).o mquartic$(OBJ_SUF).o \
	mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  qc$(OBJ_SUF).o sqfdiv$(OBJ_SUF).o \
	desc2$(OBJ_SUF).o transform$(OBJ_SUF).o minim$(OBJ_SUF).o reduce$(OBJ_SUF).o bitspace$(OBJ_SUF).o GetOpt$(OBJ_SUF).o \
	twoadic$(OBJ_SUF).o
	chmod og+r *
	$(CXX) -o mwrank version$(OBJ_SUF).o mwrank$(OBJ_SUF).o descent$(OBJ_SUF).o mrank1$(OBJ_SUF).o mrank2$(OBJ_SUF).o mequiv$(OBJ_SUF).o mquartic$(OBJ_SUF).o \
	mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  qc$(OBJ_SUF).o sqfdiv$(OBJ_SUF).o desc2$(OBJ_SUF).o \
	transform$(OBJ_SUF).o minim$(OBJ_SUF).o reduce$(OBJ_SUF).o bitspace$(OBJ_SUF).o GetOpt$(OBJ_SUF).o twoadic$(OBJ_SUF).o \
	$(LFLAGS)

tmrank: version$(OBJ_SUF).o tmrank$(OBJ_SUF).o descent$(OBJ_SUF).o mrank1$(OBJ_SUF).o mrank2$(OBJ_SUF).o mequiv$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o \
	 qc$(OBJ_SUF).o sqfdiv$(OBJ_SUF).o desc2$(OBJ_SUF).o transform$(OBJ_SUF).o minim$(OBJ_SUF).o reduce$(OBJ_SUF).o bitspace$(OBJ_SUF).o twoadic$(OBJ_SUF).o
	$(CXX) -o tmrank version$(OBJ_SUF).o tmrank$(OBJ_SUF).o descent$(OBJ_SUF).o mrank1$(OBJ_SUF).o mrank2$(OBJ_SUF).o mequiv$(OBJ_SUF).o mquartic$(OBJ_SUF).o \
	mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  qc$(OBJ_SUF).o sqfdiv$(OBJ_SUF).o desc2$(OBJ_SUF).o transform$(OBJ_SUF).o minim$(OBJ_SUF).o \
	reduce$(OBJ_SUF).o bitspace$(OBJ_SUF).o twoadic$(OBJ_SUF).o $(LFLAGS)

ratpoint: version$(OBJ_SUF).o ratpoint$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  qc$(OBJ_SUF).o
	$(CXX) -o ratpoint version$(OBJ_SUF).o ratpoint$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o \
	 qc$(OBJ_SUF).o $(LFLAGS)

ratpoint_sta: version$(OBJ_SUF).o ratpoint$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  qc$(OBJ_SUF).o
	$(CXX) -static -o ratpoint_sta version$(OBJ_SUF).o ratpoint$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o \
	 qc$(OBJ_SUF).o $(LFLAGS)

d2: d2$(OBJ_SUF).o version$(OBJ_SUF).o desc2$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o   qc$(OBJ_SUF).o reduce$(OBJ_SUF).o \
	minim$(OBJ_SUF).o transform$(OBJ_SUF).o sqfdiv$(OBJ_SUF).o
	$(CXX) -o d2 d2$(OBJ_SUF).o version$(OBJ_SUF).o desc2$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  \
	qc$(OBJ_SUF).o reduce$(OBJ_SUF).o minim$(OBJ_SUF).o transform$(OBJ_SUF).o sqfdiv$(OBJ_SUF).o $(LFLAGS)

tred: version$(OBJ_SUF).o tred$(OBJ_SUF).o reduce$(OBJ_SUF).o minim$(OBJ_SUF).o transform$(OBJ_SUF).o mquartic$(OBJ_SUF).o mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o \
	 qc$(OBJ_SUF).o
	$(CXX) -o tred version$(OBJ_SUF).o tred$(OBJ_SUF).o reduce$(OBJ_SUF).o minim$(OBJ_SUF).o transform$(OBJ_SUF).o mquartic$(OBJ_SUF).o \
	mlocsol$(OBJ_SUF).o mglobsol$(OBJ_SUF).o  qc$(OBJ_SUF).o $(LFLAGS)

tequiv: version$(OBJ_SUF).o tequiv$(OBJ_SUF).o mquartic$(OBJ_SUF).o mequiv$(OBJ_SUF).o
	$(CXX) -o tequiv version$(OBJ_SUF).o tequiv$(OBJ_SUF).o mquartic$(OBJ_SUF).o mequiv$(OBJ_SUF).o $(LFLAGS)


# DEPENDENCIES

version$(OBJ_SUF).o:  $(HEADERS) $(DOTCCS)
mwrank$(OBJ_SUF).o:  mwrank.cc descent.h mrank1.h mrank2.h mquartic.h version.h options.h
tmrank$(OBJ_SUF).o:  tmrank.cc descent.h mrank1.h mrank2.h mquartic.h version.h
descent$(OBJ_SUF).o: descent.cc descent.h mrank1.h mrank2.h mquartic.h
mrank1$(OBJ_SUF).o: mrank1.cc mrank1.h descent.h mquartic.h mequiv.h msoluble.h qc.h
mrank2$(OBJ_SUF).o: mrank2.cc mrank2.h descent.h mquartic.h msoluble.h sqfdiv.h desc2.h bitspace.h
mequiv$(OBJ_SUF).o: mequiv.cc mquartic.h mequiv.h
mquartic$(OBJ_SUF).o: mquartic.cc mquartic.h

mglobsol$(OBJ_SUF).o: mglobsol.cc mglobsol.h mquartic.h
ratpoint$(OBJ_SUF).o: ratpoint.cc qc.h mlocsol.h mglobsol.h mquartic.h qc.h version.h
desc2$(OBJ_SUF).o: desc2.cc desc2.h transform.h reduce.h mlocsol.h mglobsol.h mquartic.h
d2$(OBJ_SUF).o: d2.cc desc2.h reduce.h qc.h mlocsol.h mglobsol.h mquartic.h qc.h
reduce$(OBJ_SUF).o: reduce.cc reduce.h minim.h mquartic.h
minim$(OBJ_SUF).o: minim.cc minim.h mquartic.h
transform$(OBJ_SUF).o: transform.cc transform.h
tred$(OBJ_SUF).o: tred.cc reduce.h mquartic.h
qc$(OBJ_SUF).o: qc.cc qc.h mquartic.h
tequiv$(OBJ_SUF).o: tequiv.cc mquartic.h mequiv.h

bitspace$(OBJ_SUF).o: bitspace.h bitspace.cc

sqfdiv$(OBJ_SUF).o: sqfdiv.cc sqfdiv.h
sqftest$(OBJ_SUF).o: sqftest.cc sqfdiv.h

testarg$(OBJ_SUF).o: testarg.cc
cubic$(OBJ_SUF).o: cubic.h cubic.cc
rcubic$(OBJ_SUF).o: cubic.h rcubic.cc
lcubic$(OBJ_SUF).o: cubic.h lcubic.cc

twoadic$(OBJ_SUF).o: twoadic.cc twoadic.h mlocsol.h

GetOpt$(OBJ_SUF).o: GetOpt.cc
new$(OBJ_SUF).o: new.C new.H
	$(CXX) -g $(CFLAGS) -DHAVE_ATEXIT new.C

LIBBASE=$(RANKLIBBASE)
include ../Makefile.dynamic


