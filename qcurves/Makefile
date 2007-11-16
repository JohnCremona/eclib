# File Makefile: architecture-independent elements of qcurves Makefile
include ../Makefile

# Useful targets:
#
# install:  builds library, puts it into ../lib and includes in ../include
# tests:    builds various test programs
# progs:    builds various useful programs
#

ifeq ($(ARITH),NTL_INTS)
NTLCFLAGS = -DNTL_INTS  -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -L$(LIBDIR) -lgmp  -lpari
LIBFLAG = -ljcntlint -L$(LIBDIR)
LIBBASE = libcurvesntlint
CURVELIBFLAG = -lcurvesntlint -L.
OBJ_SUF = _ntl
EXE_SUFFIX = _ntl
else
ifeq ($(ARITH),NTL)
NTLCFLAGS = -DNTL_ALL -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -L$(LIBDIR) -lgmp  -lpari
LIBFLAG = -ljcntl -L$(LIBDIR)
LIBBASE = libcurvesntl
CURVELIBFLAG = -lcurvesntl -L.
OBJ_SUF = _n
EXE_SUFFIX = _n
else
NTLCFLAGS = 
NTLLFLAGS = 
LIBFLAG = -ljc -L$(LIBDIR)
LIBBASE = libcurves
CURVELIBFLAG = -lcurves -L.
OBJ_SUF = _
EXE_SUFFIX = _
endif
endif

CURVELIBNAME=$(LIBBASE).a

CFLAGS = -c $(OPTFLAG) -I$(INCDIR) $(NTLCFLAGS) -I$(PARIINCDIR) #-pg
LFLAGS = -lpari -L$(PARILIBDIR) $(LIBFLAG) $(NTLLFLAGS) $(IOLIB) -lm #-pg

LIBFLAGS = $(LFLAGS)

%$(OBJ_SUF).o:   %.cc
	$(CXX) $(CFLAGS) $< -o $@

OBJS = curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o points$(OBJ_SUF).o heights$(OBJ_SUF).o cperiods$(OBJ_SUF).o isogs$(OBJ_SUF).o mwprocs$(OBJ_SUF).o lambda$(OBJ_SUF).o sifter$(OBJ_SUF).o sieve_search$(OBJ_SUF).o htconst$(OBJ_SUF).o egr$(OBJ_SUF).o saturate$(OBJ_SUF).o divpol$(OBJ_SUF).o pointsmod$(OBJ_SUF).o curvemod$(OBJ_SUF).o ffmod$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o

objs: $(OBJS)

HEADERS = curve.h points.h cperiods.h isogs.h reader.h mwprocs.h lambda.h sifter.h sieve_search.h htconst.h egr.h saturate.h divpol.h pointsmod.h curvemod.h ffmod.h tlss.h elog.h
headers: $(HEADERS)

DOTCCS = curve.cc curvedata.cc curvered.cc points.cc cperiods.cc isogs.cc heights.cc mwprocs.cc lambda.cc sifter.cc indep.cc sieve_search.cc htconst.cc egr.cc saturate.cc divpol.cc pointsmod.cc curvemod.cc ffmod.cc tlss.cc elog.cc

PROGCCS =  tcurve.cc tpoints.cc theight.cc tate.cc torsion.cc conductor.cc twist.cc allisog.cc findinf.cc indep.cc tsat.cc tperiods.cc thtconst.cc thtbound.cc tegr.cc telog.cc tdivpol.cc

sources: $(HEADERS) $(DOTCCS) $(PROGCCS)

TESTS = tcurve tpoints theight thtconst tegr telog tdivpol tperiods tsat tsat2 tsat3 # thtbound
tests: $(TESTS)

PROGS = tate torsion conductor twist allisog findinf indep 
STATIC_PROGS = tate_sta torsion_sta conductor_sta twist_sta allisog_sta findinf_sta indep_sta
progs: $(PROGS)
static_progs: $(STATIC_PROGS)

lib: $(CURVELIBNAME)

$(CURVELIBNAME): $(OBJS)
	rm -f $(CURVELIBNAME)
	ar r $(CURVELIBNAME) $(OBJS)
	$(RANLIB) $(CURVELIBNAME)

install_lib: $(HEADERS) lib
	(cd $(INCDIR); rm -f $(HEADERS))
	for f in $(HEADERS); \
	 do cp $${f} $(INCDIR); \
	done
	rm -f $(LIBDIR)/$(CURVELIBNAME)
	cp $(CURVELIBNAME) $(LIBDIR)
	chmod a+rx $(LIBDIR)/$(CURVELIBNAME)

install_progs: progs #static_progs
	for f in $(PROGS); \
	do strip $${f}; \
	   cp $${f} $(PROGSDIR)/$${f}${EXE_SUFFIX}; \
	   chmod a+rx $(PROGSDIR)/$${f}${EXE_SUFFIX}; \
	   cp $${f} $(FTP_SUBDIR)/$${f}${EXE_SUFFIX}; \
	   /bin/rm -f $(FTP_SUBDIR)/$${f}${EXE_SUFFIX}.gz; \
	   gzip $(FTP_SUBDIR)/$${f}${EXE_SUFFIX}; \
	   chmod 755 $(FTP_SUBDIR)/$${f}${EXE_SUFFIX}.gz; \
	done

install: install_lib

DISTFILES = curve.h points.h reader.h curve.cc curvedata.cc curvered.cc \
            heights.cc points.cc mwprocs.h mwprocs.cc  \
	    isogs.h isogs.cc cperiods.h cperiods.cc \
            findinf.cc tate.cc conductor.cc torsion.cc twist.cc \
            allisog.cc sieve_search.h sieve_search.cc \
            indep.cc sifter.h sifter.cc htconst.h htconst.cc \
	    saturate.h saturate.cc divpol.h divpol.cc tlss.h tlss.cc \
	    pointsmod.h pointsmod.cc elog.h elog.cc \
	    curvemod.h curvemod.cc ffmod.h ffmod.cc egr.h egr.cc \
            lambda.h lambda.cc
distfiles: $(DISTFILES)

dist: distfiles
	$(CP) $(DISTFILES) $(HOME)/curves/rankdist
	chmod 644 $(HOME)/curves/rankdist/*

g0ndist: distfiles
	$(CP) $(DISTFILES) $(HOME)/curves/g0ndist
	chmod 644 $(HOME)/curves/g0ndist/*

clean:
	/bin/rm -f *~ *$(OBJ_SUF).o  $(TESTS) $(CURVELIBNAME) *_sta *.a *.so *.dylib

veryclean: clean
	$(CLEAN)
	/bin/rm -f $(PROGS)

shar: sources
	rm -f qcurves.shar qcurves.shar.Z
	shar -i files -o qcurves.shar
	compress qcurves.shar

check:  $(TESTS) $(PROGS)
	rm -f PRIMES
	./tcurve < tcurve.in > t && diff t tcurve.out 
	./tpoints > t && diff t tpoints.out 
	./theight < theight.in > t && diff t theight.out
	./thtconst < tcurve.in > t && diff t thtconst.out
	./tegr < tegr.in > t && diff t tegr.out
	./telog > t && diff t telog.out
	./tdivpol > t && diff t tdivpol.out
	./tperiods < tperiods.in >& t && diff t tperiods.out
	./tsat < tsat.in >& t && diff t tsat.out
	./tsat2 < tsat2.in >& t && diff t tsat2.out
	./tsat3 < tsat2.in >& t && diff t tsat3.out
	./conductor < conductor.in > t 2>/dev/null && diff t conductor.out
	./tate < tate.in > t 2>/dev/null && diff t tate.out
	./tate < tate.in > t 2>/dev/null && diff t tate.out
	./twist < twist.in > t 2>/dev/null && diff t twist.out
	./allisog < allisog.in > t 2>/dev/null && diff t allisog.out
	./findinf < findinf.in > t 2>/dev/null && diff t findinf.out
	./indep < indep.in > t 2>/dev/null && diff t indep.out


tcurve: tcurve$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o 
	$(CXX)  -o $(BINDIR)/tcurve tcurve$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o  $(LFLAGS)

tperiods: tperiods$(OBJ_SUF).o cperiods$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/tperiods tperiods$(OBJ_SUF).o cperiods$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o \
		curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o $(LFLAGS)

qperiods: qperiods$(OBJ_SUF).o
	$(CXX)  -o $(BINDIR)/qperiods qperiods$(OBJ_SUF).o  $(LFLAGS)

qpemu: qpemu$(OBJ_SUF).o cperiods$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o points$(OBJ_SUF).o
	$(CXX)  -o $(BINDIR)/qpemu qpemu$(OBJ_SUF).o cperiods$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o points$(OBJ_SUF).o $(LFLAGS)

tpoints: tpoints$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/tpoints tpoints$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

telog: telog$(OBJ_SUF).o elog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/telog telog$(OBJ_SUF).o elog$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

theight: theight$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/theight theight$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)


twist: twist$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/twist twist$(OBJ_SUF).o $(CURVELIBFLAG)  $(LFLAGS)

twist_sta: twist$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX) -static -o $(BINDIR)/twist_sta twist$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)


torsion: torsion$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/torsion torsion$(OBJ_SUF).o $(CURVELIBFLAG)  $(LFLAGS)

torsion_sta: torsion$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX) -static -o $(BINDIR)/torsion_sta torsion$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

tate: tate$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/tate tate$(OBJ_SUF).o  $(CURVELIBFLAG)  $(LFLAGS)

tate_sta: tate$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -static -o $(BINDIR)/tate_sta tate$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

conductor: conductor$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/conductor conductor$(OBJ_SUF).o   $(CURVELIBFLAG)  $(LFLAGS)

conductor_sta: conductor$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -static -o $(BINDIR)/conductor_sta conductor$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

csearch: csearch$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/csearch csearch$(OBJ_SUF).o   $(CURVELIBFLAG)  $(LFLAGS)

bdw2357: bdw2357$(OBJ_SUF).o   $(CURVELIBNAME)
	$(CXX) -o $(BINDIR)/bdw2357 bdw2357$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

bdw2357a: bdw2357a$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -o $(BINDIR)/bdw2357a bdw2357a$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

bdw2311: bdw2311$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -o $(BINDIR)/bdw2311 bdw2311$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

getjay: getjay$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -o $(BINDIR)/getjay getjay$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

isogsort: isogsort$(OBJ_SUF).o   $(CURVELIBNAME)
	$(CXX) -o $(BINDIR)/isogsort isogsort$(OBJ_SUF).o   $(CURVELIBFLAG) $(LFLAGS)

search: search$(OBJ_SUF).o sieve_search$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/search search$(OBJ_SUF).o sieve_search$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

findinf: findinf$(OBJ_SUF).o sieve_search$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/findinf findinf$(OBJ_SUF).o sieve_search$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

indep: indep$(OBJ_SUF).o sifter$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/indep indep$(OBJ_SUF).o sifter$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

indep_sta: indep$(OBJ_SUF).o sifter$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -static -o $(BINDIR)/indep_sta indep$(OBJ_SUF).o sifter$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

ellisog: ellisog$(OBJ_SUF).o   $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/ellisog ellisog$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

allisog: allisog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/allisog allisog$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

allisog_sta: allisog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX) -static -o $(BINDIR)/allisog_sta allisog$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

t2isog$(OBJ_SUF).o: t2isog.cc isogs.h
t2isog: t2isog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/t2isog t2isog$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

thtconst: thtconst$(OBJ_SUF).o htconst$(OBJ_SUF).o points$(OBJ_SUF).o curve$(OBJ_SUF).o cperiods$(OBJ_SUF).o curvedata$(OBJ_SUF).o	curvered$(OBJ_SUF).o egr$(OBJ_SUF).o heights$(OBJ_SUF).o mwprocs$(OBJ_SUF).o sieve_search$(OBJ_SUF).o saturate$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o divpol$(OBJ_SUF).o
	$(CXX) -o thtconst  thtconst$(OBJ_SUF).o htconst$(OBJ_SUF).o points$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o cperiods$(OBJ_SUF).o egr$(OBJ_SUF).o heights$(OBJ_SUF).o mwprocs$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o divpol$(OBJ_SUF).o sieve_search$(OBJ_SUF).o saturate$(OBJ_SUF).o $(LFLAGS)


tshift: tshift$(OBJ_SUF).o htconst$(OBJ_SUF).o points$(OBJ_SUF).o curve$(OBJ_SUF).o cperiods$(OBJ_SUF).o curvedata$(OBJ_SUF).o	curvered$(OBJ_SUF).o egr$(OBJ_SUF).o heights$(OBJ_SUF).o mwprocs$(OBJ_SUF).o sieve_search$(OBJ_SUF).o saturate$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o divpol$(OBJ_SUF).o 
	$(CXX) -o tshift  tshift$(OBJ_SUF).o htconst$(OBJ_SUF).o points$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o cperiods$(OBJ_SUF).o egr$(OBJ_SUF).o heights$(OBJ_SUF).o mwprocs$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o divpol$(OBJ_SUF).o sieve_search$(OBJ_SUF).o saturate$(OBJ_SUF).o $(LFLAGS)

htbound: htbound$(OBJ_SUF).o htconst$(OBJ_SUF).o points$(OBJ_SUF).o curve$(OBJ_SUF).o cperiods$(OBJ_SUF).o curvedata$(OBJ_SUF).o	curvered$(OBJ_SUF).o egr$(OBJ_SUF).o heights$(OBJ_SUF).o mwprocs$(OBJ_SUF).o sieve_search$(OBJ_SUF).o saturate$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o divpol$(OBJ_SUF).o 
	$(CXX) -o htbound  htbound$(OBJ_SUF).o htconst$(OBJ_SUF).o points$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o cperiods$(OBJ_SUF).o egr$(OBJ_SUF).o heights$(OBJ_SUF).o mwprocs$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o divpol$(OBJ_SUF).o sieve_search$(OBJ_SUF).o saturate$(OBJ_SUF).o $(LFLAGS)

ratpoints: ratpoints-1.4.c
	gcc ratpoints-1.4.c -Wall -o ratpoints -O2 -march=i386 -Dx86 -fomit-frame-pointer\
                            -lgmp -lgcc -lc -lm

eco_appl: eco_appl$(OBJ_SUF).o
	$(CXX)  -o eco_appl eco_appl$(OBJ_SUF).o  $(LFLAGS)

iso: iso$(OBJ_SUF).o
	$(CXX)  -o iso iso$(OBJ_SUF).o  $(LFLAGS)

# DEPENDENCIES

curve$(OBJ_SUF).o: curve.h curve.cc
curvedata$(OBJ_SUF).o: curve.h curvedata.cc
curvered$(OBJ_SUF).o: curve.h curvered.cc
tcurve$(OBJ_SUF).o: tcurve.cc curve.h
twist$(OBJ_SUF).o: twist.cc curve.h
tate$(OBJ_SUF).o: tate.cc curve.h reader.h
conductor$(OBJ_SUF).o: tate.cc curve.h reader.h
points$(OBJ_SUF).o: points.h points.cc cperiods.h curve.h
heights$(OBJ_SUF).o: heights.cc points.h curve.h
tpoints$(OBJ_SUF).o: tpoints.cc points.h curve.h
theight$(OBJ_SUF).o: theight.cc points.h curve.h
lambda$(OBJ_SUF).o:   lambda.cc points.h curve.h lambda.h
torsion$(OBJ_SUF).o: torsion.cc points.h curve.h reader.h
findinf$(OBJ_SUF).o: findinf.cc mwprocs.h points.h curve.h
mwprocs$(OBJ_SUF).o: mwprocs.cc mwprocs.h points.h curve.h sieve_search.h
cperiods$(OBJ_SUF).o: cperiods.h cperiods.cc curve.h
elog$(OBJ_SUF).o: elog.h elog.cc cperiods.h curve.h
egr$(OBJ_SUF).o: egr.cc egr.h points.h curve.h
qpemu$(OBJ_SUF).o: qpemu.cc cperiods.h points.h curve.h
crlists$(OBJ_SUF).o: crlists.h crlists.cc curve.h
isogs$(OBJ_SUF).o: isogs.h isogs.cc curve.h points.h  cperiods.h
ellisog$(OBJ_SUF).o: ellisog.cc isogs.h curve.h points.h  cperiods.h
allisog$(OBJ_SUF).o: allisog.cc isogs.h curve.h points.h  cperiods.h reader.h
tperiods$(OBJ_SUF).o: tperiods.cc cperiods.h curve.h reader.h
tesign$(OBJ_SUF).o: curve.h tesign.cc
tesign: tesign$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o 
	$(CXX) -o tesign tesign$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o $(LFLAGS)
telog$(OBJ_SUF).o: telog.cc elog.h points.h curve.h
elog$(OBJ_SUF).o: elog.h elog.cc cperiods.h curve.h
egr$(OBJ_SUF).o: egr.cc egr.h points.h curve.h

bdw2357$(OBJ_SUF).o: bdw2357.cc curve.h
bdw2357a$(OBJ_SUF).o: bdw2357a.cc curve.h
bdw2311$(OBJ_SUF).o: bdw2311.cc curve.h
getjay$(OBJ_SUF).o: getjay.cc curve.h
isogsort$(OBJ_SUF).o: isogsort.cc curve.h

indep$(OBJ_SUF).o: indep.cc sifter.h points.h curve.h
sifter$(OBJ_SUF).o: sifter.h sifter.cc points.h curve.h

htconst$(OBJ_SUF).o: htconst.h htconst.cc points.h curve.h cperiods.h elog.h egr.h
thtconst$(OBJ_SUF).o: thtconst.cc htconst.h points.h curve.h cperiods.h
thtbound$(OBJ_SUF).o: thtbound.cc htconst.h points.h curve.h cperiods.h

sieve_search$(OBJ_SUF).o: sieve_search.cc sieve_search.h

divpol$(OBJ_SUF).o: divpol.cc divpol.h points.h
saturate$(OBJ_SUF).o: saturate.cc saturate.h divpol.h points.h pointsmod.h ffmod.h tlss.h
tsat$(OBJ_SUF).o: tsat.cc saturate.h points.h pointsmod.h ffmod.h curvemod.h tlss.h
tsat: tsat$(OBJ_SUF).o saturate$(OBJ_SUF).o divpol$(OBJ_SUF).o pointsmod$(OBJ_SUF).o curvemod$(OBJ_SUF).o ffmod$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/tsat tsat$(OBJ_SUF).o saturate$(OBJ_SUF).o \
	divpol$(OBJ_SUF).o pointsmod$(OBJ_SUF).o  curvemod$(OBJ_SUF).o ffmod$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)
tsat2$(OBJ_SUF).o: tsat2.cc saturate.h divpol.h points.h pointsmod.h ffmod.h curvemod.h
tsat3$(OBJ_SUF).o: tsat3.cc saturate.h divpol.h points.h pointsmod.h ffmod.h curvemod.h htconst.h egr.h
tegr$(OBJ_SUF).o: tegr.cc egr.h saturate.h divpol.h points.h pointsmod.h ffmod.h curvemod.h htconst.h

tsat2: tsat2$(OBJ_SUF).o saturate$(OBJ_SUF).o divpol$(OBJ_SUF).o pointsmod$(OBJ_SUF).o tlss$(OBJ_SUF).o curvemod$(OBJ_SUF).o ffmod$(OBJ_SUF).o htconst$(OBJ_SUF).o elog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/tsat2 tsat2$(OBJ_SUF).o saturate$(OBJ_SUF).o divpol$(OBJ_SUF).o tlss$(OBJ_SUF).o curvemod$(OBJ_SUF).o \
	                           pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o elog$(OBJ_SUF).o	$(CURVELIBFLAG) $(LFLAGS)

tsat3: tsat3$(OBJ_SUF).o saturate$(OBJ_SUF).o divpol$(OBJ_SUF).o pointsmod$(OBJ_SUF).o tlss$(OBJ_SUF).o curvemod$(OBJ_SUF).o ffmod$(OBJ_SUF).o htconst$(OBJ_SUF).o egr$(OBJ_SUF).o elog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/tsat3 tsat3$(OBJ_SUF).o saturate$(OBJ_SUF).o divpol$(OBJ_SUF).o tlss$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o\
	                           ffmod$(OBJ_SUF).o htconst$(OBJ_SUF).o egr$(OBJ_SUF).o elog$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)


test_iso$(OBJ_SUF).o: test_iso.cc curvemod.h pointsmod.h
test_iso: test_iso$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o $(CURVELIBNAME)
	   $(CXX) -o test_iso test_iso$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

tegr: tegr$(OBJ_SUF).o egr$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o $(BINDIR)/tegr tegr$(OBJ_SUF).o egr$(OBJ_SUF).o $(CURVELIBFLAG) $(LFLAGS)

tdivpol$(OBJ_SUF).o: tdivpol.cc divpol.h points.h
tdivpol: tdivpol$(OBJ_SUF).o divpol$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o
	$(CXX)  -o $(BINDIR)/tdivpol tdivpol$(OBJ_SUF).o divpol$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o $(LFLAGS)
curvemod$(OBJ_SUF).o: curvemod.h curvemod.cc curve.h
pointsmod$(OBJ_SUF).o: pointsmod.h pointsmod.cc curve.h curvemod.h
ffmod$(OBJ_SUF).o: ffmod.h ffmod.cc pointsmod.h curve.h curvemod.h
tlss$(OBJ_SUF).o: tlss.h tlss.cc pointsmod.h curve.h ffmod.h curvemod.h


new$(OBJ_SUF).o: new.C new.H
	$(CXX) -g $(CFLAGS) -DHAVE_ATEXIT new.C
#	$(CXX) -g $(CFLAGS) new.C

lllbasis$(OBJ_SUF).o: lllbasis.cc lllbasis.h points.h curve.h

tlllbasis$(OBJ_SUF).o: tlllbasis.cc lllbasis.h
tlllbasis: tlllbasis$(OBJ_SUF).o lllbasis$(OBJ_SUF).o
	$(CXX)  -o $(BINDIR)/tlllbasis tlllbasis$(OBJ_SUF).o lllbasis$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

iso3: iso3.cc
	g++ -g -static -o iso3 iso3.cc -I/usr/local/LiDIA/include \
	-L/local/jec/LiDIA/library/.libs -lLiDIA -lgmp
iso4: iso4.cc
	g++ -g -static -o iso4 iso4.cc -I/usr/local/LiDIA/include \
	-L/local/jec/LiDIA/library/.libs -lLiDIA -lgmp
pcount: pcount.cc
	g++ -g -static -o pcount pcount.cc -I/usr/local/LiDIA/include \
	-L/local/jec/LiDIA/library/.libs -lLiDIA -lgmp

include ../Makefile.dynamic
