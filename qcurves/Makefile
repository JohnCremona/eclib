# File Makefile: architecture-independent elements of qcurves Makefile
include ../Makefile

# Useful targets:
#
# install:  builds library, puts it into ../lib and includes in ../include
# tests:    builds various test programs
# progs:    builds various useful programs
#

.PHONY = clean

ifeq ($(ARITH),NTL_INTS)
NTLCFLAGS = -DNTL_INTS  -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -L$(LIBDIR) -lgmp  -lpari
LIBFLAG = -ljcntlint -L$(LIBDIR)
LIBBASE = libcurvesntlint
CURVELIBFLAG = -L../lib -lcurvesntlint
OBJ_SUF = _ntl
EXE_SUFFIX = _ntl
else
NTLCFLAGS = -DNTL_ALL -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -L$(LIBDIR) -lgmp  -lpari
LIBFLAG = -ljcntl -L$(LIBDIR)
LIBBASE = libcurvesntl
CURVELIBFLAG = -L../lib -lcurvesntl 
OBJ_SUF = _n
EXE_SUFFIX = _n
endif

CURVELIBNAME=$(LIBBASE).a

CFLAGS = -c $(OPTFLAG) -I$(INCDIR) $(NTLCFLAGS) -I$(PARIINCDIR) #-pg
LFLAGS = $(DYN_FLAGS) $(LIBFLAG) -lpari -L$(PARILIBDIR) $(NTLLFLAGS) $(IOLIB) -lm #-pg

LIBFLAGS = $(LFLAGS)

%$(OBJ_SUF).o:   %.cc
	$(CXX) $(CFLAGS) $< -o $@

OBJS = curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o points$(OBJ_SUF).o heights$(OBJ_SUF).o cperiods$(OBJ_SUF).o isogs$(OBJ_SUF).o mwprocs$(OBJ_SUF).o lambda$(OBJ_SUF).o sifter$(OBJ_SUF).o sieve_search$(OBJ_SUF).o htconst$(OBJ_SUF).o egr$(OBJ_SUF).o saturate$(OBJ_SUF).o divpol$(OBJ_SUF).o pointsmod$(OBJ_SUF).o curvemod$(OBJ_SUF).o ffmod$(OBJ_SUF).o tlss$(OBJ_SUF).o elog$(OBJ_SUF).o getcurve$(OBJ_SUF).o 

HEADERS = curve.h points.h cperiods.h isogs.h reader.h mwprocs.h lambda.h sifter.h sieve_search.h htconst.h egr.h saturate.h divpol.h pointsmod.h curvemod.h ffmod.h tlss.h elog.h getcurve.h 

DOTCCS = curve.cc curvedata.cc curvered.cc points.cc cperiods.cc isogs.cc heights.cc mwprocs.cc lambda.cc sifter.cc indep.cc sieve_search.cc htconst.cc egr.cc saturate.cc divpol.cc pointsmod.cc curvemod.cc ffmod.cc tlss.cc elog.cc getcurve.cc 

PROGCCS =  tcurve.cc tpoints.cc theight.cc tate.cc torsion.cc conductor.cc twist.cc allisog.cc findinf.cc indep.cc tsat.cc tperiods.cc thtconst.cc thtbound.cc tegr.cc telog.cc tdivpol.cc

TESTS = tcurve tpoints theight thtconst tegr telog tdivpol tperiods tsat tsat2 tsat3 # thtbound

PROGS = tate torsion conductor twist allisog findinf indep 

STATIC_PROGS = tate_sta torsion_sta conductor_sta twist_sta allisog_sta findinf_sta indep_sta

all: $(OBJS)
headers: $(HEADERS)
objs: $(OBJS)
sources: $(HEADERS) $(DOTCCS) $(PROGCCS)
tests: $(TESTS)
progs: $(PROGS)
static_progs: $(STATIC_PROGS)
lib: $(CURVELIBNAME)

$(CURVELIBNAME): $(OBJS)
	rm -f $(CURVELIBNAME)
	ar r $(CURVELIBNAME) $(OBJS)
	$(RANLIB) $(CURVELIBNAME)

install_lib: headers lib
	(cd $(INCDIR); rm -f $(HEADERS))
	for f in $(HEADERS); \
	 do cp $${f} $(INCDIR); \
	done
	chmod a+r $(INCDIR)/*.h
	rm -f $(LIBDIR)/$(CURVELIBNAME)
	cp $(CURVELIBNAME) $(LIBDIR)
	chmod a+rx $(LIBDIR)/$(CURVELIBNAME)

install_progs: progs
	for f in $(PROGS); \
	do strip $${f}; \
	   cp $${f} $(BINDIR)/; \
	   chmod a+rx $(BINDIR)/$${f}; \
	done

install: install_lib install_progs

clean:
	rm -f *~ *$(OBJ_SUF).o  $(TESTS) $(CURVELIBNAME) *_sta *.a *.so *.dylib *.dll *.def

veryclean: clean
	rm -f $(PROGS) $(TESTS) *.exe

check:  tests progs
	rm -f PRIMES t
	./tcurve < tcurve.in > t && diff t tcurve.out 
	./tpoints > t && diff t tpoints$(OBJ_SUF).out 
	./theight < theight.in > t && diff t theight$(OBJ_SUF).out
	./thtconst < tcurve.in > t && diff t thtconst$(OBJ_SUF).out
	./tegr < tegr.in > t && diff t tegr.out
	./telog > t && diff t telog$(OBJ_SUF).out
	./tdivpol > t && diff t tdivpol.out
	./tperiods < tperiods$(OBJ_SUF).in > t 2>/dev/null && diff t tperiods$(OBJ_SUF).out
	./tsat < tsat.in > t 2>/dev/null && diff t tsat.out
	./tsat2 < tsat2.in > t 2>/dev/null && diff t tsat2.out
	./tsat3 < tsat3.in > t 2>/dev/null && diff t tsat3$(OBJ_SUF).out
	./conductor < conductor.in > t 2>/dev/null && diff t conductor.out
	./tate < tate.in > t 2>/dev/null && diff t tate.out
	./torsion < torsion$(OBJ_SUF).in > t 2>/dev/null && diff t torsion.out
	./twist < twist.in > t 2>/dev/null && diff t twist.out
	./allisog < allisog$(OBJ_SUF).in > t 2>/dev/null && diff t allisog.out
	./findinf < findinf$(OBJ_SUF).in > t 2>/dev/null && diff t findinf$(OBJ_SUF).out
	./indep < indep.in > t 2>/dev/null && diff t indep.out
	rm -f PRIMES t

tcurve: tcurve$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o 
	$(CXX)  -o ./tcurve tcurve$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o  $(LFLAGS)

tperiods: tperiods$(OBJ_SUF).o cperiods$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o getcurve$(OBJ_SUF).o
	$(CXX) -o ./tperiods tperiods$(OBJ_SUF).o cperiods$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o \
		curvered$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o ffmod$(OBJ_SUF).o getcurve$(OBJ_SUF).o $(LFLAGS)

qperiods: qperiods$(OBJ_SUF).o
	$(CXX)  -o ./qperiods qperiods$(OBJ_SUF).o  $(LFLAGS)

tpoints: tpoints$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o ./tpoints tpoints$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)

telog: telog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./telog telog$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

theight: theight$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o ./theight theight$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)

twist: twist$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./twist twist$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

twist_sta: twist$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX) -static -o ./twist_sta twist$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

torsion: torsion$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./torsion torsion$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

torsion_sta: torsion$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX) -static -o ./torsion_sta torsion$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

tate: tate$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o ./tate tate$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

tate_sta: tate$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -static -o ./tate_sta tate$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

conductor: conductor$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o ./conductor conductor$(OBJ_SUF).o $(CURVELIBNAME)  $(LFLAGS)

conductor_sta: conductor$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -static -o ./conductor_sta conductor$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

csearch: csearch$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o ./csearch csearch$(OBJ_SUF).o   $(CURVELIBNAME)  $(LFLAGS)

bdw2357: bdw2357$(OBJ_SUF).o   $(CURVELIBNAME)
	$(CXX) -o ./bdw2357 bdw2357$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)

bdw2357a: bdw2357a$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -o ./bdw2357a bdw2357a$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

bdw2311: bdw2311$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -o ./bdw2311 bdw2311$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)

getjay: getjay$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -o ./getjay getjay$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)

isogsort: isogsort$(OBJ_SUF).o   $(CURVELIBNAME)
	$(CXX) -o ./isogsort isogsort$(OBJ_SUF).o   $(CURVELIBNAME) $(LFLAGS)

search: search$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./search search$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

findinf: findinf$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./findinf findinf$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

indep: indep$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX)  -o ./indep indep$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

indep_sta: indep$(OBJ_SUF).o  $(CURVELIBNAME)
	$(CXX) -static -o ./indep_sta indep$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

ellisog: ellisog$(OBJ_SUF).o   $(CURVELIBNAME)
	$(CXX)  -o ./ellisog ellisog$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)

allisog: allisog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./allisog allisog$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

allisog_sta: allisog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX) -static -o ./allisog_sta allisog$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)

t2isog$(OBJ_SUF).o: t2isog.cc isogs.h
t2isog: t2isog$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./t2isog t2isog$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)

tsat: tsat$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./tsat tsat$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)
tsat2: tsat2$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./tsat2 tsat2$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)
tsat3: tsat3$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./tsat3 tsat3$(OBJ_SUF).o  $(CURVELIBNAME) $(LFLAGS)


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

test_iso: test_iso$(OBJ_SUF).o curvemod$(OBJ_SUF).o pointsmod$(OBJ_SUF).o $(CURVELIBNAME)
	   $(CXX) -o test_iso test_iso$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

tegr: tegr$(OBJ_SUF).o egr$(OBJ_SUF).o $(CURVELIBNAME)
	$(CXX)  -o ./tegr tegr$(OBJ_SUF).o egr$(OBJ_SUF).o $(CURVELIBNAME) $(LFLAGS)
tdivpol: tdivpol$(OBJ_SUF).o divpol$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o
	$(CXX)  -o ./tdivpol tdivpol$(OBJ_SUF).o divpol$(OBJ_SUF).o curve$(OBJ_SUF).o curvedata$(OBJ_SUF).o $(LFLAGS)


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
getcurve$(OBJ_SUF).o: getcurve.cc getcurve.h

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


tsat2$(OBJ_SUF).o: tsat2.cc saturate.h divpol.h points.h pointsmod.h ffmod.h curvemod.h
tsat3$(OBJ_SUF).o: tsat3.cc saturate.h divpol.h points.h pointsmod.h ffmod.h curvemod.h htconst.h egr.h
tegr$(OBJ_SUF).o: tegr.cc egr.h saturate.h divpol.h points.h pointsmod.h ffmod.h curvemod.h htconst.h



test_iso$(OBJ_SUF).o: test_iso.cc curvemod.h pointsmod.h

tdivpol$(OBJ_SUF).o: tdivpol.cc divpol.h points.h
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
	$(CXX)  -o ./tlllbasis tlllbasis$(OBJ_SUF).o lllbasis$(OBJ_SUF).o  $(CURVELIBFLAG) $(LFLAGS)

iso3: iso3.cc
	g++ -g -static -o iso3 iso3.cc -I/usr/local/LiDIA/include \
	-L/local/jec/LiDIA/library/.libs -lLiDIA -lgmp
iso4: iso4.cc
	g++ -g -static -o iso4 iso4.cc -I/usr/local/LiDIA/include \
	-L/local/jec/LiDIA/library/.libs -lLiDIA -lgmp
pcount: pcount.cc
	g++ -g -static -o pcount pcount.cc -I/usr/local/LiDIA/include \
	-L/local/jec/LiDIA/library/.libs -lLiDIA -lgmp

SONAME=$(LIBBASE).so.0
include ../Makefile.dynamic
