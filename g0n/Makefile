# File Makefile: architecture-independent elements of g0n Makefile

include ../Makefile

# Useful targets:
#
# tests: builds various test programs
# progs: builds various useful programs
# install: builds library, copies it and includes to ../lib, ../include

.PHONY = clean

DOT_H = periods.h cusp.h homspace.h newforms.h nflist.h moddata.h oldforms.h symb.h fixc6.h pcprocs.h

DOT_CC = h1bsd.cc newforms.cc symb.cc checkap.cc h1bsdc.cc homspace.cc nfcount.cc h1clist.cc homtest.cc nflist.cc th1nf.cc cnf.cc h1curisog.cc hpcurve.cc nftest.cc tmanin.cc curvesort.cc h1curve.cc oftest.cc cusp.cc h1degphi.cc lf1per.cc oldforms.cc fixc6.cc pcurve.cc tnf2.cc h1first.cc periods.cc vlift.cc mhcount.cc h1test.cc moddata.cc hecketest.cc modtest.cc putoldeigs.cc h1.cc nfhpcurve.cc   

HEADERS = $(DOT_H)
sources: $(DOT_H) $(DOT_CC)

TESTS = modtest homtest hecketest nftest mhcount oftest

PROGS = tmanin h1bsd h1curve pcurve nfhpcurve h1clist checkap h1bsd h1bsdcurisog h1degphi hpcurve nfcount

ifeq ($(ARITH),NTL_INTS)
LIBNAME = libjcntlint.a
NTLCFLAGS = -DNTL_INTS  -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -lgmp -lpari
LIBFLAG = -lcurvesntlint -ljcntlint -L$(LIBDIR)
CURVELIBNAME = libcurvesntlint.a
CURVELIBFLAG = -L../lib -lcurvesntlint 
RANKLIBNAME = librankntlint.a
RANKLIBFLAG = -L../lib -lrankntlint
OBJ_SUF = _ntl
EXE_SUFFIX = _ntl
G0NLIBBASE = libg0nntlint
G0NLIBFLAG = -L../lib -lg0nntlint
else
ifeq ($(ARITH),NTL)
LIBNAME = libjcntl.a
NTLCFLAGS = -DNTL_ALL  -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -lgmp -lpari
LIBFLAG = -lcurvesntl -ljcntl -L$(LIBDIR)
CURVELIBNAME = libcurvesntl.a
CURVELIBFLAG = -L../lib -lcurvesntl.
RANKLIBNAME = librankntl.a
RANKLIBFLAG = -L../lib -lrankntl
OBJ_SUF = _n
EXE_SUFFIX = _n
G0NLIBBASE = libg0nntl
G0NLIBFLAG = -L../lib -lg0nntl
else
LIBNAME = libjc.a
NTLCFLAGS = 
NTLLFLAGS = 
LIBFLAG = -L$(LIBDIR) -lcurves -ljc
CURVELIBNAME = libcurves.a
OBJ_SUF = _
EXE_SUFFIX = _
G0NLIBBASE = libg0n
G0NLIBFLAG = -L../lib -lg0n
RANKLIBNAME = librank.a
RANKLIBFLAG = -L../lib -lrank
endif
endif

G0NLIBNAME = $(G0NLIBBASE).a

OBJS = cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o newforms$(OBJ_SUF).o periods$(OBJ_SUF).o fixc6$(OBJ_SUF).o

objs: $(OBJS)

g0ndist: sources
	(cd ../procs; make g0ndist;)
	(cd ../qcurves; make g0ndist)
	$(CP) $(DOT_H) $(HOME)/curves/g0ndist
	$(CP) $(DOT_CC) $(HOME)/curves/g0ndist
	chmod 644 $(HOME)/curves/g0ndist/*

$(G0NLIBNAME): $(OBJS)
	ar r $(G0NLIBNAME) $(OBJS)
	$(RANLIB) $(G0NLIBNAME)

lib: $(G0NLIBNAME)

install_lib: $(HEADERS) $(G0NLIBNAME)
	(cd $(INCDIR); rm -f $(HEADERS))
	for f in $(HEADERS); \
	 do cp $${f} $(INCDIR); \
	done
	chmod a+r $(INCDIR)/*
	rm -f $(LIBDIR)/$(G0NLIBNAME)
	cp $(G0NLIBNAME) $(LIBDIR)
	chmod a+rx $(LIBDIR)/$(G0NLIBNAME)

install: lib install_lib


CFLAGS = -c $(OPTFLAG)  -I$(INCDIR) $(NTLCFLAGS) -I$(PARIINCDIR) -DMETHOD=2 -DUSE_XSPLIT #-DNF_DIR=\"$(NF_DIR)\"
LFLAGS = $(DYN_FLAGS) $(LIBFLAG) -lpari -L$(PARILIBDIR) $(NTLLFLAGS) $(IOLIB) -lm

LIBFLAGS = $(LFLAGS)

%$(OBJ_SUF).o:   %.cc
	$(CXX) $(CFLAGS) $< -o $@


# TARGETS


progs: $(PROGS)

tests: $(TESTS)

clean:
	/bin/rm -f *~ *$(OBJ_SUF).o *.a *.so *.dylib *.dll* *.def *.exe

veryclean: clean
	/bin/rm -f $(PROGS) $(TESTS)

check:  $(TESTS) $(PROGS)
	rm -f PRIMES t
	./modtest < modtest.in > t 2>/dev/null && diff t modtest.out 
	./homtest < homtest.in > t && diff t homtest.out
	./hecketest < hecketest.in > t 2>/dev/null && diff t hecketest.out 
	./mhcount < mhcount.in > t && diff t  mhcount.out
	rm -rf tmp_nf_dir
	mkdir tmp_nf_dir
	export NF_DIR=tmp_nf_dir && ./tmanin < tmanin.in > t 2>/dev/null && diff t tmanin.out
	export NF_DIR=tmp_nf_dir && ./nftest < nftest.in 2>/dev/null > t && diff t nftest.out
	export NF_DIR=tmp_nf_dir && ./oftest < oftest.in > t && diff t oftest.out
	export NF_DIR=tmp_nf_dir && ./nfhpcurve < nfhpcurve.in > t 2>/dev/null && diff t nfhpcurve.out
	export NF_DIR=tmp_nf_dir && ./h1clist < h1clist.in > t 2>/dev/null && diff t h1clist.out
	export NF_DIR=tmp_nf_dir && ./h1bsdcurisog < h1bsdcurisog.in > t 2>/dev/null && diff t h1bsdcurisog.out
	rm -rf tmp_nf_dir
	rm -f PRIMES t

input: input.cc
	$(CXX) input.cc -o inp

mhcount: mhcount$(OBJ_SUF).o
	$(CXX) -o mhcount mhcount$(OBJ_SUF).o   $(LFLAGS)

g0ngens: g0ngens$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CXX) -o g0ngens g0ngens$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o   $(LFLAGS)

modtest: modtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CXX) -o modtest modtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o   $(LFLAGS)

homtest: homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o
	$(CXX) -o homtest homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o  $(LFLAGS)

oftest: oftest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o homspace$(OBJ_SUF).o newforms$(OBJ_SUF).o  cusp$(OBJ_SUF).o
	$(CXX) -o oftest oftest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o homspace$(OBJ_SUF).o newforms$(OBJ_SUF).o  cusp$(OBJ_SUF).o $(LFLAGS)

hecketest: hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CXX) -o hecketest hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)



gd: gd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CXX) -o gd gd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

gd2: gd2$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CXX) -o gd2 gd2$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

tnfd: tnfd$(OBJ_SUF).o nfd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CXX) -o tnfd tnfd$(OBJ_SUF).o nfd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

nftest: nftest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o \
		newforms$(OBJ_SUF).o
	$(CXX) -o nftest nftest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		oldforms$(OBJ_SUF).o newforms$(OBJ_SUF).o  $(LFLAGS)

ecnf: ecnf$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o \
		newforms$(OBJ_SUF).o
	$(CXX) -o ecnf ecnf$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		oldforms$(OBJ_SUF).o newforms$(OBJ_SUF).o  $(LFLAGS)

nfcount: nfcount$(OBJ_SUF).o
	$(CXX) -o nfcount nfcount$(OBJ_SUF).o $(LFLAGS)

aplist: aplist$(OBJ_SUF).o
	$(CXX) -o aplist aplist$(OBJ_SUF).o $(LFLAGS)

aplist2: aplist2$(OBJ_SUF).o
	$(CXX) -o aplist2 aplist2$(OBJ_SUF).o $(LFLAGS)

findcong: findcong$(OBJ_SUF).o nflist$(OBJ_SUF).o
	$(CXX) -o findcong findcong$(OBJ_SUF).o nflist$(OBJ_SUF).o $(LFLAGS)

find1cong: find1cong$(OBJ_SUF).o nflist$(OBJ_SUF).o
	$(CXX) -o find1cong find1cong$(OBJ_SUF).o nflist$(OBJ_SUF).o $(LFLAGS)

checkap: checkap$(OBJ_SUF).o
	$(CXX) -o checkap checkap$(OBJ_SUF).o -I$(INCDIR) -L$(LIBDIR) $(LFLAGS)

tmanin: tmanin$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o tmanin tmanin$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

getoldeigs: getoldeigs$(OBJ_SUF).o
	$(CXX) -o getoldeigs getoldeigs$(OBJ_SUF).o $(LFLAGS)

putoldeigs: putoldeigs$(OBJ_SUF).o
	$(CXX) -o putoldeigs putoldeigs$(OBJ_SUF).o $(LFLAGS)


fix1290: fix1290$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o fix1290 fix1290$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1test: h1test$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1test h1test$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1first: h1first$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1first h1first$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

mh1first: h1first$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o mh1first h1first$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1: h1$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1 h1$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1bsd: h1bsd$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1bsd h1bsd$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

cong_vec: cong_vec$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o cong_vec cong_vec$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

cong_vec_h1: cong_vec_h1$(OBJ_SUF).o cong_procs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o cong_vec_h1 cong_vec_h1$(OBJ_SUF).o cong_procs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

cong_vec_x: cong_vec_x$(OBJ_SUF).o cong_procs$(OBJ_SUF).o  $(G0NLIBNAME)
	$(CXX) -o cong_vec_x cong_vec_x$(OBJ_SUF).o cong_procs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

lf1per: lf1per$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o lf1per lf1per$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

nfhpcurve: nfhpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o nfhpcurve nfhpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hpcurve: hpcurve$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o hpcurve hpcurve$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hpcurve1$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h \
homspace.h symb.h hpcurve.cc
	$(CXX) $(CFLAGS) -DSINGLE hpcurve.cc -o hpcurve1$(OBJ_SUF).o

hpcurve1: hpcurve1$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o hpcurve1 hpcurve1$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hppcurve: hppcurve$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o hppcurve hppcurve$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hpper: hpper$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o hpper hpper$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

pcurve: pcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o pcurve pcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

pcurve1$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h \
homspace.h symb.h pcurve.cc pcprocs.h
	$(CXX) $(CFLAGS) -DSINGLE pcurve.cc -o pcurve1$(OBJ_SUF).o

pcurve1: pcurve1$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o pcurve1 pcurve1$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

allpcurve: allpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o allpcurve allpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1bsdc: h1bsdc$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1bsdc h1bsdc$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

bsdtex2: bsdtex2$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o bsdtex2 bsdtex2$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1clist: h1clist$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1clist h1clist$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1curve: h1curve$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1curve h1curve$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1periods: h1periods$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1periods h1periods$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1curve1$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h \
	    periods.h h1curve.cc curvesort.cc
	$(CXX) $(CFLAGS) -DSINGLE h1curve.cc -o h1curve1$(OBJ_SUF).o

h1curve1: h1curve1$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1curve1 h1curve1$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1curisog: h1curisog$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1curisog h1curisog$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1isogbsd: h1isogbsd$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1isogbsd h1isogbsd$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

nfconvert: nfconvert$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o nfconvert nfconvert$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)


h1bsdcurisog: h1bsdcurisog$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1bsdcurisog h1bsdcurisog$(OBJ_SUF).o $(G0NLIBFLAG) $(RANKLIBFLAG) $(LFLAGS)

th1nf: th1nf$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o th1nf th1nf$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

mazurcheck$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h \
	    periods.h mazurcheck.cc curvesort.cc
	$(CXX) $(CFLAGS) mazurcheck.cc -o mazurcheck$(OBJ_SUF).o
mazurcheck: mazurcheck$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o mazurcheck mazurcheck$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)
mazurcheck2$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h \
	    periods.h mazurcheck2.cc curvesort.cc
	$(CXX) $(CFLAGS) mazurcheck2.cc -o mazurcheck2$(OBJ_SUF).o
mazurcheck2: mazurcheck2$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o mazurcheck2 mazurcheck2$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)


tnf2: tnf2$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o tnf2 tnf2$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)


h1degphi: h1degphi$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1degphi h1degphi$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1degphi_mw: h1degphi_mw$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o h1degphi_mw h1degphi_mw$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

mazurdeg: mazurdeg$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o mazurdeg mazurdeg$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

degphitab: degphitab$(OBJ_SUF).o $(G0NLIBNAME)
	$(CXX) -o degphitab degphitab$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

degphitex: degphitex$(OBJ_SUF).o
	$(CXX) -o degphitex degphitex$(OBJ_SUF).o $(LFLAGS)

test: test$(OBJ_SUF).o
	$(CXX) -o test test$(OBJ_SUF).o $(LFLAGS)

checkgens$(OBJ_SUF).o: checkgens.cc curvesort.cc
checkgens: checkgens$(OBJ_SUF).o
	   $(CXX) -o checkgens checkgens$(OBJ_SUF).o  $(LFLAGS)

# DEPENDENCIES

cusp$(OBJ_SUF).o: moddata.h symb.h cusp.h cusp.cc
homspace$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h homspace.cc
moddata$(OBJ_SUF).o: moddata.h moddata.cc
oldforms$(OBJ_SUF).o: moddata.h oldforms.h oldforms.cc
newforms$(OBJ_SUF).o: newforms.cc newforms.h moddata.h homspace.h oldforms.h 

modtest$(OBJ_SUF).o: modtest.cc moddata.h symb.h
homtest$(OBJ_SUF).o: homtest.cc moddata.h symb.h cusp.h homspace.h
hecketest$(OBJ_SUF).o: hecketest.cc moddata.h symb.h cusp.h homspace.h
wsplit$(OBJ_SUF).o: wsplit.cc moddata.h symb.h cusp.h homspace.h
gd2$(OBJ_SUF).o: gd2.cc moddata.h symb.h cusp.h homspace.h
oftest$(OBJ_SUF).o: oftest.cc moddata.h symb.h cusp.h oldforms.h
nftest$(OBJ_SUF).o: nftest.cc moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h
ecnf$(OBJ_SUF).o: ecnf.cc moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h

tmanin$(OBJ_SUF).o: tmanin.cc moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h 

th1nf$(OBJ_SUF).o: th1nf.cc moddata.h symb.h cusp.h homspace.h  oldforms.h 
tnf2$(OBJ_SUF).o: tnf2.cc moddata.h symb.h cusp.h homspace.h  oldforms.h \
         newforms.h
h1curve$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h periods.h h1curve.cc curvesort.cc newforms.h
h1curisog$(OBJ_SUF).o: h1curisog.cc moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1curve.cc curvesort.cc newforms.h
h1isogbsd$(OBJ_SUF).o: h1isogbsd.cc moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1curve.cc curvesort.cc newforms.h
h1clist$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1clist.cc curvesort.cc
h1degphi$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1degphi.cc curvesort.cc
h1degphi_mw$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1degphi.cc curvesort.cc
mazurdeg$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h mazurdeg.cc curvesort.cc
degphitab$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h degphitab.cc curvesort.cc
degphitex$(OBJ_SUF).o: degphitab.cc curvesort.cc
h1bsd$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h h1bsd.cc
h1bsdc$(OBJ_SUF).o: h1bsdc.cc moddata.h oldforms.h newforms.h symb.h cusp.h  periods.h
h1bsdc$(OBJ_SUF).o: h1bsdc.cc moddata.h oldforms.h newforms.h symb.h cusp.h periods.h curvesort.cc
lf1per$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h lf1per.cc
nfhpcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h cusp.h pcprocs.h nfhpcurve.cc
hpcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h hpcurve.cc
pcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h pcurve.cc
allpcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h allpcurve.cc

periods$(OBJ_SUF).o: periods.cc moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h periods.h fixc6.h
fixc6$(OBJ_SUF).o: fixc6.h fixc6.cc
symb$(OBJ_SUF).o: symb.h moddata.h symb.cc
getoldeigs$(OBJ_SUF).o: getoldeigs.cc
putoldeigs$(OBJ_SUF).o: putoldeigs.cc

nfd$(OBJ_SUF).o: nfd.cc nfd.h moddata.h symb.h homspace.h
tnfd$(OBJ_SUF).o: tnfd.cc  nfd.h moddata.h symb.h homspace.h

nflist$(OBJ_SUF).o: nflist.h nflist.cc
findcong$(OBJ_SUF).o: findcong.cc nflist.h
find1cong$(OBJ_SUF).o: find1cong.cc nflist.h
cong_vec$(OBJ_SUF).o: cong_vec.cc moddata.h symb.h oldforms.h homspace.h newforms.h  periods.h
cong_vec_h1$(OBJ_SUF).o: cong_vec_h1.cc cong_procs.h moddata.h symb.h oldforms.h homspace.h newforms.h  periods.h
cong_vec_x$(OBJ_SUF).o: cong_vec_x.cc cong_procs.h
cong_procs$(OBJ_SUF).o: cong_procs.h cong_procs.cc
pcurve$(OBJ_SUF).o: pcurve.cc pcprocs.h moddata.h symb.h oldforms.h homspace.h newforms.h periods.h pcprocs.h
pcprocs$(OBJ_SUF).o: pcprocs.cc moddata.h symb.h oldforms.h homspace.h newforms.h periods.h pcprocs.h

vlift$(OBJ_SUF).o: 

vlift: vlift$(OBJ_SUF).o 
	$(CXX) -o vlift vlift$(OBJ_SUF).o $(LFLAGS)

e2x: e2x.cc
	$(CXX) $(OPTFLAG) e2x.cc -o e2x $(IOLIB)

x2e: x2e.cc
	$(CXX) $(OPTFLAG) x2e.cc -o x2e $(IOLIB)

pair$(OBJ_SUF).o: pair.cc
pair: pair$(OBJ_SUF).o
	$(CXX) -static -o pair pair$(OBJ_SUF).o $(LFLAGS)
pair2$(OBJ_SUF).o: pair2.cc
pair2: pair2$(OBJ_SUF).o
	$(CXX) -static -o pair2 pair2$(OBJ_SUF).o $(LFLAGS)
pair1$(OBJ_SUF).o: pair1.cc
pair1: pair1$(OBJ_SUF).o
	$(CXX) -static -o pair1 pair1$(OBJ_SUF).o $(LFLAGS)

aphomspace$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h aphomspace.cc
aphomtest$(OBJ_SUF).o: aphomtest.cc moddata.h symb.h cusp.h homspace.h
aphomtest: homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o
	$(CXX) -o aphomtest homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o  $(LFLAGS)
aphecketest: hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CXX) -o aphecketest hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)


LIBBASE=$(G0NLIBBASE)
include ../Makefile.dynamic
