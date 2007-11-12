# File Makefile: architecture-independent elements of g0n Makefile

include ../Makefile

# Useful targets:
#
# tests: builds various test programs
#
#

DOT_H = periods.h cusp.h homspace.h newforms.h nflist.h moddata.h oldforms.h symb.h fixc6.h pcprocs.h

DOT_CC = h1bsd.cc newforms.cc symb.cc checkap.cc h1bsdc.cc homspace.cc nfcount.cc checkdegphi.cc h1clist.cc homtest.cc nflist.cc th1nf.cc cnf.cc h1curisog.cc hpcurve.cc nftest.cc tmanin.cc curvesort.cc h1curve.cc hpmcurve.cc oftest.cc tmaninx.cc cusp.cc h1degphi.cc lf1per.cc oldforms.cc tmaninxauto.cc fixc6.cc h1degphitex.cc pcurve.cc tnf2.cc fixforms.cc h1first.cc periods.cc vlift.cc g2.cc mhcount.cc pertest.cc getai.cc h1test.cc moddata.cc getoldeigs.cc hecketest.cc modtest.cc putoldeigs.cc h1.cc heegner.cc nfhpcurve.cc   

HEADERS = $(DOT_H)
sources: $(DOT_H) $(DOT_CC)

TESTS = modtest homtest hecketest nftest tmanin # oftest h1test

PROGS = tmanin h1bsd h1curve pcurve nfhpcurve h1clist checkap h1bsd h1bsdcurisog # h1degphi  h1bsdc

ifeq ($(ARITH),NTL_INTS)
LIBNAME = libjcntlint.a
LIDIACFLAGS = -DNTL_INTS  -I$(NTLINCDIR)
LIDIALFLAGS = -L$(NTLLIBDIR) -lntl -lgmp
LIBFLAG = -lcurvesntlint -ljcntlint -L$(LIBDIR)
CURVELIBNAME = libcurvesntlint.a
CURVELIBFLAG = -lcurvesntlint -L.
RANKLIBNAME = librankntlint.a
RANKLIBFLAG = -lrankntlint
OBJ_SUF = _ntl
EXE_SUFFIX = _ntl
G0NLIBNAME = libg0nntlint.a
G0NLIBFLAG = -lg0nntlint -L.
else
ifeq ($(ARITH),NTL)
LIBNAME = libjcntl.a
LIDIACFLAGS = -DNTL_ALL  -I$(NTLINCDIR)
LIDIALFLAGS = -L$(NTLLIBDIR) -lntl -lgmp
LIBFLAG = -lcurvesntl -ljcntl -L$(LIBDIR)
CURVELIBNAME = libcurvesntl.a
CURVELIBFLAG = -lcurvesntl -L.
RANKLIBNAME = librankntl.a
RANKLIBFLAG = -lrankntl
OBJ_SUF = _n
EXE_SUFFIX = _n
G0NLIBNAME = libg0nntl.a
G0NLIBFLAG = -lg0nntl -L.
else
LIBNAME = libjc.a
LIDIACFLAGS = 
LIDIALFLAGS = 
LIBFLAG = -lcurves -ljc -L$(LIBDIR)
CURVELIBNAME = libcurves.a
OBJ_SUF = _
EXE_SUFFIX = _
G0NLIBNAME = libg0n.a
G0NLIBFLAG = -lg0n -L.
RANKLIBNAME = librank.a
RANKLIBFLAG = -lrank
endif
endif

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
	rm -f $(LIBDIR)/$(G0NLIBNAME)
	cp $(G0NLIBNAME) $(LIBDIR)
	chmod a+rx $(LIBDIR)/$(G0NLIBNAME)

install: lib install_lib

CFLAGS = -c $(OPTFLAG)  -I$(INCDIR) $(LIDIACFLAGS) -I$(PARIINCDIR) -DMETHOD=2 -DUSE_XSPLIT
LFLAGS = -lpari -L$(PARILIBDIR) $(LIBFLAG) $(LIDIALFLAGS) $(IOLIB) -lm

LIBFLAGS = $(LFLAGS)

%$(OBJ_SUF).o:   %.cc
	$(CC) $(CFLAGS) $< -o $@


# TARGETS


progs: $(PROGS)

tests: $(TESTS)

clean:
	/bin/rm -f *~ *$(OBJ_SUF).o

veryclean:
	$(CLEAN)
	/bin/rm -f *~ *$(OBJ_SUF).o $(PROGS) $(TESTS)

input: input.cc
	$(CC) input.cc -o inp

mhcount: mhcount.cc
	$(CC) $(CFLAGS) mhcount.cc
	$(CC) -o mhcount mhcount.o   $(LFLAGS)

g0ngens: g0ngens$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o g0ngens g0ngens$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o   $(LFLAGS)

modtest: modtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o modtest modtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o   $(LFLAGS)

homtest: homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o
	$(CC) -o homtest homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o  $(LFLAGS)

oftest: oftest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o homspace$(OBJ_SUF).o newforms$(OBJ_SUF).o  cusp$(OBJ_SUF).o
	$(CC) -o oftest oftest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o homspace$(OBJ_SUF).o newforms$(OBJ_SUF).o  cusp$(OBJ_SUF).o $(LFLAGS)

fixforms: fixforms$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o fixforms fixforms$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

hecketest: hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o hecketest hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

wsplit: wsplit$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o wsplit wsplit$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

g2: g2$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o g2 g2$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

gd: gd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o gd gd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

gd2: gd2$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o gd2 gd2$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

tnfd: tnfd$(OBJ_SUF).o nfd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o tnfd tnfd$(OBJ_SUF).o nfd$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)

nftest: nftest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o oldforms$(OBJ_SUF).o \
		newforms$(OBJ_SUF).o
	$(CC) -o nftest nftest$(OBJ_SUF).o cusp$(OBJ_SUF).o homspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		oldforms$(OBJ_SUF).o newforms$(OBJ_SUF).o  $(LFLAGS)

nfcount: nfcount$(OBJ_SUF).o
	$(CC) -o nfcount nfcount$(OBJ_SUF).o $(LFLAGS)

aplist: aplist$(OBJ_SUF).o
	$(CC) -o aplist aplist$(OBJ_SUF).o $(LFLAGS)

aplist2: aplist2$(OBJ_SUF).o
	$(CC) -o aplist2 aplist2$(OBJ_SUF).o $(LFLAGS)

findcong: findcong$(OBJ_SUF).o nflist$(OBJ_SUF).o
	$(CC) -o findcong findcong$(OBJ_SUF).o nflist$(OBJ_SUF).o $(LFLAGS)

find1cong: find1cong$(OBJ_SUF).o nflist$(OBJ_SUF).o
	$(CC) -o find1cong find1cong$(OBJ_SUF).o nflist$(OBJ_SUF).o $(LFLAGS)

checkap: checkap$(OBJ_SUF).o
	$(CC) -o checkap checkap$(OBJ_SUF).o -I$(INCDIR) -L$(LIBDIR) $(LFLAGS)

tmanin: tmanin$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o tmanin tmanin$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

checkmanin: checkmanin$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o checkmanin checkmanin$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

mtmanin: tmanin$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o mtmanin $(G0NLIBFLAG) $(LFLAGS)

correct: correct$(OBJ_SUF).o
	 $(CC) -o correct correct$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

getoldeigs: getoldeigs$(OBJ_SUF).o
	$(CC) -o getoldeigs getoldeigs$(OBJ_SUF).o $(LFLAGS)

putoldeigs: putoldeigs$(OBJ_SUF).o
	$(CC) -o putoldeigs putoldeigs$(OBJ_SUF).o $(LFLAGS)


fix1290: fix1290$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o fix1290 fix1290$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

pmanin: pmanin$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o pmanin pmanin$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

xmanin: xmanin$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o xmanin xmanin$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

tmaninx: tmaninx$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o tmaninx tmaninx$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

cnf: cnf$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o cnf cnf$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

tmaninxauto: tmaninxauto$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o tmaninxauto tmaninxauto$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1test: h1test$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1test h1test$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1first: h1first$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1first h1first$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

mh1first: h1first$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o mh1first h1first$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1: h1$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1 h1$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1bsd: h1bsd$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1bsd h1bsd$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

cong_vec: cong_vec$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o cong_vec cong_vec$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

cong_vec_h1: cong_vec_h1$(OBJ_SUF).o cong_procs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o cong_vec_h1 cong_vec_h1$(OBJ_SUF).o cong_procs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

cong_vec_x: cong_vec_x$(OBJ_SUF).o cong_procs$(OBJ_SUF).o  $(G0NLIBNAME)
	$(CC) -o cong_vec_x cong_vec_x$(OBJ_SUF).o cong_procs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

lf1per: lf1per$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o lf1per lf1per$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

nfhpcurve: nfhpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o nfhpcurve nfhpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hpcurve: hpcurve$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o hpcurve hpcurve$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hpcurve1$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h \
homspace.h symb.h hpcurve.cc
	$(CC) $(CFLAGS) -DSINGLE hpcurve.cc -o hpcurve1$(OBJ_SUF).o

hpcurve1: hpcurve1$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o hpcurve1 hpcurve1$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hpmcurve: hpmcurve$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o hpmcurve hpmcurve$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hppcurve: hppcurve$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o hppcurve hppcurve$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

hpper: hpper$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o hpper hpper$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

pcurve: pcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o pcurve pcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

pcurve1$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h \
homspace.h symb.h pcurve.cc pcprocs.h
	$(CC) $(CFLAGS) -DSINGLE pcurve.cc -o pcurve1$(OBJ_SUF).o

pcurve1: pcurve1$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o pcurve1 pcurve1$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

allpcurve: allpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o allpcurve allpcurve$(OBJ_SUF).o pcprocs$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1bsdc: h1bsdc$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1bsdc h1bsdc$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

bsdtex2: bsdtex2$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o bsdtex2 bsdtex2$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1clist: h1clist$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1clist h1clist$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1curve: h1curve$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1curve h1curve$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1periods: h1periods$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1periods h1periods$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1curve1$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h \
	    periods.h h1curve.cc curvesort.cc
	$(CC) $(CFLAGS) -DSINGLE h1curve.cc -o h1curve1$(OBJ_SUF).o

h1curve1: h1curve1$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1curve1 h1curve1$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1curisog: h1curisog$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1curisog h1curisog$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1isogbsd: h1isogbsd$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1isogbsd h1isogbsd$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

nfconvert: nfconvert$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o nfconvert nfconvert$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)


h1bsdcurisog: h1bsdcurisog$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1bsdcurisog h1bsdcurisog$(OBJ_SUF).o $(G0NLIBFLAG) $(RANKLIBFLAG) $(LFLAGS)

heegner: heegner$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o heegner heegner$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

pertest: pertest$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o pertest pertest$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

th1nf: th1nf$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o th1nf th1nf$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

checkdegphi: checkdegphi$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o checkdegphi checkdegphi$(OBJ_SUF).o  $(G0NLIBFLAG)  $(LFLAGS)

mazurcheck$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h \
	    periods.h mazurcheck.cc curvesort.cc
	$(CC) $(CFLAGS) mazurcheck.cc -o mazurcheck$(OBJ_SUF).o
mazurcheck: mazurcheck$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o mazurcheck mazurcheck$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)
mazurcheck2$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h \
	    periods.h mazurcheck2.cc curvesort.cc
	$(CC) $(CFLAGS) mazurcheck2.cc -o mazurcheck2$(OBJ_SUF).o
mazurcheck2: mazurcheck2$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o mazurcheck2 mazurcheck2$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)


tnf2: tnf2$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o tnf2 tnf2$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)


h1degphi: h1degphi$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1degphi h1degphi$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

h1degphi_mw: h1degphi_mw$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o h1degphi_mw h1degphi_mw$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

mazurdeg: mazurdeg$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o mazurdeg mazurdeg$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

degphitab: degphitab$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o degphitab degphitab$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)

degphitex: degphitex$(OBJ_SUF).o
	$(CC) -o degphitex degphitex$(OBJ_SUF).o $(LFLAGS)

test: test$(OBJ_SUF).o
	$(CC) -o test test$(OBJ_SUF).o $(LFLAGS)

testg1: testg1$(OBJ_SUF).o $(G0NLIBNAME)
	$(CC) -o testg1 testg1$(OBJ_SUF).o $(G0NLIBFLAG)  $(LFLAGS)

checkgens$(OBJ_SUF).o: checkgens.cc curvesort.cc
checkgens: checkgens$(OBJ_SUF).o
	   $(CC) -o checkgens checkgens$(OBJ_SUF).o  $(LFLAGS)

# DEPENDENCIES

cusp$(OBJ_SUF).o: moddata.h symb.h cusp.h cusp.cc
homspace$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h homspace.cc
moddata$(OBJ_SUF).o: moddata.h moddata.cc
oldforms$(OBJ_SUF).o: moddata.h oldforms.h oldforms.cc
newforms$(OBJ_SUF).o: newforms.cc newforms.h moddata.h homspace.h oldforms.h 
manin$(OBJ_SUF).o: manin.cc manin.h moddata.h symb.h homspace.h oldforms.h newforms.h
splitter$(OBJ_SUF).o: splitter.h splitter.cc moddata.h symb.h cusp.h homspace.h
msplitter$(OBJ_SUF).o: msplitter.h splitter.h msplitter.cc moddata.h symb.h cusp.h homspace.h
getai$(OBJ_SUF).o: getai.cc getai.h

modtest$(OBJ_SUF).o: modtest.cc moddata.h symb.h
homtest$(OBJ_SUF).o: homtest.cc moddata.h symb.h cusp.h homspace.h
hecketest$(OBJ_SUF).o: hecketest.cc moddata.h symb.h cusp.h homspace.h
wsplit$(OBJ_SUF).o: wsplit.cc moddata.h symb.h cusp.h homspace.h
g2$(OBJ_SUF).o: g2.cc moddata.h symb.h cusp.h homspace.h
gd2$(OBJ_SUF).o: gd2.cc moddata.h symb.h cusp.h homspace.h
oftest$(OBJ_SUF).o: oftest.cc moddata.h symb.h cusp.h oldforms.h
nftest$(OBJ_SUF).o: nftest.cc moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h \
	
tmanin$(OBJ_SUF).o: tmanin.cc manin.h moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h 
fix1290$(OBJ_SUF).o: fix1290.cc manin.h moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h 
pmanin$(OBJ_SUF).o: pmanin.cc manin.h moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h 
tmaninx$(OBJ_SUF).o: tmaninx.cc manin.h moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h 
xmanin$(OBJ_SUF).o: xmanin.cc manin.h moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h 
tmaninxauto$(OBJ_SUF).o: tmaninxauto.cc manin.h moddata.h symb.h cusp.h homspace.h \
	oldforms.h newforms.h 
h1test$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	 periods.h h1test.cc
h1first$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	 periods.h h1first.cc newforms.h
h1$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1.cc getai.h newforms.h fixc6.h

th1nf$(OBJ_SUF).o: th1nf.cc moddata.h symb.h cusp.h homspace.h  oldforms.h 
checkdegphi$(OBJ_SUF).o: checkdegphi.cc moddata.h symb.h cusp.h homspace.h  oldforms.h 
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
pertest$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h pertest.cc
h1degphi$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1degphi.cc getai.h curvesort.cc
h1degphi_mw$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h h1degphi.cc getai.h curvesort.cc
mazurdeg$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h mazurdeg.cc getai.h curvesort.cc
degphitab$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h oldforms.h \
	periods.h degphitab.cc getai.h curvesort.cc
degphitex$(OBJ_SUF).o: degphitab.cc curvesort.cc
h1bsd$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h h1bsd.cc
h1bsdc$(OBJ_SUF).o: h1bsdc.cc moddata.h oldforms.h newforms.h symb.h cusp.h  periods.h
h1bsdc$(OBJ_SUF).o: h1bsdc.cc moddata.h oldforms.h newforms.h symb.h cusp.h periods.h curvesort.cc
lf1per$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h lf1per.cc
nfhpcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h manin.h cusp.h pcprocs.h nfhpcurve.cc
hpcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h hpcurve.cc
pcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h pcurve.cc
allpcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h allpcurve.cc
hpmcurve$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h hpmcurve.cc
heegner$(OBJ_SUF).o:  moddata.h oldforms.h newforms.h periods.h homspace.h symb.h heegner.cc

periods$(OBJ_SUF).o: periods.cc moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h periods.h fixc6.h
heegnerp$(OBJ_SUF).o: heegnerp.cc moddata.h symb.h cusp.h homspace.h oldforms.h \
	newforms.h  periods.h heegnerp.h
bsdprocs$(OBJ_SUF).o: bsdprocs.cc moddata.h symb.h cusp.h homspace.h oldforms.h \
	bsdprocs.h
fixc6$(OBJ_SUF).o: fixc6.h fixc6.cc
getai$(OBJ_SUF).o: getai.h getai.cc
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
	$(CC) -o vlift vlift$(OBJ_SUF).o $(LFLAGS)

e2x: e2x.cc
	$(CC) $(OPTFLAG) e2x.cc -o e2x $(IOLIB)

x2e: x2e.cc
	$(CC) $(OPTFLAG) x2e.cc -o x2e $(IOLIB)

pair$(OBJ_SUF).o: pair.cc
pair: pair$(OBJ_SUF).o
	$(CC) -static -o pair pair$(OBJ_SUF).o $(LFLAGS)
pair2$(OBJ_SUF).o: pair2.cc
pair2: pair2$(OBJ_SUF).o
	$(CC) -static -o pair2 pair2$(OBJ_SUF).o $(LFLAGS)
pair1$(OBJ_SUF).o: pair1.cc
pair1: pair1$(OBJ_SUF).o
	$(CC) -static -o pair1 pair1$(OBJ_SUF).o $(LFLAGS)
checknf: checknf$(OBJ_SUF).o $(G0NLIBNAME)
	 $(CC) -o checknf checknf$(OBJ_SUF).o $(G0NLIBFLAG) $(LFLAGS)
checknf$(OBJ_SUF).o: checknf.cc moddata.h symb.h cusp.h homspace.h oldforms.h newforms.h 

aphomspace$(OBJ_SUF).o: moddata.h symb.h cusp.h homspace.h aphomspace.cc
aphomtest$(OBJ_SUF).o: aphomtest.cc moddata.h symb.h cusp.h homspace.h
aphomtest: homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o
	$(CC) -o aphomtest homtest$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o  $(LFLAGS)
aphecketest: hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o
	$(CC) -o aphecketest hecketest$(OBJ_SUF).o cusp$(OBJ_SUF).o aphomspace$(OBJ_SUF).o symb$(OBJ_SUF).o moddata$(OBJ_SUF).o \
		 $(LFLAGS)
