# File Makefile: architecture-independent elements of procs Makefile
include ../Makefile

# Useful targets:
#
# install:  builds library, puts it into ../lib and includes in ../include
# tests:    builds various test programs
#

.PHONY = clean

ifeq ($(ARITH),NTL_INTS)
LIBBASE = libjcntlint
NTLCFLAGS = -DNTL_INTS -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -lgmp -lpari
OBJ_SUF = _ntl
else
ifeq ($(ARITH),NTL)
LIBBASE = libjcntl
NTLCFLAGS = -DNTL_ALL -I$(NTLINCDIR)
NTLLFLAGS = -L$(NTLLIBDIR) -lntl -lgmp -lpari
OBJ_SUF = _n
else
LIBBASE = libjc.a
NTLCFLAGS =
NTLLFLAGS =
OBJ_SUF = ._
endif
endif

LIBNAME=$(LIBBASE).a

CFLAGS = -c $(OPTFLAG) $(NTLCFLAGS) -I$(PARIINCDIR)
LFLAGS = -lpari -L$(PARILIBDIR) $(NTLLFLAGS) -lm $(IOLIB) 


LIBFLAGS = $(LFLAGS)

HEADERS = interface.h templates.h arith.h xmod.h marith.h gpslave.h compproc.h vec.h vector.h mat.h matrix.h sub.h subspace.h rat.h kbessel.h svec.h svector.h smat.h smatrix.h smat_elim.h smatrix_elim.h mvector.h mmatrix.h msubspace.h method.h splitbase.h xsplit.h conic.h legendre.h quadratic.h unimod.h illl.h hilbert.h timer.h cubic.h gf.h polys.h realroots.h parifact.h

DOTCCS = interface.cc unimod.cc arith.cc marith.cc gpslave.cc compproc.cc vec.cc vector.cc mat.cc matrix.cc sub.cc subspace.cc kbessel.cc mvector.cc mmatrix.cc msubspace.cc vectest.cc mvectest.cc mattest.cc mmattest.cc space.cc mspace.cc comptest.cc rattest.cc ptest.cc mptest.cc tbessel.cc svec.cc svector.cc svectest.cc smat.cc smatrix.cc smattest.cc random.cc smat_elim.cc smatrix_elim.cc xsplit.cc conic.cc tconic.cc legendre.cc quadratic.cc illl.cc hilbert.cc thilbert.cc tleg.cc timer.cc cubic.cc polys.cc realroots.cc parifact.cc rcubic.cc lcubic.cc

ARITH_OBJS  = arith$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o conic$(OBJ_SUF).o legendre$(OBJ_SUF).o unimod$(OBJ_SUF).o quadratic$(OBJ_SUF).o cubic$(OBJ_SUF).o illl$(OBJ_SUF).o hilbert$(OBJ_SUF).o parifact$(OBJ_SUF).o
TRANS_OBJS =  compproc$(OBJ_SUF).o kbessel$(OBJ_SUF).o
#XSPLIT_OBJS = xsplit0$(OBJ_SUF).o xsplit1$(OBJ_SUF).o xsplit2$(OBJ_SUF).o xsplit5$(OBJ_SUF).o xsplit6$(OBJ_SUF).o
XSPLIT_OBJS = xsplit2$(OBJ_SUF).o # temporary
LINALG_OBJS = vector$(OBJ_SUF).o matrix$(OBJ_SUF).o subspace$(OBJ_SUF).o
MLINALG_OBJS = mvector$(OBJ_SUF).o mmatrix$(OBJ_SUF).o msubspace$(OBJ_SUF).o
SMAT_OBJS = smatrix$(OBJ_SUF).o smatrix_elim$(OBJ_SUF).o
SVEC_OBJS = svector$(OBJ_SUF).o

OBJS = interface$(OBJ_SUF).o $(ARITH_OBJS) $(TRANS_OBJS) $(LINALG_OBJS) $(MLINALG_OBJS) $(XSPLIT_OBJS) timer$(OBJ_SUF).o $(SVEC_OBJS) $(SMAT_OBJS) polys$(OBJ_SUF).o realroots$(OBJ_SUF).o 

TESTS = vectest1 vectest2 mattest1 mattest2 space1 space2 smattest comptest rattest ptest mptest tbessel mvectest mmattest mspace tconic tleg thilbert rcubic lcubic svectest1 svectest2

all: $(OBJS)

tests: $(TESTS)
progs: $(TESTS)

headers: $(HEADERS)
objs: $(OBJS)
sources: $(HEADERS) $(DOTCCS)

lib: headers objs
	ar r $(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBNAME)

install: sources
	make lib 
	chmod a+r *
	(cd $(INCDIR); rm -f $(HEADERS))
	$(CP) $(HEADERS) $(INCDIR)
	chmod a+r $(INCDIR)/*.h
	rm -f $(LIBDIR)/$(LIBNAME)
	ar r $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)
	chmod a+rx $(LIBDIR)/$(LIBNAME)

DISTFILES = arith.h xmod.h arith.cc marith.h marith.cc compproc.h compproc.cc conic.cc conic.h tconic.cc legendre.h legendre.cc tleg.cc interface.h templates.h timer.h timer.cc matrix.h mvector.h mat.h vec.h vector.h vector.cc vec.cc matrix.cc mat.cc mvector.cc quadratic.h quadratic.cc illl.h illl.cc mmatrix.h mmatrix.cc unimod.h hilbert.h hilbert.cc unimod.cc interface.cc gpslave.h gpslave.cc cubic.h cubic.cc sub.h sub.cc subspace.h subspace.cc parifact.h parifact.cc gf.h polys.h polys.cc realroots.h realroots.cc
distfiles:  $(DISTFILES)

dist: distfiles
	$(CP) $(DISTFILES) $(HOME)/curves/rankdist
	chmod 644 $(HOME)/curves/rankdist/*

g0ndist: sources
	 $(CP) $(HEADERS) $(HOME)/curves/g0ndist
	 $(CP) $(DOTCCS) $(HOME)/curves/g0ndist
	 chmod 644 $(HOME)/curves/g0ndist/*

clean:
	/bin/rm -f *$(OBJ_SUF).o *~ *.a *.so *.dylib

veryclean: clean
	(cd $(BINDIR); /bin/rm -f $(TESTS) $(PROGS))

shar: sources
	rm -f procs.shar procs.shar.Z procs.shar.gz
	shar -i files -o procs.shar
	gzip procs.shar

tar: sources
	rm -f procs.tar
	tar -cvf procs.tar -I files
	gzip procs.tar

check: $(TESTS)
	rm -f PRIMES
	./vectest1 < vectest.in >  t && diff t vectest.out
	./vectest2 < vectest.in >  t && diff t vectest.out
	./mattest1 < mattest.in >  t && diff t mattest.out
	./mattest2 < mattest.in >  t && diff t mattest.out
	./space1 < space.in > t && diff t space.out
	./space2 < space.in > t && diff t space.out
	./svectest1 < svectest.in > t && diff t svectest.out
	./svectest2 < svectest.in > t && diff t svectest.out
	./smattest < smattest.in > t && diff t smattest.out
	./mvectest < mvectest.in > t && diff t mvectest.out
	./mmattest < mmattest.in > t && diff t mmattest.out
	./mspace < mspace.in > t && diff t mspace.out
	./rattest < rattest.in > t && diff t rattest.out
	./ptest < ptest.in > t && diff t ptest.out
	./mptest < mptest.in > t && diff t mptest.out
	./tconic < tconic.in > t && diff t tconic.out
	./tleg < tleg.in > t && diff t tleg.out
	./thilbert < thilbert.in > t && diff t thilbert.out
	./rcubic < rcubic.in > t && diff t rcubic.out
	./lcubic < lcubic.in > t && diff t lcubic.out

vectest1: vectest.cc vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=1 vectest.cc -o vectest1$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/vectest1 vectest1$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

vectest2: vectest.cc vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=2 vectest.cc -o vectest2$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/vectest2 vectest2$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

mattest1: mattest.cc matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o timer$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=1 mattest.cc -o mattest1$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/mattest1 mattest1$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o timer$(OBJ_SUF).o $(LFLAGS)

mattest2: mattest.cc matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o timer$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=2 mattest.cc -o mattest2$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/mattest2 mattest2$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o timer$(OBJ_SUF).o $(LFLAGS)

space1: space.cc subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=1 space.cc -o space1$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/space1 space1$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o \
	 interface$(OBJ_SUF).o $(LFLAGS)

space2: space.cc subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=2 space.cc -o space2$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/space2 space2$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o \
	interface$(OBJ_SUF).o $(LFLAGS)

mvectest: mvectest$(OBJ_SUF).o mvector$(OBJ_SUF).o vector$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/mvectest mvectest$(OBJ_SUF).o mvector$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o  \
	compproc$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

svectest1: svectest.cc svector$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=1 svectest.cc -o svectest1$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/svectest1 svectest1$(OBJ_SUF).o svector$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

svectest2: svectest.cc svector$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) $(CFLAGS) -DSCALAR_OPTION=2 svectest.cc -o svectest2$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/svectest2 svectest2$(OBJ_SUF).o svector$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

#svectest: svectest$(OBJ_SUF).o svec$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o
#	$(CXX) -o $(BINDIR)/svectest svectest$(OBJ_SUF).o svec$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o $(LFLAGS)

lintest: lintest$(OBJ_SUF).o linsyst$(OBJ_SUF).o linear$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/lintest lintest$(OBJ_SUF).o linsyst$(OBJ_SUF).o linear$(OBJ_SUF).o $(LFLAGS)

mmattest: mmattest$(OBJ_SUF).o mmatrix$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o mvector$(OBJ_SUF).o arith$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/mmattest mmattest$(OBJ_SUF).o  mmatrix$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o \
	mvector$(OBJ_SUF).o arith$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

smattest: smattest$(OBJ_SUF).o smatrix_elim$(OBJ_SUF).o svector$(OBJ_SUF).o smatrix$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o 
	  $(CXX) -Wall -g -o $(BINDIR)/smattest smattest$(OBJ_SUF).o smatrix_elim$(OBJ_SUF).o svector$(OBJ_SUF).o smatrix$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o $(LFLAGS)

smat2test: smat2test$(OBJ_SUF).o smatrix2_elim$(OBJ_SUF).o svector$(OBJ_SUF).o smatrix2$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o 
	  $(CXX) -Wall -g -o $(BINDIR)/smat2test smat2test$(OBJ_SUF).o smatrix2_elim$(OBJ_SUF).o svector$(OBJ_SUF).o smatrix2$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o $(LFLAGS)

matxtest: matxtest$(OBJ_SUF).o  matx$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/matxtest matxtest$(OBJ_SUF).o  matx$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o $(LFLAGS)

mspace: mspace$(OBJ_SUF).o msubspace$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o compproc$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o 
	$(CXX) -o $(BINDIR)/mspace mspace$(OBJ_SUF).o msubspace$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o \
	 parifact$(OBJ_SUF).o compproc$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

periods: periods$(OBJ_SUF).o parifact$(OBJ_SUF).o compproc$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/periods periods$(OBJ_SUF).o compproc$(OBJ_SUF).o $(LFLAGS)

comptest: comptest$(OBJ_SUF).o compproc$(OBJ_SUF).o arith$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/comptest comptest$(OBJ_SUF).o compproc$(OBJ_SUF).o arith$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

rattest: rattest$(OBJ_SUF).o arith$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/rattest rattest$(OBJ_SUF).o  arith$(OBJ_SUF).o $(LFLAGS)

ptest: ptest$(OBJ_SUF).o arith$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/ptest ptest$(OBJ_SUF).o arith$(OBJ_SUF).o $(LFLAGS) 

mptest: mptest$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o parifact$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/mptest mptest$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)  parifact$(OBJ_SUF).o

tbessel: tbessel$(OBJ_SUF).o kbessel$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/tbessel tbessel$(OBJ_SUF).o kbessel$(OBJ_SUF).o $(LFLAGS)


read_rat: read_rat$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/read_rat read_rat$(OBJ_SUF).o $(LFLAGS)

divisors: divisors$(OBJ_SUF).o  marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o parifact$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/divisors divisors$(OBJ_SUF).o  marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o parifact$(OBJ_SUF).o $(LFLAGS)
divisors$(OBJ_SUF).o: divisors.cc

fact1: fact1.cc
	$(CXX) fact1.cc -o fact1 -O3 -DGCC_INLINE -lpari-2.2
fact2: fact2$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o parifact$(OBJ_SUF).o
	$(CXX) fact2$(OBJ_SUF).o  parifact$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o interface$(OBJ_SUF).o -o fact2 $(LFLAGS) -lpari-2.2
fact2$(OBJ_SUF).o: fact2.cc marith.h interface.h
	 $(CXX) $(CFLAGS) fact2.cc -o fact2$(OBJ_SUF).o -DGCC_INLINE 

parifact$(OBJ_SUF).o: parifact.cc parifact.h

fact3: fact3$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/fact3 fact3$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS) 


CC.o: CC.cc RRextra.h

CCt: CCt.o CC.o RRextra.o
	$(CC) -o CCt CCt.o CC.o RRextra.o $(LFLAGS)

RRextra.o: RRextra.h RRextra.cc
RRt.o: RRt.cc RRextra.h
	$(CC) -c -g $(OPTFLAG) -I$(NTLINCDIR) -I$(NTLINCDIR) RRt.cc 
RRt: RRt.o RRextra.o
	$(CC) -o RRt RRt.o RRextra.o -L$(NTLLIBDIR) -lLiDIA $(LFLAGS)


gpslave$(OBJ_SUF).o: gpslave.cc interface.h templates.h marith.h arith.h xmod.h \
  gpslave.h

test: test$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o test test$(OBJ_SUF).o subspace$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

fifo: fifo$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/fifo fifo$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS) 


tconic: tconic$(OBJ_SUF).o conic$(OBJ_SUF).o legendre$(OBJ_SUF).o quadratic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o illl$(OBJ_SUF).o unimod$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/tconic tconic$(OBJ_SUF).o conic$(OBJ_SUF).o legendre$(OBJ_SUF).o quadratic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o illl$(OBJ_SUF).o unimod$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

tcon2: tcon2$(OBJ_SUF).o conic$(OBJ_SUF).o legendre$(OBJ_SUF).o quadratic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o illl$(OBJ_SUF).o unimod$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/tcon2 tcon2$(OBJ_SUF).o conic$(OBJ_SUF).o legendre$(OBJ_SUF).o quadratic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o mmatrix$(OBJ_SUF).o mvector$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o illl$(OBJ_SUF).o unimod$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

tleg: tleg$(OBJ_SUF).o legendre$(OBJ_SUF).o conic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o mvector$(OBJ_SUF).o mmatrix$(OBJ_SUF).o vector$(OBJ_SUF).o matrix$(OBJ_SUF).o illl$(OBJ_SUF).o quadratic$(OBJ_SUF).o unimod$(OBJ_SUF).o interface$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/tleg tleg$(OBJ_SUF).o legendre$(OBJ_SUF).o conic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o  mvector$(OBJ_SUF).o mmatrix$(OBJ_SUF).o vector$(OBJ_SUF).o matrix$(OBJ_SUF).o illl$(OBJ_SUF).o quadratic$(OBJ_SUF).o unimod$(OBJ_SUF).o interface$(OBJ_SUF).o $(LFLAGS)

thilbert: thilbert$(OBJ_SUF).o hilbert$(OBJ_SUF).o  quadratic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o interface$(OBJ_SUF).o\
	compproc$(OBJ_SUF).o mvector$(OBJ_SUF).o vector$(OBJ_SUF).o conic$(OBJ_SUF).o legendre$(OBJ_SUF).o illl$(OBJ_SUF).o arith$(OBJ_SUF).o unimod$(OBJ_SUF).o
	$(CXX) -o thilbert thilbert$(OBJ_SUF).o hilbert$(OBJ_SUF).o quadratic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o interface$(OBJ_SUF).o \
	compproc$(OBJ_SUF).o mvector$(OBJ_SUF).o vector$(OBJ_SUF).o conic$(OBJ_SUF).o legendre$(OBJ_SUF).o illl$(OBJ_SUF).o arith$(OBJ_SUF).o unimod$(OBJ_SUF).o $(LFLAGS)

auto: auto$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o
	$(CXX) -o $(BINDIR)/auto auto$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o $(LFLAGS)

auto$(OBJ_SUF).o: auto.cc matrix.h vector.h arith.h xmod.h

test$(OBJ_SUF).o: test.cc
tconic$(OBJ_SUF).o: tconic.cc conic.h legendre.h quadratic.h unimod.h marith.h arith.h xmod.h interface.h templates.h
tcon2$(OBJ_SUF).o: tcon2.cc conic.h legendre.h quadratic.h unimod.h marith.h arith.h xmod.h interface.h templates.h
conic$(OBJ_SUF).o: conic.cc marith.h interface.h templates.h arith.h xmod.h \
  quadratic.h conic.h legendre.h
hilbert$(OBJ_SUF).o: hilbert.cc marith.h interface.h templates.h arith.h xmod.h \
  quadratic.h hilbert.h
thilbert$(OBJ_SUF).o: hilbert.h thilbert.cc

tleg$(OBJ_SUF).o: tleg.cc legendre.h conic.h quadratic.h unimod.h mmatrix.h
legendre$(OBJ_SUF).o: legendre.cc mmatrix.h matrix.h mat.h vec.h interface.h \
  templates.h arith.h xmod.h mvector.h marith.h unimod.h quadratic.h \
  conic.h legendre.h illl.h
tbessel$(OBJ_SUF).o: tbessel.cc kbessel.h
kbessel$(OBJ_SUF).o: kbessel.cc kbessel.h
nagtest$(OBJ_SUF).o: nagtest.cc
rattest$(OBJ_SUF).o: rattest.cc rat.h

interface$(OBJ_SUF).o: interface.cc interface.h templates.h
unimod$(OBJ_SUF).o: unimod.cc interface.h templates.h marith.h arith.h xmod.h \
  unimod.h
arith$(OBJ_SUF).o: arith.cc arith.h interface.h templates.h xmod.h
marith$(OBJ_SUF).o: marith.cc interface.h templates.h marith.h arith.h xmod.h \
  gpslave.h
polys$(OBJ_SUF).o: polys.cc polys.h gf.h
realroots$(OBJ_SUF).o: realroots.h realroots.cc

ptest$(OBJ_SUF).o: ptest.cc arith.h xmod.h interface.h templates.h
space$(OBJ_SUF).o: space.cc subspace.h matrix.h vector.h arith.h xmod.h
comptest$(OBJ_SUF).o: interface.h templates.h comptest.cc compproc.h marith.h arith.h xmod.h
compproc$(OBJ_SUF).o: compproc.cc compproc.h interface.h templates.h
linear$(OBJ_SUF).o: linear.cc linear.h
linsyst$(OBJ_SUF).o: linsyst.cc linsyst.h linear.h
lintest$(OBJ_SUF).o: lintest.cc linsyst.h linear.h
mptest$(OBJ_SUF).o: interface.h templates.h marith.h mptest.cc arith.h xmod.h
mvectest$(OBJ_SUF).o: interface.h templates.h mvectest.cc mvector.h marith.h arith.h xmod.h
mvector$(OBJ_SUF).o: mvector.cc mvector.h vector.h vec.h interface.h templates.h \
  arith.h xmod.h marith.h
mmatrix$(OBJ_SUF).o: mmatrix.cc mmatrix.h matrix.h mat.h vec.h interface.h \
  templates.h arith.h xmod.h mvector.h marith.h
mmattest$(OBJ_SUF).o: interface.h templates.h mmattest.cc mmatrix.h mvector.h vector.h marith.h
msubspace$(OBJ_SUF).o: msubspace.cc msubspace.h mmatrix.h matrix.h mat.h vec.h \
  interface.h templates.h arith.h xmod.h mvector.h marith.h
mspace$(OBJ_SUF).o: interface.h templates.h mspace.cc msubspace.h mmatrix.h mvector.h marith.h

### svec is not used now
svec$(OBJ_SUF).o: svec.cc vector.cc  svec.h vector.h arith.h xmod.h
svectest$(OBJ_SUF).o: svectest.cc svec.h svec.cc vector.h vector.cc arith.h xmod.h
### svec is not used now

smatrix$(OBJ_SUF).o: smat.cc smat.h smatrix.cc smatrix.h vector.h vec.h matrix.h mat.h sub.h subspace.h arith.h xmod.h interface.h templates.h

smatrix2$(OBJ_SUF).o: smat2.cc smat2.h smatrix.cc smatrix.h smat.h svector.h svec.h vector.h vec.h matrix.h mat.h sub.h subspace.h arith.h xmod.h interface.h templates.h

svector(OBJ_SUF).o: svec.cc svec.h svector.cc svector.h vector.h vec.h matrix.h mat.h sub.h subspace.h arith.h xmod.h interface.h templates.h

smatrix_elim$(OBJ_SUF).o: smatrix_elim.cc smatrix_elim.h smat_elim.h smat.h sub.h \
  mat.h vec.h interface.h templates.h arith.h xmod.h smat_elim.cc

smatrix2_elim$(OBJ_SUF).o: smatrix2_elim.cc smatrix2_elim.h \
  smat2_elim.cc smat2_elim.h smatrix2.h smat2.h svector.h svec.h sub.h subspace.h \
  mat.h vec.h interface.h templates.h arith.h xmod.h smat_elim.cc

smattest$(OBJ_SUF).o: smattest.cc random.cc smat.h smatrix.h smat_elim.h smatrix_elim.h sub.h subspace.h mat.h matrix.h vec.h vector.h arith.h xmod.h subspace.h interface.h templates.h

smat2test$(OBJ_SUF).o: smat2test.cc random.cc smat2.h smatrix2.h smat2_elim.h smatrix2_elim.h svec.h svector.h sub.h subspace.h mat.h matrix.h vec.h vector.h arith.h xmod.h subspace.h interface.h templates.h

xsplit: $(XSPLIT_OBJS)

xsplit0$(OBJ_SUF).o: xsplit.cc xsplit.h method.h arith.h interface.h templates.h \
  xmod.h vector.h vec.h matrix.h mat.h subspace.h sub.h splitbase.h \
  smatrix_elim.h smat_elim.h smat.h
	$(CXX) $(CFLAGS) -DMETHOD=0 xsplit.cc -o xsplit0$(OBJ_SUF).o

xsplit1$(OBJ_SUF).o: xsplit.cc xsplit.h method.h arith.h interface.h templates.h \
  xmod.h vector.h vec.h matrix.h mat.h subspace.h sub.h splitbase.h \
  smatrix_elim.h smat_elim.h smat.h
	$(CXX) $(CFLAGS) -DMETHOD=1 xsplit.cc -o xsplit1$(OBJ_SUF).o

xsplit2$(OBJ_SUF).o: xsplit.cc xsplit.h method.h arith.h interface.h templates.h \
  xmod.h vector.h vec.h matrix.h mat.h subspace.h sub.h splitbase.h \
  smatrix_elim.h smat_elim.h smat.h
	$(CXX) $(CFLAGS) -DMETHOD=2 xsplit.cc -o xsplit2$(OBJ_SUF).o

xsplit3$(OBJ_SUF).o: xsplit.cc xsplit.h method.h arith.h interface.h templates.h \
  xmod.h vector.h vec.h matrix.h mat.h subspace.h sub.h marith.h \
  mvector.h mmatrix.h msubspace.h splitbase.h smatrix_elim.h smat_elim.h \
  smat.h
	$(CXX) $(CFLAGS) -DMETHOD=3 xsplit.cc -o xsplit3$(OBJ_SUF).o

xsplit4$(OBJ_SUF).o: xsplit.cc xsplit.h method.h arith.h interface.h templates.h \
  xmod.h vector.h vec.h matrix.h mat.h subspace.h sub.h marith.h \
  mvector.h mmatrix.h msubspace.h splitbase.h smatrix_elim.h smat_elim.h \
  smat.h
	$(CXX) $(CFLAGS) -DMETHOD=4 xsplit.cc -o xsplit4$(OBJ_SUF).o

xsplit5$(OBJ_SUF).o: xsplit.cc xsplit.h method.h arith.h interface.h templates.h \
  xmod.h vector.h vec.h matrix.h mat.h subspace.h sub.h splitbase.h \
  smatrix_elim.h smat_elim.h smat.h
	$(CXX) $(CFLAGS) -DMETHOD=5 xsplit.cc -o xsplit5$(OBJ_SUF).o

xsplit6$(OBJ_SUF).o: xsplit.cc xsplit.h method.h arith.h interface.h templates.h \
  xmod.h vector.h vec.h matrix.h mat.h subspace.h sub.h splitbase.h \
  smatrix_elim.h smat_elim.h smat.h
	$(CXX) $(CFLAGS) -DMETHOD=6 xsplit.cc -o xsplit6$(OBJ_SUF).o

linalg: $(LINALG_OBJS)
vector$(OBJ_SUF).o: vector.cc vec.h interface.h templates.h arith.h xmod.h vec.cc
svector$(OBJ_SUF).o: svector.cc svector.h svec.cc svec.h interface.h templates.h arith.h xmod.h vec.cc vec.h
matrix$(OBJ_SUF).o: matrix.cc mat.h vec.h interface.h templates.h arith.h xmod.h \
  mat.cc
subspace$(OBJ_SUF).o: subspace.cc sub.h mat.h vec.h interface.h templates.h arith.h \
  xmod.h sub.cc


#template for compiling *.cc -> *$(OBJ_SUF).o :

%$(OBJ_SUF).o:   %.cc
	$(CXX) $(CFLAGS) $< -o $@

rcubic: rcubic$(OBJ_SUF).o cubic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o realroots$(OBJ_SUF).o
	$(CXX) -o rcubic rcubic$(OBJ_SUF).o cubic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o realroots$(OBJ_SUF).o $(LFLAGS)

lcubic: lcubic$(OBJ_SUF).o cubic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o realroots$(OBJ_SUF).o
	$(CXX) -o lcubic lcubic$(OBJ_SUF).o cubic$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o parifact$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o interface$(OBJ_SUF).o realroots$(OBJ_SUF).o $(LFLAGS)

quadratic$(OBJ_SUF).o: quadratic.cc marith.h interface.h templates.h arith.h xmod.h \
  unimod.h quadratic.h
cubic$(OBJ_SUF).o: cubic.cc marith.h interface.h templates.h arith.h xmod.h \
  unimod.h cubic.h
rcubic$(OBJ_SUF).o: cubic.h rcubic.cc
lcubic$(OBJ_SUF).o: cubic.h lcubic.cc

testlat$(OBJ_SUF).o: testlat.cc
testlat: testlat$(OBJ_SUF).o
	 $(CXX) -o testlat testlat$(OBJ_SUF).o $(LFLAGS)

illl$(OBJ_SUF).o: illl.cc mvector.h vector.h vec.h interface.h templates.h arith.h \
  xmod.h marith.h illl.h
#ifdef NTL_INTS
#	$(CXX) -c -g $(NTLCFLAGS) illl.cc
#endif
# since with -O gcc-2.95 gives an internal compiler error!
# -- doesn't any more

tilll$(OBJ_SUF).o: illl.h tilll.cc mmatrix.h mvector.h marith.h
tilll: tilll$(OBJ_SUF).o illl$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o mvector$(OBJ_SUF).o mmatrix$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o \
                               arith$(OBJ_SUF).o compproc$(OBJ_SUF).o 
	$(CXX) -o tilll tilll$(OBJ_SUF).o illl$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o mvector$(OBJ_SUF).o mmatrix$(OBJ_SUF).o \
                       matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o $(LFLAGS)
leg_lll$(OBJ_SUF).o: illl.h tilll.cc mmatrix.h mvector.h marith.h legendre.h conic.h
leg_lll: leg_lll$(OBJ_SUF).o illl$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o mvector$(OBJ_SUF).o mmatrix$(OBJ_SUF).o matrix$(OBJ_SUF).o vector$(OBJ_SUF).o \
                                   arith$(OBJ_SUF).o compproc$(OBJ_SUF).o legendre$(OBJ_SUF).o conic$(OBJ_SUF).o
	$(CXX) -o leg_lll leg_lll$(OBJ_SUF).o illl$(OBJ_SUF).o marith$(OBJ_SUF).o gpslave$(OBJ_SUF).o mvector$(OBJ_SUF).o mmatrix$(OBJ_SUF).o \
           matrix$(OBJ_SUF).o vector$(OBJ_SUF).o arith$(OBJ_SUF).o compproc$(OBJ_SUF).o legendre$(OBJ_SUF).o conic$(OBJ_SUF).o $(LFLAGS)

timer$(OBJ_SUF).o: timer.cc interface.h templates.h timer.h

new.o: new.C new.H
	$(CXX) -g $(CFLAGS) -DHAVE_ATEXIT new.C


include ../Makefile.dynamic

