## Process this file with automake to produce Makefile.in

AM_CPPFLAGS     = -I../include $(NTL_CFLAGS) $(PARI_CFLAGS)
AM_LDFLAGS      = -ljc -L../lib $(NTL_LIBS) $(PARI_LIBS)

PROCS_TESTS = smattest comptest rattest bigrattest ptest mptest tbessel mvectest mmattest mspace tconic tleg thilbert rcubic lcubic tp2points tilll vectest1 vectest2 mattest1 mattest2 space1 space2 svectest1 svectest2

procs_tests: $(PROCS_TESTS)

ptest_SOURCES=ptest.cc
mptest_SOURCES=mptest.cc
smattest_SOURCES=smattest.cc
comptest_SOURCES=comptest.cc
rattest_SOURCES=rattest.cc
bigrattest_SOURCES=bigrattest.cc
mvectest_SOURCES=mvectest.cc
mmattest_SOURCES=mmattest.cc
mspace_SOURCES=mspace.cc
tbessel_SOURCES=tbessel.cc
tconic_SOURCES=tconic.cc
tleg_SOURCES=tleg.cc
thilbert_SOURCES=thilbert.cc
rcubic_SOURCES=rcubic.cc
lcubic_SOURCES=lcubic.cc
tp2points_SOURCES=tp2points.cc
tilll_SOURCES=tilll.cc
vectest1_SOURCES=vectest.cc
vectest1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
vectest2_SOURCES=vectest.cc
vectest2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2
mattest1_SOURCES=mattest.cc
mattest1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
mattest2_SOURCES=mattest.cc
mattest2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2
space1_SOURCES=space.cc
space1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
space2_SOURCES=space.cc
space2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2
svectest1_SOURCES=svectest.cc
svectest1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
svectest2_SOURCES=svectest.cc
svectest2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2

check_procs: procs_tests
	./vectest1 < in/vectest.in >  t && diff t out/vectest.out
	./vectest2 < in/vectest.in >  t && diff t out/vectest.out
	./mattest1 < in/mattest.in >  t && diff t out/mattest.out
	./mattest2 < in/mattest.in >  t && diff t out/mattest.out
	./space1 < in/space.in > t && diff t out/space.out
	./space2 < in/space.in > t && diff t out/space.out
	./svectest1 < in/svectest.in > t && diff t out/svectest.out
	./svectest2 < in/svectest.in > t && diff t out/svectest.out
	./smattest < in/smattest.in > t && diff t out/smattest.out
	./mvectest < in/mvectest.in > t && diff t out/mvectest.out
	./mmattest < in/mmattest.in > t && diff t out/mmattest.out
	./mspace < in/mspace.in > t && diff t out/mspace.out
	./rattest < in/rattest.in > t && diff t out/rattest.out
	./bigrattest < in/bigrattest.in > t && diff t out/bigrattest.out
	./ptest < in/ptest.in > t && diff t out/ptest.out
	./mptest < in/mptest.in > t && diff t out/mptest.out
	./tconic < in/tconic.in > t && diff t out/tconic.out
	./tleg < in/tleg.in > t && diff t out/tleg.out
	./thilbert < in/thilbert.in > t && diff t out/thilbert.out
	./rcubic < in/rcubic.in > t && diff t out/rcubic$(OBJ_SUF).out
	./lcubic < in/lcubic.in > t && diff t out/lcubic.out
	./tp2points < in/tp2points.in > t && diff t out/tp2points$(OBJ_SUF).out
	rm -f PRIMES t 1

QCURVES_TESTS = tcurve tpoints theight thtconst tegr telog tdivpol tperiods tsat tsat2 tsat3 conductor tate torsion twist allisog findinf indep
qcurves_tests: $(QCURVES_TESTS)

tcurve_SOURCES=tcurve.cc 
tpoints_SOURCES=tpoints.cc 
theight_SOURCES=theight.cc 
thtconst_SOURCES=thtconst.cc 
tegr_SOURCES=tegr.cc 
telog_SOURCES=telog.cc 
tdivpol_SOURCES=tdivpol.cc 
tperiods_SOURCES=tperiods.cc 
tsat_SOURCES=tsat.cc 
tsat2_SOURCES=tsat2.cc 
tsat3_SOURCES=tsat3.cc 
conductor_SOURCES=conductor.cc 
tate_SOURCES=tate.cc 
torsion_SOURCES=torsion.cc 
twist_SOURCES=twist.cc 
allisog_SOURCES=allisog.cc 
findinf_SOURCES=findinf.cc 
indep_SOURCES=indep.cc

check_qcurves: qcurves_tests
	./tcurve < in/tcurve.in > t && diff t out/tcurve.out
	./tpoints > t && diff t out/tpoints.out
	./theight < in/theight.in > t && diff t out/theight.out
	./thtconst < in/tcurve.in > t && diff t out/thtconst.out
	./tegr < in/tegr.in > t && diff t out/tegr.out
	./telog > t && diff t out/telog.out
	./tdivpol > t && diff t out/tdivpol.out
	./tperiods < in/tperiods.in > t 2>/dev/null && diff t out/tperiods.out
	./tsat < in/tsat.in > t 2>/dev/null && diff t out/tsat.out
	./tsat2 < in/tsat2.in > t 2>/dev/null && diff t out/tsat2.out
	./tsat3 < in/tsat3.in > t 2>/dev/null && diff t out/tsat3.out
	./conductor < in/conductor.in > t 2>/dev/null && diff t out/conductor.out
	./tate < in/tate.in > t 2>/dev/null && diff t out/tate.out
	./torsion < in/torsion.in > t 2>/dev/null && diff t out/torsion.out
	./twist < in/twist.in > t 2>/dev/null && diff t out/twist.out
	./allisog < in/allisog.in > t 2>/dev/null && diff t out/allisog.out
	./findinf < in/findinf.in > t 2>/dev/null && diff t out/findinf.out
	./indep < in/indep.in > t 2>/dev/null && diff t out/indep.out
	rm -f PRIMES t 1

QRANK_TESTS = tmrank tequiv tred mwrank ratpoint
qrank_tests: $(QRANK_TESTS)

tmrank_SOURCES=tmrank.cc
tequiv_SOURCES=tequiv.cc
tred_SOURCES=tred.cc
mwrank_SOURCES=mwrank.cc
ratpoint_SOURCES=ratpoint.cc

check_qrank: qrank_tests
	./tequiv < in/tequiv.in > t && diff t out/tequiv.out
	./tred < in/tred.in > t 2>/dev/null && diff t out/tred.out
	./ratpoint < in/ratpoint.in > t 2>/dev/null && diff t out/ratpoint.out
	./mwrank < in/mwrank.in > t 2>/dev/null && diff t out/mwrank.out
	rm -f PRIMES t 1

G0N_TESTS = modtest homtest hecketest nftest mhcount oftest tmanin h1bsd h1curve pcurve nfhpcurve nfhpmcurve h1clist checkap h1bsdcurisog h1degphi hpcurve nfcount ecnf h1first moreap
g0n_tests: $(G0N_TESTS)

modtest_SOURCES=modtest.cc
homtest_SOURCES=homtest.cc
hecketest_SOURCES=hecketest.cc
nftest_SOURCES=nftest.cc
mhcount_SOURCES=mhcount.cc
oftest_SOURCES=oftest.cc
tmanin_SOURCES=tmanin.cc
h1bsd_SOURCES=h1bsd.cc
h1curve_SOURCES=h1curve.cc
pcurve_SOURCES=pcurve.cc
nfhpcurve_SOURCES=nfhpcurve.cc
nfhpmcurve_SOURCES=nfhpmcurve.cc
h1clist_SOURCES=h1clist.cc
checkap_SOURCES=checkap.cc
h1bsdcurisog_SOURCES=h1bsdcurisog.cc
h1degphi_SOURCES=h1degphi.cc
hpcurve_SOURCES=hpcurve.cc
nfcount_SOURCES=nfcount.cc
ecnf_SOURCES=ecnf.cc
h1first_SOURCES=h1first.cc
moreap_SOURCES=moreap.cc

# This is the name of a directory which will be used for newform
# input/output in the "make check" tests and removed afterwards.

export NF_DIR:=nftmp

check_g0n:  g0n_tests
	./modtest < in/modtest.in > t 2>/dev/null && diff t out/modtest.out
	./homtest < in/homtest.in > t && diff t out/homtest.out
	./hecketest < in/hecketest.in > t 2>/dev/null && diff t out/hecketest.out
	./mhcount < in/mhcount.in > t && diff t out/mhcount.out
	rm -rf $(NF_DIR)
	mkdir $(NF_DIR)
	./tmanin < in/tmanin.in > t 2>/dev/null && diff t out/tmanin.out
	./nftest < in/nftest.in 2>/dev/null > t && diff t out/nftest.out
	./h1first < in/h1first.in 2>/dev/null > t && diff t out/h1first.out
	./ecnf < in/ecnf.in 2>/dev/null > t && diff t out/ecnf.out
	./oftest < in/oftest.in > t && diff t out/oftest.out
	./nfhpcurve < in/nfhpcurve.in > t 2>/dev/null && diff t out/nfhpcurve.out
	./nfhpmcurve < in/nfhpmcurve.in > t 2>/dev/null && diff t out/nfhpmcurve.out
	./h1clist < in/h1clist.in > t 2>/dev/null && diff t out/h1clist.out
	./h1bsdcurisog < in/h1bsdcurisog.in > t 2>/dev/null && diff t out/h1bsdcurisog.out
	rm -rf $(NF_DIR)
	rm -f PRIMES t

tests: procs_tests qcurves_tests qrank_tests g0n_tests
all: tests
check: check_procs check_qcurves check_qrank check_g0n


bin_PROGRAMS=$(PROCS_TESTS) $(QCURVES_TESTS) $(QRANK_TESTS) $(G0N_TESTS)


ACLOCAL_AMFLAGS = -I m4
