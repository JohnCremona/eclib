## Process this file with automake to produce Makefile.in

LDFLAGS_BOOST = $(BOOST_LDFLAGS) $(BOOST_ASIO_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
BOOST_LIBS = $(BOOST_ASIO_LIB) $(BOOST_THREAD_LIB) 

AM_CPPFLAGS     = -I$(top_srcdir)/libsrc $(FLINT_CFLAGS) $(NTL_CFLAGS) $(PARI_CFLAGS) $(BOOST_CPPFLAGS) -DNTL_ALL
AM_LDFLAGS         = $(FLINT_LDFLAGS) $(NTL_LDFLAGS) $(PARI_LDFLAGS) $(BOOST_LDFLAGS)
LIBS            = $(top_builddir)/libsrc/libec.la $(FLINT_LIBS) $(PARI_LIBS) $(NTL_LIBS) $(BOOST_LIBS)

check_run = echo "Testing $${prog}..."; ./$${prog}$(EXEEXT) < $(test_input_dir)/$${prog}.in > $${prog}.testout 2>/dev/null && diff $${prog}.testout $(test_output_dir)/$${prog}.out || exit $$?

test_input_dir = $(srcdir)/in
test_output_dir = $(srcdir)/out

################ PROCS #############################################

PROCS_TESTS = smattest comptest rattest bigrattest ptest mptest tbessel mvectest mmattest mspace thilbert tp2points tilll vectest1 vectest2 mattest1 mattest2 space1 space2 svectest1 svectest2 tcon2 tlegcert timing threading
PROCS_TESTS_SUFFIXED = $(PROCS_TESTS:=$(EXEEXT))
procs_tests: $(PROCS_TESTS_SUFFIXED)

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
thilbert_SOURCES=thilbert.cc
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
tcon2_SOURCES=tcon2.cc
tlegcert_SOURCES=tlegcert.cc
timing_SOURCES=timing.cc
threading_SOURCES=threading.cc
threading_LDFLAGS=${AM_LDFLAGS} -lboost_timer

check_procs: procs_tests
	@echo Running procs checks...
	for prog in $(PROCS_TESTS); do $(check_run); done

################ QCURVES #############################################

QCURVES_TESTS = tcurve tpoints theight thtconst tegr telog tdivpol tperiods tsat tsat2 tsat3 conductor tate torsion twist allisog checkgens
QCURVES_TESTS_SUFFIXED = $(QCURVES_TESTS:=$(EXEEXT))
qcurves_tests: $(QCURVES_TESTS_SUFFIXED)

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
checkgens_SOURCES=checkgens.cc

check_qcurves: qcurves_tests
	@echo Running qcurves checks...
	cp $(test_input_dir)/allgens.268605 .
	for prog in $(QCURVES_TESTS); do $(check_run); done
	rm -f ./allgens.268605

################ QRANK #############################################

QRANK_TESTS = tmrank tequiv d2
QRANK_TESTS_SUFFIXED = $(QRANK_TESTS:=$(EXEEXT))
qrank_tests: $(QRANK_TESTS_SUFFIXED)

tmrank_SOURCES=tmrank.cc
tequiv_SOURCES=tequiv.cc
d2_SOURCES=d2.cc

check_qrank: qrank_tests
	@echo Running qrank checks...
	for prog in $(QRANK_TESTS); do $(check_run); done
	rm -f PRIMES 1

################ G0N #############################################

G0N_TESTS = modtest homtest hecketest mhcount tmanin nftest oftest tnfd 
G0N_TESTS_SUFFIXED = $(G0N_TESTS:=$(EXEEXT))
g0n_tests: $(G0N_TESTS_SUFFIXED)

# NB The order of the above does matter: tmanin must be run before
# nftest and oftest

modtest_SOURCES=modtest.cc
homtest_SOURCES=homtest.cc
hecketest_SOURCES=hecketest.cc
nftest_SOURCES=nftest.cc
mhcount_SOURCES=mhcount.cc
oftest_SOURCES=oftest.cc
tmanin_SOURCES=tmanin.cc
tnfd_SOURCES=tnfd.cc

# This is the name of a directory which will be used for newform
# input/output in the "make check" tests and removed afterwards.

export NF_DIR:=nftmp

check_g0n:  g0n_tests
	@echo Running g0n checks...
	rm -rf $(NF_DIR)
	mkdir $(NF_DIR)
	for prog in $(G0N_TESTS); do $(check_run); done
	rm -rf $(NF_DIR)

check: check_procs check_qcurves check_qrank check_g0n
	rm -f *.testout PRIMES 1


check_PROGRAMS=$(PROCS_TESTS) $(QCURVES_TESTS) $(QRANK_TESTS) $(G0N_TESTS)

# We list here the directories in and out which contain the test input
# and expected output files for the tests run by "make check".

EXTRA_DIST = in out

ACLOCAL_AMFLAGS = -I m4
