## Process this file with automake to produce Makefile.in

if NO_MPFP
test_input_dir = $(srcdir)/in_no_ntl
test_output_dir = $(srcdir)/out_no_ntl
MPFP_SWITCH="-DNO_MPFP"
else
test_input_dir = $(srcdir)/in_ntl
test_output_dir = $(srcdir)/out_ntl
MPFP_SWITCH=
endif

LDFLAGS_BOOST = $(BOOST_LDFLAGS) $(BOOST_ASIO_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
BOOST_LIBS = $(BOOST_ASIO_LIB) $(BOOST_THREAD_LIB) 

AM_CPPFLAGS     = -I$(top_srcdir)/libsrc $(FLINT_CFLAGS) $(NTL_CFLAGS) $(PARI_CFLAGS) $(BOOST_CPPFLAGS) $(MPFP_SWITCH) $(PTHREAD_CFLAGS)
AM_LDFLAGS         = $(FLINT_LDFLAGS) $(NTL_LDFLAGS) $(PARI_LDFLAGS) $(LDFLAGS_BOOST) $(PTHREAD_CFLAGS) $(PTHREAD_LDFLAGS)
LIBS            = $(top_builddir)/libsrc/libec.la $(FLINT_LIBS) $(PARI_LIBS) $(NTL_LIBS) $(BOOST_LIBS) $(PTHREAD_LIBS)

check_run = echo -n "Testing $${prog}..."; ./$${prog}$(EXEEXT) < $(test_input_dir)/$${prog}.in > $${prog}.testout 2>/dev/null && echo "$${prog} completed" && diff $${prog}.testout $(test_output_dir)/$${prog}.out || exit $$?

################ PROCS #############################################

PROCS_TESTS = smattest1 smattest2 smattest3 comptest rattest bigrattest ptest mptest tbessel mvectest mmattest mspace thilbert tp2points tilll vectest1 vectest2 vectest3 mattest1 mattest2 mattest3 space1 space2 space3 svectest1 svectest2 svectest3 tcon2 tlegcert timing threading tversion
PROCS_TESTS_SUFFIXED = $(PROCS_TESTS:=$(EXEEXT))
procs_tests: $(PROCS_TESTS_SUFFIXED)

ptest_SOURCES=ptest.cc
mptest_SOURCES=mptest.cc
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
vectest3_SOURCES=vectest.cc
vectest3_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=3
mattest1_SOURCES=mattest.cc
mattest1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
mattest2_SOURCES=mattest.cc
mattest2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2
mattest3_SOURCES=mattest.cc
mattest3_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=3
space1_SOURCES=space.cc
space1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
space2_SOURCES=space.cc
space2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2
space3_SOURCES=space.cc
space3_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=3
svectest1_SOURCES=svectest.cc
svectest1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
svectest2_SOURCES=svectest.cc
svectest2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2
svectest3_SOURCES=svectest.cc
svectest3_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=3
smattest1_SOURCES=smattest.cc
smattest1_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=1
smattest2_SOURCES=smattest.cc
smattest2_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=2
smattest3_SOURCES=smattest.cc
smattest3_CPPFLAGS=${AM_CPPFLAGS} -DSCALAR_OPTION=3
tcon2_SOURCES=tcon2.cc
tlegcert_SOURCES=tlegcert.cc
timing_SOURCES=timing.cc
threading_SOURCES=threading.cc
tversion_SOURCES=tversion.cc

check_procs: procs_tests
	@echo Running procs checks...
	for prog in $(PROCS_TESTS); do $(check_run); done

################ QCURVES #############################################

QCURVES_TESTS = tcurve tpoints theight thtconst tlatconst tegr telog tdivpol tperiods tsatbnd tsat tsat2 tsat3 conductor tate torsion twist allisog checkgens tegros
QCURVES_TESTS_SUFFIXED = $(QCURVES_TESTS:=$(EXEEXT))
qcurves_tests: $(QCURVES_TESTS_SUFFIXED)

tcurve_SOURCES=tcurve.cc
tpoints_SOURCES=tpoints.cc
theight_SOURCES=theight.cc
thtconst_SOURCES=thtconst.cc
tlatconst_SOURCES=tlatconst.cc
tegr_SOURCES=tegr.cc
telog_SOURCES=telog.cc
tdivpol_SOURCES=tdivpol.cc
tperiods_SOURCES=tperiods.cc
tsatbnd_SOURCES=tsatbnd.cc
tsat_SOURCES=tsat.cc
tsat2_SOURCES=tsat2.cc
tsat3_SOURCES=tsat3.cc
conductor_SOURCES=conductor.cc
tate_SOURCES=tate.cc
torsion_SOURCES=torsion.cc
twist_SOURCES=twist.cc
allisog_SOURCES=allisog.cc
checkgens_SOURCES=checkgens.cc
tegros_SOURCES=tegros.cc

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
export SNF_DIR:=snftmp
export TCURVE_DIR:=tcurves

check_g0n: g0n_tests
	@echo Running g0n checks...
	rm -rf $(NF_DIR)
	rm -rf $(SNF_DIR)
	rm -rf $(TCURVE_DIR)
	mkdir $(NF_DIR)
	mkdir $(SNF_DIR)
	mkdir $(TCURVE_DIR)
	for prog in $(G0N_TESTS); do $(check_run); done
	rm -rf $(NF_DIR)
	rm -rf $(SNF_DIR)
	rm -rf $(TCURVE_DIR)

check:
	$(MAKE) check_procs check_qcurves check_qrank check_g0n
	rm -f *.testout PRIMES 1


check_PROGRAMS=$(PROCS_TESTS) $(QCURVES_TESTS) $(QRANK_TESTS) $(G0N_TESTS)

# We list here the directories in and out which contain the test input
# and expected output files for the tests run by "make check".

EXTRA_DIST = in_ntl in_no_ntl out_ntl out_no_ntl

ACLOCAL_AMFLAGS = -I m4
