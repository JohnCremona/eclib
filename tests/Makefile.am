## Process this file with automake to produce Makefile.in

AM_CPPFLAGS     = -I$(top_srcdir)/libsrc $(NTL_CFLAGS) $(PARI_CFLAGS) $(GMP_CFLAGS) -DNTL_ALL
LDADD           = $(top_builddir)/libsrc/libec.la $(PARI_LIBS) $(NTL_LIBS) $(GMP_LIBS)

test_input_dir = $(srcdir)/in
test_output_dir = $(srcdir)/out

################ PROCS #############################################

PROCS_TESTS = smattest comptest rattest bigrattest ptest mptest tbessel mvectest mmattest mspace thilbert tp2points tilll vectest1 vectest2 mattest1 mattest2 space1 space2 svectest1 svectest2 tcon2 tlegcert
procs_tests: $(PROCS_TESTS)

# NB tcon2 and tlegcert are not currently tested in check_procs

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

check_procs: procs_tests
	@echo Running procs checks...
	./vectest1 < $(test_input_dir)/vectest.in >  vectest1.testout && diff vectest1.testout $(test_output_dir)/vectest.out
	./vectest2 < $(test_input_dir)/vectest.in >  vectest2.testout && diff vectest2.testout $(test_output_dir)/vectest.out
	./mattest1 < $(test_input_dir)/mattest.in >  mattest1.testout && diff mattest1.testout $(test_output_dir)/mattest.out
	./mattest2 < $(test_input_dir)/mattest.in >  mattest2.testout && diff mattest2.testout $(test_output_dir)/mattest.out
	./space1 < $(test_input_dir)/space.in > space1.testout && diff space1.testout $(test_output_dir)/space.out
	./space2 < $(test_input_dir)/space.in > space2.testout && diff space2.testout $(test_output_dir)/space.out
	./svectest1 < $(test_input_dir)/svectest.in > svectest1.testout && diff svectest1.testout $(test_output_dir)/svectest.out
	./svectest2 < $(test_input_dir)/svectest.in > svectest2.testout && diff svectest2.testout $(test_output_dir)/svectest.out
	./smattest < $(test_input_dir)/smattest.in > smattest.testout && diff smattest.testout $(test_output_dir)/smattest.out
	./mvectest < $(test_input_dir)/mvectest.in > mvectest.testout && diff mvectest.testout $(test_output_dir)/mvectest.out
	./mmattest < $(test_input_dir)/mmattest.in > mmattest.testout && diff mmattest.testout $(test_output_dir)/mmattest.out
	./mspace < $(test_input_dir)/mspace.in > mspace.testout && diff mspace.testout $(test_output_dir)/mspace.out
	./rattest < $(test_input_dir)/rattest.in > rattest.testout && diff rattest.testout $(test_output_dir)/rattest.out
	./bigrattest < $(test_input_dir)/bigrattest.in > bigrattest.testout && diff bigrattest.testout $(test_output_dir)/bigrattest.out
	./ptest < $(test_input_dir)/ptest.in > ptest.testout && diff ptest.testout $(test_output_dir)/ptest.out
	./mptest < $(test_input_dir)/mptest.in > mptest.testout && diff mptest.testout $(test_output_dir)/mptest.out
	./thilbert < $(test_input_dir)/thilbert.in > thilbert.testout && diff thilbert.testout $(test_output_dir)/thilbert.out
	./tp2points < $(test_input_dir)/tp2points.in > tp2points.testout && diff tp2points.testout $(test_output_dir)/tp2points$(OBJ_SUF).out

################ QCURVES #############################################

QCURVES_TESTS = tcurve tpoints theight thtconst tegr telog tdivpol tperiods tsat tsat2 tsat3 conductor tate torsion twist allisog
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

check_qcurves: qcurves_tests
	@echo Running qcurves checks...
	./tcurve < $(test_input_dir)/tcurve.in > tcurve.testout && diff tcurve.testout $(test_output_dir)/tcurve.out
	./tpoints > tpoints.testout && diff tpoints.testout $(test_output_dir)/tpoints.out
	./theight < $(test_input_dir)/theight.in > theight.testout && diff theight.testout $(test_output_dir)/theight.out
	./thtconst < $(test_input_dir)/tcurve.in > thtconst.testout && diff thtconst.testout $(test_output_dir)/thtconst.out
	./tegr < $(test_input_dir)/tegr.in > tegr.testout && diff tegr.testout $(test_output_dir)/tegr.out
	./telog > telog.testout && diff telog.testout $(test_output_dir)/telog.out
	./tdivpol > tdivpol.testout && diff tdivpol.testout $(test_output_dir)/tdivpol.out
	./tperiods < $(test_input_dir)/tperiods.in > tperiods.testout 2>/dev/null && diff tperiods.testout $(test_output_dir)/tperiods.out
	./tsat < $(test_input_dir)/tsat.in > tsat.testout 2>/dev/null && diff tsat.testout $(test_output_dir)/tsat.out
	./tsat2 < $(test_input_dir)/tsat2.in > tsat2.testout 2>/dev/null && diff tsat2.testout $(test_output_dir)/tsat2.out
	./tsat3 < $(test_input_dir)/tsat3.in > tsat3.testout 2>/dev/null && diff tsat3.testout $(test_output_dir)/tsat3.out
	./conductor < $(test_input_dir)/conductor.in > conductor.testout 2>/dev/null && diff conductor.testout $(test_output_dir)/conductor.out
	./tate < $(test_input_dir)/tate.in > tate.testout 2>/dev/null && diff tate.testout $(test_output_dir)/tate.out
	./torsion < $(test_input_dir)/torsion.in > torsion.testout 2>/dev/null && diff torsion.testout $(test_output_dir)/torsion.out
	./twist < $(test_input_dir)/twist.in > twist.testout 2>/dev/null && diff twist.testout $(test_output_dir)/twist.out
	./allisog < $(test_input_dir)/allisog.in > allisog.testout 2>/dev/null && diff allisog.testout $(test_output_dir)/allisog.out

################ QRANK #############################################

QRANK_TESTS = tmrank tequiv d2
qrank_tests: $(QRANK_TESTS)
# NB d2 is not currently tested in check_qrank

tmrank_SOURCES=tmrank.cc
tequiv_SOURCES=tequiv.cc

check_qrank: qrank_tests
	@echo Running qrank checks...
	./tequiv < $(test_input_dir)/tequiv.in > tequiv.testout && diff tequiv.testout $(test_output_dir)/tequiv.out
	./tmrank < $(test_input_dir)/tmrank-short.in > tmrank.testout 2>/dev/null && diff tmrank.testout $(test_output_dir)/tmrank-short.out
	rm -f PRIMES 1

################ G0N #############################################

G0N_TESTS = modtest homtest hecketest nftest mhcount oftest tmanin h1degphi
g0n_tests: $(G0N_TESTS)

# NB aplist, checkgens and tnfd are not currently tested in check_g0n

modtest_SOURCES=modtest.cc
homtest_SOURCES=homtest.cc
hecketest_SOURCES=hecketest.cc
nftest_SOURCES=nftest.cc
mhcount_SOURCES=mhcount.cc
oftest_SOURCES=oftest.cc
tmanin_SOURCES=tmanin.cc
h1degphi_SOURCES=h1degphi.cc

# This is the name of a directory which will be used for newform
# input/output in the "make check" tests and removed afterwards.

export NF_DIR:=nftmp

check_g0n:  g0n_tests
	@echo Running g0n checks...
	./modtest < $(test_input_dir)/modtest.in > modtest.testout 2>/dev/null && diff modtest.testout $(test_output_dir)/modtest.out
	./homtest < $(test_input_dir)/homtest.in > homtest.testout && diff homtest.testout $(test_output_dir)/homtest.out
	./hecketest < $(test_input_dir)/hecketest.in > hecketest.testout 2>/dev/null && diff hecketest.testout $(test_output_dir)/hecketest.out
	./mhcount < $(test_input_dir)/mhcount.in > mhcount.testout && diff mhcount.testout $(test_output_dir)/mhcount.out
	rm -rf $(NF_DIR)
	mkdir $(NF_DIR)
	./tmanin < $(test_input_dir)/tmanin.in > tmanin.testout 2>/dev/null && diff tmanin.testout $(test_output_dir)/tmanin.out
	./nftest < $(test_input_dir)/nftest.in 2>/dev/null > nftest.testout && diff nftest.testout $(test_output_dir)/nftest.out
	./oftest < $(test_input_dir)/oftest.in > oftest.testout && diff oftest.testout $(test_output_dir)/oftest.out
	rm -rf $(NF_DIR)

tests: procs_tests qcurves_tests qrank_tests g0n_tests
all: tests
check: check_procs check_qcurves check_qrank check_g0n
	rm -f *.testout PRIMES 1

extra_progs = aplist checkgens tnfd d2 tcon2 tlegcert
aplist_SOURCES=aplist.cc
checkgens_SOURCES=checkgens.cc
tnfd_SOURCES=tnfd.cc
d2_SOURCES=d2.cc
tcon2_SOURCES=tcon2.cc
tlegcert_SOURCES=tlegcert.cc

noinst_PROGRAMS=$(PROCS_TESTS) $(QCURVES_TESTS) $(QRANK_TESTS) $(G0N_TESTS) $(extra_progs)

# We list here the directories in and out which contain the test input
# and expected output files for the tests run by "make check".

EXTRA_DIST = in out

ACLOCAL_AMFLAGS = -I m4
