## Process this file with automake to produce Makefile.in

AM_CPPFLAGS     = -I$(top_srcdir)/libsrc $(NTL_CFLAGS) $(PARI_CFLAGS)
AM_LDFLAGS      = -ljc -L$(top_builddir)/libsrc $(NTL_LIBS) $(PARI_LIBS)

test_input_dir = $(srcdir)/in
test_output_dir = $(srcdir)/out

################ PROCS (conics & cubics programs) #################

PROCS_PROGS =  solve_conic solve_legendre reduce_cubics list_cubics
procs_progs: $(PROCS_PROGS)

solve_conic_SOURCES=solve_conic.cc
solve_legendre_SOURCES=solve_legendre.cc
reduce_cubics_SOURCES=reduce_cubics.cc
list_cubics_SOURCES=list_cubics.cc

check_procs: procs_progs
	@echo Checking conics and cubics programs...
	./solve_conic < $(test_input_dir)/solve_conic.in > solve_conic.testout && diff solve_conic.testout $(test_output_dir)/solve_conic.out
	./solve_legendre < $(test_input_dir)/solve_legendre.in > solve_legendre.testout && diff solve_legendre.testout $(test_output_dir)/solve_legendre.out
	./reduce_cubics < $(test_input_dir)/reduce_cubics.in > reduce_cubics.testout && diff reduce_cubics.testout $(test_output_dir)/reduce_cubics$(OBJ_SUF).out
	./list_cubics < $(test_input_dir)/list_cubics.in > list_cubics.testout && diff list_cubics.testout $(test_output_dir)/list_cubics.out

################ QCURVES (elliptic curves programs, excluding two-descent) #################

QCURVES_PROGS = point_search indep_test
qcurves_progs: $(QCURVES_TESTS)

point_search_SOURCES=point_search.cc
indep_test_SOURCES=indep_test.cc

check_qcurves: qcurves_progs
	@echo Checking qcurves programs...
	./point_search < $(test_input_dir)/point_search.in > point_search.testout 2>/dev/null && diff point_search.testout $(test_output_dir)/point_search.out
	./indep_test < $(test_input_dir)/indep_test.in > indep_test.testout 2>/dev/null && diff indep_test.testout $(test_output_dir)/indep_test.out

################ QCURVES (elliptic curves two-descent programs) #################

QRANK_PROGS = mwrank reduce_quartics quartic_points
qrank_progs: $(QRANK_PROGS)

mwrank_SOURCES=mwrank.cc
reduce_quartics_SOURCES=reduce_quartics.cc
quartic_points_SOURCES=quartic_points.cc

check_qrank: qrank_progs
	@echo Checking qrank programs...
	./mwrank < $(test_input_dir)/mwrank.in > mwrank.testout 2>/dev/null && diff mwrank.testout $(test_output_dir)/mwrank.out
	./reduce_quartics < $(test_input_dir)/reduce_quartics.in > reduce_quartics.testout 2>/dev/null && diff reduce_quartics.testout $(test_output_dir)/reduce_quartics.out
	./quartic_points < $(test_input_dir)/quartic_points.in > quartic_points.testout 2>/dev/null && diff quartic_points.testout $(test_output_dir)/quartic_points.out
	rm -f PRIMES 1

################ G0N (modular symbols programs) #################

G0N_PROGS = h1bsd h1curve pcurve nfhpcurve nfhpmcurve h1clist checkap h1bsdcurisog nfcount ecnf h1first moreap qexp
g0n_progs: $(G0N_PROGS)

# NB aplist, checkgens and tnfd are not currently tested in check_g0n

h1bsd_SOURCES=h1bsd.cc
h1curve_SOURCES=h1curve.cc
pcurve_SOURCES=pcurve.cc
nfhpcurve_SOURCES=nfhpcurve.cc
nfhpmcurve_SOURCES=nfhpmcurve.cc
h1clist_SOURCES=h1clist.cc
checkap_SOURCES=checkap.cc
h1bsdcurisog_SOURCES=h1bsdcurisog.cc
nfcount_SOURCES=nfcount.cc
ecnf_SOURCES=ecnf.cc
h1first_SOURCES=h1first.cc
moreap_SOURCES=moreap.cc
qexp_SOURCES=qexp.cc

# This is the name of a directory which will be used for newform
# input/output in the "make check" tests and removed afterwards.

export NF_DIR:=nftmp

check_g0n:  g0n_progs
	@echo Checking g0n programs...
	rm -rf $(NF_DIR)
	mkdir $(NF_DIR)
	./ecnf < $(test_input_dir)/ecnf.in 2>/dev/null > ecnf.testout && diff ecnf.testout $(test_output_dir)/ecnf.out
	./nfhpcurve < $(test_input_dir)/nfhpcurve.in > nfhpcurve.testout 2>/dev/null && diff nfhpcurve.testout $(test_output_dir)/nfhpcurve.out
	./nfhpmcurve < $(test_input_dir)/nfhpmcurve.in > nfhpmcurve.testout 2>/dev/null && diff nfhpmcurve.testout $(test_output_dir)/nfhpmcurve.out
	./h1first < $(test_input_dir)/h1first.in 2>/dev/null > h1first.testout && diff h1first.testout $(test_output_dir)/h1first.out
	./h1clist < $(test_input_dir)/h1clist.in > h1clist.testout 2>/dev/null && diff h1clist.testout $(test_output_dir)/h1clist.out
	./h1bsdcurisog < $(test_input_dir)/h1bsdcurisog.in > h1bsdcurisog.testout 2>/dev/null && diff h1bsdcurisog.testout $(test_output_dir)/h1bsdcurisog.out
	echo 0 25 1 `seq 50 60` 0 | ./nfhpcurve > /dev/null 2>/dev/null
	./qexp < $(test_input_dir)/qexp.in > qexp.testout 2>/dev/null && diff qexp.testout $(test_output_dir)/qexp.out
	rm -rf $(NF_DIR)

#################################################################################

progs: procs_progs qcurves_progs qrank_progs g0n_progs
all: progs
check: check_procs check_qcurves check_qrank check_g0n
	rm -f *.testout PRIMES 1

if ALLPROGS
bin_PROGRAMS=$(PROCS_PROGS) $(QCURVES_PROGS) $(QRANK_PROGS) $(G0N_PROGS)
else
bin_PROGRAMS=mwrank
endif

# We list here the directories in and out which contain the test input
# and expected output files for the tests run by "make check".

EXTRA_DIST = in out

ACLOCAL_AMFLAGS = -I m4
