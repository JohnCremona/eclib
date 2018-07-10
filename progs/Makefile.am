## Process this file with automake to produce Makefile.in

prog_input_dir = $(srcdir)/in
if NTL_ALL
NTL_SWITCH = -D NTL_ALL
prog_output_dir = $(srcdir)/out_ntl
else
prog_output_dir = $(srcdir)/out_no_ntl
endif

LDFLAGS_BOOST = $(BOOST_LDFLAGS) $(BOOST_ASIO_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
BOOST_LIBS = $(BOOST_ASIO_LIB) $(BOOST_THREAD_LIB) 

AM_CPPFLAGS     = -I$(top_srcdir)/libsrc $(FLINT_CFLAGS) $(NTL_CFLAGS) $(PARI_CFLAGS) $(BOOST_CPPFLAGS) $(NTL_SWITCH)
AM_LDFLAGS         = $(FLINT_LDFLAGS) $(NTL_LDFLAGS) $(PARI_LDFLAGS) $(BOOST_LDFLAGS)
LIBS            = $(top_builddir)/libsrc/libec.la $(FLINT_LIBS) $(PARI_LIBS) $(NTL_LIBS) $(BOOST_LIBS)

check_run = echo "Testing $${prog}..."; ./$${prog}$(EXEEXT) < $(prog_input_dir)/$${prog}.in > $${prog}.testout 2>/dev/null && diff $${prog}.testout $(prog_output_dir)/$${prog}.out || exit $$?


################ PROCS (conics & cubics programs) #################

PROCS_PROGS = solve_conic solve_legendre reduce_cubics list_cubics cubics
PROCS_PROGS_SUFFIXED = $(PROCS_PROGS:=$(EXEEXT))
procs_progs: $(PROCS_PROGS_SUFFIXED)

solve_conic_SOURCES=solve_conic.cc
solve_legendre_SOURCES=solve_legendre.cc
reduce_cubics_SOURCES=reduce_cubics.cc
list_cubics_SOURCES=list_cubics.cc
cubics_SOURCES=cubics.cc

check_procs: procs_progs
	@echo Checking conics and cubics programs...
	for prog in $(PROCS_PROGS); do $(check_run); done

################ QCURVES (elliptic curves programs, excluding two-descent) #################

QCURVES_PROGS = point_search indep_test
QCURVES_PROGS_SUFFIXED = $(QCURVES_PROGS:=$(EXEEXT))
qcurves_progs: $(QCURVES_PROGS_SUFFIXED)

point_search_SOURCES=point_search.cc
indep_test_SOURCES=indep_test.cc

check_qcurves: qcurves_progs
	@echo Checking qcurves programs...
	for prog in $(QCURVES_PROGS); do $(check_run); done

################ QCURVES (elliptic curves two-descent programs) #################

# When allprogs is disabled the only binary to be built and installed
# is mwrank, but we want "make check" to build and test all.

QRANK_PROGS = mwrank reduce_quartics quartic_points
QRANK_CHECK_PROGS = reduce_quartics quartic_points
QRANK_PROGS_SUFFIXED = $(QRANK_PROGS:=$(EXEEXT))
QRANK_CHECK_PROGS_SUFFIXED = $(QRANK_CHECK_PROGS:=$(EXEEXT))
qrank_progs: $(QRANK_PROGS_SUFFIXED)

mwrank_SOURCES=mwrank.cc
reduce_quartics_SOURCES=reduce_quartics.cc
quartic_points_SOURCES=quartic_points.cc

check_qrank: qrank_progs
	@echo Checking qrank programs...
	for prog in $(QRANK_PROGS); do $(check_run); done
	rm -f PRIMES 1

################ G0N (modular symbols programs) #################

G0N_PROGS = ecnf nfhpcurve nfhpmcurve h1first h1clist h1bsdcurisog qexp h1bsd h1curve pcurve checkap moreap nfcount aplist
G0N_PROGS_SUFFIXED = $(G0N_PROGS:=$(EXEEXT))
g0n_progs: $(G0N_PROGS_SUFFIXED)

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
aplist_SOURCES=aplist.cc

# This is the name of a directory which will be used for newform
# input/output in the "make check" tests and removed afterwards.
# This does not work without the export!

export NF_DIR:=nftmp
export SNF_DIR:=snftmp
export TCURVE_DIR:=tcurves

check_g0n: g0n_progs
	@echo Checking g0n programs...
	rm -rf $(NF_DIR)
	rm -rf $(SNF_DIR)
	rm -rf $(TCURVE_DIR)
	mkdir $(NF_DIR)
	mkdir $(SNF_DIR)
	mkdir $(TCURVE_DIR)
	for prog in $(G0N_PROGS); do $(check_run); done
	rm -rf $(NF_DIR)
	rm -rf $(SNF_DIR)
	rm -rf $(TCURVE_DIR)

#################################################################################

check:
	$(MAKE) check_procs check_qcurves check_qrank check_g0n
	rm -f *.testout PRIMES 1

# When allprogs is disabled the only binary to be built and installed
# is mwrank, but we want "make check" to build and test all.

if ALLPROGS
bin_PROGRAMS=$(PROCS_PROGS) $(QCURVES_PROGS) $(QRANK_PROGS) $(G0N_PROGS)
check_PROGRAMS=
else
bin_PROGRAMS=mwrank$(EXEEXT)
check_PROGRAMS=$(PROCS_PROGS) $(QCURVES_PROGS) $(QRANK_CHECK_PROGS) $(G0N_PROGS)
endif

# We list here the directories in and out which contain the test input
# and expected output files for the tests run by "make check".

EXTRA_DIST = in out_ntl out_no_ntl

ACLOCAL_AMFLAGS = -I m4
