## Process this file with automake to produce Makefile.in

if NO_MPFP
MPFP_SWITCH="-DNO_MPFP"
else
MPFP_SWITCH=
endif

LDFLAGS_BOOST = $(BOOST_LDFLAGS) $(BOOST_ASIO_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
BOOST_LIBS = $(BOOST_ASIO_LIB) $(BOOST_THREAD_LIB)

AM_CPPFLAGS = $(FLINT_CFLAGS) $(NTL_CFLAGS) $(PARI_CFLAGS) $(BOOST_CPPFLAGS) $(MPFP_SWITCH) $(PTHREAD_CFLAGS)
AM_LDFLAGS  = $(FLINT_LDFLAGS) $(NTL_LDFLAGS) $(PARI_LDFLAGS) $(LDFLAGS_BOOST) $(PTHREAD_CFLAGS) $(PTHREAD_LDFLAGS)
LIBS        = $(FLINT_LIBS) $(NTL_LIBS) $(PARI_LIBS) $(BOOST_LIBS) $(PTHREAD_LIBS)

lib_LTLIBRARIES = libec.la

PROCS_DOTHS = eclib/interface.h eclib/templates.h eclib/pari_init.h eclib/int.h eclib/convert.h eclib/arith.h eclib/xmod.h eclib/modulus.h eclib/marith.h eclib/bigcomplex.h eclib/compproc.h eclib/vector.h eclib/matrix.h eclib/subspace.h eclib/rat.h eclib/bigrat.h eclib/kbessel.h eclib/svector.h eclib/smatrix.h eclib/ssubspace.h eclib/smatrix_elim.h eclib/splitbase.h eclib/xsplit.h eclib/conic.h eclib/legendre.h eclib/quadratic.h eclib/unimod.h eclib/illl.h eclib/hilbert.h eclib/timer.h eclib/cubic.h eclib/gf.h eclib/polys.h eclib/realroots.h eclib/pari_interface.h eclib/p2points.h eclib/xsplit_data.h eclib/threadpool.h eclib/logger.h eclib/linalg.h eclib/flinterface.h

QCURVES_DOTHS = eclib/curve.h eclib/points.h eclib/cperiods.h eclib/isogs.h eclib/reader.h eclib/mwprocs.h eclib/lambda.h eclib/sifter.h eclib/sieve_search.h eclib/htconst.h eclib/egr.h eclib/saturate.h eclib/divpol.h eclib/pointsmod.h eclib/curvemod.h eclib/ffmod.h eclib/tlss.h eclib/elog.h eclib/getcurve.h eclib/egros.h

QRANK_DOTHS = eclib/mrank2.h eclib/mequiv.h eclib/mlocsol.h eclib/mglobsol.h eclib/msoluble.h eclib/mquartic.h eclib/mrank1.h eclib/qc.h eclib/sqfdiv.h eclib/version.h eclib/minim.h eclib/reduce.h eclib/transform.h eclib/desc2.h eclib/bitspace.h eclib/options.h eclib/GetOpt.h eclib/twoadic.h eclib/descent.h

G0N_DOTHS = eclib/periods.h eclib/cusp.h eclib/homspace.h eclib/newforms.h eclib/moddata.h eclib/oldforms.h eclib/symb.h eclib/fixc6.h eclib/pcprocs.h eclib/nfd.h eclib/curvesort.h

DOTHS = $(PROCS_DOTHS) $(QCURVES_DOTHS) $(QRANK_DOTHS) $(G0N_DOTHS)

PROCS_DOTCCS = interface.cc int.cc pari_init.cc convert.cc unimod.cc modulus.cc arith.cc marith.cc compproc.cc vector.cc matrix.cc subspace.cc kbessel.cc svector.cc smatrix.cc ssubspace.cc smatrix_elim.cc xsplit.cc conic.cc legendre.cc quadratic.cc illl.cc hilbert.cc timer.cc cubic.cc polys.cc realroots.cc pari_interface.cc p2points.cc gf.cc xsplit_data.cc threadpool.cc logger.cc

QCURVES_DOTCCS = curve.cc curvedata.cc curvered.cc points.cc cperiods.cc isogs.cc heights.cc mwprocs.cc lambda.cc sifter.cc sieve_search.cc htconst.cc egr.cc saturate.cc divpol.cc pointsmod.cc curvemod.cc ffmod.cc tlss.cc elog.cc getcurve.cc egros.cc

QRANK_DOTCCS = mequiv.cc  mrank1.cc mlocsol.cc mglobsol.cc mquartic.cc mrank2.cc qc.cc sqfdiv.cc version.cc minim.cc reduce.cc transform.cc desc2.cc bitspace.cc GetOpt.cc twoadic.cc descent.cc

G0N_DOTCCS = newforms.cc symb.cc homspace.cc cusp.cc oldforms.cc fixc6.cc periods.cc moddata.cc pcprocs.cc nfd.cc curvesort.cc

DOTCCS = $(PROCS_DOTCCS) $(QCURVES_DOTCCS) $(QRANK_DOTCCS) $(G0N_DOTCCS)

# NB Do NOT include the following file(s) above, since they are not
# compiled by themselves but only included in other files which are
# compiled.  But they do need to be in the distribution made by "make
# dist" so we must list them here:

EXTRA_DIST = random.cc

libec_la_SOURCES = $(DOTCCS)
libec_la_LDFLAGS = -no-undefined $(AM_LDFLAGS) -version-info @LT_CURRENT@:@LT_REVISION@:@LT_AGE@
nobase_include_HEADERS = $(DOTHS)
ACLOCAL_AMFLAGS = -I m4
