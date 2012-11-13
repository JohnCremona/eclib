## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(FLINT_CFLAGS) $(NTL_CFLAGS) $(PARI_CFLAGS) $(GMP_CFLAGS) -DNTL_ALL -DUSE_PARI_FACTORING -DNEW_OP_ORDER
LIBS        = $(FLINT_LIBS) $(NTL_LIBS) $(PARI_LIBS) $(GMP_LIBS)

lib_LTLIBRARIES = libec.la

PROCS_DOTHS = eclib/interface.h eclib/templates.h eclib/arith.h eclib/xmod.h eclib/marith.h eclib/gpslave.h eclib/compproc.h eclib/vec.h eclib/vector.h eclib/mat.h eclib/matrix.h eclib/sub.h eclib/subspace.h eclib/rat.h eclib/bigrat.h eclib/kbessel.h eclib/svec.h eclib/svector.h eclib/smat.h eclib/smatrix.h eclib/smat_elim.h eclib/smatrix_elim.h eclib/mvector.h eclib/mmatrix.h eclib/msubspace.h eclib/method.h eclib/splitbase.h eclib/xsplit.h eclib/conic.h eclib/legendre.h eclib/quadratic.h eclib/unimod.h eclib/illl.h eclib/hilbert.h eclib/timer.h eclib/cubic.h eclib/gf.h eclib/polys.h eclib/realroots.h eclib/parifact.h eclib/p2points.h

QCURVES_DOTHS = eclib/curve.h eclib/points.h eclib/cperiods.h eclib/isogs.h eclib/reader.h eclib/mwprocs.h eclib/lambda.h eclib/sifter.h eclib/sieve_search.h eclib/htconst.h eclib/egr.h eclib/saturate.h eclib/divpol.h eclib/pointsmod.h eclib/curvemod.h eclib/ffmod.h eclib/tlss.h eclib/elog.h eclib/getcurve.h

QRANK_DOTHS = eclib/mrank2.h eclib/mequiv.h eclib/mlocsol.h eclib/mglobsol.h eclib/msoluble.h eclib/mquartic.h eclib/mrank1.h eclib/qc.h eclib/sqfdiv.h eclib/version.h eclib/minim.h eclib/reduce.h eclib/transform.h eclib/desc2.h eclib/bitspace.h eclib/options.h eclib/GetOpt.h eclib/twoadic.h eclib/descent.h

G0N_DOTHS = eclib/periods.h eclib/cusp.h eclib/homspace.h eclib/newforms.h eclib/moddata.h eclib/oldforms.h eclib/symb.h eclib/fixc6.h eclib/pcprocs.h eclib/nfd.h eclib/curvesort.h

DOTHS = $(PROCS_DOTHS) $(QCURVES_DOTHS) $(QRANK_DOTHS) $(G0N_DOTHS)

PROCS_DOTCCS = interface.cc unimod.cc arith.cc marith.cc gpslave.cc compproc.cc vector.cc matrix.cc subspace.cc kbessel.cc mvector.cc mmatrix.cc msubspace.cc svector.cc smatrix.cc smatrix_elim.cc xsplit.cc conic.cc legendre.cc quadratic.cc illl.cc hilbert.cc timer.cc cubic.cc polys.cc realroots.cc parifact.cc p2points.cc gf.cc

QCURVES_DOTCCS = curve.cc curvedata.cc curvered.cc points.cc cperiods.cc isogs.cc heights.cc mwprocs.cc lambda.cc sifter.cc sieve_search.cc htconst.cc egr.cc saturate.cc divpol.cc pointsmod.cc curvemod.cc ffmod.cc tlss.cc elog.cc getcurve.cc

QRANK_DOTCCS = mequiv.cc  mrank1.cc mlocsol.cc mglobsol.cc mquartic.cc mrank2.cc qc.cc sqfdiv.cc version.cc minim.cc reduce.cc transform.cc desc2.cc bitspace.cc GetOpt.cc twoadic.cc descent.cc

G0N_DOTCCS = newforms.cc symb.cc homspace.cc cusp.cc oldforms.cc fixc6.cc periods.cc moddata.cc pcprocs.cc nfd.cc curvesort.cc

DOTCCS = $(PROCS_DOTCCS) $(QCURVES_DOTCCS) $(QRANK_DOTCCS) $(G0N_DOTCCS)

# NB Do NOT include the following files above, since they are not
# compiled by themselves but only included in other files which are
# compiled.  But they do need to be in the distribution made by "make
# dist" so we must list them here:

EXTRA_DIST = vec.cc mat.cc sub.cc smat.cc svec.cc smat_elim.cc random.cc

libec_la_SOURCES = $(DOTCCS)
libec_la_LDFLAGS = -no-undefined $(AM_LDFLAGS)
nobase_include_HEADERS = $(DOTHS)
ACLOCAL_AMFLAGS = -I m4
