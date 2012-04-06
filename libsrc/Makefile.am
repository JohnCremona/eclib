## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libjc.la

PROCS_DOTHS = interface.h templates.h arith.h xmod.h marith.h gpslave.h compproc.h vec.h vector.h mat.h matrix.h sub.h subspace.h rat.h bigrat.h kbessel.h svec.h svector.h smat.h smatrix.h smat_elim.h smatrix_elim.h mvector.h mmatrix.h msubspace.h method.h splitbase.h xsplit.h conic.h legendre.h quadratic.h unimod.h illl.h hilbert.h timer.h cubic.h gf.h polys.h realroots.h parifact.h p2points.h

QCURVES_DOTHS = curve.h points.h cperiods.h isogs.h reader.h mwprocs.h lambda.h sifter.h sieve_search.h htconst.h egr.h saturate.h divpol.h pointsmod.h curvemod.h ffmod.h tlss.h elog.h getcurve.h

QRANK_DOTHS = mrank2.h mequiv.h mlocsol.h mglobsol.h msoluble.h mquartic.h mrank1.h qc.h sqfdiv.h version.h minim.h reduce.h transform.h desc2.h bitspace.h options.h GetOpt.h twoadic.h descent.h

G0N_DOTHS = periods.h cusp.h homspace.h newforms.h moddata.h oldforms.h symb.h fixc6.h pcprocs.h nfd.h

DOTHS = $(PROCS_DOTHS) $(QCURVES_DOTHS) $(QRANK_DOTHS) $(G0N_DOTHS)

PROCS_DOTCCS = interface.cc unimod.cc arith.cc marith.cc gpslave.cc compproc.cc vector.cc matrix.cc subspace.cc kbessel.cc mvector.cc mmatrix.cc msubspace.cc svector.cc smatrix.cc smatrix_elim.cc xsplit.cc conic.cc legendre.cc quadratic.cc illl.cc hilbert.cc timer.cc cubic.cc polys.cc realroots.cc parifact.cc p2points.cc gf.cc

QCURVES_DOTCCS = curve.cc curvedata.cc curvered.cc points.cc cperiods.cc isogs.cc heights.cc mwprocs.cc lambda.cc sifter.cc sieve_search.cc htconst.cc egr.cc saturate.cc divpol.cc pointsmod.cc curvemod.cc ffmod.cc tlss.cc elog.cc getcurve.cc

QRANK_DOTCCS = mequiv.cc  mrank1.cc mlocsol.cc mglobsol.cc mquartic.cc mrank2.cc qc.cc sqfdiv.cc version.cc minim.cc reduce.cc transform.cc desc2.cc bitspace.cc GetOpt.cc twoadic.cc descent.cc

G0N_DOTCCS = newforms.cc symb.cc homspace.cc cusp.cc oldforms.cc fixc6.cc periods.cc moddata.cc pcprocs.cc nfd.cc

DOTCCS = $(PROCS_DOTCCS) $(QCURVES_DOTCCS) $(QRANK_DOTCCS) $(G0N_DOTCCS)

# NB Do NOT include the following files above, since they are not
# compiled by themselves but only included in other files which are
# compiled.  But they do need to be in the distribution made by "make
# dist" so we must list them here:

EXTRA_DIST = vec.cc mat.cc sub.cc smat.cc svec.cc smat_elim.cc random.cc

libjc_la_SOURCES = $(DOTCCS)
include_HEADERS = $(DOTHS)
ACLOCAL_AMFLAGS = -I m4
