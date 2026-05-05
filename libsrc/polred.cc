// FILE POLRED.CC: implementation of functions for reducing ZZX polynomials (polredabs) via libpari
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
//
// This file is part of the eclib package.
//
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
//
//////////////////////////////////////////////////////////////////////////

#include "eclib/polred.h"
#include "eclib/pari_init.h"
#include "eclib/convert.h"
#include <assert.h>

using PARI::degpol;
using PARI::t_POL;
using PARI::content0;
using PARI::gdiv;
using PARI::gen_1;
using PARI::RgX_renormalize_lg;
using PARI::polredbest;
using PARI::polredabs0;
using PARI::nf_ORIG;
using PARI::lift;
using PARI::ZX_equal;
using PARI::poleval;
using PARI::nfisincl0;
using PARI::nfinit0;

//#define DEBUG_POLY

// convert a t_POL (with Z or Q coefficients) to a ZZX with common denominator
ZZX t_POL_to_ZZX(GEN P, ZZ& d)
{
#ifdef DEBUG_POLY
  pari_printf(" Converting t_POL %Ps to ZZX\n", P);
#endif
  int deg = degpol(P);
  GEN cont = content0(P, gen_1);
#ifdef DEBUG_POLY
  pari_printf(" Content = %Ps\n", cont);
#endif
  GEN P1 = gdiv(P,cont); // should be integral and primitive
#ifdef DEBUG_POLY
  pari_printf(" P1 = %Ps and content  = %Ps\n", P1, cont);
#endif
  ZZ n = PARI_to_NTL(PARI::numerator(cont,gen_1)); // numerator of content
  d = PARI_to_NTL(PARI::denominator(cont,gen_1));  // denominator of content
  ZZX f;
  for (int i=0; i<=deg; i++)
    SetCoeff(f, i, n * PARI_to_NTL(gel(P1, i+2)));
#ifdef DEBUG_POLY
  if (d==1)
    cout << " Result is " << str(f) << endl;
  else
    cout << " Result is (" << str(f) << ") / " << d << endl;
#endif
  return f;
}

GEN ZZX_to_t_POL(const ZZX& f)
{
#ifdef DEBUG_POLY
  cout << " Converting  ZZX " << str(f) << " to t_POL" << endl;
#endif
  int d = deg(f);
  GEN P = cgetg(d+3, t_POL);
  P[1] = evalvarn(0); // set variable to #0, i.e. 'x'
  for (int i=0; i<=d; i++)
    gel(P,i+2) = NTL_to_PARI(coeff(f,i));
  P = RgX_renormalize_lg(P,d+3);
#ifdef DEBUG_POLY
  pari_printf(" Result is %Ps\n", P);
#endif
  return P;
}

//#define DEBUG_POLRED

// polredabs (if canonical) or polredbest of an *irreducible*
// polynomial in Z[X]

// (1) return monic integral g defining the same field as f

ZZX polred(const ZZX& f, int canonical)
{
  if (!IsIrreducible(f))
    {
      cerr << "polred() called with f = " << str(f) << " which is reducible" << endl;
      return f;
    }
  pari_sp av = avma;
  GEN G = ZZX_to_t_POL(f);
  int nsteps = 0;
  if (canonical)
    G = polredabs0(G, 0);
  else
    {
#ifdef DEBUG_POLRED
      pari_printf("applying polredbest to %Ps\n", G);
#endif
      GEN G1 = polredbest(G, 0);
      nsteps++;
#ifdef DEBUG_POLRED
      pari_printf("polredbest step %d gives %Ps\n", nsteps, G1);
#endif
      while (!ZX_equal(G,G1))
        {
          G = G1;
          G1 = polredbest(G, 0);
          nsteps++;
#ifdef DEBUG_POLRED
          pari_printf("polredbest step %d gives %Ps\n", nsteps, G1);
#endif
        }
    }
#ifdef DEBUG_POLRED
  if (canonical)
    pari_printf("polredabs0(f, 0) returns %Ps\n", G);
  else
    pari_printf("repeating polredbest stabilises after %d steps at %Ps\n", nsteps, G);
#endif
  ZZ d;
  ZZX g = t_POL_to_ZZX(G, d); // d will be 1
  assert (d==1);
  avma = av;
  return g;
}

// (2) also sets h such that a=h(b) (so f(h(b))=0) where f(a)=g(b)=0

ZZX polred(const ZZX& f, ZZX& h, ZZ& d, int canonical)
{
  if (!IsIrreducible(f))
    {
      cerr << "polred() called with f = " << str(f) << " which is reducible" << endl;
      return f;
    }
  pari_sp av = avma;
  GEN F = ZZX_to_t_POL(f), A, G;
  int nsteps = 0;
  if (canonical)
    {
      GEN G_H = polredabs0(F, nf_ORIG);
#ifdef DEBUG_POLRED
      pari_printf("polredabs0(f,nf_ORIG) returns [G, H] = %Ps\n", G_H);
#endif
      G = gel(G_H,1);
      A = lift(gel(G_H,2));
#ifdef DEBUG_POLRED
      pari_printf("  G = %Ps\n  A = %Ps\n", G, A);
#endif
    }
  else
    {
#ifdef DEBUG_POLRED
      pari_printf("applying polredbest to %Ps\n", F);
#endif
      GEN G_H = polredbest(F, nf_ORIG);
      nsteps++;
      G = F;
      GEN G1 = gel(G_H,1);
      A = lift(gel(G_H,2));
#ifdef DEBUG_POLRED
      pari_printf("After 1 step of polredbest,   G1 = %Ps,  A = %Ps\n", G1, A);
#endif
      while (!ZX_equal(G,G1))
        {
          G = G1;
          G_H = polredbest(G, nf_ORIG);
          nsteps++;
          G1 = gel(G_H,1);
          if (!ZX_equal(G,G1))
            A = lift(poleval(A, gel(G_H,2)));
#ifdef DEBUG_POLRED
          pari_printf("After %d steps of polredbest,  G1 = %Ps,  A = %Ps\n", nsteps, G1,  A);
#endif
        }
    }

#ifdef DEBUG_POLRED
  if (canonical)
    pari_printf("polredabs0(f,nf_ORIG) returns G = %Ps, A = %Ps\n", G, A);
  else
    pari_printf("polredbest(f,nf_ORIG) returns G = %Ps, A = %Ps\n", G, A);
#endif
  ZZX g = t_POL_to_ZZX(G, d); // this d=1
  assert (d==1);

  // Now H may not be unique (if the field has nontrivial
  // automorphisms).  To be deterministic, we find all possibles, sort
  // and return the first one.

  GEN Hlist = nfisincl0(F, G, 0);
  int nH = lg(Hlist)-1;
#ifdef DEBUG_POLRED
  cout << nH << " possible images: ";
  pari_printf("%Ps\n", Hlist);
#endif
  if (nH==1) // only one choice
    {
      // NB H need not have integer coefficients!
      h = t_POL_to_ZZX(A, d); // this d may be >1
      avma = av;
      return g;
    }
  vector<ZZX> hlist(nH);
  vector<ZZ> dlist(nH);
  for (int i=0; i<nH; i++)
    hlist[i] = t_POL_to_ZZX(gel(Hlist,i+1), dlist[i]);

#ifdef DEBUG_POLRED
  cout << "After conversion:\nh-polys: " << hlist << "\ndenoms: " << dlist << endl;
#endif

  // Now we sort the h's and pick the first
  auto hmin = std::min_element(hlist.begin(), hlist.end(), poly_cmp);
  h = *hmin;
  d = dlist[std::distance(hlist.begin(), hmin)];
#ifdef DEBUG_POLRED
  cout << "Smallest is h =  " << h << ", d = " << d << endl;
#endif

  avma = av;
  return g;
}

//#define DEBUG_NF_INIT

void nfinit(const ZZX& f, ZZ& ind, vector<Qvec>& zbasis, mat_m& bcm)
{
  pari_sp av = avma;
  int d = deg(f);
#ifdef DEBUG_NF_INIT
  cout << "In nfinit(" << str(f) << ")" << endl;
#endif
  GEN F = ZZX_to_t_POL(f);
  GEN nf = nfinit0(F, 0, DEFAULTPREC);
#ifdef DEBUG_NF_INIT
  pari_printf(" - PARI::nfinit0() returns%Ps\n", nf);
#endif

  ind = PARI_to_NTL(gel(nf, 4));
#ifdef DEBUG_NF_INIT
  cout << " - index of equation order in maximal order = ind = " << ind << endl;
#endif

  GEN zb = gel(nf, 7); // list of d integer polys mod F
  // pari_printf(" - zb = %Ps\n", zb);
  vec_m co(d);
  zbasis.clear();
  for (int i=0; i<d; i++)
    {
      GEN zbi = lift(gel(zb, i+1)); // integer poly
      // pari_printf(" - zbi = %Ps\n", zbi);
      int di = degpol(zbi);
      for (int j=1; j<=d; j++)
        co[j] = (j<=di+1? PARI_to_NTL(gel(zbi, j+1)) : ZZ(0)); // coeff of x^(j-1)
      Qvec bi(co,ind);
      zbasis.push_back(bi);
#ifdef DEBUG_NF_INIT
      cout << i << ": " << bi << endl;
#endif
    }
#ifdef DEBUG_NF_INIT
  cout << " - coefficients of integral basis w.r.t. power basis:\n";
  for (auto bi: zbasis)
    cout << bi << endl;
#endif

  GEN zbi = gel(nf, 8); // dxd matrix of integers
  bcm.init(d,d);
  for (int i=1; i<=d; i++)
    for (int j=1; j<=d; j++)
      bcm.set(i,j, PARI_to_NTL(gcoeff(zbi, i,j)));
#ifdef DEBUG_NF_INIT
  cout << " - coefficients of power basis w.r.t. integral basis:\n";
  cout << bcm << endl;
#endif

  avma = av;
  return;
}
