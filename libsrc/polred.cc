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
  GEN G = ZZX_to_t_POL(f), A;
  int nsteps = 0;
  if (canonical)
    {
      GEN G_H = polredabs0(G, nf_ORIG);
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
      pari_printf("applying polredbest to %Ps\n", G);
#endif
      GEN G_H = polredbest(G, nf_ORIG);
      nsteps++;
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
  // NB H need not have integer coefficients!
  h = t_POL_to_ZZX(A, d); // this d may be >1
  avma = av;
  return g;
}

