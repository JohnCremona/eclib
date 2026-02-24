// FILE POLRED.CC: implementation of functions for reducing ZZX polynomials (polredabs) via libpari

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
using PARI::polredabs;
using PARI::polredabs0;
using PARI::nf_ORIG;
using PARI::lift;

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
  ZZ num = PARI_to_NTL(PARI::numerator(cont,gen_1)); // numerator of content
  d = PARI_to_NTL(PARI::denominator(cont,gen_1));    // denominator of content
  ZZX f;
  for (int i=0; i<=deg; i++)
    SetCoeff(f, i, num * PARI_to_NTL(gel(P1, i+2)));
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

// polredabs of an *irreducible* polynomial in Z[X]
// (1) return monic integral g defining the same field as f
ZZX polredabs(const ZZX& f)
{
  if (!IsIrreducible(f))
    {
      cerr << "polredabs() called with f = " << str(f) << " which is reducible" << endl;
      return f;
    }
  pari_sp av = avma;
  GEN G = polredabs(ZZX_to_t_POL(f));
#ifdef DEBUG_POLY
  pari_printf("polredabs(f) returns %Ps\n", G);
#endif
  ZZ d;
  ZZX g = t_POL_to_ZZX(G, d); // d will be 1
  assert (d==1);
  avma = av;
  return g;
}

// (2) also sets h such that a=h(b) (so f(h(b))=0) where f(a)=g(b)=0

ZZX polredabs(const ZZX& f, ZZX& h, ZZ& d)
{
  if (!IsIrreducible(f))
    {
      cerr << "polredabs() called with f = " << str(f) << " which is reducible" << endl;
      return f;
    }
  pari_sp av = avma;
  GEN G_H = polredabs0(ZZX_to_t_POL(f), nf_ORIG);
#ifdef DEBUG_POLY
  pari_printf("polredabs0(f,nf_ORIG) returns [G, H] = %Ps\n", G_H);
#endif
  GEN G = gel(G_H,1);
#ifdef DEBUG_POLY
  pari_printf("  G = %Ps\n", G);
#endif
  GEN H = lift(gel(G_H,2));
#ifdef DEBUG_POLY
  pari_printf("  H = %Ps\n", H);
#endif
  ZZX g = t_POL_to_ZZX(G, d); // this d=1
  assert (d==1);
  // NB H need not have integer coefficients!
  h = t_POL_to_ZZX(H, d); // thid d may be >1
  avma = av;
  return g;
}

