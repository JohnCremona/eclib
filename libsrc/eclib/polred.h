// FILE POLRED.H: declaration of functions for reducing ZZX polynomials (polredabs) via libpari

#if     !defined(_POLRED_H)
#define _POLRED_H      1       //flags that this file has been included

#include <assert.h>
#include "templates.h"
#include "int.h"
#include "polys.h"
#include "pari_init.h"
#undef recip // pariold.h #defines recip = serreverse

// convert a t_POL (with Z or Q coefficients) to a ZZX with common denominator
ZZX t_POL_to_ZZX(GEN P, ZZ& d);

// conversion from ZZX to t_POL (with Z coefficients)
GEN ZZX_to_t_POL(const ZZX& f);

// polredabs of an *irreducible* polynomial in Z[X]
// (1) return monic integral g defining the same field as f
ZZX polredabs(const ZZX& f);
// (2) also sets h such that a=h(b)/d (so f(h(b)/d)=0) where f(a)=g(b)=0
ZZX polredabs(const ZZX& f, ZZX& h, ZZ& d);
#endif
