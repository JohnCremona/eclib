// FILE POLRED.H: declaration of functions for reducing ZZX polynomials (polredabs) via libpari
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

#if     !defined(_POLRED_H)
#define _POLRED_H      1       //flags that this file has been included

#include <assert.h>
#include "templates.h"
#include "polys.h"
#include "qvecmat.h"
#include "pari_init.h"
#undef recip // pariold.h #defines recip = serreverse

// convert a t_POL (with Z or Q coefficients) to a ZZX with common denominator
ZZX t_POL_to_ZZX(GEN P, ZZ& d);

// conversion from ZZX to t_POL (with Z coefficients)
GEN ZZX_to_t_POL(const ZZX& f);

// polredabs of an *irreducible* polynomial in Z[X]

// These use pari's polredabs0() function with flag=0 or flag=nf_ORIG
// respectively when canonical=1: "Expensive but provide canonical
// representatives", so slow if degree large (say 20 or more). If
// canonical=0 (default), use polredbest() instead: "Fast, generally
// best. Possibly smaller discriminant (!)"

// (1) return monic integral g defining the same field as f
ZZX polred(const ZZX& f, int canonical=0);
// (2) also sets h such that a=h(b)/d (so f(h(b)/d)=0) where f(a)=g(b)=0
ZZX polred(const ZZX& f, ZZX& h, ZZ& d, int canonical=0);

// interface to nfinit: given f (monic irreducible) defining a field
// K=Q(t), returns d = [O_K:Z[t]], zbasis of OK (as rational
// coordinate vectors w.r.t. t-power basis), integral base-change
// matrix (integral basis coordinates of powers of t)

void nfinit(const ZZX& f, ZZ& ind, vector<Qvec>& zbasis, mat_m& bcm);

#endif
