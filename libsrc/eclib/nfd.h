// File NFD.H: class for newforms of any dimension
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
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

#ifndef _ECLIB_NFD_H
#define _ECLIB_NFD_H      1

#include "homspace.h"

class nfd {
private:
  mat_m tp0,  // a defining matrix
    V, W, Winv, WinvV, Winv_scaled;
  mat projcoord;
  scalar coord_fac;
  bigint Wdetnum, Wdetdenom;
  vector<bigint> minpol;  // min poly of alpha, field generator
public:
  vector<bigint> Hscales;
  vector<bigint> Sscales;
  subspace_m S;   // the basis
  long N;         // the level
  homspace* H1; // the ambient modular symbol space
  mat K;    // basis of ker(delta) on h1
  vec Kcol; // one column of K
  int rk;   // #rows of K
  int dimH, dimS; // dimensions of H1 and S
  bigint dH, dS, dHS; // denominators of H1, S (relative) and S (absolute)
  nfd(void) {;}
  nfd(homspace* h1, int one_p, int w_split, int mult_one, int verbose=0);
  void display(void) const;
  mat_m oldheckeop(long p);
  mat_m heckeop(long p);
  vec_m ap(long p);
};

#endif
