// sqfdiv.h : declaration of class sqfdiv for managing square-free divisors
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
 
#ifndef _ECLIB_SQFDIV_H
#define _ECLIB_SQFDIV_H 1       //flags that this file has been included

#include "marith.h"

class sqfdiv {
  vector<ZZ>* primebase;  // includes all relevant primes
  ZZ d;                // product of current subset
  long np;                  // number in current subset
  int positive;             // flag for sign
  long factor;              // counts log_2 of saving index since initialisation
  vector<ZZ> subgp;       // subgp factored out (complete)
  vector<ZZ> gens;        // generators of latter
  long nsub, maxnsub, ngens, maxngens;       // current, max number in subgp
  vector<long> pivs;
public:
  sqfdiv(const ZZ& dd, int posd, vector<ZZ>* plist);
  void usediv(const ZZ& e);  
  vector<ZZ> getdivs() const;
  vector<ZZ> getsupp(int bothsigns=0) const;
  vector<ZZ> getsubgp() {return vector<ZZ>(subgp.begin(),subgp.begin()+nsub);}
  long getfactor() {return factor;}
  void display();
};

ZZ sqfred(const ZZ& a, const vector<ZZ>& plist);

inline ZZ sqfred(const ZZ& a) { return sqfred(a,pdivs(a));}

ZZ sqfmul(const ZZ& a, const ZZ& b);

ZZ makenum(const vector<ZZ>& supp, long mask);

long makeindex(const vector<ZZ>& supp, const ZZ& n, ZZ& n0);

// support(n) is like pdivs(n) but includes -1 always 
// (except for n=0, but it should never be called with 0)

vector<ZZ> support(const ZZ& n);

#endif
