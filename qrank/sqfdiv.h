// sqfdiv.h : declaration of class sqfdiv for managing square-free divisors
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 
class sqfdiv {
  vector<bigint>* primebase;  // includes all relevant primes
  bigint d;                // product of current subset
  long np;                  // number in current subset
  int positive;             // flag for sign
  long factor;              // counts log_2 of saving index since initialisation
  vector<bigint> subgp;       // subgp factored out (complete)
  vector<bigint> gens;        // generators of latter
  long nsub, maxnsub, ngens, maxngens;       // current, max number in subgp
  vector<long> pivs;
public:
  sqfdiv(const bigint& dd, int posd, vector<bigint>* plist);
  void usediv(const bigint& e);  
  vector<bigint> getdivs() const;
  vector<bigint> getsupp(int bothsigns=0) const;
  vector<bigint> getsubgp() {return vector<bigint>(subgp.begin(),subgp.begin()+nsub);}
  long getfactor() {return factor;}
  void display();
};

bigint sqfred(const bigint& a, const vector<bigint>& plist);

inline bigint sqfred(const bigint& a) { return sqfred(a,pdivs(a));}

bigint sqfmul(const bigint& a, const bigint& b);

bigint makenum(const vector<bigint>& supp, long mask);

long makeindex(const vector<bigint>& supp, const bigint& n, bigint& n0);

// support(n) is like pdivs(n) but includes -1 always 
// (except for n=0, but it should never be called with 0)

vector<bigint> support(const bigint& n);
