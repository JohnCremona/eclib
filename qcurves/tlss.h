// tlss.h: definition of class TLSS for sieving E(Q)/pE(Q) at one prime q
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
 

// NB: TLSS = Tate--Lichtenbaum--Samir-Siksek: we use a simple
// discrete log a la Siksek when the p-torsion in E(F_q) is cyclic,
// else use the Tate-Lichtenbaum pairing

// allow for multiple includes
#ifndef _TLSS_
#define _TLSS_


class TLSS {
private: 
  int p;                // the prime to saturate at.
  int rank;             // p-rank of E mod q (0,1,2)
  bigint q;             // the modulus
  bigint q1p;           //  = (q-1)/p;

  galois_field Fq;          // F_q
  vector<gf_element> mu_p;  // all p'th roots mod q
  curvemodqbasis Emodq;     // E over F_q (including its structure)
  vector<pointmodq> Pi;    // basis for p-torsion of E(F_q) (length = rank = 0,1,2)
  vector<ffmodq> TLpolys;   // in the function field Fq(E)
  int verbose;
  void init_tlpolys(void);

public:
  TLSS(void) :Emodq() {;}
  void assign(const curvemodqbasis& E) {Emodq=E; Fq=get_field(Emodq); q=Fq.characteristic(); }
  void init(int pp, int verb=0);
  void init(int pp, const vector<bigint>& pdivpol, int verb=0);
  ~TLSS() {; }
  
  // apply map to P, result in (ntp*)[0..p-1]:
  vector<int> map1point(const Point& P) const;
  // apply map to all P in Plist, result is a (ntp*#Plist) matrix:
  mat map_points(const vector<Point>& Plist) const;
  // give the current p-rank
  int get_rank() const {return rank;}

};


#endif
