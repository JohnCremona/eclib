// mglobsol.h: declaration of class quartic_sieve and functions for quartic solubility testing
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
 
#include "mquartic.h"
#include "sieve_search.h"

// Function for naive search, no sieving:

int ratpoint(const quartic& g, const bigint& min, const bigint& max, bigint& xx, bigint&yy, bigint& zz);

// class for fancier sieve-assisted search

class quartic_sieve : public point_processor {
private:
  quartic *g;
  bigint a,b,c,d,e,roota,roote;
  bigint pu,pv,pw;   // coords of point found
  int verbose, easy, use_stoll;
  long ulim;
  int num_aux;
  long* auxs;
  int** xgood_mod_aux;
  int** squares;
  long* umod;
  long nwprimes; long* wprimes;  long* uprimes;
  long npoints, maxnpoints;
  int process(const bigint& x, const bigint& y, const bigint& z)
  {pu=x; pv=y; pw=z; npoints++; 
  //cout<<"[x,y,z]=["<<x<<","<<y<<","<<z<<"]\n"; 
  return (npoints>=maxnpoints);
  }
 // (x,y,z) as returned by ms's sieve; the point is (x/z,y/z^2)
public:
  quartic_sieve(void) {;}
  quartic_sieve(quartic * gg, int moduli_option=2, int verb=0); 
  ~quartic_sieve();
  long search(double h_lim, long maxnpts=1, int posxonly=0);
  long stoll_search(double h_lim, int posxonly=0);
  long search_range(int lower, bigfloat lower_bound, 
		   int upper, bigfloat upper_bound, int posxonly=0);
  void getpoint(bigint& x, bigint& y, bigint& z) {x=pu;y=pv;z=pw;}
};

