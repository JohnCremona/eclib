// mwprocs.h: definition of class mw for Mordell-Weil basis
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
 
#include "interface.h"
#include "compproc.h"

#include "matrix.h"
#include "subspace.h"

#include "points.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "ffmod.h"
#include "divpol.h"
#include "tlss.h"
#include "elog.h"
#include "saturate.h"
#include "sieve_search.h"

// allow for multiple includes
#ifndef _MWPROCS_
#define _MWPROCS_

void ratapprox(bigfloat x, bigint& a, bigint& b);
void ratapprox(bigfloat x, long& a, long& b);

// global to this file -- size of height_pairing matrix, equal
// to maximum rank of any curve likely to study

const int MAXRANK = 30;
const int MAXSATPRIME = 20;  // default saturation limit

class mw : public point_processor {
private: 
  Curvedata *E;
  vector<Point> basis;
  int rank, maxrank;
  bigfloat *height_pairs;
  bigfloat reg;
  int verbose, process_points;
  bigfloat& mat_entry(int i, int j);
  bigint a1,a2,a3,a4,a6;
  int iso;
  saturator satsieve;
public:
  mw(Curvedata*, int verb=0, int pp=1, int maxr=999);
  ~mw();

 // processing of new points, with saturation at primes up to sat 
 // (default MAXSATPRIME,  none if sat==0)
  int process(const bigint& x, const bigint& y, const bigint& z);
  int process(const bigint& x, const bigint& y, const bigint& z, int sat);
 // as returned by ms's sieve; the point is (x/z,y/z^(3/2)) and z is square
  int process(const Point& P, int sat=MAXSATPRIME);
  int process(const vector<Point>& Plist, int sat=MAXSATPRIME);
  // saturate the current basis:
  int saturate(bigint& index, vector<long>& unsat, long sat_bd=-1, int odd_primes_only=0);
  void search(bigfloat h_lim, int moduli_option=0, int verb=0);
  void search_range(bigfloat xmin, bigfloat xmax, bigfloat h_lim, 
		    int moduli_option=2, int verb=0);
  bigfloat regulator(void) {return reg;}
  vector<Point> getbasis() {vector<Point> b(basis.begin(),basis.begin()+rank); return b;}
  int getrank() {return rank;}
};

inline bigfloat& mw::mat_entry(int i, int j)
{
  return *(height_pairs + (i*MAXRANK) + j);
}

class sieve {
private:
  Curvedata *E;
  bigint a1,a2,a3,a4,a6;
  bigint d1,d2,d3,d4,d6,c2,c3,c4,c6;
  long a,c;
  mw * mwbasis;
  int verbose, posdisc, firstnl;
  bigfloat xmin,x1,x2,x3;
  int num_aux;
  long* auxs;
  int** xgood_mod_aux;
  int** x1good_mod_aux;
  int** squares;
  long* amod;
  long *modhits;
  long npoints, ascore, amodc, alim, clim0, clim1, clim2, clim;
  int* cflag;  int use_cflag;
  void a_search(const long& amin, const long& amax);
  void a_simple_search(const long& amin, const long& amax);
public:
  sieve(void) {;}
  sieve(Curvedata * EE, mw* mwb, int moduli_option, int verb=0); 
  ~sieve();
  void search(bigfloat h_lim);
  void search_range(bigfloat xmin, bigfloat xmax, bigfloat h_lim);
  void stats();   // report sieving statistics
};

int order_real_roots(vector<double>& bnd, vector<bigcomplex> roots);
//checks (and returns) how many roots are actually real, and puts those in 
//bnd, in increasing order, by calling set_the_bound
int set_the_bounds(vector<double>& bnd, bigfloat x0, bigfloat x1, bigfloat x2);
//This transforms (if possible) x0, x1 and x1 into double;  the search 
//should be made on [x0,x1]U[x2,infty] so if x1 or x2 overflows, the search 
//is on [x0,infty].  The function returns 3 in the first case, 1 in the second.
//If x0 overflows, it returns 0.  A warning is printed out.

#endif
