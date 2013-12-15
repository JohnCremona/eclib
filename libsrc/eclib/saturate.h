// saturate.h: declaration of class saturator for sieving E(Q)/pE(Q)
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
 

// This is used for proving that points are p-saturated

// allow for multiple includes
#ifndef _SATURATE_
#define _SATURATE_

#ifndef MAX_REPEAT
#define MAX_REPEAT 10
#endif

// automatic saturation will only saturate at primes less than this,
// unless overridden by the sat_bd parameter:
const long SAT_MAX_PRIME = 100000;

class saturator {
private: 
  Curvedata *E;         // the curve 
  vector<Point> Plist;  // the points
  vector<Point> Plistp;  // the p-cotorsion
  vector<Point> Plistx;  // the points plus p-cotorsion
  vector<Point> AllTorsion; // all torsion on E
  vector<bigint> pdivpol; // coefficients of p-division poly (not always used)
  int rank;             // = #Plistx
  bigint disc;          // discriminant of E
  int p;                // current prime to saturate at
  int log_index;         // current points have index p^log_index in original
  primevar qvar;          // loops over possible sieving primes q
  vector<curvemodqbasis> Eqlist;   // E mod q for q=3,5,... (good reduction)
  vector<curvemodqbasis>::iterator Eqptr;
  int newq;               // =1 iff we are using q not yet cahced

  mat TLimage; 
  int TLrank, stuck_counter, verbose, use_div_pols;

  // apply TL map (mod current q) to P, result in (ntp*)[0..p-1]:
  vector<int> TLmap1(const Point& P);
  // apply TL map (mod current q) to Plistx, result is a (ntp*rank) matrix:
  mat TLmap();
  //

public:
  saturator(Curvedata* EE, int verb=0)
    :E(EE), verbose(verb)
    {
#ifdef NTL_INTS
      use_div_pols=0;//1;
#else
      use_div_pols=0;
#endif
      disc = getdiscr(*E);
      AllTorsion = torsion_points(*EE);
    }
  ~saturator() {; }
  
  // initialize point list
  void set_points(const vector<Point>& PP) {Plist = PP;}
  // find next usable q and use it
  void nextq();
  // keep on using q until saturated or stuck:
  int test_saturation(int pp, int ms=MAX_REPEAT);
  // try harder with same p (no initialization, but same point set)
  int test_saturation_extra(int pp, int ms=MAX_REPEAT);
  // enlarge basis if dim(kernel)>0:
  int enlarge();
  // repeat testing saturation and enlarging until done:
  // returns log_p of index
  int do_saturation(int pp, int maxntries=10);
  // As above but for all primes up to p, returns index
  int do_saturation_upto(int maxp, int maxntries=10);
  // As above but for all primes in plist, returns success flag, sets
  // index and unsat = list of primes in plist at which
  // saturation failed
  int do_saturation(vector<int> plist, 
		    bigint& index, vector<int>& unsat, int maxntries=10);
  int do_saturation(vector<long> plist, 
		    bigint& index, vector<long>& unsat, int maxntries=10);
  // auto-saturate, after finding an upper bound on saturation index
  int saturate(vector<long>& unsat, bigint& index, long sat_bd=-1, int egr=1, int maxntries=10, int odd_primes_only=0);

  // replace the generating points & reset matrices and ranks
  // (then can use test_saturation_extra())
  // (for re-saturating at a prime p already used)
  void reset_points(const vector<Point>& PP);

  // get p-rank of input + p-cotorsion
  long getprank() const {return rank;};
  // get rank of current image
  long getTLrank() const {return TLrank;};
  // get current q:
  long get_q() const {return (long)qvar;};
  // get current generators (not including p-cotorsion)
  vector<Point> getgens() const {return Plist;};
  // get current generators (including p-cotorsion)
  vector<Point> getxgens() const {return Plistx;};
  // get # steps since last rank increase
  long stuckfor() const {return stuck_counter;};
  // get a nonzero kernel vector (if any)
  vec kernel_vector();
  // test if saturated: 
  int is_saturated() {return rank==TLrank;}
};

// auxiliary functions, not part of the saturator class:

// Use this where Torsion already computed:
vector<Point> pCoTorsion(const vector<Point>& AllTorsion, int p);

// Use this for a one-off computation where Torsion not yet computed:
inline vector<Point> pCoTorsion(Curvedata* EE, int p)
{
  return pCoTorsion(torsion_points(*EE),p);
}

// Saturate the given points automatically (using egr strategy if set)
// Returns success flag and index, and (if success==0) a list of
// primes at which saturation failed.  The array points will be
// changed (if necessary) into a basis for the saturation.
// If sat_bd>-1 it is used to truncate the upper bound

int saturate_points(Curvedata& C, vector<Point>& points, 
		    bigint& index, vector<long>& unsat,  
		    long sat_bd=-1, int egr=1, int verbose=0);


// Bound for the index of saturation for the given set of points If
// egr is set it determines the egr subgroup of the group the points
// generate and only searches for points with egr, This might be faster
// in some cases...
bigint index_bound(Curvedata* C,  vector<Point>& points, 
		   int egr=1, int verbose=0);

// Tamagawa primes: primes dividing any Tamagawa number
vector<long> tamagawa_primes(const Curvedata& C);

#endif
