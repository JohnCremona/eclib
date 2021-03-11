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
#ifndef _ECLIB_SATURATE_H
#define _ECLIB_SATURATE_H

#include <eclib/matrix.h>
#include <eclib/ffmod.h>


#ifndef MAX_REPEAT
#define MAX_REPEAT 10
#endif

// automatic saturation will only saturate at primes less than this,
// unless overridden by the sat_bd parameter:
const long SAT_MAX_PRIME = 100000;


// Bound for the index of saturation for the given set of points If
// egr is set it determines the egr subgroup of the group the points
// generate and only searches for points with egr, This might be
// faster in some cases.

// The points parameter is not const, since the heights of the points
// may be computed and set.

// We can recover the curve from the points, unless the list is empty
// in which case the bound is 1.

long index_bound(vector<Point>& points, int egr=1, int verbose=0);

// Tamagawa primes: primes dividing any Tamagawa number
vector<long> tamagawa_primes(const Curvedata& C);


class saturator {
private: 
  Curvedata *E;         // the curve 
  vector<Point> Plist;  // the points
  vector<Point> Plistp;  // the p-cotorsion
  vector<Point> Plistx;  // the points plus p-cotorsion
  vector<Point> AllTorsion; // all torsion on E
  ZPoly pdivpol;            // p-division poly (not always used)
  long the_index_bound;     // set initially, but may get reduced
  vector<long> tam_primes;  // primes dividing any Tamagawa index
  int rank;             // = #Plistx
  bigint disc;          // discriminant of E
  int p;                // current prime to saturate at
  int log_index;         // current points have index p^log_index in original
  primevar qvar;          // loops over possible sieving primes q
  map<bigint, curvemodqbasis> Emodq;
  map<bigint, bigint> Emodq_order;
  int newq;                // =1 iff we are using q not yet cached
  map<bigint,int> q_tally;   // key=q, value = count of number of times q used
  bigint maxq;               // largest q used
  int maxp;                  // p for which largest q used

  mat_l TLimage;
  int TLrank, stuck_counter, verbose, use_div_pols;

  // apply TL map (mod current q) to P, result in (ntp*)[0..p-1]:
  vector<int> TLmap1(const Point& P);
  // apply TL map (mod current q) to Plistx, result is a (ntp*rank) matrix:
  mat_l TLmap();
  //

public:
  saturator(Curvedata* EE, int verb=0)
    :E(EE), verbose(verb)
    {
      use_div_pols=0;
      disc = getdiscr(*E);
      AllTorsion = torsion_points(*EE);
      tam_primes = tamagawa_primes(*E);
      maxq = 0;
      the_index_bound = 0; // means not set
    }
  ~saturator() {; }

  // initialize point list
  void set_points(const vector<Point>& PP) {
    Plist = PP;
    the_index_bound = 0;
  }

  // initialize index bound
  void set_index_bound(int egr=1);

  // return current index bound (compute if necessary)
  long get_index_bound(int egr=1);

  // test whether p is greater than the saturation index and not a Tamagawa prime
  int trivially_saturated(long p);

  // find next usable q and use it
  void nextq();

  // test p-saturation by using q until saturated (return 1) or stuck
  // (return 0):
  int test_saturation(int pp, int ms=MAX_REPEAT);

  // try harder with same p (no initialization, but same point set)
  int test_saturation_extra(int pp, int ms=MAX_REPEAT);

  // enlarge basis if dim(kernel)>0:
  int enlarge();

  // repeat testing p-saturation and enlarging until done: returns
  // log_p of index, or -1 if failed
  int do_saturation(int pp, int maxntries=10);

  // As above but for all primes in plist, returns success flag, sets
  // index and unsat = list of primes in plist at which saturation
  // failed = empty iff success flag=0
  int do_saturation(vector<int> plist,
		    long& index, vector<int>& unsat, int maxntries=10);
  int do_saturation(vector<long> plist,
                    long& index, vector<long>& unsat, int maxntries=10);

  // Either (when sat_bd=-1) auto-saturate, finding an upper bound ib
  // on saturation index, or (when sat_bd>0) saturate only for primes
  // up to min(ib,sat_bd).  If sat_low_bd>0 then we do no p-saturation
  // for p<sat_low_bd: e.g. set to 3 if we already know 2-saturaton.
  int saturate(vector<long>& unsat, long& index,
               long sat_bd=-1, long sat_low_bd=2,
               int egr=1, int maxntries=10);

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
  vec_l kernel_vector();
  // test if saturated:
  int is_saturated() {return rank==TLrank;}
  void show_q_tally(); // display list of q used with multiplicity.
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
// If sat_bd>-1 it is used to truncate the upper bound.
// No p-saturation is done for p<sat_low_bd.

int saturate_points(Curvedata& C, vector<Point>& points,
		    long& index, vector<long>& unsat,
		    long sat_bd=-1, long sat_low_bd=2,
                    int egr=1, int verbose=0);

#endif
