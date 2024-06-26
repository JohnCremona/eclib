// htconst.h:  declarations of functions for height bounds
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
 
// allow for multiple includes
#ifndef _ECLIB_HTCONST_H_
#define _ECLIB_HTCONST_H_

#include <eclib/elog.h>
#include <eclib/egr.h>
#include <eclib/sieve_search.h>

double silverman_bound(const Curvedata& CD);
double cps_bound(const Curvedata& CD);
inline double height_constant(const Curvedata& CD)
{
  //  return silverman_bound(CD);
  double b1 = silverman_bound(CD), b2 = cps_bound(CD);
  return min(b1,b2);
}
double egr_height_constant(const Curvedata& CD);

bigfloat lower_height_bound_search(const Curvedata& CD, const bigfloat& reg);
bigfloat Gamma_n(long n); // Gamma(n) = (n-1)!
bigfloat Gamma_n_plus_half(long n); // Gamma(n+1/2) = (2n)!sqrt(pi) / (4^n*n!)

bigfloat lattice_const(int r);

// returns lower bound for height of non-torsion points, following
// Cremona & Siksek in ANTS7
bigfloat lower_height_bound(const Curvedata& CD, int egr);

// Class to find point of minimal height by searching.

// If egr_flag is set it ignores points which do not have everywhere
// good reduction (at all finite primes)

class point_min_height_finder : public point_processor {
  Curvedata *E;
  ComponentGroups CG;  // used if egr_flag to test for egr
  bigint a1,a2,a3,a4,a6; 
  vector<bigint> c;
  int iso, egr_flag, verbose;
  bigfloat min_ht;
  Point Pmin;
  vector<Point> all_points;  //store points found
 public:
  explicit point_min_height_finder(Curvedata* EE, int egr=0, int verb=0);
  ~point_min_height_finder() {};
  int process(const bigint& x, const bigint& y, const bigint& z) override;
  void search(bigfloat h_lim);
  bigfloat get_min_ht() const {return min_ht;}
  Point get_min_ht_point() const {return Pmin;}
  vector<Point> points() const {return all_points;}
};

// class Interval represents a closed interval [lh,rh] where either
// empty=1; or empty=0 and lh <= rh; flags rhinf, lhinf denote
// rh=infty and lh=-infty resp.

class Interval {
  bigfloat lh, rh;
  bool empty, lhinf, rhinf;
public:
  Interval()     : empty(0), lhinf(1), rhinf(1) {show(1);}
  Interval(const bigfloat& a, const bigfloat& b)
    :lh(a), rh(b), empty(a>b), lhinf(0), rhinf(0)  {show(2);}
  explicit Interval(const bigfloat& a)
    :lh(a), empty(0), lhinf(0), rhinf(1)  {show(3);}
  Interval(const Interval& I)
    :lh(I.lh), rh(I.rh), empty(I.empty), lhinf(I.lhinf), rhinf(I.rhinf)  {show(4);}
  void operator=(const Interval& I) 
  {lh=I.lh; rh=I.rh; empty=I.empty; lhinf=I.lhinf; rhinf=I.rhinf;show(5);}
  friend ostream& operator<< (ostream&s, const Interval&);
  //  void show(int i) {cout<<i<<": "<<(*this)<<endl;}
  void show(int i) {;}
  bool is_empty() const {return empty;}
  void intersect(const Interval& I);
  friend Interval intersect(const Interval& I, const Interval& J);
};

inline Interval intersect(const Interval& I, const Interval& J)
{
  Interval K(I);  K.intersect(J); return K;
}

vector<Interval> intersect(const vector<Interval>& L1, const vector<Interval>& L2);
//void output(const vector<Interval>& L);

// class Interval01 represents a closed subinterval [lh,rh] of
// [0,1], where either empty=1; or empty=0 and lh <= rh.

class Interval01 {
  bigfloat lh, rh;
  bool empty;
public:
  Interval01()     : lh(to_bigfloat(0)), rh(to_bigfloat(1)), empty(0)  {;}
  Interval01(const bigfloat& a, const bigfloat& b)
    :lh(a), rh(b), empty(a>b)  {;}
  friend ostream& operator<< (ostream&s, const Interval01&);
  //  void show(int i) {cout<<i<<": "<<(*this)<<endl;}
  void show(int i) {;}
  bool is_empty() const {return empty;}
  void intersect(const Interval01& I);
  friend Interval01 intersect(const Interval01& I, const Interval01& J);
  friend Interval01 operator/(const Interval01& I, const long n);
  friend Interval01 operator+(const Interval01& I, const bigfloat& shift);
};

inline Interval01 intersect(const Interval01& I, const Interval01& J)
{
  Interval01 K(I);  K.intersect(J); return K;
}

vector<Interval01> intersect(const vector<Interval01>& L1, 
			     const vector<Interval01>& L2);

// Class to compute lower bound for height of non-torsion points of good
// reduction, following Cremona & Siksek in ANTS7

class CurveHeightConst : public CurveRed, Cperiods {
  bigfloat c;    // archimidean constribution
  bigfloat e3;   // largest (or only) real root
  bigfloat lower, upper;
  int n_max;
  long e_p(long p); // fetch /compute exponent at p
  map<long,long> ann;  // stored exponents of E0/E1 for first few primes
  bigfloat D(long n);  // fetch /compute denomContrib at n
  map<long,bigfloat> DE; // stored D_E(n) values
  bigfloat Bnmu(long n, const bigfloat& mu) // = B_n(mu)
  { return exp(n*n*mu+c-D(n));  }
  int test_target(const bigfloat& target, long k);
  vector<Interval> canonicalHeightInterval(const bigfloat& target, long k);
  vector<Interval01> canonicalHeightInterval01(const bigfloat& target, long k);
  void compute_phase1();
  void compute_phase2();
  vector<Interval> solveLEQ(long n, const bigfloat& B); 
  vector<Interval> solveGEQ(long n, const bigfloat& B); 
  vector<Interval01> solveLEQ01(long n, const bigfloat& B); 
  vector<Interval01> solveGEQ01(long n, const bigfloat& B); 
  vector<bigfloat> solveEllNPower(long n, const bigfloat& x1); 
  vector<bigfloat> solveEllNPower01(long n, const bigfloat& x1); 
  bigfloat psi(const bigfloat& x);
  vector<bigfloat> ordinates(const bigfloat& x);
  vector<bigcomplex> ztopoint(const bigcomplex& z)
  {return ellztopoint(z, bigcomplex(I2bigfloat(a1)), bigcomplex(I2bigfloat(a2)), bigcomplex(I2bigfloat(a3)));}
  bigcomplex pointtoz(const bigfloat& x, const bigfloat& y)
  {return ellpointtoz(*this, *this, x, y);}
public:
  explicit CurveHeightConst(CurveRed& CR);
  void compute() {compute_phase1(); compute_phase2(); }
  bigfloat get_value() const {return lower;} 
};

#endif
