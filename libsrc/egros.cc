// egros.cc:   Implementation of functions for elliptic curves with good reduction outside S
//////////////////////////////////////////////////////////////////////////
//
// Copyright 2022 John Cremona
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

#include <eclib/egros.h>

// Test whether a curve with good reduction outside S and this j-invariant could exist
// (using criteria from Cremona-Lingham)
int is_j_possible(const bigrational& j, const vector<bigint>& S)
{
  static const bigint three(3);
  bigint nj(num(j)), dj(den(j));
  bigint mj = nj-1728*dj;
  if (is_zero(mj)) // j==1728: Cremona-Lingham Prop. 4.1
    return 1;
  if (is_zero(nj)) // j==0: Cremona-Lingham Prop. 4.1
    return std::find(S.begin(), S.end(), three) != S.end();
  if (!is_S_integral(j, S))
    return 0;
  return // Cremona-Lingham Prop. 3.2
    is_nth_power(prime_to_S_part(nj, S), 3)
    &&
    is_nth_power(prime_to_S_part(mj, S), 2);
}

// Return integers representing QQ(S,n)
vector<bigint> twist_factors(const vector<bigint>& S, int n)
// only intended for n=2,4,6
{
  static const bigint one(1);
  vector<bigint> wlist = {one,-one};
  for (auto p: S)
    {
      vector<bigint> ppowers = {one};
      for (int i=1; i<n; i++)
        ppowers.push_back(ppowers[i-1]*p);
      wlist = multiply_lists(wlist, ppowers);
    }
  return wlist;
}

// Return list of curves with good reduction outside S and j=1728
// using Cremona-Lingham Prop.4.2 and the remark following
vector<CurveRed> egros_from_j_1728(const vector<bigint>& S)
{
  static const bigint zero(0);
  static const bigint two(2);
  vector<CurveRed> Elist;
  int no2 = std::find(S.begin(), S.end(), two) == S.end();
  vector<bigint> wlist = twist_factors(S, 4);
  for (auto w: wlist)
    {
      if (no2) w *= 4;
      Curve E(zero,zero,zero,w,zero);
      Curvedata Emin(E, 1);
      CurveRed Ered(Emin);
      if (Ered.has_good_reduction_outside_S(S))
        Elist.push_back(Ered);
    }
  return Elist;
}

// Return list of curves with good reduction outside S and j=1728
// using Cremona-Lingham Prop.4.1 and the remark following
vector<CurveRed> egros_from_j_0(const vector<bigint>& S)
{
  static const bigint zero(0);
  static const bigint two(2);
  static const bigint three(3);
  vector<CurveRed> Elist;
  int no3 = std::find(S.begin(), S.end(), three) == S.end();
  if (no3)
    return Elist;
  int no2 = std::find(S.begin(), S.end(), two) == S.end();
  vector<bigint> wlist = twist_factors(S, 6);
  for (auto w: wlist)
    {
      if (no2) w *= 16;
      Curve E(zero,zero,zero,zero,w);
      Curvedata Emin(E, 1);
      CurveRed Ered(Emin);
      if (Ered.has_good_reduction_outside_S(S))
        Elist.push_back(Ered);
    }
  return Elist;
}

vector<CurveRed> egros_from_j(const bigrational& j, const vector<bigint>& S)
{
  static const bigint zero(0);
  static const bigint two(2);
  static const bigint three(3);

  vector<CurveRed> Elist;

  // Return empty list if necessary conditions fail:
  if (!is_j_possible(j, S))
    return Elist;

  bigint n = num(j);
  bigint m = n-1728*den(j);

 // Call special function if j=1728:
  if (is_zero(m))
    return egros_from_j_1728(S);

  // Call special function if j=0:
  if (is_zero(n))
    return egros_from_j_0(S);

  // Now j is not 0 or 1728, we take quadratic twists of a base curve:

  vector<bigint> Sx = S;
  vector<bigint> Sy = pdivs(n*m*(n-m));
  vector<bigint> extra_primes;
  for (auto p: Sy)
    {
      if (std::find(Sx.begin(), Sx.end(), p) == Sx.end())
        {
          Sx.push_back(p);
          extra_primes.push_back(p);
        }
    }
  vector<bigint> wlist = twist_factors(Sx, 2);

  bigint a4 = -3*n*m;
  bigint a6 = -2*n*m*m;  // the base curve [0,0,0,a4,a6] has the right j-invariant

  // We'll test twists of [0,0,0,a4,a6], whose discriminant is
  // 1728n^2m^3(n-m). For primes p>3 not in S we already have
  // ord_p(n)=0(3), ord_p(m)=0(2) and ord_p(n-m)=ord_p(denom(j))=0, so
  // ord_p(disc)=0(6).  For there to be any good twists we want
  // ord_p(disc)=0(12).  The twist by w is [0,0,0,w^2*a4,w^3*a6] which
  // has disc w^6 times the that of the base curve.  So for the
  // 'extra' primes p (not in S) with ord_p(n)=3(6) we must have ord_p(w) odd,
  // while for p with ord_p(m)=2(4) we must have ord_p(w) odd.

  vector<bigint> a4a6primes;
  for (auto p: extra_primes)
    {
      if ((p==two) || (p==three))
        continue;
      if ((val(p,n)%6==3) || (val(p,m)%4==2))
        a4a6primes.push_back(p);
    }
  // cout << "extra_primes = "<<a4a6primes<<endl;
  // cout << "a4a6primes =   "<<a4a6primes<<endl;

  int no2 = std::find(S.begin(), S.end(), two) == S.end();

  for (auto w: wlist)
    {
      for (auto p: a4a6primes)
        if(val(p,w)%2==0) continue;

      if (no2)
        w *= 16;
      bigint w2 = w*w;
      bigint w3 = w*w2;
      Curve E(zero,zero,zero,w2*a4,w3*a6);
      Curvedata Emin(E, 1);
      CurveRed Ered(Emin);
      if (Ered.has_good_reduction_outside_S(S))
        Elist.push_back(Ered);
    }
  return Elist;
}

int conductor_exponent_bound(const bigint& p)
{
  static const bigint two(2);
  static const bigint three(3);
  return (p==two? 8 : (p==three? 5 : 2));
}

int is_N_possible_helper(const bigint& N, const vector<bigint>& support)
{
  for (auto p: support)
    {
      int np = val(p,N);
      if ((np < 2) or (np > conductor_exponent_bound(p))) return 0;
    }
  return 1;
}

// Test whether N is a possible conductor for j=0: 3|N, no p||N and
// usual bounds on ord_p(N)

int is_N_possible_j_0(const bigint& N, const vector<bigint>& support)
{
  static const bigint three(3);
  return div(three,N) and is_N_possible_helper(N,support);
}

// Test whether N is a possible conductor for j=1728: 2|N, no p||N and
// usual bounds on ord_p(N)

int is_N_possible_j_1728(const bigint& N, const vector<bigint>& support)
{
  static const bigint two(2);
  return div(two,N) and is_N_possible_helper(N,support);
}
