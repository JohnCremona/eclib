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

// Test whether a curve with good reduction outside S and this j-invariant could exist
// (using criteria from Cremona-Lingham)
int is_j_possible(const bigrational& j, const vector<bigint>& S)
{
  static const bigrational j1728(1728);
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
      vector<bigint> pwlist;
      pwlist.reserve(wlist.size()*n);
      for (auto w: wlist)
        for (auto pp: ppowers)
          pwlist.push_back(pp*w);
      wlist = pwlist;
    }
  return wlist;
}

// Return list of curves with good reduction outside S and j=1728
// using Cremona-Lingham Prop.4.2 and the remark following
vector<CurveRed> egros_from_j_1728(const vector<bigint>& S)
{
  vector<CurveRed> Elist;
  int no2 = std::find(S.begin(), S.end(), BIGINT(2)) == S.end();
  vector<bigint> wlist = twist_factors(S, 4);
  for (auto w: wlist)
    {
      if (no2) w *= 4;
      Curve E(0,0,0,w,0);
      CurveData Emin(E, 1);
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
  vector<CurveRed> Elist;
  int no3 = std::find(S.begin(), S.end(), BIGINT(3)) == S.end();
  if (no3)
    return Elist;
  int no2 = std::find(S.begin(), S.end(), BIGINT(2)) == S.end();
  vector<bigint> wlist = twist_factors(S, 6);
  for (auto w: wlist)
    {
      if (no2) w *= 16;
      Curve E(0,0,0,0,w);
      CurveData Emin(E, 1);
      CurveRed Ered(Emin);
      if (Ered.has_good_reduction_outside_S(S))
        Elist.push_back(Ered);
    }
  return Elist;
}

vector<CurveRed> egros_from_j(const bigrational& j, const vector<bigint>& S)
{
  vector<CurveRed> Elist;

  // Return empty list if necessary conditions fail:
  if (!is_j_possible(j, S))
    return Elist;
  // Call special function if j=1728:
  if (j==1728)
    return egros_from_j_1728(S);
  // Call special function if j=0:
  if (is_zero(j))
    return egros_from_j_0(S);

  // Now j is not 0 or 1728, we take quadratic twists of a base curve:

  vector<bigint> wlist = twist_factors(S, 2);

  bigint n = num(j);
  bigint m = n-1728*den(j);
  a4 = -3*n*m;
  a6 = -2*n*m**2;  // the base curve [0,0,0,a4,a6] has the right j-invariant
  for (auto w: wlist)
    {
      if (no2)
        w *= 16;
      bigint w2 = w*w;
      bigint w3 = w*w2;
      Curve E(0,0,0,w2*a4,w3*a6);
      CurveData Emin(E, 1);
      CurveRed Ered(Emin);
      if (Ered.has_good_reduction_outside_S(S))
        Elist.push_back(Ered);
    }
}

// Test whether N is a possible conductor for j=0: 3|N, no p||N and
// usual bounds on ord_p(N)
int is_N_possible_j_0(const bigint& N, const vector<bigint>& support)
{
  if (!div(three,N))
    return 0;
  for (auto p: support)
    {
      int np = val(p,N);
      int okp = (np>=2) && (np<= (p==two? 8 : (p==three? 5 : 2)));
      if (!okp)
        return 0;
    }
  return 1;
}

// Test whether N is a possible conductor for j=1728: 2|N, no p||N and
// usual bounds on ord_p(N)
int is_N_possible_j_1728(const bigint& N, const vector<bigint>& support)
{
  if (!div(two,N))
    return 0;
  for (auto p: support)
    {
      int np = val(p,N);
      int okp = (np>=2) && (np<= (p==two? 8 : (p==three? 5 : 2)));
      if (!okp)
        return 0;
    }
  return 1;
}
