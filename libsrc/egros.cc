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
 
int is_j_possible(const bigrational& j, const vector<bigint>& S)
{
  if (j==1728)
    return 1;
  if (is_zero(j))
    return std::find(S.begin(), S.end(), BIGINT(3)) != S.end();
  if (!is_S_integral(j, S))
    return 0;
  return is_nth_power(prime_to_S_part(num(j), S), 3)
    &&
    is_nth_power(prime_to_S_part(num(j)-1728*den(j), S), 2);
}

vector<bigint> twist_factors(const vector<bigint>& S, int n)
// only intended for n=2,4,6
{
  vector<bigint> wlist = {1,-1};
  for (auto pi=S.begin(); pi!=S.end(); ++pi)
    {
      bigint p = *pi;
      vector<bigint> ppowers = {BIGINT(1)};
      for (int i=1; i<n; i++)
        ppowers.push_back(ppowers[i-1]*p);
      vector<bigint> pwlist;
      pwlist.reserve(wlist.size()*n);
      for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
        for (auto pp = ppowers.begin(); pp!=ppowers.end(); ++pp)
          pwlist.push_back((*pp)*(*wi));
      wlist = pwlist;
    }
  return wlist;
}

vector<CurveRed> egros_from_j_1728(const vector<bigint>& S)
{
  vector<CurveRed> Elist;
  int no2 = std::find(S.begin(), S.end(), BIGINT(2)) == S.end();
  vector<bigint> wlist = twist_factors(S, 4);
  for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
    {
      bigint w = *wi;
      if (no2) w *= 4;
      Curve E(0,0,0,w,0);
      CurveData Emin(E, 1);
      CurveRed Ered(Emin);
      if (Ered.has_good_reduction_outside_S(S))
        Elist.push_back(Ered);
    }
  return Elist;
}

vector<CurveRed> egros_from_j_0(const vector<bigint>& S)
{
  vector<CurveRed> Elist;
  int no3 = std::find(S.begin(), S.end(), BIGINT(3)) == S.end();
  if (no3)
    return Elist;
  int no2 = std::find(S.begin(), S.end(), BIGINT(2)) == S.end();
  vector<bigint> wlist = twist_factors(S, 6);
  for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
    {
      bigint w = *wi;
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
  if (!is_j_possible(j, S))
    return Elist;
  if (j==1728)
    return egros_from_j_1728(S);
  if (is_zero(j))
    return egros_from_j_0(S);

  vector<bigint> wlist = twist_factors(S, 2);

  bigint n = num(j);
  bigint m = n-1728*den(j);
  a4 = -3*n*m;
  a6 = -2*n*m**2;
  for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
    {
      bigint w = *wi;
      bigint w2 = w*w;
      bigint w3 = w*w2;
      if (no2)
        w *= 16;
      Curve E(0,0,0,w2*a4,w3*a6);
      CurveData Emin(E, 1);
      CurveRed Ered(Emin);
      if (Ered.has_good_reduction_outside_S(S))
        Elist.push_back(Ered);
    }
}

