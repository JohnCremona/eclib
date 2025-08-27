// FILE CUSP.CC
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

#include "eclib/cusp.h"

// This function tests cusps for Gamma_0(N)-equivalence, unless
// plusflag is +1 in which case it tests for equivalence under
// <Gamma_0(N),-I>

int cusplist::cuspeq(const rational& c1, const rational& c2, int plusflag) const
{
  //  cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2<<endl;
  if (c1==c2) return 1;
  long p1 = num(c1), p2 = num(c2), q1 = den(c1), q2 = den(c2);
  if ((N->gcd(q1))!=(N->gcd(q2))) return 0;
  long s1,r1,s2,r2;
  bezout(p1,q1,s1,r1);  s1*=q2;
  bezout(p2,q2,s2,r2);  s2*=q1;
  long q3 = N->gcd(q1*q2);
  int ans = ((s1-s2)%q3==0);       // 1 iff [c1]=[c2]
  //  cout << "ans = "<<ans<<endl;
  if (ans || (plusflag!=+1)) return ans;
  ans = ((s1+s2)%q3==0);         // 1 iff [c1]=[-c2]
  //  cout << "ans = "<<ans<<endl;
  return ans;
}

long cusplist::index(const rational& c)
{  // adds c to list if not there already, and return index (offset by 1)
  auto it = std::find_if(list.begin(), list.end(),
                         [c, this] (const rational& d) { return cuspeq(c,d, N->plusflag);});
  if (it!=list.end())
    return 1 + std::distance(list.begin(), it);  // note offset
  list.push_back(c);
  //   cout<<"Adding c="<<c<<" as cusp number "<<number<<endl;
  return list.size();
}

long cusplist::index_1(const rational& c)
{ // adds c to list if not there already, and return index (offset by 1)
  // For use with minus space; only one of [c],[-c] is stored and the
  // index returned is negative if [-c] is the one listed and 0 if
  // [c]=[-c] (which are not listed)
  if (cuspeq(c,-c,0))
    return 0;

  int j = 0;
  for (unsigned int i=0; i<list.size() && j==0; i++)
    {
      if (cuspeq(c,list[i], 0))
        j = (i+1);  // note offset
      else
        if (cuspeq(-c,list[i], 0))
          j = -int(i+1);
    }
  if (j)
    return j;
  list.push_back(c);
  return list.size();
}

long cusplist::index_2(const rational& c)
{ // adds c to list if not there already, and return index (offset by 1)
  // For use with minus space; only store [c] if [c]=[-c]
  if (!cuspeq(c,-c,0)) {return 0;}
  for (unsigned int i=0; i<list.size(); i++)
    if (cuspeq(c,list[i], 0))
      return (i+1);  // note offset
  list.push_back(c);
  return list.size();
}
