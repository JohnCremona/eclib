// hilbert.cc: implementation of Hilbert symbol functions
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
 
#include <eclib/hilbert.h>
#include <eclib/marith.h>

//#define DEBUG_HILBERT

// In all the functions below, the value of the Hilbert symbol is 0 or
// 1 (as an int) rather than +1 or -1, for efficiency;  

int local_hilbert(const bigint& a, const bigint& b, const bigint& p)
{
  static const bigint zero(0);
  static const bigint  two(2);
  long alpha, beta;
  bigint u,v;
  int ans;

  if(is_zero(a)) {cout<<"Error in local_hilbert(): a==0\n"; return -1;}
  if(is_zero(b)) {cout<<"Error in local_hilbert(): b==0\n"; return -1;}

  if(is_zero(p)||is_negative(p)) // p=0 or p=-1 mean the infinite prime
    {
      if(is_positive(a)) return 0;
      if(is_positive(b)) return 0;
      return 1;
    }

  u=a; alpha = divide_out(u,p)%2;  // so a=u*p^alpha *square
  v=b; beta  = divide_out(v,p)%2;  // so b=v*p^beta  *square
  
  if(p==two)
    {
      //      ans = eps4(u)&eps4(v);
      ans = ((u+1)%4==0);
      if(ans) ans = ((v+1)%4==0);
      if(alpha) if(omega8(v)) ans=!ans;
      if(beta)  if(omega8(u)) ans=!ans;
      return ans;
    }

  // now p is odd

  ans = alpha&beta;
  if(ans) ans = ((p+1)%4==0);
  if(alpha) if(legendre(v,p)==-1) ans=!ans;
  if(beta)  if(legendre(u,p)==-1) ans=!ans;
  return ans;
}

int global_hilbert(const bigint& a, const bigint& b, const vector<bigint>& plist, bigint& badp)
{
#ifdef DEBUG_HILBERT
  cout<<"In global_hilbert("<<a<<","<<b<<"), plist = "<<plist<<endl;
#endif
  if(local_hilbert(a,b,0))
    {
      badp=0;
      return 1;
    }
#ifdef DEBUG_HILBERT
  cout<<"Passed local condition at infinity..."<<endl;
#endif
  return std::any_of(plist.begin(), plist.end(),
                     [a,b,&badp] (const bigint& p) {badp=p;return local_hilbert(a,b,p);});
}

int global_hilbert(const bigint& a, const bigint& b, bigint& badp)
{
  vector<bigint> plist=vector_union(pdivs(a),pdivs(b));
  return global_hilbert(a,b,plist,badp);
}

int global_hilbert(const quadratic& q, const bigint& d, bigint& badp)
{
  bigint D = q.disc();
  vector<bigint> plist = vector_union(pdivs(D),pdivs(d));
  return global_hilbert(q[0]*d,D,plist,badp);
}

int global_hilbert(const quadratic& q, const bigint& d, const vector<bigint>& plist, bigint& badp)
{
  return global_hilbert(q[0]*d,q.disc(),plist,badp);
}
