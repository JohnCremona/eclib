// hilbert.h: declarations of Hilbert symbol functions
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
 
#ifndef _HILBERT_H
#define _HILBERT_H      1
                           //flags that this file has been included

// In all the functions below, the value of the Hilbert symbol is 0 or
// 1 (as an int) rather than +1 or -1, for efficiency;  

inline int eps4(const bigint& u) // u must be odd; m1pow(eps4(u))=chi4(u)
{
  return (u+1)%4==0;  // so 1 mod 4 gives 0, 
                      //    3 mod 4 gives 1
}

inline int omega8(const bigint& u) // u must be odd; m1pow(omega8(u))=chi2(u)
{
  return ((u-3)%8==0)||((u+3)%8==0);  // so 1,7 mod 8 give 0, 
                                      //    3,5 mod 8 give 1
}

// Use p=0 for the infinite prime

int local_hilbert(const bigint& a, const bigint& b, const bigint& p);

inline int local_hilbert(const bigint& a, const bigint& b, const long& p)
{
  return local_hilbert(a,b,BIGINT(p));
}

// Returns 0 if soluble at all primes in list (or all dividing a*b, or
// all dividing disc(q)*d) and at infinity; otherwise returns 1 and
// puts the first insoluble prime into badp

int global_hilbert(const bigint& a, const bigint& b, const vector<bigint>& plist, bigint& badp);

int global_hilbert(const bigint& a, const bigint& b, bigint& badp);

int global_hilbert(const quadratic& q, const bigint& d, const vector<bigint>& plist, bigint& badp);

int global_hilbert(const quadratic& q, const bigint& d, bigint& badp);

#endif
