// bitspace.h: declaration of class bitspace for handling F_2 spaces
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
 
#if     !defined(_ECLIB_BITSPACE_H)
#define _ECLIB_BITSPACE_H      1       //flags that this file has been included

#include "templates.h"

class bitspace {
private: 
  long maxdim;
  long dim;
  vector<long> pivs;          // holds the position of the ith pivot
  vector<unsigned long> gens; // holds the ith basis element
  unsigned long bitmask;      // holds the bits of the pivs
public:
  explicit bitspace(long d);
  unsigned long getbitmask() {return bitmask;}
  long reduce(unsigned long& v, long start=0) const; 
  // reduces v mod this, returns minimal i such that the reduced v has 
  // ith bit set, or -1 if v reduces to 0.  Assumes already reduced for i<start
  int mask(unsigned long i) {return (i&bitmask)!=0;}
  void augment(unsigned long v, long piv);
  // uses reduced v to augment this, given piv as suitable pivot position in v
  int augment(unsigned long v)
  // uses v to augment this unless it it dependent, returns 1 if new
    {
      long j = reduce(v);
      if(j<0) return 0; 
      augment(v,j); return 1;
    }
};

inline int testbit(long a, long i) {return (a& (1<<i));}
inline int setbit( long& a, long i) {return (a|=(1<<i));}

// return the dot product (0/1) of a and b bitwise, with 0<=a,b<2^r
int dotbits(long a, long b, int r);
// return list of bits of a
vector<int> bits(long a, int r);
// recover a from its bit vector of length r
long from_bits(const vector<int>& aa, int r);
inline long from_bits(const vector<int>& aa) {return from_bits(aa, aa.size());}

// return a basis for the orthogonal complement of a<2^r (viewed as a bit vector of length r)
vector<long> dotperp(long a, int r);
// return a basis for the orthogonal complement of the span of a in alist (viewed as bit vectors of length r)
vector<long> dotperp(const vector<long>& alist, int r);

#endif
