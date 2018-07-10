// bitspace.h: declaration of class bitspace for handling F_2 spaces
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
 
#if     !defined(_ECLIB_BITSPACE_H)
#define _ECLIB_BITSPACE_H      1       //flags that this file has been included

class bitspace {
private: 
  long maxdim;
  long dim;
  long * pivs; // holds the position of the ith pivot
  unsigned long * gens; // holds the ith basis element
  unsigned long bitmask;   // holds the bits of the pivs
public:
  bitspace(long d);
  ~bitspace();
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

#endif
