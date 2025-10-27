// bitspace.cc: implementation of class bitspace for handling F_2 spaces
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
 
#include "eclib/interface.h"
#include "eclib/bitspace.h"

using NTL::bit;

bitspace::bitspace(long d)
{
  if(d<0)
    {
      cout<<"Error in bitspace constructor with negative dimension "<<d
	  <<"! replacing with 0\n";
      d=0;
    }
  if(d>NTL_BITS_PER_LONG)
    {
      cout<<"Error in bitspace constructor with dimension "<<d
	  <<">" <<NTL_BITS_PER_LONG<<"! replacing with "
          <<NTL_BITS_PER_LONG<<"\n";
      d=NTL_BITS_PER_LONG;
    }
  maxdim=d;
  pivs.resize(maxdim);
  gens.resize(maxdim);
  dim=0;
  bitmask=0;
}

long bitspace::reduce(unsigned long& v, long start) const
{
  long j;
  for(j=start; j<dim; j++) if(testbit(v,pivs[j])) {v=v^gens[j];}
  for(j=maxdim-1; j>=0; j--) if(testbit(v,j)) {return j;}
  return -1;
}

void bitspace::augment(unsigned long v, long piv)
{
  gens[dim]=v;
  pivs[dim]=piv;
  bitmask |= (1<<piv);
  dim++;
}

// return the dot product (0/1) of a and b bitwise, with 0<=a,b<2^r
int dotbits(long a, long b, int r)
{
  int x=0;
  for (int i=0; i<r; i++)
    x ^= (bit(a,i)*bit(b,i));
  return x;
}

// return list of bits of a
vector<int> bits(long a, int r)
{
  vector<int> ans(r);
  for(int i=0; i<r; i++)
    ans[i] = bit(a,i);
  // assert (a==from_bits(ans,r));
  return ans;
}

// recover a from its bit vector of length r
long from_bits(const vector<int>& aa, int r)
{
  long a=0;
  for(int i=0; i<r; i++)
    if (aa[i])
      a |= (1<<i); // sets the i'th bit to 1
  // assert (aa==bits(a,r));
  return a;
}
