// bitspace.cc: implementation of class bitspace for handling F_2 spaces
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
 
#include <iostream>
#include <eclib/interface.h>
#include <eclib/bitspace.h>



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
  pivs = new long[maxdim];
  gens = new unsigned long[maxdim];
  dim=0;
  bitmask=0;
}

bitspace::~bitspace()
{
  delete[]pivs;
  delete[]gens;
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

