// bitspace.cc: implementation of class bitspace for handling F_2 spaces
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
 
#include <iostream>
#include "bitspace.h"

using namespace std;

bitspace::bitspace(long d)
{
  if(d<0)
    {
      cout<<"Error in bitspace constructor with negative dimension "<<d
	  <<"! replacing with 0\n";
      d=0;
    }
  if(d>32)
    {
      cout<<"Error in bitspace constructor with dimension "<<d
	  <<">32! replacing with 32\n";
      d=32;
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

