// file fixc6.cc: implementation of fixc6 class
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

// the constructor initializes the static data member with all the
// fixed values

#include "marith.h"
#include <map>
#include "fixc6.h"

fixc6::fixc6()
{
#ifdef MPFP // Multi-Precision Floating Point
  return;
#else

  long n=1; int i; bigint c4,c6;
  ifstream datafile("fixc6.data");
  while(n) 
    {
      datafile>>n>>i>>c6;
      if(n!=0) fixc6table[pair<long,int>(n,i)] = c6;
    }
  datafile.close();
  datafile.open("fixc4.data");
  n=1;
  while(n) 
    {
      datafile>>n>>i>>c4;
      if(n!=0) fixc4table[pair<long,int>(n,i)] = c4;
    }
  datafile.close();

#endif // MPFP
} // end of constructor

void fixc6::operator()(long N, int i, bigint& c4, bigint& c6)
{
  pair<long,int> key(N,i+1);
  map<pair<long,int>,bigint>::const_iterator j = fixc6table.find(key); 
  if(j!=fixc6table.end()) c6=j->second;
  j = fixc4table.find(key); 
  if(j!=fixc4table.end()) c4=j->second;
  return; 
}

map< pair<long,int>, bigint > fixc6::fixc6table;
map< pair<long,int>, bigint > fixc6::fixc4table;
fixc6 c4c6fixer;  // the one and only instance of the class

