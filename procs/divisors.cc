// divisors.cc: divisors of a positive integer

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
 //

#include "marith.h"
#define MAXPRIME 1000000

int main(int argc, char *argv[])
{
  initprimes("PRIMES",0);
  //  cout<<argv[1]<<endl;
  bigint n = BIGINT((const char*)argv[1]);
//  cout<<n<<" : "<<flush;
  vector<bigint> dlist = posdivs(n);
  //  cout<<dlist<<endl;
  copy(dlist.begin(),dlist.end(), ostream_iterator<bigint>(cout, " "));
  cout<<endl;
 }
