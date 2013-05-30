//tcon2.cc: conic test program
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
 
#include <eclib/marith.h>
#include <eclib/quadratic.h>
#include <eclib/conic.h>

#define VERBOSITY 0
#ifndef CONIC_METHOD
#define CONIC_METHOD 4
#endif
//#define TEST_PARAM

int main()
{
  initprimes("PRIMES",VERBOSITY);

  bigint a,b,c,d,x0,y0,z0,disc;
  bigint q[3];

  cin >> a >> b >> c >> d;
  q[0]=a; q[1]=b; q[2]=c;
  int res = solve_conic(q,d,x0,y0,z0,CONIC_METHOD);
  if(!res) {x0=y0=z0=0;}
  cout<<x0<<" "<<y0<<" "<<z0<<endl;
}



