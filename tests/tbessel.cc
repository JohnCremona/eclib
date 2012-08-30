// tbessel.cc: test program for Bessel functions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include <iostream>
#include <eclib/kbessel.h>

using namespace std;

int main()
{
  cout.precision(15);
  double x;
  int debug=0;
//  cout << "Debug? "; cin>>debug;
  while(cout<<"Enter x: ", cin>>x, x!=0.0)
    {
      cout << "x = " << x;
      cout << "\tK0(x) = " << kbessel(0,x,debug) << flush;
      cout << "\tK1(x) = " << kbessel(1,x,debug) << endl;
    }      
}
