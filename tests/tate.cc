// tate.cc: program to call Tate's algorithm and display curve details
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
 
#include <eclib/curve.h>
#include <eclib/getcurve.h>

void line(int w) {int wid=w; while(wid--)cout << "-";cout<<endl;}

int main(void)
{
  initprimes("PRIMES",0);
  line(80);
        
  int verb=1;
  bigint v;
  vector<bigrational> ai(5);

  while (getcurve(ai,verb))
    {
      Curvedata D(ai,v); 
      cout<<"Integral model: "<<(Curve)D<<endl;
      CurveRed C(D);  // minimises here
      cout << "Curve ";
      C.display(cout);
      line(80);
    }
}
          

