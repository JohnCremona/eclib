// tate.cc: program to call Tate's algorithm and display curve details
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
 
#include "curve.h"
#include "reader.h"

void line(int w) {int wid=w; while(wid--)cout << "-";cout<<endl;}

int main(void)
{
  initprimes(string("PRIMES").c_str(),0);
  line(80);
        
  CurveReader in;
  Curve E;

  while (in>>E)
    {
      Curvedata D(E,1);  // minimises here
      cout<<"Standard minimal model: "<<(Curve)D<<endl;
      CurveRed C(D);
      cout << endl << "Curve ";
      C.display(cout);
      line(80);
    }
}
          

