// twist.cc: program to compute & display quadratic twists
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
 
#include <eclib/curve.h>
#include <eclib/getcurve.h>

int main(){

  initprimes("PRIMES",0);
  int verbose=0;
  vector<bigrational> ai(5);

while (getcurve(ai,verbose))
{
  bigint c4,c6,twist2;
  bigint v;
  long twist=1;

  Curvedata D(Curvedata(ai,v),1);
  Curvedata E = D;
  CurveRed CR = CurveRed(E);
  cout << "\nCurve is:  " << endl;  CR.display(cout);

  while(1)
    {
      cout << "\nEnter a twist value: ";
      cout << "\n(0 to set a new `original' curve, 1 to twist original,\n";
      cout << "or any other integer to twist immediately preceding output)\n ";
      cin>>ws;  if(cin.eof()) {cout<<endl; exit(0);}
      cin >> twist;
      if(twist==0) break;
      if(twist==1) { E=D;}
      else
	{
	  E.getci(c4, c6);
	  twist *= 4;        //corrects for twist not congruent -1 mod 4
	  twist2 = twist*twist;
	  c4 *= twist2;
	  c6 *= twist*twist2;
	  E = Curvedata(Curve(c4, c6), 1); //minimal form
	  CR=CurveRed(E);
	  cout << "\n\nE * " << twist/4 << " is:  \n"; CR.display(cout);
	  cout << endl;
	}
    }
}


} //ends main
