// FILE ECNF.CC: program for newform construction from an elliptic curve
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
//
//
#include <eclib/interface.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/cperiods.h>
#include <eclib/newforms.h>
#include <eclib/curve.h>
#include <eclib/getcurve.h>

// defining to be 1 this causes display_modular_symbol_map() to be
// called with the check option on, which takes a lot longer but
// checks that the algebraic modular symbols agree with numerical
// ones, even in their sign.
#define DEBUG 0

int main(void)
{
  int verbose=0;
  vector<bigrational> ai(5);
  bigint v;

  // Read in curves, minimise and construct CurveRed (needed for
  // conductor and Traces of Frobenius etc.)
  while (getcurve(ai,verbose))
    {
      Curvedata CD(ai,v);
      CurveRed CR(CD);
      bigint N = getconductor(CR);
      int n = I2int(N);
      cout << ">>> Level = conductor = " << n << " <<<" << endl;
      cout << "Minimal curve = " << (Curve)(CR) << endl;
      cout<<endl;

  // Construct newforms class
  int sign=1;

  cout<<"Enter sign (1,-1,0 for both):"; cin>>sign;

  newforms nf(n, DEFAULT_MODULUS, verbose);

  // Create the newform from the curve (first create the homspace,
  // then split off the eigenspace)
  nf.createfromcurve(sign,CR);

  // Display newform info
  cout << "Newform information:"<<endl;
  nf.display();
  cout<<endl;

  // Display modular symbol info
  cout << "Modular symbol map (";
  if (sign!=-1) cout << "+";
  if (sign==0) cout << ",";
  if (sign!=1) cout << "-";
  cout << ")" << endl;
  nf.display_modular_symbol_map(DEBUG);

  // Compute more modular symbols as prompted:

  rational r; long nu,de;
  cout << "\nComputation of further modular symbols" << endl << endl;
  cout << "Base point? (enter 0 for 0, or 1 for oo) ";
  int base_at_infinity;
  cin >> base_at_infinity;
  string base = (base_at_infinity?"oo":"0");
  cout << "Values of {" << base << ",r} for rational r:"<<endl;
  while(1)
    {
      cout<<"Enter numerator and denominator of r: ";
      cin>>ws;
      if(cin.eof()) {cout<<endl; break;}
      cin >> nu >> de; r=rational(nu,de);
      if((nu==0)&&(de==0)) {cout<<endl; break;}
      if(sign==+1)
        cout<<"{"<<base<<","<<r<<"} -> "<< nf.plus_modular_symbol(r, 0, base_at_infinity)<<endl;
      if(sign==-1)
        cout<<"{"<<base<<","<<r<<"} -> "<< nf.minus_modular_symbol(r, 0, base_at_infinity)<<endl;
      if(sign==0)
	{
	  pair<rational,rational> s = nf.full_modular_symbol(r, 0, base_at_infinity);
	  cout<<"{"<<base<<","<<r<<"} -> ("<< s.first << "," << s.second << ")" <<endl;
	}
    }

  cout << "All modular symbols with bounded denominator" << endl << endl;
  int dmax;
  cout << "Enter maximum denominator (0 for none): "; cin >> dmax;
  for(de=1; de<=dmax; de++)
    for(nu=0; nu<de; nu++)
      {
	if(gcd(nu,de)==1)
	  {
	    r=rational(nu,de);
            if(sign==+1)
              cout<<"{"<<base<<","<<r<<"} -> "<< nf.plus_modular_symbol(r, 0, base_at_infinity)<<endl;
            if(sign==-1)
              cout<<"{"<<base<<","<<r<<"} -> "<< nf.minus_modular_symbol(r, 0, base_at_infinity)<<endl;
            if(sign==0)
              {
                pair<rational,rational> s = nf.full_modular_symbol(r, 0, base_at_infinity);
                cout<<"{"<<base<<","<<r<<"} -> ("<< s.first << "," << s.second << ")" <<endl;
              }
	  }
      }
    }   // end of curve input loop
}       // end of main()
