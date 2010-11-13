// FILE ECNF.CC: program for newform construction from an elliptic curve
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
//
//
#include "interface.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"
#include "newforms.h"
#include "curve.h"

int main(void)
{
  int verbose=0;

  // Read in the curve, minimise and construct CurveRed (needed for
  // conductor and Traces of Frobenius etc.)
  Curve C;
  cout << "Enter curve: "; cin >> C;
  Curvedata CD(C,1); // minimise
  CurveRed CR(CD);
  bigint N = getconductor(CR);
  int n = I2int(N);
  cout << ">>> Level = conductor = " << n << " <<<" << endl;
  cout << "Minimal curve = " << (Curve)(CR) << endl;
  cout<<endl;

  // Construct newforms class (this does little work)
  int sign=1;

  cout<<"Enter sign (1,-1,0 for both):"; cin>>sign;

  newforms nf(n,verbose);

  // Create the newform from the curve (first create the homspace,
  // then split off the eigenspace)
  nf.createfromcurve(sign,CR);

  // Display newform info
  cout << "Newform information:"<<endl;
  nf.display();
  cout<<endl;

  // Display modular symbol info
  cout << "Modular symbol map:"<<endl;
  nf.display_modular_symbol_map();

  // Compute more modular symbols as prompted:

  rational r; long nu,de;
  cout << "Computation of further modular symbols {0,r} for rational r:"<<endl;
  while(1)
    {
      cout<<"Enter numerator and denominator of r: "; 
      cin>>ws;  
      if(cin.eof()) {cout<<endl; break;}
      cin >> nu >> de; r=rational(nu,de);
      if((nu==0)&&(de==0)) {cout<<endl; break;}
      if(sign==+1)
        cout<<"{0,"<<r<<"} -> "<< nf.plus_modular_symbol(r)<<endl;
      if(sign==-1)
        cout<<"{0,"<<r<<"} -> "<< nf.minus_modular_symbol(r)<<endl;
      if(sign==0)
	{
	  pair<rational,rational> s = nf.full_modular_symbol(r);
	  cout<<"{0,"<<r<<"} -> ("<< s.first << "," << s.second << ")" <<endl;
	}
    }

  for(de=1; de<20; de++)
    for(nu=0; nu<de; nu++)
      {
	if(gcd(nu,de)==1)
	  {
	    r=rational(nu,de);
            if(sign==+1)
              cout<<"{0,"<<r<<"} -> "<< nf.plus_modular_symbol(r)<<endl;
            if(sign==-1)
              cout<<"{0,"<<r<<"} -> "<< nf.minus_modular_symbol(r)<<endl;
            if(sign==0)
              {
                pair<rational,rational> s = nf.full_modular_symbol(r);
                cout<<"{0,"<<r<<"} -> ("<< s.first << "," << s.second << ")" <<endl;
              }
	  }
      }
}       // end of main()
