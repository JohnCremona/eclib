// getcurve.cc: implementation of function getcurve() for curve input
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
#include "getcurve.h"

int getcurve(Curvedata& CD, Curvedata& CD_orig, 
	     bigint& u, bigint& r, bigint& s, bigint& t, 
	     int& change, int verb)
{
  Curve C0;
  int sing=1;
  while(sing)
    {
      if(verb) cout  << "Enter curve: ";
      cin >> C0;
      if(C0.isnull()) return 0;  // quitting condition
      CD_orig = Curvedata(C0,0); // DON'T change coords
      CD = CD_orig.minimalize(u,r,s,t);
      if (verb) cout << endl;
      sing=CD.isnull();
      if(sing) // input curve was singular, non-null
	{
	  if(verb) cout<<"Curve "<<C0<<" is singular\n";
	}
    }
  cout << "Curve "<< (Curve)C0<<" :\t";
  change=(Curve)CD!=C0;
  if(change&&verb)
    {
      cout<<"Working with minimal curve "<<(Curve)CD<<"\n";
      cout<<"\t[u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]\n";
    }
  return 1;
}
