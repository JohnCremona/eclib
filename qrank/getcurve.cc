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

int getcurve(Curvedata& CD, int verb)
{
  Curve C0;
  if(verb) cout  << "Enter curve: ";
  cin>>ws;  if(cin.eof()) return 0; // quit if EOF reached
  cin >> C0;
  if (verb) cout << endl;
  if(C0.isnull()) return 0;  // quit if null curve entered
  cout << "Curve "<< C0<<" :\t";
  CD = Curvedata(C0,0);      // DON'T change coords
  if(CD.isnull()) // input curve was singular, non-null
    {
      cout<<" singular"<<endl;
      return 0;
    }
  return 1;
}
