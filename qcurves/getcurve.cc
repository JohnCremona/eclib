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
  CD = Curvedata(C0,0);      // DON'T change coords
  if(CD.isnull()) // input curve was singular, non-null
    {
      cout<<C0<<" is singular"<<endl;
      return 0;
    }
  return 1;
}

//#define DEBUG_Q_INPUT

int getcurve(vector<bigrational>& ai, int verb)
{
  // read the coefficients, either as "a1 a2 a3 a4 a6" or as
  // "[a1,a2,a3,a4,a6]" using code essentially the same as in
  // ../qcurves/curve.cc
  ai.resize(5);
  if(verb) cerr  << "Enter curve: ";
  cin>>ws;  if(cin.eof()) return 0; // quit if EOF reached
  char c;
  cin.get(c);
#ifdef DEBUG_Q_INPUT
  cout<<"First char read = "<<c<<"\n";
#endif
  if(c=='[')
    {
#ifdef DEBUG_Q_INPUT
      cout<<"Reading [a1,a2,a3,a4,a6]...\n";
#endif
      cin >> ai[0] >> c;
      if(c!=',')
	{
	  cout << "syntax error on curve input" << endl;
	  abort();
	}
      cin >> ai[1] >> c; 
      if(c!=',')
	{
	  cout << "syntax error on curve input" << endl;
	  abort();
	}
      cin >> ai[2] >> c; 
      if(c!=',')
	{
	  cout << "syntax error on curve input" << endl;
	  abort();
	}
      cin >> ai[3] >> c; 
      if(c!=',')
	{
	  cout << "syntax error on curve input" << endl;
	  abort();
	}
      cin >> ai[4] >> c; 
      if(c!=']')
	{
	  cout << "syntax error on curve input" << endl;
	  abort();
	}
#ifdef DEBUG_Q_INPUT
      cout<<ai<<endl;
#endif
    }
  else
    {
#ifdef DEBUG_Q_INPUT
      cout<<"Reading a1 a2 a3 a4 a6 ...\n";
#endif
      cin.unget();
      cin >> ai[0] >> ai[1] >> ai[2] >> ai[3] >> ai[4];
#ifdef DEBUG_Q_INPUT
      cout<<ai<<endl;
#endif
    }
  // test for null curve input
  if((num(ai[0])==0)&&(num(ai[1])==0)&&(num(ai[2])==0)&&(num(ai[3])==0)&&(num(ai[4])==0))  return 0;  // quit if null curve entered
  return 1;
}
