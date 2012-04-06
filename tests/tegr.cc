// tegr.cc -- test for finding egr subgroup from a set of points
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
//

#include "matrix.h"
#include "cperiods.h"
#include "points.h"
#include "sieve_search.h"
#include "cperiods.h"
#include "elog.h"
#include "egr.h"
#include "htconst.h"

void codeletter(int i, char* code, int width=0);

int main()
{
  //  set_precision(string("Enter number of decimal places").c_str());
  initprimes(string("PRIMES").c_str(),0);
  int j, npts;

  long N, nclass, ncurve;
  char code[20];
  Curve E;

  while(1) {
    cin >> N; if(N==0) exit(0);
  cin >> nclass >> ncurve;
  cin >> E;
  Curvedata C(E);
  codeletter((nclass-1),code);
  cout<<endl;
  cout<<"==============================================================="<<endl;
  cout<<endl;
  cout << N<<code<<ncurve<<" = "<< E << endl;
  Point P(C);
  cin >> npts;
  vector<Point> points; points.reserve(npts);
  j=0; 
  while(j<npts)
    { 
      cin >> P;
      if ( P.isvalid() ) {points.push_back(P); j++;}
      else {cout<<"point "<<P<<" not on curve.\n\n"; }
    }
  cout<<npts<<" points entered:"<<points<<endl;
  //  bigfloat reg = regulator(points);
  //  cout<<"Regulator = "<<reg<<endl;

  bigint tam = Tamagawa_exponent(C);
  cout<<"Tamagawa exponent = "<<tam<<endl;

  vector<Point> egr_points=points;
  bigint egri = egr_index(egr_points);
  cout<<"Index of egr subgroup = "<<egri<<endl;
  }
}

void codeletter(int i, char* code, int width)
{
  int n=width;    // pads string to this width with blanks
  code[n]='\0';
  while (n) code[--n]=' ';

  int nc = i%26;
  char c = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[nc];
  n = 1 + (i-nc)/26;
  if(width==0) code[n]='\0';
  while (n) code[--n]=c;
}

//end of file tegr.cc
