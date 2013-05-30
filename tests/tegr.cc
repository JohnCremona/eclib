// tegr.cc -- test for finding egr subgroup from a set of points
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

#include <eclib/matrix.h>
#include <eclib/cperiods.h>
#include <eclib/points.h>
#include <eclib/sieve_search.h>
#include <eclib/cperiods.h>
#include <eclib/elog.h>
#include <eclib/egr.h>
#include <eclib/htconst.h>
#include <eclib/curvesort.h>

int main()
{
  //  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);
  int j, npts;

  long N, nclass, ncurve;
  Curve E;

  while(1) {
    cin >> N; if(N==0) exit(0);
  cin >> nclass >> ncurve;
  cin >> E;
  Curvedata C(E);
  cout<<endl;
  cout<<"==============================================================="<<endl;
  cout<<endl;
  cout << N<<codeletter(nclass-1)<<ncurve<<" = "<< E << endl;
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

//end of file tegr.cc
