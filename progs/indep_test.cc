// indep_test.cc: program to test input points for (in)dependence 
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
 
// NB This is based on 2-saturation, trying to show points are
// independent in E(Q)/2E(Q), and is now largely obsolete, being
// superceded by general saturation algorithms

#include <eclib/points.h>
#include <eclib/sifter.h>

int main()
{
  set_precision(100);
  initprimes("PRIMES",0);

  int verbose = 1;
  cout<<"verbose (0/1)? ";
  cin>>ws;  if(cin.eof()) {cout<<endl; exit(0);}
  cin >>verbose;

    while (1)
    {
      cout<<"\nInput a curve: ";
      cin>>ws;  if(cin.eof()) {cout<<endl; exit(0);}
      Curve E;
      cin >> E;
      if ( E.isnull() ) exit(0);
      Curvedata C(E, 0);
      cout << "Curve " << (Curve)C << endl;
      Point P(C);
      cout<<"enter number of points: ";
      cin>>ws;  if(cin.eof()) {cout<<endl; exit(0);}
      int npts;
      cin >> npts;
      vector<Point> points; points.reserve(npts);
      int j=0;
      while(j<npts)
	{
	  cout<<"\n  enter point "<<(j+1)<<" : ";
	  cin>>ws;  if(cin.eof()) {cout<<endl; exit(0);}
	  cin >> P;
	  if ( P.isvalid() ) {points.push_back(P); j++;}
	  else {cout<<"point "<<P<<" not on curve.\n\n"; }
	}
      if(verbose) cout<<npts<<" points entered.\n";

      long naux=npts+10;
      cout << "Enter number of primes to use: "; 
      cin>>ws;  if(cin.eof()) {cerr<<endl; exit(0);}
      cin>>naux;

      sifter box(&C, naux, verbose);
      box.process(points);
      long rank = box.getrank();
      if(rank==npts)
	cout<<"Points are all independent, their rank is "<<rank<<endl;
      else
	cout<<"Points may be dependent, rank is at least "<<rank<<endl;
    }
}


//end of file indep.cc





