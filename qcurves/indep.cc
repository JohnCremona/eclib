// indep.cc: program to test input points for (in)dependence 
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
 
// NB This is based on 2-saturation, trying to show points are
// independent in E(Q)/2E(Q), and is now largely obsolete, being
// superceded by general saturation algorithms

#include "points.h"
#include "sifter.h"

int main()
{
  set_precision(30);
  initprimes("PRIMES",0);

  int verbose = 1;
  long rank, npts, j;
  cout<<"verbose (0/1)? ";             cin >>verbose;
  Curve E;

    while (1)
    {
      cout<<"\nInput a curve: ";      cin >> E;
      if ( E.isnull() ) break;
      Curvedata C(E);
      cout << "Curve " << (Curve)C << endl;
      Point P(C);
      cout<<"enter number of points: ";      cin >> npts;
      vector<Point> points; points.reserve(npts);
      j=0; 
      while(j<npts)
	{ 
	  cout<<"\n  enter point "<<(j+1)<<" : ";
	  cin >> P;
	  if ( P.isvalid() ) {points.push_back(P); j++;}
	  else {cout<<"point "<<P<<" not on curve.\n\n"; }
	}
      if(verbose) cout<<npts<<" points entered.\n";

      long naux=npts+10;
      cout << "Enter number of primes to use: "; cin>>naux;
      
      sifter box(&C, naux, verbose);
      box.process(points);
      rank = box.getrank();
      if(rank==npts)
	cout<<"Points are all independent, their rank is "<<rank<<endl;
      else
	cout<<"Points may be dependent, rank is at least "<<rank<<endl;
    }
}


//end of file indep.cc





