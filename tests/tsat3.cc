// tsat3.cc -- test for saturate.h/cc reading from gens files directly
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

#include <eclib/interface.h>
#include <eclib/method.h>
#include <eclib/curve.h>
#include <eclib/points.h>
#include <eclib/cperiods.h>
#include <eclib/polys.h>
#include <eclib/curvemod.h>
#include <eclib/pointsmod.h>
#include <eclib/saturate.h>
#include <eclib/elog.h>
#include <eclib/sieve_search.h>
#include <eclib/mwprocs.h>
#include <eclib/curvesort.h>

#define USE_EGR
#define SAT_BND -1 // saturation bound, or 1 for automatic

#ifdef USE_EGR
const int use_egr=1;
#else
const int use_egr=0;
#endif

#define INPUT_CLASS_IS_LETTER // we only use letters now!

// Utility function for parsing input of lists of integers such as [], [2], [2,2] 
vector<int> input_list(istream & is);

int main()
{
#ifdef MPFP
  //  set_precision("Enter precision in bits");
  set_precision(100);
#endif
  initprimes("PRIMES",0);
  int verbose = 0;
  cerr<<"verbose (0/1)? ";             cin >>verbose;
  int j, npts;

  long N, ncurve, nclass;
  string code;
  Curve E;

  while(!(feof(stdin))) {

  // Input the curve's ID and the curve:

  cin >> N; if(N==0) break;
#ifdef INPUT_CLASS_IS_LETTER
  cin >> code;
#else
  cin >> nclass; 
  code = codeletter(nclass-1);
#endif
  cin >> ncurve;
  cin >> E;
  Curvedata C(E, 0);
  cout<<endl;
  cout << N<<code<<ncurve<<" = "<< E << endl;

  // Input the rank:

  Point P(C);
  cin >> npts;
  cout<<"rank =  "<<npts<<endl;

  // Input the torsion structure:
  vector<int> torsion_group = input_list(cin);
  int trank = torsion_group.size();
  cout<<"torsion group structure =  "<<torsion_group<<" (torsion rank "<<trank<<")"<<endl;

  // Input the points of infinite order (if any):

  vector<Point> points; points.reserve(npts);
  j=0; 
  while(j<npts)
    { 
      cin >> P;
      if ( !P.isvalid() ) 
	{
	  cout<<"point "<<P<<" not on curve.\n"<<endl; 
	}
      points.push_back(P); 
      j++;
    }
  cout<<"Input non-torsion points: "<<points<<endl;

  // Input the points of finite order (if any):

  vector<Point> tpoints; points.reserve(trank);
  j=0; 
  while(j<trank)
    { 
      cin >> P;
      if ( !P.isvalid() ) 
	{
	  cout<<"point "<<P<<" not on curve.\n"<<endl; 
	}
      tpoints.push_back(P); 
      j++;
    }
  cout<<"Input torsion points: "<<tpoints<<endl;

  long index = 1;
  vector<long> unsatprimes;
  int success = 1;

  if (npts)
    {
      success = saturate_points(C, points, index, unsatprimes, SAT_BND, 2, use_egr, (verbose));
    }

  if(success)
    {
      cout<<"Saturation complete --";
      if(index==1) 
	cout<<"input points were saturated"<<endl;
      else
	{
	  cout<<"input points had index "<<index
	      <<" in their saturation."<<endl;
	  cout<<"Basis for saturation:\t"<<points<<endl;
	}
    }
  else
    cout<<"Saturation failed at "<<unsatprimes<<endl;

  for (int i=0; i<npts; i++)
    { 
      Point P = points[i];
      cout << "Generator "<<(i+1)<<" is "<<P<<"; ";
      cout << "height "<<height(P);
      if(!P.isvalid()) cout<<" --warning: NOT on curve!";
      cout<<endl;
    }
  cout<<endl;
  cout << "Regulator = "<<regulator(points)<<endl<<endl;

  // Finally output a line similar to the input line:

#ifdef INPUT_CLASS_IS_LETTER
  cout<<N<<"\t"<<code<<"\t"<<ncurve<<"\t"<<E<<"\t"<<npts;
#else
  cout<<N<<"\t"<<nclass<<"\t"<<ncurve<<"\t"<<E<<"\t"<<npts;
#endif
  for(j=0; j<npts; j++) cout<<"\t"<<points[j];
  cout<<endl;
  cout<<"=============================================================="<<endl;
  }
}

// Utility function for parsing input of lists of integers such as [], [2], [2,2] 
vector<int> input_list(istream & is)
{
  char c; int a; vector<int> ai;
  is>>c;  // swallow first [
  is>>ws>>c;
  if (c==']') return ai;
  is.unget();
  while (c!=']')
    {
      is >> a >> c; // c is a comma or ]
      ai.push_back(a);
   }
  return ai; 
}

//end of file tsat3.cc

