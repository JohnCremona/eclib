// tsat3.cc -- test for saturate.h/cc reading from gens files directly
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

#include "interface.h"
#include "matrix.h"
#include "curve.h"
#include "points.h"
#include "cperiods.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "saturate.h"
#include "elog.h"
#include "sieve_search.h"
#include "mwprocs.h"
#include "curvesort.cc"

#define USE_EGR
#define SAT_BND 100000 // saturation bound:  use -1 for global default

#ifdef USE_EGR
const int use_egr=1;
#else
const int use_egr=0;
#endif

#define INPUT_CLASS_IS_LETTER // we only use letters now!

//void codeletter(int i, char* code, int width=0);
// Utility function for parsing input of lists of integers such as [], [2], [2,2] 
vector<int> input_list(istream & is);

int main()
{
  //  set_precision(string("Enter number of decimal places").c_str());
  set_precision(100);
  initprimes(string("PRIMES").c_str(),0);
  int verbose = 0;
  //  cout<<"verbose (0/1)? ";             cin >>verbose;
  int j, npts;

  long N, ncurve;
  char code[20];
  Curve E;

  while(!(feof(stdin))) {

  // Input the curve's ID and the curve:

  cin >> N; if(N==0) break;
#ifdef INPUT_CLASS_IS_LETTER
  cin >> code;
#else
  long nclass;
  cin >> nclass; 
  codeletter((nclass-1),code);
#endif
  cin >> ncurve;
  cin >> E;
  Curvedata C(E);
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
	  cout<<"point "<<P<<" not on curve.\n\n"; 
	  abort();
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
	  cout<<"point "<<P<<" not on curve.\n\n"; 
	  abort();
	}
      tpoints.push_back(P); 
      j++;
    }
  cout<<"Input torsion points: "<<tpoints<<endl;

  bigint index;
  vector<long> unsatprimes;
  int success = saturate_points(C, points, index, unsatprimes, SAT_BND, use_egr, (verbose));

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




#if(0)
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
#endif

//end of file tsat3.cc

