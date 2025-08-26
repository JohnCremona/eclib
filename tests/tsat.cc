// tsat.cc -- test for saturate.h/cc
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
//

#include <eclib/interface.h>
#include <eclib/curve.h>
#include <eclib/getcurve.h>
#include <eclib/points.h>
#include <eclib/cperiods.h>
#include <eclib/polys.h>
#include <eclib/curvemod.h>
#include <eclib/pointsmod.h>
#include <eclib/saturate.h>
#include <eclib/elog.h>
#include <eclib/sieve_search.h>
#include <eclib/mwprocs.h>


int randint(int top)
{
  int ans=1+(int) (double(top)*rand()/(RAND_MAX+1.0));
  ans-=(top/2);
  return ans;
}

int main()
{
  //  set_precision("Enter precision in bits");
  set_precision(100);
  initprimes("PRIMES",0);
  int verbose = 1;
  cerr<<"verbose (0/1)? ";             cin >>verbose;

  // Curve E;
  // cout<<"\nInput a curve: ";      cin >> E;
  // Curvedata C(E);
  // cout << "Curve " << E << endl;

  Curvedata C;
  while (getcurve(C, verbose))
    {

      cout << "======================================================\n\n";
      cout << "E = " << (Curve)C <<endl;
      saturator sieve(&C,1, verbose);

  int j=0, npts;
  Point P(C);
  cerr<<"enter number of points: ";      cin >> npts;
  vector<Point> points; points.reserve(npts);

  while(j<npts)
    { 
      cerr<<"\n  enter point "<<(j+1)<<" : ";
      cin >> P;
      if ( P.isvalid() ) {points.push_back(P); j++;}
      else {cerr<<"point "<<P<<" not on curve.\n\n"; }
    }
  cout<<npts<<" points entered.\n";

  int pmin, pmax;

  cerr<<"minimum prime p to saturate at (or 0)? ", cin>>pmin;
  cerr<<"maximum prime p to saturate at (or -1 for automatic)? ", cin>>pmax;
  if (pmax>0)
    {
      pmax = NextPrime(pmax);
      cout<<"\nSaturating at primes from " <<pmin <<" up to "<<pmax<<endl;
    }
  else
    {
      cout<<"\nSaturating at all primes";
      if (pmin>2) cout << " from " <<pmin;
      cout<<endl;
    }

  sieve.set_points(points);
  cout<<"Original generators:\n"<<points<<endl;
  //  bigfloat reg = regulator(points);
  //  cout<<"Regulator = "<<reg<<endl;

  long index;
  vector<long> unsat;
  int ok = sieve.saturate(unsat, index, pmax, pmin, 10);

  cout<<"Finished saturation" << endl;
  if (ok || pmax>0)
    {
      cout << "Saturation was successful: ";
        if(index>1)
          {
            cout<<" index gain = "<<index<<"."<<endl;
            vector<Point> newpoints = sieve.getgens();
            cout<<"New generators:\n"<<newpoints<<endl;
            //  bigfloat newreg = regulator(newpoints);
            //  cout<<"New regulator = "<<newreg<<endl;
            //  cout<<"Ratio = "<<(reg/newreg)<<endl;
          }
        else
          {
            cout<<"points were saturated."<<endl;
          }
    }
  else
    {
      cout << "Saturation not successful: ";
      if (unsat.size()==0)
        {
          cout << " the index bound was too large.\n";
        }
      else
        {
          cout << " for the following primes p, we both failed to prove p-saturation ";
          cout << "and failed to enlarge by index p" << endl;
          cout << unsat << endl;
        }
    }
  if (verbose)
    sieve.show_q_tally();
    }
}

//end of file tsat.cc





