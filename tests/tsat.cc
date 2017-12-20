// tsat.cc -- test for saturate.h/cc
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
#include <eclib/matrix.h>
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


int randint(int top)
{
  int ans=1+(int) (double(top)*rand()/(RAND_MAX+1.0));
  ans-=(top/2);
  return ans;
}

int main()
{
  //  set_precision("Enter precision in bits");
  set_precision(600);
  initprimes("PRIMES",0);
  int verbose = 1;
  cout<<"verbose (0/1)? ";             cin >>verbose;
  int j, npts;

  /* test curve
  Curve E(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-7),BIGINT(6));

  Curvedata C(E);
  cout << "Curve " << E << endl;
  
  vector<Point> PP;
  Point P1(C,2,0,1);  
  Point P2(C,-1,3,1); 
  Point P3(C,4,6,1);  
  PP.push_back(P1); 
  PP.push_back(3*P2+P1); 
  PP.push_back(P3);
  cout << "Points " << PP << endl;
  */

  Curve E;
  cout<<"\nInput a curve: ";      cin >> E;
  Curvedata C(E);
  cout << "Curve " << E << endl;
  saturator sieve(&C,verbose);

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
  cout<<npts<<" points entered.\n";

  int pmax;

  cout<<"prime p to saturate at? ", cin>>pmax;
  pmax = NextPrime(pmax);
  cout<<"\nSaturating at prime "<<pmax<<endl;

  //  if(npts==1) points[0]=randint(10)*points[0];
  /*
  if(npts>1) 
    {
      int a=randint(10),b=randint(10),c=randint(10),d=randint(10);
      cout<<"a,b,c,d="<<a<<","<<b<<","<<c<<","<<d<<endl;
      cout<<"det = "<<(a*d-b*c)<<endl;
      Point Q = a*points[0]+b*points[1];
      points[0]=c*points[0]+d*points[1];
      points[1]=Q;
    }
  */
  sieve.set_points(points);
  cout<<"Original generators:\n"<<points<<endl;
  //  bigfloat reg = regulator(points);
  //  cout<<"Regulator = "<<reg<<endl;

  /*
  int index = sieve.do_saturation_upto(pmax);
  cout<<"Finished p-saturation for p up to "<<pmax;
  */
  int log_index = sieve.do_saturation(pmax);
  cout<<"Finished p-saturation for p =  "<<pmax;
  if(log_index>0) 
    {
      cout<<", index gain = "<<pmax<<"^"<<log_index<<endl;
      vector<Point> newpoints = sieve.getgens();
      cout<<"New generators:\n"<<newpoints<<endl;
    //  bigfloat newreg = regulator(newpoints);
    //  cout<<"New regulator = "<<newreg<<endl;
    //  cout<<"Ratio = "<<(reg/newreg)<<endl;
    }
  else
    {
      cout<<", points were saturated"<<endl;
    }
}

//end of file tsat.cc





