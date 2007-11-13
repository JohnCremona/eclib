// tsat.cc -- test for saturate.h/cc
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

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)
#define NextPrime next_prime
#endif


int randint(int top)
{
  int ans=1+(int) (double(top)*rand()/(RAND_MAX+1.0));
  ans-=(top/2);
  return ans;
}

int main()
{
  //  set_precision("Enter number of decimal places");
  set_precision(200);
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
  int index = sieve.do_saturation(pmax);
  cout<<"Finished p-saturation for p =  "<<pmax;
  if(index>1) 
    {
      cout<<", index gain = "<<index<<endl;
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





