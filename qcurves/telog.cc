// telog.cc -- to test elog.h/cc
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

#include "curve.h"
#include "points.h"
#include "cperiods.h"
#include "elog.h"

int test1(Curvedata& CD, Cperiods& per, const Point& P)
{
  bigcomplex z=elliptic_logarithm(CD,per,P);
  cout<<"Elliptic log of P is "<<z<<endl;
  Point Q=ellztopoint(CD,per,z,getZ(P));
  cout<<"Reconstructed P = "<<Q<<endl;
  return (P==Q);
}

void test2(Curvedata& CD, Cperiods& per, const Point& P, int m)
{
  vector<Point> Qlist = division_points(CD, per, m*P, m);
  cout<<"(m*"<<P<<")/m = "<<Qlist<<endl;
  if(Qlist.size()==0) return;
  cout<<"Checking..."<<endl;
  unsigned int i; for(i=0; i<Qlist.size(); i++) 
    {
      Point Q = Qlist[i];
      if(m*Q==m*P) cout<<Q<<" OK"<<endl;
      else  cout<<Q<<" wrong"<<endl;
    }
}

// torsion test:
void test3(Curvedata& CD, Cperiods& per, int m)
{
  cout<<"torsion point test (m="<<m<<")"<<endl;
  vector<Point> Qlist = torsion_points(CD, per, m);
  cout<<"m-torsion points: "<<Qlist<<endl;
  if(Qlist.size()==0) return;
  cout<<"Checking..."<<endl;
  unsigned int i; for(i=0; i<Qlist.size(); i++) 
    {
      Point Q = Qlist[i];
      if((m*Q).iszero()) cout<<Q<<" OK"<<endl;
      else  cout<<Q<<" wrong"<<endl;
    }
}

int main(){
  //  set_precision(string("Enter number of decimal places").c_str());
  set_precision(50);
  initprimes(string("PRIMES").c_str(),0);

  // a curve with rank 3 and D>0:

  Curve c(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-7),BIGINT(6));
  Curvedata cd(c,1);

  cout << "Testing some points:\n";
  Point P0(cd, BIGINT(0),BIGINT(2)) ;
  Point P1(cd, BIGINT(1),BIGINT(0)) ;
  Point P2(cd, BIGINT(2),BIGINT(0)) ;

  cout << "The points are P0 = " << P0 << 
    ", P1 = " << P1 << ", and P2 = " << P2 << endl ;

  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  if (!P1.isvalid()) cout << "P1 is not on the curve!\n";
  if (!P2.isvalid()) cout << "P2 is not on the curve!\n";

  /*
  // a curve with rank 3 and D<0:

  Curve c(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-1),BIGINT(6));
  Curvedata cd(c,1);

  cout << "Testing some points:\n";
  Point P0(cd, BIGINT(0),BIGINT(2)) ;
  Point P1(cd, BIGINT(1),BIGINT(2)) ;
  Point P2(cd, BIGINT(12),BIGINT(41)) ;

  cout << "The points are P0 = " << P0 << 
    ", P1 = " << P1 << ", and P2 = " << P2 << endl ;

  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  if (!P1.isvalid()) cout << "P1 is not on the curve!\n";
  if (!P2.isvalid()) cout << "P2 is not on the curve!\n";
  */
  /*
  // a curve (7998K1) with rank 1 and 5-torsion:

  Curve c(BIGINT(1),BIGINT(0),BIGINT(0),BIGINT(5355560),BIGINT(7740216896));
  Curvedata cd(c,1);

  Point P0(cd, BIGINT(-248),BIGINT(80104)) ;

  cout << "The point is P0 = " << P0 << endl ;

  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  */

  cout<<"Curve "<<c<<endl;
  Cperiods cp(cd);
  cout<<"Periods: "<<cp<<endl<<endl;

  if(test1(cd,cp,P0)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;

  if(test1(cd,cp,P1)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;
  if(test1(cd,cp,P2)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;

  if(test1(cd,cp,2*P0)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;
  if(test1(cd,cp,2*P1)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;
  if(test1(cd,cp,2*P2)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;

  if(test1(cd,cp,P0+P1)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;
  if(test1(cd,cp,P1-P2)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;
  if(test1(cd,cp,P2-P0)) cout<<"OK!\n"<<endl;
  else cout<<"WRONG!\n"<<endl;


  test2(cd,cp,P0,2);cout<<endl;
  test2(cd,cp,P0,3);cout<<endl;
  test2(cd,cp,P0,5);cout<<endl;

  cout << "================================ "<<endl;

  Curve c3(BIGINT(1),BIGINT(1),BIGINT(0),BIGINT(-202),BIGINT(1025));
  Curvedata cd3(c3,1);
  cout << "Curve "<<c3<<endl;

  Cperiods cp3(cd3);
  cout << "Periods: "<<cp3<<endl<<endl;

  Point Q(cd3,BIGINT(-8),BIGINT(51));
  cout << "The point Q is = " << Q << endl ;
  test1(cd3,cp3,Q);
  Point P3(cd3, BIGINT(8),BIGINT(-3)) ;
  cout << "The point P3 is = " << P3 << endl ;

  test1(cd3,cp3,P3);
  test2(cd3,cp3,P3,3);cout<<endl;

  cout << "================================ "<<endl;
  /*
  test3(cd,cp,2);cout<<endl;
  test3(cd,cp,3);cout<<endl;
  test3(cd,cp,5);cout<<endl;
  */
  return 0;
} //ends main

//ends file telog.cc
