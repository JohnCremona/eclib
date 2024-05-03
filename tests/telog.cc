// telog.cc -- to test elog.h/cc
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

#include <eclib/curve.h>
#include <eclib/points.h>
#include <eclib/cperiods.h>
#include <eclib/elog.h>

int test1(Curvedata& CD, Cperiods& per, const Point& P)
{
  bigcomplex z=elliptic_logarithm(CD,per,P);
  cout<<"Elliptic log of P is "<<z<<endl;
  Point Q=ellztopoint(CD,per,z,P.getZ());
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
      if((m*Q).is_zero()) cout<<Q<<" OK"<<endl;
      else  cout<<Q<<" wrong"<<endl;
    }
}

int main(){
#ifdef MPFP
  //  set_precision("Enter precision in bits");
  set_precision(175);
  long original_output_precision = RR::OutputPrecision();
  RR::SetOutputPrecision(original_output_precision-3);
#endif
  initprimes("PRIMES",0);

  // a curve with rank 3 and D>0:

  Curve c(0,0,1,-7,6);
  Curvedata cd(c,1);

  cout << "Testing some points:\n";
  Point P0(cd, bigint(0),bigint(2)) ;
  Point P1(cd, bigint(1),bigint(0)) ;
  Point P2(cd, bigint(2),bigint(0)) ;

  cout << "The points are P0 = " << P0 << 
    ", P1 = " << P1 << ", and P2 = " << P2 << endl ;

  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  if (!P1.isvalid()) cout << "P1 is not on the curve!\n";
  if (!P2.isvalid()) cout << "P2 is not on the curve!\n";

  /*
  // a curve with rank 3 and D<0:

  Curve c(0,0,1,-1,6);
  Curvedata cd(c,1);

  cout << "Testing some points:\n";
  Point P0(cd, bigint(0),bigint(2)) ;
  Point P1(cd, bigint(1),bigint(2)) ;
  Point P2(cd, bigint(12),bigint(41)) ;

  cout << "The points are P0 = " << P0 << 
    ", P1 = " << P1 << ", and P2 = " << P2 << endl ;

  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  if (!P1.isvalid()) cout << "P1 is not on the curve!\n";
  if (!P2.isvalid()) cout << "P2 is not on the curve!\n";
  */
  /*
  // a curve (7998K1) with rank 1 and 5-torsion:

  Curve c(1,0,0,5355560,7740216896);
  Curvedata cd(c,1);

  Point P0(cd, bigint(-248),bigint(80104)) ;

  cout << "The point is P0 = " << P0 << endl ;

  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  */

  cout<<"Curve "<<c<<endl;
  Cperiods cp(cd);
  //cout<<"Periods: "<<cp<<endl;
  cout<<endl;

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

  Curve c3(1,1,0,-202,1025);
  Curvedata cd3(c3,1);
  cout << "Curve "<<c3<<endl;

  Cperiods cp3(cd3);
  //cout << "Periods: "<<cp3<<endl;
  cout<<endl;

  Point Q(cd3,bigint(-8),bigint(51));
  cout << "The point Q is = " << Q << endl ;
  test1(cd3,cp3,Q);
  Point P3(cd3, bigint(8),bigint(-3)) ;
  cout << "The point P3 is = " << P3 << endl ;

  test1(cd3,cp3,P3);
  test2(cd3,cp3,P3,3);cout<<endl;

  cout << "================================ "<<endl;
  /*
  test3(cd,cp,2);cout<<endl;
  test3(cd,cp,3);cout<<endl;
  test3(cd,cp,5);cout<<endl;
  */

#ifdef MPFP
  RR::SetOutputPrecision(original_output_precision);
#endif

  return 0;
} //ends main

//ends file telog.cc
