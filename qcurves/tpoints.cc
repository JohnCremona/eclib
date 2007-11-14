// tpoints.cc -- to test points.h/cc
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

#include "points.h"

int main(){
  //  set_precision("Enter number of decimal places");
  set_precision(30);
  initprimes("PRIMES",1);

  Curve c(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-7),BIGINT(6));
  Curvedata cd(c,1);

  cout << "Testing some points:\n";
  Point P0(cd, BIGINT(0),BIGINT(2)) ;
  Point P1(cd, BIGINT(1),BIGINT(0)) ;
  Point P2(cd, BIGINT(2),BIGINT(0)) ;
  cout << "The points are P0 = " << P0 << 
    ", P1 = " << P1 << ", and P2 = " << P2 << endl ;
  cout << "Now in Pari format:\n";
  cout << "The points are P0 = "; output_pari(cout,P0); 
  cout << ", P1 = ";              output_pari(cout,P1); 
  cout << ", and P2 = ";          output_pari(cout,P2); cout << endl ;
  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  if (!P1.isvalid()) cout << "P1 is not on the curve!\n";
  if (!P2.isvalid()) cout << "P2 is not on the curve!\n";
  cout << "Their negatives are -P0 = " << -P0 << 
    ", -P1 = " << -P1 << ", and -P2 = " << -P2 << endl ;
        
  cout << "Computing their heights:\n";
  bigfloat ht0 = height(P0) ;
  bigfloat ht1 = height(P1) ;
  bigfloat ht2 = height(P2) ;

  cout << "Heights are " << ht0 << ", " << ht1 << ", and " << ht2 << endl ;
        
  Point origin(cd);
  cout << "The origin is " << origin << endl;
  cout << "Now some additions etc,:\n";
  Point sum = P0 + P1 ;
  cout << "P0 + P1 = " << sum << endl ;
  sum = P0 - P1 ;
  cout << "P0 - P1 = " << sum << endl ;
  sum -= P2 ;
  cout << "P0 - P1 - P2 = " << sum << endl ;
  sum = P0.twice();
  cout << "P0.twice() = " << sum << endl;
  sum = P0 + P0;
  cout << "P0 + P0 = " << sum << endl;
  sum = 3*P0;
  cout << "3*P0 = " << sum << endl;
  sum = P0 - P0;
  cout << "P0 - P0 = " << sum << endl;
  sum = P0 + 3 * P1 - P2 ;
  cout << "P0 +3 P1 - P2 = " << sum << endl ;
  sum = 2*P0 + 2* P1 + P2 ;
  cout << "2P0 +2 P1 + P2 = " << sum << endl ;

  cout << "Now we try a systematic exploration" << endl ;
  int loop=1;
  //  cout << "input 0 to skip: ";  cin >> loop;
  if (loop)
    {
      Point Q0,Q1,Q2; int i0,i1,i2;
      for(i0 = -3, Q0=-3*(P0+P1+P2) ; i0 < 4 ; i0++, Q0=Q0+P0){
	for(i1 = -3, Q1=Q0 ; i1 < 4 ; i1++, Q1=Q1+P1){
	  for(i2 = -3, Q2=Q1 ; i2 < 4 ; i2++, Q2=Q2+P2){
	    cout << i0 << ", " << i1 <<", " << i2<< ": " << Q2 << endl ;
	  }
	}
      }
    }

  sum = P0 -P1 -P2 ;
  cout << "P0 -P1 -P2 = " << sum << endl ;
  cout << "Height of     P0-P1-P2 = " << flush;
  bigfloat htsum = height(sum) ;
  cout << htsum << endl ;
        
  cout <<"3 (P0 -P1 -P2) = " << flush;  
  Point triplesum = 3 * sum ;
  cout << triplesum << endl ;

  cout << "Height of 3*(P0-P1-P2) = " << flush;
  bigfloat ht3sum = height(triplesum) ;
  cout << ht3sum << endl ;
        
  cout << "The quotient is " << ht3sum/htsum << endl ;
        
  /*
  vector<Point> pointlist(3); // won't initialize
  pointlist[0] = P0; pointlist[1] = P1; pointlist[2] = P2 ;
  bigfloat reg = regulator(pointlist) ;

  cout << "The regulator of P0, P1, P2 is " << reg << endl ;      
  */
  return 0;
} //ends main

//ends file tpoints.cc
