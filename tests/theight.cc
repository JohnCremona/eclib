// THEIGHT.CC -- test of height procedures
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

#include <eclib/points.h>

int main(){
#ifdef MPFP
  set_precision(350);
#endif
  initprimes("PRIMES",1);
  Curve E;
  bigint u, r, s, t;

  while (cerr<<"Input a curve: ", cin>>E, !E.isnull())
  {
    Curvedata C(E, 0);
    cout<<"Curve "<< E <<endl;
    vector<bigint> badp = getbad_primes(C);

    Curvedata Emin = C.minimalize(u,r,s,t);
    int is_min = is_one(abs(u)); // then E was already minimal, no adjustments needed
    if(!is_min)
      {
        cout<<"Input curve is not minimal at "<<pdivs(u)<<endl;
        cout<<"Heights will be computed after mapping points to the minimal model "<<(Curve)Emin<<endl;
      }
    Point P(C);

    while
      (
       cerr<<"Input a point on curve, [0] to finish:\n",
       cin >> P,
       !(!cin)&& P.isvalid()
       )
        {
          cout << "Point " << P;
          int ord = order(P);
          if(ord>0)cout<< " has order " << ord;
          else     cout<< " has infinite order";
          cout << endl;

          Point Pmin = (is_min? P: transform(P, &Emin, u, r, s, t));

          cout << "Local heights:\n";
          bigfloat gh = to_bigfloat(0);
	  bigint d = gcd(Pmin.getZ(),Pmin.getX());
          vector<bigint> pdivsz=pdivs(d);
	  for ( const auto& q : pdivsz)
            {
	      cout << q << ":\t\t"  << flush;
	      bigfloat ph = pheight(Pmin,q);
	      gh+=ph;
	      cout << ph << endl;
            }
	  cout << "Sum so far =\t" << gh << endl;
	  bigfloat lxd = 2*log(I2bigfloat(d));
	  cout << "log(den(x(P))) = " << lxd << endl;
	  for ( const auto& pr : badp)
            {
	      if(div(pr,d)) continue;
              cout << pr << ":\t\t"  << flush;
              bigfloat  ph = pheight(Pmin,pr);
              gh+=ph;
              cout << ph << endl;
            }
	  cout << "Sum so far =\t" << gh << endl;
          bigfloat rh = realheight(Pmin);
          cout << "R:\t\t" << rh << endl;
	  gh += rh;
          cout << "\nSum of local heights: " << gh << endl;

          //
          // N.B. The call to height() calls realheight() and pheight() again,
          //      since as height() was not called earlier, P's height field was 
          //      not set.
          //
          // N.B. Computing the global height automatically handles
          // the move to a minimal model, so height(P) and
          // height(Pmin) are the same
          cout<<"global height of "<<P<<" is "<<height(P)<<endl;
          //          cout<<"global height of "<<Pmin<<" is "<<height(Pmin)<<endl;

        } // end point loop
  }       // end curve loop
  //  exit(0);
// Some tests borrowed from OM's test program:
  E = Curve(0, 0, 1, -7, 6);
  Curvedata C(E, 0);
  cout << "\n\nTesting points on the curve " << E << endl;
  Point P0(C, bigint(0),bigint(2));
  Point P1(C, bigint(1),bigint(0));
  Point P2(C, bigint(2),bigint(0));
  cout << "The points are P0 = " << P0 << 
    ", P1 = " << P1 << ", and P2 = " << P2 << endl;
  if (!P0.isvalid()) cout << "P0 is not on the curve!\n";
  if (!P1.isvalid()) cout << "P1 is not on the curve!\n";
  if (!P2.isvalid()) cout << "P2 is not on the curve!\n";
  cout << "Their negatives are -P0 = " << -P0 << 
    ", -P1 = " << -P1 << ", and -P2 = " << -P2 << endl;

  cout << "Computing their heights:\n";
  bigfloat ht0 = height(P0);
  bigfloat ht1 = height(P1);
  bigfloat ht2 = height(P2);

  cout << "Heights are " << ht0 << ", " << ht1 << ", and " << ht2 << endl;
        
  Point origin(C);
  cout << "The origin is " << origin << endl;
  cout << "Now some additions etc,:\n";
  Point sum = P0 + P1;
  cout << "P0 + P1 = " << sum << endl;
  sum = P0 - P1;
  cout << "P0 - P1 = " << sum << endl;
  sum -= P2;
  cout << "P0 - P1 - P2 = " << sum << endl;
  sum = P0.twice();
  cout << "P0.twice() = " << sum << endl;
  sum = P0 + P0;
  cout << "P0 + P0 = " << sum << endl;
  sum = 3*P0;
  cout << "3*P0 = " << sum << endl;
  sum = P0 - P0;
  cout << "P0 - P0 = " << sum << endl;
  sum = P0 + 3 * P1 - P2;
  cout << "P0 +3 P1 - P2 = " << sum << endl;
  sum = 2*P0 + 2* P1 + P2;
  cout << "2P0 +2 P1 + P2 = " << sum << endl;
  /*
  cout << "Now we try a systematic exploration" << endl;
  Point Q; int i0,i1,i2;
  int k=3;
  for(i0 = -k; i0 <= k; i0++){
    for(i1 = -k; i1 <= k; i1++){
      for(i2 = -k; i2 <= k; i2++){
	Q = i0*P0+i1*P1+i2*P2;
        cout << i0 << ", " << i1 <<", " << i2<< ": " << Q << endl;
      }
    }
  }
  */
  sum = P0 -P1 -P2;
  bigfloat htsum = height(sum);
  cout << "P0 -P1 -P2 = " << sum << "\n\t with height " << htsum<<endl;

  Point doublesum = 2 * sum;
  bigfloat ht2sum = height(doublesum);
  cout << "2 (P0 -P1 -P2) = " << doublesum << "\n\t with height " << ht2sum << endl;
  cout << "The quotient is " << ht2sum/htsum << endl;

  Point triplesum = 3 * sum;
  bigfloat ht3sum = height(triplesum);
  cout << "3 (P0 -P1 -P2) = " << triplesum << "\n\t with height " << ht3sum << endl;
  cout << "The quotient is " << ht3sum/htsum << endl;
        
  vector<Point> pointlist(3);
  pointlist[0] = P0; pointlist[1] = P1; pointlist[2] = P2;
//cout << "Making a PointArray out of P0, P1, P2.  It is " << pointlist << endl;
//cout << "Calling regulator" << endl;
  bigfloat reg = regulator(pointlist);
//cout << "Back from calling regulator" << endl;
  cout << "The regulator of P0, P1, P2 is " << reg << endl;      
  return 0;
  }         // end main()

