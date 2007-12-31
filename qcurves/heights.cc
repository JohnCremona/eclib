// heights.cc: implementation of height functions declared in points.h
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
 
#include "points.h"  // which includes curve.h

bigfloat height(Point& P)
{
  // WARNING -- no check made of validity of point on curve
  bigfloat zero(to_bigfloat(0));
  if ( P.iszero() ) return zero;
  if (P.height >= zero) return P.height;  // already calculated it
  if (order(P) > 0) {P.height = zero; return zero; } // zero height if torsion
// N.B. So if we ever ask a point its height it will compute its order.
// otherwise need to calculate it

// Add local heights at finite primes dividing discr(E) OR denom(P).
// The components for primes dividing denom(P) add to log(denom(x(P)));
//   since P=(XZ:Y:Z^3), denom(P)=Z=gcd(XZ,Z^3), called "zroot" here,
//   and so the contribution is log(denom(x(P))) = 2*log(zroot).
//   This avoids factorizing the denominator.

  const bigint& zroot = gcd(getX(P),getZ(P));   // = cube root of Z
  bigfloat ans = realheight(P);
  ans += 2*log(I2bigfloat(zroot));

  const vector<bigint>& bad_p = getbad_primes( *(P.E) );
  vector<bigint>::const_iterator pr = bad_p.begin();
  while(pr!=bad_p.end())
    {
      const bigint& p = *pr++;
      if(ndiv(p,zroot)) ans += pheight(P,p);
    }
  P.height = ans;
  return ans;
}

//#define DEBUG_HEIGHT

bigfloat pheight(const Point& P, const bigint& pr)
{
#ifdef DEBUG_HEIGHT
cout<<"In pheight with P = "<<P<<" and pr = "<<pr<<endl;
cout<<I2double(pr)<<"\n";
cout<<"(as a bigfloat, pr = "<<I2bigfloat(pr)<<")"<<endl;
#endif
  bigint a1,a2,a3,a4,a6,b2,b4,b6,b8,c4,c6,discr;
  P.E->getai(a1,a2,a3,a4,a6);
  P.E->getbi(b2,b4,b6,b8);
  P.E->getci(c4,c6);
  discr = getdiscr(*(P.E));
  long n = val(pr, discr);
#ifdef DEBUG_HEIGHT
cout<<"n = val(pr, discr) = " << n << endl;
#endif
  bigint x,y,z;
  P.getcoordinates(x,y,z);
  const bigint& zroot = gcd(x,z); // = cube root of z
  long vpz = 3*val(pr,zroot);
#ifdef DEBUG_HEIGHT
cout<<"vpz = val(pr, z) = " << vpz << endl;
#endif
  const bigint& x2 = x*x;
  const bigint& z2 = z*z;
  const bigint& xz = x*z;
  const bigint& yz = y*z;
  long a = val(pr, 3*x2 + 2*a2*xz + a4*z2 - a1*yz) - 2*vpz;
  long b = val(pr, 2*y + a1*x + a3*z) - vpz;
  long c = val(pr, 3*x2*x2 + b2*x2*xz + 3*b4*x2*z2 + 3*b6*xz*z2 + b8*z2*z2)
          -4*vpz;
#ifdef DEBUG_HEIGHT
cout<<"a = " << a << endl;
cout<<"b = " << b << endl;
cout<<"c = " << c << endl;
#endif
// some obvious changes enable calculation of lambda as a rational
// some improvements can be made if this is never to be done
// eg in the above, no need to work with projective coords, just use real x/z
  bigfloat halfn = to_bigfloat(n); halfn /= to_bigfloat(2);
  bigfloat lambda;
  
  if ( (a<=0) || (b<=0) ) 
    {
      lambda = vpz - val(pr,x);
      if(lambda<0) lambda=0;
    }
  else if ( ndiv(pr, c4) )
    {
     bigfloat m = to_bigfloat(b);
     if(halfn<m) m=halfn;  // m = min( b , halfn );
     lambda = (m*(m-n)) / n; 
    }
  else if ( c>=(3*b) ) 
     lambda = (-2*b) / to_bigfloat(3);
  else 
     lambda = -c / to_bigfloat(4);
  
  bigfloat ans = lambda * log( I2bigfloat(pr) );
#ifdef DEBUG_HEIGHT
cout<<"...returning lambda = " << lambda << ", pheight = "<<ans<<endl;
#endif
  return ans;
}

#undef DEBUG_HEIGHT
//#define DEBUG_HEIGHT

bigfloat realheight(const Point& P)
{
  bigfloat x,y;
  P.getrealcoordinates(x,y);
#ifdef DEBUG_HEIGHT
  cout<<"Computing real height of P = " << P <<", x(P) = "<<x<<endl;
#endif
  return realheight(x,P.E);
}

bigfloat realheight(const bigfloat& x, const Curvedata* E)
{
  bigint bb2,bb4,bb6,bb8;
  E->getbi(bb2,bb4,bb6,bb8);
  bigfloat b2 = I2bigfloat(bb2), b4 = I2bigfloat(bb4), 
           b6 = I2bigfloat(bb6), b8 = I2bigfloat(bb8);
  bigfloat b2dash = b2 - 12;
  bigfloat b4dash = b4 - b2 + 6;
  bigfloat b6dash = b6 - 2*b4 + b2 - 4;
  bigfloat b8dash = b8 - 3*b6 + 3*b4 - b2 + 3;
#ifdef DEBUG_HEIGHT
  cout<<"b2, b4, b6, b8 = "<<b2<<", "<<b4<<", "<<b6<<", "<<b8<<"\n";
  cout<<"b2dash, b4dash, b6dash, b8dash = "<<b2dash<<", "<<b4dash<<", "<<b6dash<<", "<<b8dash<<"\n";
#endif

  bigfloat t, w, z, zw;
  bigfloat H = to_bigfloat(4); 
            // max(4.0, max(abs(b2), max(2*abs(b4), max(2*abs(b6), abs(b8)))));
  t=abs(b2);   if(t>H) H=t;
  t=2*abs(b4); if(t>H) H=t;
  t=2*abs(b6); if(t>H) H=t;
  t=abs(b8);   if(t>H) H=t;
  long precision = decimal_precision();
#ifdef DEBUG_HEIGHT
  cout<<"decimal precision = "<<precision<<endl;
#endif
  long nlim=I2long(Iround(ceil( (5.0/3.0)*precision + 0.5 + 0.75*log( 7.0 + (4.0/3.0)*log(H) ))));
#ifdef DEBUG_HEIGHT
  cout<<"H = "<<H<<"; log(H) = "<<log(H)<<"; using "<<nlim<<" terms in the sum.\n";
#endif

  long beta;
  if ( abs(x) < 0.5 ) {t = 1 / (x + 1); beta = 0; }
  else {t = 1 / x; beta = 1; }
  bigfloat mu = -log( abs(t) ), dmu;  
  bigfloat f = to_bigfloat(1);
#ifdef DEBUG_HEIGHT
  cout<<"initial mu = "<<mu<<"\n";
#endif

  for (long n = 0; n <= nlim; n++)
  {
    f /= 4;
    if ( beta )
      {w = (((b6*t + 2*b4)*t + b2)*t + 4)*t;
       z = 1 - t*t*(b4 + t*(2*b6 + t*b8));
       zw = z + w;
      }
    else
      {w = (((b6dash*t + 2*b4dash)*t + b2dash)*t + 4)*t;
       z = 1 - t*t*(b4dash + t*(2*b6dash + t*b8dash));
       zw = z - w;
      }
    if ( abs(w) <= 2*abs(z) )
      {
	dmu=f*log(abs(z));
	mu += dmu; 
	t = w/z; 
      }
    else
      {
	dmu = f*log(abs(zw));
	mu += dmu; 
	t = w/zw;
#ifdef DEBUG_HEIGHT
	cout<<"switching...\n";
#endif
	beta = ! beta; 
      }
#ifdef DEBUG_HEIGHT
    cout<<"n="<<n<<": z, dmu, mu = "<<"\t"<<z<<"\n\t\t\t"<<dmu<<"\n\t\t\t"<<mu<<"\n";
#endif
  }
#ifdef DEBUG_HEIGHT
cout << "returning real height = " << mu << endl;
#endif

  return mu;
}
#undef DEBUG_HEIGHT

bigfloat height_pairing(Point& P, Point& Q) 
{
  // we avoid doing any real work, especially addition of points,
  // if we can.
  if(P.iszero() || Q.iszero())    return to_bigfloat(0);
  else if(P == Q)    return height(P) ;
  else
    {
      bigfloat hP = height(P);
      bigfloat hQ = height(Q);
      Point PQ = P+Q;
      bigfloat hPQ = height(PQ);
      bigfloat ans =  (hPQ - hP - hQ)/2;
      return ans;
    }
}

//#define DEBUG_REG
// regulator of a list of n points
bigfloat regulator(vector<Point>& P)   // nb not const; sets heights when found
{
#ifdef DEBUG_REG
  cout<<"In regulator with PointArray = " << P << endl;
#endif
  int i,j,k,d;
  int n = P.size();
  if( n <= 0) return to_bigfloat(1);
  if(n == 1) return height(P[0]) ;
  if(n == 2 )
    {
     bigfloat pair00 = height(P[0]) ;
     bigfloat pair01 = height_pairing(P[0], P[1]);  // nb this will set height
     bigfloat pair11 = height(P[1]);                // of P[1]; is efficient
     bigfloat reg = pair00 * pair11 - pair01 * pair01;
     return  reg; 
    }
  if (n == 3)
    {
#ifdef DEBUG_REG
      cout<<"n=3, computing height pairing matrix..."<<flush;
#endif
     bigfloat pair[3][3] ;
     for (i = 0; i < 3; i++)
       {pair[i][i] = height(P[i]) ;
        for (j = i + 1; j < 3; j++)
          {pair[i][j] = pair[j][i] = height_pairing(P[i], P[j]) ; }
       }
#ifdef DEBUG_REG
      cout<<"done.  Matrix = " << endl;
      for (i = 0; i < 3; i++)
        {for (j = 0; j < 3; j++) cout << pair[i][j] << "\t";
         cout << endl;
       }
#endif
     bigfloat reg = (pair[0][0] * ( pair[1][1] * pair[2][2] - pair[1][2] * pair[1][2] )
                 - pair[0][1] * ( pair[0][1] * pair[2][2] - pair[1][2] * pair[0][2] )
                 + pair[0][2] * ( pair[0][1] * pair[1][2] - pair[1][1] * pair[0][2] )
                   );
#ifdef DEBUG_REG
      cout<<"regulator = " << reg << endl;
#endif
     return reg;
    }
  if (n == 4)
    {
     bigfloat pair[4][4] ;
     for (i = 0; i < 4; i++)
       {pair[i][i] = height(P[i]) ;
        for (j = i + 1; j < 4; j++)
          {pair[i][j] = height_pairing(P[i], P[j]) ; }
       }
     //
     // the following explicit expression is courtesy of Maple and AB
     // (Maple uses **2, we converted these into explicit squares)
     //
     // it purports to be the expression for the determinant of our symmetric
     // pairing matrix
     //
     bigfloat reg = (
       ((2 * pair[1][2] * pair[3][3]-2* pair[1][3] *pair[2][3])*pair[0][1]+
        (-pair[1][1]*pair[3][3]+pair[1][3]*pair[1][3])*pair[0][2])
         * pair[0][2]+pair[0][0]
         *(pair[1][1]*pair[2][2]*pair[3][3]
         - pair[1][1]*pair[2][3]*pair[2][3]-pair[1][2]*pair[1][2]*pair[3][3]
         + 2*pair[1][2]*pair[1][3]*pair[2][3]-pair[1][3]*pair[1][3]*pair[2][2])
         + (-pair[2][2]*pair[3][3]+pair[2][3]*pair[2][3])*pair[0][1]*pair[0][1]
           + ((2*pair[1][3]*pair[2][2]-2*pair[1][2]*pair[2][3])*pair[0][1]
              + (2*pair[1][1]*pair[2][3]-2*pair[1][3]*pair[1][2])*pair[0][2]
              + (-pair[1][1]*pair[2][2]+pair[1][2]*pair[1][2])*pair[0][3])
             *pair[0][3]
             );
     return reg;
    }
  if ( n <= 50)
    {bigfloat pair[50][50] ;
     // initialize the matrix of pairings
     for (i = 0; i < n; i++)
       {pair[i][i] = height(P[i]) ;
        for (j = i + 1; j < n; j++)
          {pair[j][i] = pair[i][j] = height_pairing(P[i], P[j]) ; }
       }
     // Gaussian elimination 
     // for the first n - 1 rows
     for (j = 0 ; j < n - 1; j ++)
       {// use row j to pivot with
         bigfloat pivot = pair[j][j] ;
         // kill off rows below row j
         for (i = j + 1; i < n ; i ++)
           {bigfloat multiplier = pair[i][j] / pivot ;
            // subtract multiplier * row j from row i 
            //noting the beginning of row i is already zeroed
            for (k = j ; k < n; k++)
              {pair[i][k] -= multiplier * pair[j][k] ; }
           }
       }
     // now reg is the product of the diagonal entries
     bigfloat reg = to_bigfloat(1) ;
     for (d = 0; d < n; d++) reg *= pair[d][d] ;
     return reg ;
    }
  // else
  //  n> 50 not yet (could fold into last case)
  cout << "## If you really want the regulator of more than 50 points,\n";
  cout << "then edit heights.cc youself!" << endl;
  abort();
  return to_bigfloat(1) ;
}

// end of HEIGHTS.CC
