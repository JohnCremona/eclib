// heights.cc: implementation of height functions declared in points.h
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
 
#include <eclib/points.h>  // which includes curve.h
#ifdef MPFP  // use NTL to compute the determinant
#include <NTL/mat_RR.h>
#else
#define MAX_RANK_REG 50 // cannot ask for regulator of more than 50 points.
#endif

//#define DEBUG_HEIGHT

bigfloat height(Point& P)
{
#ifdef DEBUG_HEIGHT
  cout<<"Computing height of P = "<<P<<endl;
  cout<<"(P.E = "<<P.E<<")\n";
  cout<<"(E = "<<(Curve)*(P.E)<<")\n";
  cout<<"(current height attribute is "<<P.height<<")\n";
#endif
  if (!(P.isvalid()))
    {
      cerr<<"Run-time error: point "<<P<<" is not valid on its curve "<<(Curve)(P.getcurve())<<endl;
      exit(1);
    }
  bigfloat zero(to_bigfloat(0));
  if (P.height >= zero) return P.height;  // already calculated it
  if (P.is_zero())  {P.height = zero; return zero; } // zero height if torsion
  if (order(P) > 0) {P.height = zero; return zero; } // zero height if torsion

  // N.B. So if we ever ask a point its height it will compute its order.
  // otherwise need to calculate it

  // The local height at p will only be correctly computed by
  // pheight() if the curve is minimal at p

  Curvedata* E = P.E;
  vector<bigint> bad_p = getbad_primes(*E);
  Curvedata Emin;
  Point Pmin(Emin); // assigns Pmin.E to a pointer to Emin
  // NB the is_minimal function returns 0 when minimization has not
  // been done; the curve may still be minimal
  if (!is_minimal(*E))
    {
      bigint u, r, s, t;
      Emin = E->minimalize(u,r,s,t);
      Pmin = transform(P, &Emin, u, r, s, t);
      bad_p = getbad_primes(Emin);
    }
  else
    {
      Pmin = P;
    }
  // Add local heights at finite primes dividing discr(E) OR denom(P).
  // The components for primes dividing denom(P) add to log(denom(x(P)));
  //   since P=(XZ:Y:Z^3), denom(P)=Z=gcd(XZ,Z^3), called "zroot" here,
  //   and so the contribution is log(denom(x(P))) = 2*log(zroot).
  //   This avoids factorizing the denominator.

  const bigint& zroot = gcd(Pmin.getX(),Pmin.getZ());   // = cube root of Z
  bigfloat h = realheight(Pmin);
#ifdef DEBUG_HEIGHT
  cout<<" - real height = "<<h<<"\n";
#endif

  // contribution from primes dividing the denoinator:
  h += 2*log(I2bigfloat(zroot));
#ifdef DEBUG_HEIGHT
  cout<<" - after adding log(denom), height = "<<h<<"\n";
#endif

#ifdef DEBUG_HEIGHT
  cout<<" - E (min) = "<<(Curve)Emin<<" with bad primes "<<bad_p<<"\n";
#endif

  for ( const auto& p : bad_p)
    {
#ifdef DEBUG_HEIGHT
      cout<<" - bad prime p = "<<p<<"\n";
#endif
      // we already have included the local height at p for primes p
      // dividing the denominator
      if(ndiv(p,zroot))
        {
          bigfloat pht = pheight(Pmin,p);
#ifdef DEBUG_HEIGHT
          cout<<" - local height at p is "<<pht<<"\n";
#endif
            h += pht;
        }
#ifdef DEBUG_HEIGHT
      else
        cout<<" - p divides denominator (zroot="<<zroot<<"), so ignoring\n";
#endif
    }
  P.height = h;
#ifdef DEBUG_HEIGHT
  cout << "height(P) returns "<<h<<endl;
#endif
  return h;
}

#undef DEBUG_HEIGHT

bigfloat pheight(const Point& P, const bigint& pr)
// NB The local height at p will only be correctly computed by
// pheight() if the curve is minimal at p
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
  
  bigfloat h = lambda * log( I2bigfloat(pr) );
#ifdef DEBUG_HEIGHT
cout<<"...returning lambda = " << lambda << ", pheight = "<<h<<endl;
#endif
  return h;
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

#ifdef MPFP // Multi-Precision Floating Point
  long original_prec, new_prec;
  original_prec = bit_precision();
  new_prec = original_prec + 100;
  //  cout<<"Setting bit precision to "<<new_prec<<endl;
  set_bit_precision(new_prec); // does not change output precision
#endif

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

  // NB We use decimal precision here since the formula for nlim (from
  // Silverman) is given in terms of decimal places required.
  long precision = decimal_precision();
#ifdef DEBUG_HEIGHT
  cout<<"decimal precision = "<<precision<<endl;
#endif
  long nlim=I2long(Iround(ceil( (5.0/3.0)*precision + 0.5 + 0.75*log( 7.0 + (4.0/3.0)*log(H) ))));
  nlim *=2;
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


#ifdef MPFP // Multi-Precision Floating Point
//  cout<<"Setting bit precision back to "<<original_prec<<endl;
 set_bit_precision(original_prec); // does not change output precision
#endif
  return mu;
}
#undef DEBUG_HEIGHT

bigfloat height_pairing(Point& P, Point& Q) 
{
  // we avoid doing any real work, especially addition of points,
  // if we can.
  if(P.is_zero() || Q.is_zero())    return to_bigfloat(0);
  else if(P == Q)    return height(P) ;
  else
    {
      bigfloat hP = height(P);
      bigfloat hQ = height(Q);
      Point PQ = P + Q;
      return (height(PQ) - hP - hQ)/2;
    }
}

//#define DEBUG_REG
// regulator of a list of n points
bigfloat regulator(vector<Point>& P)   // nb not const; sets heights when found
{
#ifdef DEBUG_REG
  cout<<"In regulator with PointArray = " << P << endl;
#endif
  int n = P.size();
  if( n <= 0) return to_bigfloat(1);
  if(n == 1) return height(P[0]) ;
  if(n == 2 )
    {
     bigfloat pair00 = height(P[0]);
     bigfloat pair11 = height(P[1]);
     Point Q = P[0] + P[1];
     bigfloat h = height(Q);
     bigfloat pair01 = (h-pair00-pair11)/2;
     return pair00 * pair11 - pair01 * pair01;
    }
#ifdef MPFP  // use NTL to compute the determinant
  // initialize the matrix of pairings
  mat_RR height_matrix;
  height_matrix.SetDims(n,n);
  for (int i = 0; i < n; i++)
    {
      height_matrix[i][i] = height(P[i]) ;
    }
  for (int i = 0; i < n; i++)
    {
      for (int j = i + 1; j < n; j++)
          {
            Point Q = P[i] + P[j];
            bigfloat h = (height(Q) - height_matrix[i][i] - height_matrix[j][j])/2 ;
            height_matrix[j][i] = h;
            height_matrix[i][j] = h;
          }
    }
  return determinant(height_matrix);

#else // use a naive determinant method

  if (n == 3)
    {
     bigfloat pair[3][3] ;
     for (int i = 0; i < 3; i++)
       {pair[i][i] = height(P[i]) ;
        for (int j = i + 1; j < 3; j++)
          {pair[i][j] = pair[j][i] = height_pairing(P[i], P[j]) ; }
       }
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
     for (int i = 0; i < 4; i++)
       {pair[i][i] = height(P[i]) ;
        for (int j = i + 1; j < 4; j++)
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
  if ( n > MAX_RANK_REG)
    {
      //  n> 50 not yet (could fold into last case)
      cout << "## Assuming that the regulator of more than "<<MAX_RANK_REG<<" points is 0" << endl;
      return to_bigfloat(0) ;
    }
  bigfloat pair[MAX_RANK_REG][MAX_RANK_REG] ;
  // initialize the matrix of pairings
  for (int i = 0; i < n; i++)
    {
      pair[i][i] = height(P[i]) ;
      for (int j = i + 1; j < n; j++)
        {
          pair[j][i] = pair[i][j] = height_pairing(P[i], P[j]) ;
        }
    }
  // Gaussian elimination
  // for the first n - 1 rows
  for (int j = 0 ; j < n - 1; j ++)
    {// use row j to pivot with
      bigfloat pivot = pair[j][j] ;
      // kill off rows below row j
      for (int i = j + 1; i < n ; i ++)
        {bigfloat multiplier = pair[i][j] / pivot ;
          // subtract multiplier * row j from row i
          //noting the beginning of row i is already zeroed
          for (int k = j ; k < n; k++)
            {
              pair[i][k] -= multiplier * pair[j][k] ;
            }
        }
    }
  // now reg is the product of the diagonal entries
  bigfloat reg = to_bigfloat(1) ;
  for (int d = 0; d < n; d++) reg *= pair[d][d] ;
  return reg ;
#endif
}

// end of HEIGHTS.CC
