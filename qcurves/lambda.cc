// lambda.h   Declarations of functions which compute Silverman's
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
 //            finite set Lambda_bad for a curve


// N.B. (1) Uses my height normalization, double S's.
// (3) Uses the local height normalization WITHOUT the log|Delta|
// (2) Intended for use in computing Heegner points (not yet implemented)


#include "points.h"
#include "lambda.h"

#define MAX_NUM_LAMBDA 1000

vector<bigfloat> lambda_bad_1(const bigint& p, long kcode, long npd, long& nlambda)
{
  bigfloat logp = log(I2bigfloat(p)), n=to_bigfloat(npd);

  if((kcode%10)==0) // Type I_m
    {
      long i, m = kcode/10;
      nlambda=1+(m/2);
      vector<bigfloat> ans; ans.reserve(nlambda);
      for(i=0; i<nlambda; i++) ans.push_back( ( (i*i)/n  - i ) * logp);
      return ans;
    }

  if((kcode%10)==1) // Type I*m
    {
      bigfloat m = to_bigfloat(kcode-1)/10;
      nlambda = 3;
      vector<bigfloat> ans; ans.reserve(nlambda);
      ans.push_back(to_bigfloat(0));
      ans.push_back(-logp);
      ans.push_back(-(1 +m/4)*logp);
      return ans;
    }

  if((kcode==2)||(kcode==7)||(p>3))
    {
      nlambda = 1;
      vector<bigfloat> ans(1,to_bigfloat(0));
      return ans;
    }

  nlambda = 2;
  vector<bigfloat> ans; ans.reserve(nlambda);
  ans.push_back(to_bigfloat(0));
  long nn = (kcode<5? kcode: kcode+3);
  ans.push_back( -(nn*logp)/6 );
  return ans;
}

vector<bigfloat> lambda_bad(const CurveRed& C, long& nlambda, int verbose)
{
  vector<bigfloat> ans;
  nlambda = 1;
  ans.push_back(to_bigfloat(0));
  bigint discr = getdiscr(C);
  vector<bigint> plist = getbad_primes(C);
  long i, j, nl, nnl;
  vector<bigint>::const_iterator pr=plist.begin();
  while(pr!=plist.end())
    {
      bigint p = *pr++;
      if (ndiv(p*p,discr))
	{
	  if(verbose)
	    cout<<"Lambda_bad("<<p<<") has only one element, 0.\n";
	  continue;
	}
      // else do some real work
      long kcode = getKodaira_code(C,p).code;
      long npd   = getord_p_discr(C,p);
      vector<bigfloat> list = lambda_bad_1(p,kcode,npd,nl);
      if(verbose)
	{
	  cout << "Lambda_bad("<<p<<") has " << nl << " element(s): ";
	  cout << list << endl;
	}
      nnl = nlambda*nl;
      ans.reserve(nnl);
      for(i=0; i<nlambda; i++)
	for(j=0; j<nl; j++)
	  ans.push_back(ans[i]+list[j]);
      nlambda = nnl;
    }  // end of loop on p
  return ans;
}

int make_point_from_x(Curvedata* CD, const bigint& xa, const bigint& xd, Point* P)
{
  bigint a(xa),b,c,d(xd);
  if(d<0) {a=-a; d=-d;}
  if(!isqrt(d,c)) return 0;
  bigint b2,b4,b6,b8;
  CD->getbi(b2,b4,b6,b8);
  bigint d2=d*d;
  bigint d3=d2*d;
  bigint e,e2 = ((4*a+b2*d)*a + 2*b4*d2)*a + b6*d3;
  if(!isqrt(e2,e)) return 0;
  bigint a1,a2,a3,a4,a6;
  CD->getai(a1,a2,a3,a4,a6);
  b = (e-a1*a*c-a3*c*d)/2;
  *P = Point(CD,a*c,b,c*d);
  return 1;
}

int make_point_from_x(Curvedata* CD, const bigfloat& x, long maxdd, Point* P)
{
  bigint a,b,c,d;
//cout<<"In ratapprox2 with x = " << x << endl;
  bigint x0, x1, x2, y0, y1, y2;
  bigfloat xx, diff, xc;
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
  diff = 1;
  bigint maxdenom = pow(BIGINT(10),maxdd);
  while ( !is_approx_zero(diff) && (y2<maxdenom))
    { c = Iround( xx ); xc=I2bigfloat(c);
      x2 = x0 + c*x1; x0 = x1; x1 = x2;
      y2 = y0 + c*y1; y0 = y1; y1 = y2;
      if(make_point_from_x(CD,x2,y2,P)) return 1;
      diff = abs( x - (I2bigfloat(x2)/I2bigfloat(y2)) );
//cout<<"x2,y2,diff = " << x2 << ", " << y2 << ", " << diff << endl;
      if ( is_approx_zero(abs(xx - xc)) ) diff = 0;
      else xx = to_bigfloat(1)/(xx - xc);
    }
  a = x2; d = y2;
  if ( d < 0 )
    {::negate(a); ::negate(d); }

  if(!isqrt(d,c)) return 0;
  bigint b2,b4,b6,b8;
  CD->getbi(b2,b4,b6,b8);
  bigint d2=d*d;
  bigint d3=d2*d;
  bigint e,e2 = ((4*a+b2*d)*a + 2*b4*d2)*a + b6*d3;
  if(!isqrt(e2,e)) return 0;
  bigint a1,a2,a3,a4,a6;
  CD->getai(a1,a2,a3,a4,a6);
  b = (e-a1*a*c-a3*c*d)/2;
  *P = Point(CD,a*c,b,c*d);
  return 1;
}

int make_point_from_x_and_ht(Curvedata* CD, vector<bigfloat> lambdas, const bigfloat& xp, const bigfloat& ht, Point* P)
{
  bigfloat rh = realheight(xp,CD);
  vector<bigfloat>::const_iterator lam = lambdas.begin();
  while(lam!=lambdas.end())
    {
      bigfloat logd = (ht-rh-(*lam++))/2;
      bigfloat approxd = exp(logd);
      bigint xa, xd2, xd = Iround(approxd);
      if(xd>0) 
	{
	  xd2 = xd*xd;
	  xa = Iround(xp*I2bigfloat(xd2));
	  if(make_point_from_x(CD,xa,xd2,P)) return 1;
	}

      bigint xdx; long id, delta=10;
      for(id=1;id<=delta;id++)
	{
	  xdx=xd+delta;
	  if(xdx>0) 
	    {
	      xd2 = xd*xd;
	      xa = Iround(xp*I2bigfloat(xd2));
	      if(make_point_from_x(CD,xa,xd2,P)) return 1;
	    }
	  xdx=xd-delta;
	  if(xdx>0) 
	    {
	      xd2 = xd*xd;
	      xa = Iround(xp*I2bigfloat(xd2));
	      if(make_point_from_x(CD,xa,xd2,P)) return 1;
	    }
	}
    }
  return 0;
}
