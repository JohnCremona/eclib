// realroots.cc: implementation of funtions for real roots of polynomials
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
 
#include "compproc.h"
#include "realroots.h"

bigfloat safe_sqrt(const bigfloat& x)
{
  static bigfloat zero=to_bigfloat(0);
  if(x>zero) return sqrt(x);
  return zero;
}


bigfloat cube_root(const bigfloat& x)
{
  if(is_zero(x)) return x;
  if(x<0) return -exp(log(-x)/3);
  return exp(log(x)/3);
}
/*
bigfloat cube_root(const bigfloat& x)
{
  static bigfloat third = to_bigfloat(1)/to_bigfloat(3);
  if(x<0) return -pow(-x, third);
  else    return  pow( x, third);
}
*/

// coeff contains deg+1 reals starting with the leading coefficient
// which must be nonzero
// 
// we assume the roots are distinct

//#define DEBUG_REALROOTS

vector<bigfloat> realroots( const vector<bigfloat>& coeff )
{
#ifdef DEBUG_REALROOTS
 cout<<"In realroots with coeff = "<<coeff<<endl;
#endif
  // trim leading zeros:
  vector<bigfloat> tcoeff;
  unsigned int i=0; while(is_zero(coeff[i])) i++;
  while(i<coeff.size()) tcoeff.push_back(coeff[i++]);
#ifdef DEBUG_REALROOTS
  cout<<"realroots: trimmed coeffs = "<<tcoeff<<endl;
#endif

  long deg = tcoeff.size()-1;
#ifdef DEBUG_REALROOTS
  cout<<"deg = "<<deg<<endl;
#endif
  vector<bigfloat> ans;
  bigfloat a = tcoeff[0];
#ifdef DEBUG_REALROOTS
  cout<<"a = "<<a<<" (should be nonzero!)"<<endl;
#endif
  if(deg==0) return ans;
  if(deg==1) 
    {
      ans.push_back(-tcoeff[1]/a); 
      return ans;
    }
  if(deg==2)
    {
      bigfloat b=tcoeff[1], c=tcoeff[2]; 
      bigfloat disc = b*b-4*a*c;
      if(disc<0) return ans;
      bigfloat r1 = (safe_sqrt(disc)-b)/(2*a);
      ans.push_back(r1);
      ans.push_back(-(b/a)-r1);
      return ans;
    }
  if(deg==3)
    {
#ifdef DEBUG_REALROOTS
      cout<<"degree 3 case "<<endl;
#endif
      bigfloat b=tcoeff[1], c=tcoeff[2], d=tcoeff[3]; 
      bigfloat P = b*b-3*a*c;
      bigfloat Q = b*c-9*a*d;
      bigfloat R = c*c-3*b*d;
      bigfloat D3 = Q*Q-4*P*R;  // =-3*disc
      if(D3>0) // one real root
	{
#ifdef DEBUG_REALROOTS
	  cout<<"one real root "<<endl;
#endif
	  if(is_zero(P)) 
	    {
#ifdef DEBUG_REALROOTS
	      cout<<"Case P=0 (P="<<P<<")"<<endl;
	      cout<<"a = "<<a<<endl;
	      cout<<"Q = "<<Q<<endl;
	      cout<<"3*a*Q = "<<3*a*Q<<endl;
#endif
	      bigfloat eta3  = d/a - (c*R)/(3*a*Q);
#ifdef DEBUG_REALROOTS
	      cout<<"eta3 = "<<eta3<<endl;
#endif
	      bigfloat eta   = cube_root(eta3);
#ifdef DEBUG_REALROOTS
	      cout<<"eta = "<<eta<<endl;
#endif
	      bigfloat alpha = -eta - R;
#ifdef DEBUG_REALROOTS
	      cout<<"root = "<<alpha<<endl;
#endif
	      ans.push_back(alpha);
	      return ans;
	    }
	  else
	    {
	      bigfloat U =  2*b*b*b + 27*a*a*d - 9*a*b*c;
#ifdef DEBUG_REALROOTS
	      cout<<"Case P!=0 (P="<<P<<")"<<endl;
	      cout<<"U = "<<U<<endl;
#endif

	      bigfloat delta = safe_sqrt(D3);
#ifdef DEBUG_REALROOTS
	      cout<<"delta = "<<delta<<endl;
#endif
	      bigfloat gamma1 = (-Q+delta)/(2*P);  // roots of Hessian
	      bigfloat gamma2 = (-Q-delta)/(2*P);  //
#ifdef DEBUG_REALROOTS
	      cout<<"gamma1 = "<<gamma1<<endl;
	      cout<<"gamma2 = "<<gamma2<<endl;
#endif

	      bigfloat lambda=U-3*a*delta;  // lambda+mu=2*U, lambda*mu=4*P^3
	      bigfloat mu    =2*U-lambda;
#ifdef DEBUG_REALROOTS
	      cout<<"lambda = "<<lambda<<endl;
	      cout<<"mu     = "<<mu<<endl;
#endif
	      lambda=cube_root(lambda);	      mu=cube_root(mu);
	      bigfloat alpha = (lambda*gamma1-mu*gamma2)/(lambda-mu);
#ifdef DEBUG_REALROOTS
	      cout<<"root = "<<alpha<<endl;
#endif
	      ans.push_back(alpha);
	      return ans;
	    }
	}
      else // all roots real
	{
	  vector<bigcomplex> croots = solvecubic(b/a,c/a,d/a);
	  for(int i=0; i<3; i++) ans.push_back(real(croots[i]));
	  return ans;
	}
    }
  if(deg==4)
    {
    bigfloat b=tcoeff[1], c=tcoeff[2], d=tcoeff[3], e=tcoeff[4]; 
    bigfloat ii = 12*a*e - 3*b*d + c*c;
    bigfloat jj =  (72*a*e + 9*b*d - 2*c*c) * c - 27*(a*d*d + b*b*e);
    bigfloat disc = 4*ii*ii*ii-jj*jj;
    bigfloat H = 8*a*c - 3*b*b, R = b*b*b + 8*d*a*a - 4*a*b*c;
    bigfloat Q = H*H-16*a*a*ii;  // = 3*Q really
    int type, nrr;
    if(disc<0) 
      {type=3; nrr=2;}       // 2 real roots
    else 
      {
	if((H<0)&&(Q>0)) 
	  {type=2; nrr=4;}   // 4 real roots
	else 
	  {type=1; nrr=0;}   // 0 real roots
      }
    bigcomplex c1(to_bigfloat(0)), c2(-3*ii), c3(jj);
    vector<bigcomplex> cphi = solvecubic( c1, c2, c3);
    vector<bigcomplex> roots(4);
    bigfloat a4=4*a;
    bigfloat oneover4a = to_bigfloat(1)/a4;
    
    if(type==1) return ans;

    if(type<3)  // all roots are real
      {
	// all the phi are real;  order them so that a*phi[i] decreases
	bigfloat phi1 = real(cphi[0]);
	bigfloat phi2 = real(cphi[1]);
	bigfloat phi3 = real(cphi[2]);
	if(a>0)      orderreal(phi1,phi2,phi3); 
	else         orderreal(phi3,phi2,phi1); 

	bigfloat r1 = safe_sqrt((a4*phi1-H)/to_bigfloat(3));
	bigfloat r2 = safe_sqrt((a4*phi2-H)/to_bigfloat(3));
	bigfloat r3 = safe_sqrt((a4*phi3-H)/to_bigfloat(3));
	if(R<0)  r3 = -r3;
	ans.push_back(( r1 + r2 - r3 - b) * oneover4a);
	ans.push_back(( r1 - r2 + r3 - b) * oneover4a);
	ans.push_back((-r1 + r2 + r3 - b) * oneover4a);
	ans.push_back((-r1 - r2 - r3 - b) * oneover4a);
	// Those are all real and in descending order of size
	return ans;
      }
    bigfloat realphi;       // will hold the real root, which will be cphi[2]
    if (is_real(cphi[1])) 
      {
	realphi=real(cphi[1]);
	cphi[1]=cphi[2];
	cphi[2]=realphi;
      }
    else 
      if (is_real(cphi[2])) 
	{
	  realphi=real(cphi[2]);
	}
      else 
	{
	  realphi=real(cphi[0]);
	  cphi[0]=cphi[2];
	  cphi[2]=realphi;
	}
    bigcomplex r1 = sqrt((a4*cphi[0]-H)/to_bigfloat(3));
    bigfloat r3   = safe_sqrt((a4*realphi-H)/to_bigfloat(3));
    if(R<0)  r3 = -r3;
    ans.push_back((( 2*real(r1) - r3 - b)) * oneover4a);
    ans.push_back(((-2*real(r1) - r3 - b)) * oneover4a);
    return ans;
    }
  return ans; // not implemented for degree>4
}

// As above but only root in the interval [-1,1]

vector<bigfloat> realroots11( const vector<bigfloat>& coeff )
{
#ifdef DEBUG_REALROOTS
  cout<<"In realroots11 with coeff = "<<coeff<<endl;
#endif
  vector<bigfloat> ans0 = realroots(coeff);
#ifdef DEBUG_REALROOTS
  cout<<"realroots11: ans0 = "<<ans0<<endl;
#endif
  vector<bigfloat> ans;
  for(unsigned int i=0; i<ans0.size(); i++)
    {
      bigfloat x = ans0[i];
      if((x<=1)&&(x>=-1)) ans.push_back(x);
    }
#ifdef DEBUG_REALROOTS
  cout<<"realroots11: ans = "<<ans<<endl;
#endif
  return ans;
}
