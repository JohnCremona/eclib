// file pcprocs.cc: implementation of functions used to compute periods
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

#include "marith.h"
#include "compproc.h"
#include "moddata.h"
#include "symb.h"
#include "oldforms.h"
#include "homspace.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "newforms.h"
#include "periods.h"
#include "pcprocs.h"

void ratapprox(bigfloat x, long& a, long& b)
{
  long c, x0, x1, x2, y0, y1, y2;
  bigfloat xx, diff;
  bigfloat eps = to_bigfloat(1.0e-6);
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
//cout<<"ratapprox("<<x<<"): \n";
  diff = 1;
  while ( diff > eps )
    { c = longify( xx + ( (xx>0) ? 0.5 : -0.5 ) ); // ie round(xx)
      x2 = x0 + c*x1; x0 = x1; x1 = x2;
      y2 = y0 + c*y1; y0 = y1; y1 = y2;
      diff = abs( x - (to_bigfloat(x2)/to_bigfloat(y2)) );
//cout<<"a = "<<x2<<", b = "<<y2<<", diff = "<<diff<<endl;
      if ( abs(xx - c) < eps ) diff = 0;
      else xx = 1/(xx - c);
    }
  a = x2; b = y2;
  if ( b < 0 ) {a=-a; b=-b; }
  if ( x < 0 ) {a=-a;}
}

// Function used to test whether a denominator found by ratapprox() is
// "trustworthy": always with multiprecision, but only if below
// a fixed bound otherwise.
#ifdef MPFP // Multi-Precision Floating Point
inline int trust_denom(long d) { return (d<1000);}
#else
inline int trust_denom(long d) { return (d<251);}
#endif

// Given a newform (the i'th in newforms) at level n with a real
// period x0; Finds a matrix [a0,b0;Nc0,d0] whose integral is
// dotplus*x0+dotminus*y0*i N.B. The value of x0 may be changed, if we
// come across some period whose real part isnot an integer multiple
// of the original x0.

// We compute periods of lots of matrices in Gamma_0(N), over all
// symbols {0,b/d} for d<dmax.  We are looking for (i) a symbol whose
// imaginary part is nonzero; (ii) a symbol whose real and imaginary
// parts are both non-zero, which we store.

int newforms::find_matrix(long i, long dmax, int&rp_known, bigfloat&x0, bigfloat&y0)
{
  int have_both=0, have_ip=0;
  int have_rp = get_real_period(i,x0,verbose);
  rp_known = have_rp;
  // have_rp is set if we know a real period; rp_known is set if we
  // know that it is the real period (strictly, the least real part of
  // a period), which is usually the case here, the exception being
  // when sfe=-1 and n is square.

  // the following are to allow us a choice for the real period; 
  // we'll use the value which has greater precision
  long lplus = nflist[i].lplus;
  if(nflist[i].dp0!=0) lplus=0;
  int rp_fixed = !have_rp;

  periods_direct integrator(this,&(nflist[i]));
  long nrx=1, nry=1, drx=1, dry=1;
  long nrx0=1, nry0=1, drx0=1, dry0=1;
  long a, b, c, d;
  long& dotplus=(nflist[i].dotplus);
  long& dotminus=(nflist[i].dotminus);
  long dotplus0=1, dotminus0=1;
  
  for(d=2; (d<dmax)||(!have_both); d++)
    {
      if(gcd(modulus,d)!=1) continue; long d2=d/2;
      for(b=1; (b<=d2); b++)
	{
	  if(bezout(d,modulus*b,a,c)!=1) continue;
	  c=-c;
	  if(verbose>1) 
	    cout << "Matrix ("<<a<<","<<b<<";"<<modulus*c<<","<<d<<"):\n";
	  integrator.compute(a,b,c,d);
	  bigfloat x = abs(integrator.rper());
	  if(have_rp)
	   {
	     bigfloat ratio=x/x0;
	     ratapprox(ratio, nrx, drx);
	     if(verbose>1) 
	       cout<<"real part = " << x << ", x/x0 = " << ratio   
		   << " =~= "<<nrx<<"/"<<drx<<"\n";
	     if(rp_known)
	       if(drx>1)
		 {
		   cout<<"******************************real part of period not an multiple of x0?\n";
		   if(verbose<=1) 
		     cout<<"real part = " << x << ", x/x0 = " << ratio   
			 << " =~= "<<nrx<<"/"<<drx<<"\n";		   
		   if(drx>10)
		     {
		       drx=1;
		       nrx=I2long(Iround(ratio));
		       cout << "Using rounded value nrx=" << nrx <<endl;
		     }
		 }
	     if(trust_denom(drx)) dotplus0=lcm(dotplus0,drx);

// fix the value of x0 if the current x is more accurate:
	     if(!rp_fixed&&(nrx>0))
	       {
		 if(d<lplus) 
		   {
		     x0=x*to_bigfloat(drx)/to_bigfloat(nrx);
		     rp_fixed=1;
		     if(verbose>1)
		       cout<<"replacing original x0 by "<<x0
			   <<" (which is more accurate since "<<d<<" < "<<lplus
			   <<"), to make the preceding ratio exact.\n";
		   }
		 rp_fixed=1;  // d will not get any smaller...
	       }
	   }
	  else
	    if(x>0.001)
	      {
		x0=x; have_rp=1; nrx=drx=1;
		if(verbose>1) cout<<"real period = " << x0 << "\n";
	      }
	  bigfloat y = abs(integrator.iper());
	  if(have_ip)
	   {
	     bigfloat ratio=y/y0;
	     ratapprox(ratio, nry, dry);
	     if(trust_denom(dry)) dotminus0=lcm(dotminus0,dry);
	     if(verbose>1) 
	       cout<<"imag part = " << y << ", y/y0 = " << ratio   
		   << " =~= "<<nry<<"/"<<dry<<"\n";
	   }
	  else
	    if(y>0.001)
	      {
		y0=y; have_ip=1; nry=dry=1;
		if(verbose>1) cout<<"imag period = " << y0 << "\n";
	      }
	  if(!have_both && (x>0.001) && (y>0.001) && trust_denom(dry))
	    {
	      have_both=1;
	      if(trust_denom(drx)) {nrx0=nrx; drx0=drx;} 
	      nry0=nry; dry0=dry;
	      //	      cout<<"nrx0="<<nrx0<<endl;	      
	      //	      cout<<"drx0="<<drx0<<endl;	      
	      nflist[i].a=a;
	      nflist[i].b=b;
	      nflist[i].c=c;
	      nflist[i].d=d;
	   }
       } // end of b loop
     } // end of d loop
  x0/=to_bigfloat(dotplus0);
  y0/=to_bigfloat(dotminus0);
  dotplus =(dotplus0 *nrx0)/drx0;
  dotminus=(dotminus0*nry0)/dry0;
  if(verbose>1){
    cout<<"dotplus0 ="<<dotplus0<<endl;
    cout<<"dotminus0="<<dotminus0<<endl;
    cout<<"dotplus  ="<<dotplus<<endl;
    cout<<"dotminus ="<<dotminus<<endl;
  }
  if(!have_both) {a=d=1; b=c=0; dotplus=dotminus=0;}
  return have_both;
}

int get_curve(long n, long fac, long maxnx, long maxny,
	      const bigfloat& x0, const bigfloat& y0, 
	      long& nx, long& ny, int& type, int detail)
{
  static bigfloat zero=to_bigfloat(0);
  long fac6=(odd(fac)?fac:2*fac);
  if(detail&&(fac>1)) cout<<"c6 factor " << fac6 << endl;
	 
  bigcomplex w1, w2; bigcomplex c4, c6;
  bigfloat x1=x0, y1=y0;
  for(nx=1; (nx<=maxnx); nx++)
    {x1=x0/to_bigfloat(nx);
    for(ny=1; (ny<=maxny); ny++)
      {
	y1=y0/to_bigfloat(ny);
	for(type=1; (type<=2); type++) 
	  {
	    if(type==2){w1=bigcomplex(x1,zero); w2=bigcomplex(zero,y1);}
	    else {w1=bigcomplex(2*x1,zero); w2=bigcomplex(x1,y1);}
	    bigcomplex tau=normalize(w1,w2);
	    getc4c6(w1,w2,c4,c6);
	    bigint ic4 = fac*Iround(real(c4)/fac);
	    bigint ic6 = fac6*Iround(real(c6)/fac6);
	    int close=abs(I2bigfloat(ic4)-real(c4))<0.00001;
	    int validc4c6 = valid_invariants(ic4,ic6);
	    if((validc4c6||close)&&detail)
	      {
		cout << "type = " << type << ", nx = " << nx << ", ny = " << ny << "\n"; 
		cout << "w1 = " << w1 << ", w2 = " << w2 << "\n";
		cout << "c4 = " << real(c4) << ", c6 = " << real(c6) << "\n";
		cout << "ic4 = " << ic4 << ", ic6 = " << ic6 << "\n";
	      }
//	    if(validc4c6&&close)
	    if(validc4c6)
	      {
		Curve C(ic4,ic6);
		Curvedata CD(C,1);
		CurveRed CR(CD);
		bigint cond = getconductor(CR);
		if(cond==n)
		  {
		    if(detail)cout<<"Curve is ";
		    cout << (Curve)CD << "  N = " << cond << "  ";
// Check periods were correct:
  unsigned int disagree=0;
  Cperiods cpC(CD);
  bigcomplex wRC, wRIC;
  cpC.getwRI(wRC, wRIC);
  int Ctype = get_lattice_type(cpC);
  bigfloat x1C, y1C;
  if(Ctype==1) {x1C=real(wRC)/to_bigfloat(2); y1C=imag(wRIC);}
  else {x1C=real(wRC); y1C=imag(wRIC);}
  if(type!=Ctype)
    {
      disagree|=1;
      cout<<"Period lattice type of constructed curve does not match"<<endl;
      cout<<"Lattice type of C: "<<Ctype<<endl;
      cout<<"Guessed type:      "<<type<<endl;
    }
  if((abs((x1-x1C)/x1C)>0.001))
    {
      disagree|=2;
      cout<<"Real period of constructed curve does not match that"
	  <<" of the newform"<<endl;
      cout<<"Real period of C: "<<x1C<<endl;
      cout<<"Real period of f: "<<x1<<endl;
      cout<<"Ratio = "<<(x1/x1C)<<endl;
    }
  if((abs((y1-y1C)/y1C)>0.001))
    {
      disagree|=4;
      cout<<"Imag period of constructed curve does not match"<<endl;
      cout<<"Imag part of second period of C: "<<y1C<<endl;
      cout<<"Guessed imag part for f:         "<<y1<<endl;
      cout<<"Ratio of imaginary parts = "<<(y1/y1C)<<endl;
    }
  if(disagree)
    {
      if(disagree&1)
	{
	  cout<<"Changing type to "<<Ctype<<endl;
	  type=Ctype;
	}
      if(disagree&2)
	{
	  cout<<"Changing real scaling from "<<nx;
	  nx = I2long(Iround(x0/x1C));
	  cout<<" to "<<(x0/x1C)<<" = "<<nx<<endl;
	}
      if(disagree&4)
	{
	  cout<<"Changing imag scaling from "<<ny;
	  ny = I2long(Iround(y0/y1C));
	  cout<<" to "<<(y0/y1C)<<" = "<<ny<<endl;
	}
    }
		    return 1;
		  }     // end of if(cond==n)
		else
		{
		    if(detail) cout<<"c4,c6 valid but conductor wrong, continuing..."<<endl;
			}
	      }         // end of if(validc4c6)	   
	  }             // end of type loop
      }                 // end of ny loop
    }                   // end of nx loop
  return 0;
}

int newforms::find_lminus(long i, long lmax, const bigfloat& y1)
{
  long nry, dry, ell;
  lfchi lx(this,&(nflist[i]));
  long mm=0;
  for(primevar l; ((l<lmax)||(lmax==0))&&(mm==0); l++)
    {
      ell = l; 
      if(ell%4!=3) continue;
      if(legendre(-modulus,ell)!=nflist[i].sfe) continue;  // skip this l
      lx.compute(ell);
      bigfloat y = abs(lx.scaled_value());
      if(verbose>1) cout<<"L(f,"<<ell<<",1) = "<<y<<"\n";
      if(y>0.001)
	{
	  nflist[i].lminus=ell;
	  bigfloat ratio = y/y1;
	  if(verbose>1) cout<<"ratio = "<<ratio<<endl;
	  ratapprox(ratio, nry, dry);
	  mm=nry;
	  if(dry!=1)
	    {
	      cout << "******************************L(f,"<<ell<<")/ip = "
		   <<ratio
		   <<" is not integral! (denom = "<<dry<<")"<<endl;
	      if(dry>10)
		{
		  mm=I2long(Iround(ratio));
		  cout << "Using rounded value mminus=" << mm <<endl;
		}
	    }
	  if(verbose>1) 
	    cout << "lminus = "<<ell<< "\tmminus = " << mm << "\n";
	  nflist[i].mminus=mm;
	  return 1;
	}
    } // end of primes loop
  return 0;
}

