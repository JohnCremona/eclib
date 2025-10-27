// file pcprocs.cc: implementation of functions used to compute periods
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

#include "eclib/periods.h"
#include "eclib/pcprocs.h"

// Function used to test whether a denominator found by ratapprox() is
// "trustworthy": always with multiprecision, but only if below
// a fixed bound otherwise.
#ifdef MPFP // Multi-Precision Floating Point
//inline int trust_denom(long d) { return (d<501);} // 1000);}
inline int trust_denom(long d) { return (d<1201);} // 1000);}
#else
inline int trust_denom(long d) { return (d<251);}
#endif

// Given a newform (the i'th in newforms) at level n with a real
// period x0; Finds a matrix [a0,b0;Nc0,d0] whose integral is
// dotplus*x0+dotminus*y0*i N.B. The value of x0 may be changed, if we
// come across some period whose real part is not an integer multiple
// of the original x0.

// We compute periods of lots of matrices in Gamma_0(N), over all
// symbols {0,b/d} for d<dmax.  We are looking for (i) a symbol whose
// imaginary part is nonzero; (ii) a symbol whose real and imaginary
// parts are both non-zero, which we store.

int newforms::find_matrix(long i, long dmax, int&rp_known, bigfloat&x0, bigfloat&y0)
{
  int have_both=0;
  int have_rp = get_real_period(i,x0,verbose);
  // NB The code below relies on x0 being *positive*
  int have_ip = 0;
  rp_known = have_rp;
  // have_rp is set if we know a real period; rp_known is set if we
  // know that it is the real period (strictly, the least real part of
  // a period), which should always be the case.

  // The code below allows for the situation where we do not know the
  // correct scaling factor for the real period, so that x0 will be an
  // unknown integer multiple of the minimal real period; in that
  // sitiuation we would have have_rp=1 but rp_known=0.

  // the following are to allow us a choice for the real period;
  // we'll use the value which has greater precision
  long lplus = nflist[i].lplus;
  if(nflist[i].dp0!=0) lplus=0;
  int rp_fixed = !have_rp;

  long nrx=1, nry=1, drx=1, dry=1;
  long nrx0=1, nry0=1, drx0=1, dry0=1;
  long a, b, b1, c, d;
  long& dotplus=(nflist[i].dotplus);
  long& dotminus=(nflist[i].dotminus);
  long nf_b = nflist[i].b;
  long nf_d = nflist[i].d;
  long dotplus0=1, dotminus0=1;
  periods_direct integrator(this,&(nflist[i]));

  for(d=nf_d; (d<dmax)||(!have_both); d++)
    {
      if(gcd(N,d)!=1) continue;
      long d2=d/2;
      b1 = 1;
      if(d==nf_d) b1=nf_b;
      //      for(b=b1; (b<=d2)&&(!have_both); b++)
      for(b=b1; (b<=d2); b++)
	{
	  if(bezout(d,N*b,a,c)!=1) continue;
	  c=-c;
	  if(verbose>1) 
	    cout << "Matrix ("<<a<<","<<b<<";"<<N*c<<","<<d<<"):\n";
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
             //             cout<<"dry = "<<dry<<", trusted? "<<trust_denom(dry)<<endl;
	     if(trust_denom(dry)) dotminus0=lcm(dotminus0,dry);
	     if(verbose>1) 
	       cout<<"imag part = " << y << ", y/y0 = " << ratio   
		   << " =~= "<<nry<<"/"<<dry<<"\n"<<"dotminus0 updated to "<<dotminus0<<endl;
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
             //              cout<<"nrx0="<<nrx0<<endl;
             //              cout<<"drx0="<<drx0<<endl;
	      nflist[i].a=a;
	      nflist[i].b=b;
	      nflist[i].c=c;
	      nflist[i].d=d;
	   }
       } // end of b loop
     } // end of d loop
  x0/=to_bigfloat(dotplus0);
  if (x0<0)
    {
      x0 *= -1;
      dotplus0 *= -1;
    }
  y0/=to_bigfloat(dotminus0);
  if (y0<0)
    {
      y0 *= -1;
      dotminus0 *= -1;
    }
  dotplus =(dotplus0 *nrx0)/drx0;
  dotminus=(dotminus0*nry0)/dry0;
  if(verbose>1){
    cout<<"dotplus0 ="<<dotplus0<<endl;
    cout<<"dotminus0="<<dotminus0<<endl;
    cout<<"dotplus  ="<<dotplus<<endl;
    cout<<"dotminus ="<<dotminus<<endl;
  }
  // if(!have_both) {a=d=1; b=c=0; dotplus=dotminus=0;}
  return have_both;
}

// Given a newform (the i'th in newforms) at level n ,with known data
// including a matrix [a0,b0;Nc0,d0] whose integral is
// dotplus*x0+dotminus*y0*i (where x0 and y0 are real and imaginary
// periods).  Computes both x0 and y0.  rp_known, ip_known are success
// flags.

int newforms::get_both_periods(long i, bigfloat&x0, bigfloat&y0) const
{
  x0=y0=to_bigfloat(0);
  if(nflist[i].a==0)
    {
      //      cout<<"Cannot compute get_periods(): matrix not known."<<<endl;
      return 0;
    }
  periods_direct integrator(this,&(nflist[i]));
  integrator.compute(nflist[i].a,nflist[i].b,nflist[i].c,nflist[i].d);
  int dot = nflist[i].dotplus;
  if (dot) // else we're in a minus space
    {
      x0 = integrator.rper();
      x0 /= to_bigfloat(dot);
    }
  dot = nflist[i].dotminus;
  if (dot) // else we're in a plus space
    {
      y0 = integrator.iper();
      y0 /= to_bigfloat(dot);
    }
  return 1;
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
	    normalize(w1,w2);
	    getc4c6(w1,w2,c4,c6);
            if (detail>1)
              {
                cout<<"w1 = "<<w1<<endl;
                cout<<"w2 = "<<w2<<endl;
                cout<<"type = "<<type<<endl;
                cout<<"c4 = "<<c4<<endl;
                cout<<"c6 = "<<c6<<endl;
              }
	    ZZ ic4 = fac*Iround(real(c4)/fac);
	    ZZ ic6 = fac6*Iround(real(c6)/fac6);
            if (detail>1)
              {
                cout<<"ic4 = "<<ic4<<endl;
                cout<<"ic6 = "<<ic6<<endl;
              }
	    int validc4c6 = 0;
            validc4c6 = valid_invariants(ic4,ic6);
            if((validc4c6)&&detail)
	      {
		cout << "type = " << type << ", nx = " << nx << ", ny = " << ny << "\n"; 
                //                cout << "x = " << x1 << ", y = " << y1 << "\n";
		cout << "w1 = " << w1 << ", w2 = " << w2 << "\n";
		cout << "c4 = " << real(c4) << ", c6 = " << real(c6) << "\n";
		cout << "ic4 = " << ic4 << ", ic6 = " << ic6 << "\n";
	      }
            else
              {
                if (detail>1) cout << "Invalid c-invariants" << endl;
              }

	    if(validc4c6)
	      {
		Curve C(ic4,ic6);
		Curvedata CD(C,1);
		CurveRed CR(CD);
		ZZ cond = getconductor(CR);
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
  long nry, dry;
  lfchi lx(this,&(nflist[i]));
  long mm=0;
  for(primevar l; ((l<lmax)||(lmax==0))&&(mm==0); l++)
    {
      long ell = l;
      if(ell%4!=3) continue;
      if(legendre(-N,ell)!=nflist[i].sfe) continue;  // skip this l
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
              if (verbose>1)
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
	    }
	  if(verbose>1) 
	    cout << "lminus = "<<ell<< "\tmminus = " << mm << "\n";
	  nflist[i].mminus=mm;
	  return 1;
	}
    } // end of primes loop
  return 0;
}

