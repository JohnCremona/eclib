// mlocsol.cc: implementation of functions for local solubility of quartics 
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
 
//#define DEBUG_NEW_LOCSOL

#include "eclib/mlocsol.h"
#include "eclib/hilbert.h"

using NTL::vec_pair_ZZ_pX_long;
using NTL::pair_ZZ_pX_long;
using NTL::ZZ_pX;

int psquare(const ZZ& aa, const ZZ& p);  
/* tests if aa is a p-adic square */

int lemma6(const ZZ& a,const ZZ& b,const ZZ& c,const ZZ& d,const ZZ& e, 
	   const ZZ& p, int nu,const ZZ& x);

/* returns -1 for insoluble, 0 for undecided, +1 for soluble --- odd p */

int lemma7(const ZZ& a,const ZZ& b,const ZZ& c,const ZZ& d,const ZZ& e, 
	   const ZZ& p, int nu, const ZZ& x);

/* returns -1 for insoluble, 0 for undecided, +1 for soluble --- p=2 */

int zpsoluble(const ZZ& a,const ZZ& b,const ZZ& c,const ZZ& d,const ZZ& e, 
	      const ZZ& p, int pzp);

/* Checks for solublility in Zp, or in pZp if "pzp" is 1 */

int psquare(const ZZ& aa, const ZZ& p)  /* tests if aa is a p-adic square */
{
  if (is_zero(aa)) return 1;
  long v = val(p,aa);
  if (odd(v)) return 0;
  ZZ a(aa); while(v--) a/=p;
  if(p==2) return posmod(a,8)==1;
  else     return legendre(a,p) == 1 ;
}  /* of psquare */
 

int lemma6(const ZZ& a,const ZZ& b,const ZZ& c,const ZZ& d,const ZZ& e, 
	   const ZZ& p, int nu, const ZZ& x)
// returns -1 for insoluble, 0 for undecided, +1 for soluble --- odd p 
{
   ZZ gx = (((a*x+b)*x+c)*x+d)*x+e;
   if (psquare(gx,p)) return +1;
   ZZ gdashx = ((4*a*x+3*b)*x+2*c)*x+d;
   long lambda = val(p,gx);
   if(is_zero(gdashx))
// then effectively mu = infinity
     {
       if (lambda >= 2*nu) return 0;
       return -1;
     }
// now gdashx !=0:     
   long mu = val(p,gdashx);
   if ((lambda-mu >= nu) && (nu >  mu)) return +1;
   if ((lambda >= 2*nu)  && (mu >= nu)) return 0;
   return -1;
}  /* end of lemma6 */

int lemma7(const ZZ& a,const ZZ& b,const ZZ& c,const ZZ& d,const ZZ& e, 
	   const ZZ& p, int nu, const ZZ& x)
// returns -1 for insoluble, 0 for undecided, +1 for soluble --- p=2 
{
   ZZ gx = (((a*x+b)*x+c)*x+d)*x+e;
   if (psquare(gx,p)) return +1;
   ZZ gdashx = ((4*a*x+3*b)*x+2*c)*x+d;
   long lambda = val(p,gx);
   ZZ oddgx = gx;
   if (oddgx==0) oddgx= 1;
      else while (even(oddgx)) oddgx /= 2;
   int odd4 = (posmod(oddgx,4)==1);
   if(is_zero(gdashx)) // mu = infinity
     {
       if (lambda >= 2*nu) return 0;
       if ((lambda == 2*nu-2) && odd4) return 0;
       return -1;
     }
   long mu = val(p,gdashx); // Now gdashx is nonzero
   if ((lambda-mu >= nu) && (nu >  mu)) return +1;
   if ((nu > mu)  && (lambda==mu+nu-1) && even(lambda)) return +1;
   if ((nu > mu)  && (lambda==mu+nu-2) && even(lambda) && odd4) return +1;
   if ((mu >= nu) && (lambda >= 2*nu)) return 0;
   if ((mu >= nu) && (lambda == 2*nu-2) && odd4) return 0;
   return -1;
}  /* end of lemma7 */

int zpsol(const ZZ& a,const ZZ& b,const ZZ& c,const ZZ& d,const ZZ& e, 
	  const ZZ& p, const ZZ& x0, long nu)
// Checks for solublility in Zp with x=x0 (mod p^nu)
// Fully recursive (depth-first) version
{
  // cout << "In zpsol with (p,nu) = (" << p << "," << nu << "), x0="<<x0<<"\n";

 long result =  (p==2) ? lemma7(a,b,c,d,e,p,nu, x0)
                      : lemma6(a,b,c,d,e,p,nu, x0);
 if(result==+1) return 1;
 if(result==-1) return 0;
//else result==0, so refine to look modulo p^(nu+1):
// But before we do a depth-first recursive search, we try all residues 
 ZZ i, x=x0, pnu=pow(p,nu);
 if(nu==0)
   {
     for(i=0; i<p; ++i, x+=pnu)
       {
	 if(p==2) if(lemma7(a,b,c,d,e,p,nu+1, x)==+1) return 1;
	 if(lemma6(a,b,c,d,e,p,nu+1, x)==+1) return 1;
       }
   }

 for(i=0; i<p; ++i, x+=pnu)
  {
    if(zpsol(a,b,c,d,e,p,x,nu+1)) return 1;
  }
 return 0;
}

int qpsoluble(const quartic& g, const ZZ& p)
{ 
  static const ZZ zero(0);
  ZZ a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  if (zpsol(a,b,c,d,e,p,zero,0)) return 1;
  else return zpsol(e,d,c,b,a,p,zero,1);
} /* end of qpsoluble */

int qpsoluble(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d, 
	      const ZZ& e, const ZZ& p)
{ 
  static const ZZ zero(0);
  if (zpsol(a,b,c,d,e,p,zero,0)) return 1;
  else return zpsol(e,d,c,b,a,p,zero,1);
} /* end of qpsoluble */

int qpsoluble(const ZZ& a, const ZZ& c, const ZZ& e, const ZZ& p)
{ 
  static const ZZ zero(0);
  if (zpsol(a,zero,c,zero,e,p,zero,0)) return 1;
  else return zpsol(e,zero,c,zero,a,p,zero,1);
} /* end of qpsoluble */

int Rsoluble(const quartic& g)
{
  return ((g.gettype()>1)|| is_positive(g.geta()));
}

int Rsoluble(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,  
	      const ZZ& e)
{
  quartic g(a,b,c,d,e);
  return Rsoluble(g);
}


//#define CHECK_LOC_SOL 1

int locallysoluble(const quartic& g, const vector<ZZ>& plist, ZZ& badp)
{
  ZZ a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  return locallysoluble(a,b,c,d,e,plist,badp);
}

int locallysoluble(const ZZ& a, const ZZ& c, const ZZ& e, 
		   const vector<ZZ>& plist, ZZ& badp)
{
  static const ZZ zero(0);
  ZZ d = c*c-4*a*e;
  if(global_hilbert(a,d,plist,badp)) return 0;
  return locallysoluble(a,zero,c,zero,e,plist,badp);
}

int locallysoluble(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,
		   const ZZ& e, const vector<ZZ>& plist, ZZ& badp)
{
  // First check R-solubility
  if (!Rsoluble(a,b,c,d,e))
    {
      badp = ZZ(0);
      return 0;
    }
  if(is_zero(b)&&is_zero(d)) // do a quick Hilbert check:
    if(global_hilbert(a,c*c-4*a*e,plist,badp))
      return 0;

  return std::all_of(plist.begin(), plist.end(),
                     [&badp, a,b,c,d,e] (const ZZ& p) {badp=p; return new_qpsoluble(a,b,c,d,e,p,0);});
}  /* end of locallysoluble */

/* Samir Siksek's Local Solubility Test for odd p */

// The following is the cross-over between the B&SD method & the SS method
// i.e. only uses SS method for p greater than this

//#define CROSS_OVER_P 10000  // for testing only!
#define CROSS_OVER_P 1000
//#define CROSS_OVER_P 100
//#define CROSS_OVER_P 7 // for testing only!

int new_qpsoluble(const quartic& g, const ZZ& p, int verbose)
{
  ZZ a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  return new_qpsoluble(a,b,c,d,e,p,verbose);
}

int new_qpsoluble_ace(const ZZ& a, const ZZ& c, const ZZ& e, 
		      const ZZ& p, int verbose)
{
  ZZ b; b=0;
  return new_qpsoluble(a,b,c,b,e,p,verbose);
}

int new_qpsoluble(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d, 
			const ZZ& e, const ZZ& p, int verbose)
{
  int verb;
#ifdef DEBUG_NEW_LOCSOL
  verb=1;
#else
  verb=verbose;
#endif
  if(p<CROSS_OVER_P) 
    {
      if(verb) 
	{
	  cout<<"new_qpsoluble with p<"<<CROSS_OVER_P;
	  cout<<" passing to old qpsoluble.\n";
	}
      return qpsoluble(a,b,c,d,e,p);
    }
  if(verb) 
    {
      cout<<"Using new_qpsoluble with p = "<<p<<endl;
    }
  if (new_zpsol(a,b,c,d,e,p,verbose)) return 1;
  else return new_zpsol(e,d,c,b,a,p,verbose);
} /* end of new_qpsoluble */

int new_zpsol(const ZZ& a,const ZZ& b,const ZZ& c,const ZZ& d,const ZZ& e, const ZZ& p, int verbose)
{
  vector<ZZ> coeff = {a, b, c, d, e};
  return local_sol(p,coeff,verbose);
}

/* Samir's Local Solubility Test for odd p */

vec_pair_ZZ_pX_long fact_c(vector<ZZ> c, int verbose=0)
{
  ZZ_pX f;
  for (long i=0; i<5; i++)
    SetCoeff(f,i,to_ZZ_p(c[i]));
  if(verbose) cout<<"Factorizing "<<f<<" after making monic: ";
  MakeMonic(f);
  if(verbose) cout<<f<<endl;
  verbose=0;
  return  CanZass(f, verbose);
}

int local_sol(const ZZ& p, vector<ZZ> c, int verbose)
{
  ZZ_p::init(p);
  if (verbose)
    {  cout << "---------------------------------------------\n";
       cout << "LOCAL_SOL \n";
       cout << c[4] << " " << c[3] << " " << c[2] << " " << c[1] << " ";
       cout << c[0] << "      p=" << p << endl;
    }

  ZZ r[2],t;
  long fl,i;
  ZZ p2=sqr(p);
  pair_ZZ_pX_long term;
  ZZ_pX F;

  int zeromodp=1;
  for (i=0; (i<5) && zeromodp; i++) {zeromodp=div(p,c[i]);}
  if (zeromodp)
    { // Case II (Case I is below)
      if (verbose) cout << "f is 0 mod p: Case II" << endl;
      zeromodp=1;
      for (i=0; (i<5) && zeromodp; i++) { zeromodp=div(p2,c[i]); }
      vector<ZZ> dd(5);
      if (zeromodp)
         { for (i=0; i<5; i++)   { dd[i]=c[i]/p2; }
	   if (verbose) cout << "f is 0 mod p^2, recursing" << endl;
           return local_sol(p,dd,verbose);
         }
      for (i=0; i<5; i++) { dd[i]=c[i]/p; }
      vec_pair_ZZ_pX_long fact_f=fact_c(dd);
      // Is there a non-repeated root
      if (verbose) cout << "Factorization of f/p = "<<fact_f << endl;
      for (i=0; i<fact_f.length(); i++)
        { term=fact_f[i];
	  F=term.a;
	  long e=term.b;
          if ((deg(F)==1) && (e==1))
	    { if (verbose)
	      cout << "Non-Repeated Root " << -ConstTerm(F)
		   <<", returning 1"<<endl;
              return 1;
	    }
        }

      /* Go thru each repeated root and make the
         required transformation
      */
      vector<ZZ> d(5);
      fl=0;
      for (i=0; i<fact_f.length() && (!fl); i++)
        { term=fact_f[i];
          F=term.a;
          long e=term.b;
          if ((deg(F)==1) && (e!=1))
            { ZZ rt=-rep(ConstTerm(F));
              if (verbose)
		cout << "Repeated Root=" << rt << ", recursing"<<endl;
 // Using f(pX+r)/p^2
              d[4]=dd[4]*p2*p;
              d[3]=p2*(dd[3]+4*dd[4]*rt);
              d[2]=p*(dd[2]+6*dd[4]*rt*rt+3*dd[3]*rt);
              d[1]=(2*dd[2]*rt+4*dd[4]*rt*rt*rt+3*dd[3]*rt*rt+dd[1]);
              d[0]=((((dd[4]*rt+dd[3])*rt+dd[2])*rt+dd[1])*rt+dd[0])/p;
              fl=local_sol(p,d,verbose);
            }
        }
      return fl;
    }
// CASE I
  if (verbose) cout << "f is not 0 mod p: Case I" << endl;
  ZZ unit;
  zeromodp=1;
  for (i=4; i>=0 && zeromodp; i--)
     { unit=c[i];
       zeromodp=div(p,unit);
     }
  // If leading non-zero term is a square return 1
  if (legendre(unit,p)==1)      { return 1; }
  // If f is a constant mod p and constant not a square return 0
  if (i==-1) { return 0; }

  // Factorize f

  vec_pair_ZZ_pX_long fact_f=fact_c(c);
  long nc=fact_f.length();
  // Check if of the form unit*g^2
  if (verbose) cout << "Factorization of f = "<< fact_f << endl;
  for (i=0; i<nc; i++)
    { if ((fact_f[i].b)%2!=0) { return 1; } }
  if (verbose) cout << "Roots mod p: ";
  // Compute roots mod p
  long num_r=0;
  for (i=0; i<nc; i++)
     { F=fact_f[i].a;
       if (deg(F)==1)
         { r[num_r]=-rep(ConstTerm(F));
           if (verbose) cout << r[num_r] << "\t";
           num_r=num_r+1;
         }
     }
  if (num_r==0) 
    { 
      if(verbose) cout<<"none, returning 0"<<endl;
      return 0; 
    }
  if(verbose) cout<<endl;

  // Now need to determine g s.t. f=unit*g^2; here g is called F
  if (nc==1)
    { F=fact_f[0].a;
      if ((fact_f[0].b)%4==0) { F=F*F; }
    }
  else
    { 
      F=fact_f[0].a*fact_f[1].a;
    }
  if(verbose)
    cout<<"f = unit*g^2 where unit="<<unit<<", g="<<F<<endl;
  ZZ te,g0,g1,g2;
  g2=rep(coeff(F,0)); // lift to ZZ
  g1=rep(coeff(F,1));
  g0=0;
  if (deg(F)==2) { g0=1; }
  if (verbose) 
    cout << "g = " << g0 << " " << g1 << " " << g2 << endl;
  // Now determine h
  ZZ h0,h1,h2,h3,h4;
  h4=(c[4]-unit*g0*g0)/p;
  h3=(c[3]-2*unit*g1*g0)/p;
  h2=(c[2]-2*unit*g2*g0-unit*g1*g1)/p;
  h1=(c[1]-2*unit*g1*g2)/p;
  h0=(c[0]-unit*g2*g2)/p;
  if (verbose)
    { cout << "h =" << h4 << " " << h3 << " " 
	   << h2 << " " << h1 << " " << h0 << endl;
    }
  /* For each root which is also a root of h
     transform the equations and call again
  */
  vector<ZZ> d(5);
  fl=0;
  for (i=0; i<num_r && !fl; i++)
    { te=(((h4*r[i]+h3)*r[i]+h2)*r[i])%p;
      te=((te+h1)*r[i]+h0)%p;
      if (is_zero(te))
        { if (verbose)
        	  cout << "Using " << r[i] << " and recursing"<<endl;
          d[4]=c[4]*p2;
          d[3]=p*(c[3]+4*c[4]*r[i]);
          d[2]=(c[2]+6*c[4]*r[i]*r[i]+3*c[3]*r[i]);
          d[1]=(2*c[2]*r[i]+4*c[4]*r[i]*r[i]*r[i]+3*c[3]*r[i]*r[i]+c[1])/p;
          d[0]=((((c[4]*r[i]+c[3])*r[i]+c[2])*r[i]+c[1])*r[i]+c[0])/(p2);
          fl=local_sol(p,d,verbose);
        }
    }
  return fl;
}


/* END OF FILE MLOCSOL.CC */
