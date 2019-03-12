// mlocsol.cc: implementation of functions for local solubility of quartics
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

#include <NTL/ZZ_pXFactoring.h>

#include <eclib/mlocsol.h>
#include <eclib/hilbert.h>

int psquare(const bigint& aa, const bigint& p);
/* tests if aa is a p-adic square */

int lemma6(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e,
	   const bigint& p, int nu, bigint& x);

/* returns -1 for insoluble, 0 for undecided, +1 for soluble --- odd p */

int lemma7(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e,
	   const bigint& p, int nu, bigint& x);

/* returns -1 for insoluble, 0 for undecided, +1 for soluble --- p=2 */

int zpsol(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e,
	  const bigint& p, bigint& x0, long nu, int s_=0, const quadratic& q=q0_);


int psquare(const bigint& aa, const bigint& p)  /* tests if aa is a p-adic square */
{
  if (is_zero(aa)) return 1;
  long v = val(p,aa);
  if (odd(v)) return 0;
  bigint a(aa); while(v--) a/=p;
  if(p==2) return posmod(a,8)==1;
  else     return legendre(a,p) == 1 ;
}  /* of psquare */


int lemma6(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e,
	   const bigint& p, int nu, bigint& x)
// returns -1 for insoluble, 0 for undecided, +1 for soluble --- odd p
{
   bigint gx = (((a*x+b)*x+c)*x+d)*x+e;
   if (psquare(gx,p)) {return +1;}
   bigint gdashx = ((4*a*x+3*b)*x+2*c)*x+d;
   long lambda = val(p,gx);
   if(is_zero(gdashx))
// then effectively mu = infinity
     {
       if (lambda >= 2*nu) return 0;
       return -1;
     }
// now gdashx !=0:
   long mu = val(p,gdashx);
   if ((lambda-mu >= nu) && (nu >  mu)) return 1;
   if ((lambda >= 2*nu)  && (mu >= nu)) return 0;
   return -1;
}  /* end of lemma6 */

int lemma7(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e,
	   const bigint& p, int nu, bigint& x)
// returns -1 for insoluble, 0 for undecided, +1 for soluble --- p=2
{
	if (nu==0) return 0;
   bigint gx = (((a*x+b)*x+c)*x+d)*x+e;
   if (psquare(gx,p)) {return +1;}
   bigint gdashx = ((4*a*x+3*b)*x+2*c)*x+d;
   long lambda = val(p,gx);
   bigint oddgx = gx;
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
   if ((lambda-mu >= nu) && (nu >  mu)) return 1;
   if ((nu > mu)  && (lambda==mu+nu-1) && even(lambda)) return 1;
   if ((nu > mu)  && (lambda==mu+nu-2) && even(lambda) && odd4) return 1;
   if ((mu >= nu) && (lambda >= 2*nu)) return 0;
   if ((mu >= nu) && (lambda == 2*nu-2) && odd4) return 0;
   return -1;
}  /* end of lemma7 */


int zpsol(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e,
	  const bigint& p, bigint& x0, long nu, int s_, const quadratic& q)
// Checks for solublility in Zp with x=x0 (mod p^nu);
// If it's and s_!=0, then (x0,1) is solution with val(p,q(x0,1))<=nu if s_=1 or
//                         (1,x0) is solution with val(p,q(1,x0))<=nu if s_=2;
// Fully recursive (depth-first) version
{
  // cout << "In zpsol with (p,nu) = (" << p << "," << nu << "), x0="<<x0<<"\n";
 int quit=0;
 long result =  (p==2) ? lemma7(a,b,c,d,e,p,nu, x0): lemma6(a,b,c,d,e,p,nu, x0);

 if(result==+1)
	{
		if (nu<20) return -1; //Avoid infinity loop case q and g has common solution;
		if (s_==1)
		{
			if (val(p,q(x0,BIGINT(1)))<=nu) return 1;////////////////////////////val(q)<nu
			result=0;
			quit=1;
		}
		else if (s_==2)
		{
			if (val(p,q(BIGINT(1),x0))<=nu) return 1;///////////////////////////
			result=0;
			quit=1;
		}
		else
		{
			return 1;
		}
	}
 if(result==-1) return 0;
//else result==0, so refine to look modulo p^(nu+1):
// But before we do a depth-first recursive search, we try all residues

 if (quit)
 {//x0 is solution but val(q(x0))>nu
	bigint i, xx, pnu=pow(p,nu), xxp, pnup;
	pnup=pnu*p;
	xx=x0+pnu;
	xxp=x0+pnup;
	for (i=1; i<p; i++, xx+=pnu, xxp+=pnup)
	{
		if ((zpsol(a,b,c,d,e,p,xx,nu+1,s_,q)))
		{
			x0 = xx;
			return 1;
		}
		if (zpsol(a,b,c,d,e,p,xxp,nu+2,s_,q))
		{
			x0 = xxp;
			return 1;
		}
	}
	return 0;
  }

 bigint i, pnu=pow(p,nu), x=x0;
 for(i=0; i<p; ++i, x+=pnu)
  {
      if(zpsol(a,b,c,d,e,p,x,nu+1,s_,q))
      {
		  x0 = x;
		  return 1;
	  }
  }
 return 0;
}

int qpsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
	      const bigint& e, const bigint& p, int s_, vector< vector<bigint> >& xplist, const quadratic& q)
{
//  static const bigint zero = BIGINT(0);
  if (s_) {s_=1;}
  bigint xp=BIGINT(0);
  if (zpsol(a,b,c,d,e,p,xp,0,s_,q))
   {
     if (s_)
     {
       vector<bigint> vp{xp, BIGINT(1)};
       xplist.push_back(vp);
     }
	  return 1;
   }
  else
  {
	  xp=BIGINT(0);
	  if (s_) {s_=2;}
	  if (zpsol(e,d,c,b,a,p,xp,1,s_,q))
	  {
      if (s_)
      {
        vector<bigint> vp{BIGINT(1), xp};
	      xplist.push_back(vp);
      }
	    return 1;
	  }
	  return 0;
  }
} /* end of qpsoluble */

int qpsoluble(const quartic& g, const bigint& p, int s_, vector< vector<bigint> >& xplist, const quadratic& q)
{
  bigint a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  return qpsoluble(a,b,c,d,e,p,s_,xplist,q);
} /* end of qpsoluble */

int qpsoluble(const bigint& a, const bigint& c, const bigint& e, const bigint& p, int s_, vector< vector<bigint> >& xplist, const quadratic& q)
{
  static const bigint zero = BIGINT(0);
  return qpsoluble(a,zero,c,zero,e,p,s_,xplist,q);
} /* end of qpsoluble */

int Rsoluble(const quartic& g, int s_, vector< vector<bigint> >& xplist, const quadratic& q)
{
  int R_soluble=((g.gettype()>1)||(g.geta()>BIGINT(0)));
  if (!s_)
  {
    return R_soluble;
  }
  if (R_soluble)
  {//Need at least 3 solutions so that q!=0 in one of these solutions;
	  int i;
	  vector<bigfloat> realroots_;
    vector<bigcomplex> complexroots = g.getroots();
	  for (i=0; i<4; i++)
	  {
		  if (is_real(complexroots[i]))
		  {
			  realroots_.push_back(real(complexroots[i]));
		  }
	  }
	  sort(realroots_.begin(), realroots_.end());
    vector<bigint> vp;
	  if (realroots_.size()>0)
	  {
		  if (g.geta() > 0)
		  {
        if (q(BIGINT(realroots_[0]-1),BIGINT(1))!=BIGINT(0)) {vp={BIGINT(realroots_[0]-1),BIGINT(1)};}
        else if (q(BIGINT(realroots_[0]-2),BIGINT(1))!=BIGINT(0)) {vp={BIGINT(realroots_[0]-2),BIGINT(1)};}
		    else {vp={BIGINT(realroots_[0]-3),BIGINT(1)};}
			  xplist.push_back(vp);
		  }
		  else
		  {
			  bigint xp_0, yp_0;
			  ratapprox((realroots_[0] + realroots_[1])/2, xp_0, yp_0);
        if (q(xp_0,yp_0)!=0) {vp={xp_0,yp_0};}
        else
        {
          ratapprox((realroots_[0]*2)/3 + realroots_[1]/3, xp_0, yp_0);
          if (q(xp_0,yp_0)!=0) {vp={xp_0,yp_0};}
          else
          {
            ratapprox(realroots_[0]/3 + (realroots_[1]*2)/3, xp_0, yp_0);
            vp={xp_0,yp_0};
          }
        }
			  xplist.push_back(vp);
		  }
	  }
	  else
	  {
      if (q(BIGINT(1), BIGINT(0))!=BIGINT(0)) {vp={BIGINT(1), BIGINT(0)};}
      else if (q(BIGINT(0), BIGINT(1))!=BIGINT(0)) {vp={BIGINT(0), BIGINT(1)};}
		  else {vp={BIGINT(1), BIGINT(1)};}
		  xplist.push_back(vp);
	  }
	  return 1;
  }
  return 0;
}

int Rsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
	      const bigint& e, int s_, vector< vector<bigint> >& xplist, const quadratic& q)
{
  quartic g(a,b,c,d,e);
  return Rsoluble(g,s_,xplist,q);
}


//#define CHECK_LOC_SOL 1

int locallysoluble(const quartic& g, const vector<bigint>& plist, bigint& badp, int s_, vector< vector<bigint> >& xplist, const quadratic& q)
{
  bigint a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  return locallysoluble(a,b,c,d,e,plist,badp, s_, xplist, q);
}

int locallysoluble(const bigint& a, const bigint& c, const bigint& e,
		   const vector<bigint>& plist, bigint& badp, int s_, vector< vector<bigint> >& xplist, const quadratic& q)
{
  static const bigint zero = BIGINT(0);
  bigint d = c*c-4*a*e;
//   cout<<"In locallysoluble(a,c,e) with plist = "<<plist<<endl;
//   cout<<"(a,c,e)=("<<a<<","<<c<<","<<e<<")\n";
//   cout<<"d="<<d<<endl;
//   int h = global_hilbert(a,d,plist,badp);
//   cout<<"global_hilbert() returns "<<h<<endl;
//   if(h) {cout<<"badp="<<badp<<endl; return 0;}
  if(global_hilbert(a,d,plist,badp)) return 0;
  return locallysoluble(a,zero,c,zero,e,plist,badp, s_, xplist, q);
}

int locallysoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d, 
                    const bigint& e, const vector<bigint>& plist, bigint& badp, int s_, 
                    vector< vector<bigint> >& xplist, const quadratic& q)
{
  //  cout<<"("<<a<<","<<b<<","<<c<<","<<d<<","<<e<<")\n";
  xplist=v00_;
  // First check R-solubility
  if (!Rsoluble(a,b,c,d,e,s_,xplist,q))
    {
      badp = BIGINT(0);
      return 0;
    }
  int sol = 1;

  if(is_zero(b)&&is_zero(d)) // do a quick Hilbert check:
    {
      bigint d = c*c-4*a*e;
      if(global_hilbert(a,d,plist,badp)) return 0;
    }
  for (vector<bigint>::const_iterator p = plist.begin(); (p!=plist.end())&&sol; p++)
   {
     badp = *p;
     if (s_) {sol = qpsoluble(a,b,c,d,e,badp,s_,xplist,q);}
     else {sol = new_qpsoluble(a,b,c,d,e,badp,0);}
     //Need get solutions(xplist) for new_qpsoluble case;
#ifdef CHECK_LOC_SOL
     cout<<"Checking solubility with old method (p="<<p<<")...\t";
     if(sol != qpsoluble(a,b,c,d,e,badp))
       {
	 cout<<"\nProblem with local solubility of ("<<a<<","<<b<<","<<c<<","<<d<<","
	     <<e<<") mod "<<badp<<"!\n";
	 cout<<"New method returns "<<sol<<", old method disagrees!\n";
	 sol=!sol;
       }
     else cout << "...OK.\n";
#endif
   }
   return sol;
}  /* end of locallysoluble */


/* Samir Siksek's Local Solubility Test for odd p */

// The following is the cross-over between the B&SD method & the SS method
// i.e. only uses SS method for p greater than this

//#define CROSS_OVER_P 10000  // for testing only!
#define CROSS_OVER_P 1000
//#define CROSS_OVER_P 100
//#define CROSS_OVER_P 7 // for testing only!

int new_qpsoluble(const quartic& g, const bigint& p, int verbose)
{
  bigint a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  return new_qpsoluble(a,b,c,d,e,p,verbose);
}

int new_qpsoluble_ace(const bigint& a, const bigint& c, const bigint& e,
		      const bigint& p, int verbose)
{
  bigint b; b=0;
  return new_qpsoluble(a,b,c,b,e,p,verbose);
}

int new_qpsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
			const bigint& e, const bigint& p, int verbose)
{
  int verb=verbose;
#ifdef DEBUG_NEW_LOCSOL
  verb=1;
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

int new_zpsol(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e, const bigint& p, int verbose)
{
  bigint *coeff =new bigint[5];
  coeff[0]=a;  coeff[1]=b;  coeff[2]=c;  coeff[3]=d;  coeff[4]=e;
  int res=local_sol(p,coeff,verbose);
  delete[] coeff;
  return res;
}

/* Samir's Local Solubility Test for odd p */

#define Factorization vec_pair_ZZ_pX_long
#define Term pair_ZZ_pX_long
#define Poly ZZ_pX

Factorization fact_c(bigint *c, int verbose=0)
{
  Poly f; ZZ_p ci;
  for (long i=0; i<5; i++)
    { ci=to_ZZ_p(c[i]); SetCoeff(f,i,ci); }
  if(verbose) cout<<"Factorizing "<<f<<" after making monic: ";
  MakeMonic(f);
  if(verbose) cout<<f<<endl;
  verbose=0;
  return  CanZass(f, verbose);
}

int local_sol(const bigint& p, bigint *c, int verbose)
{
  ZZ_p::init(p);
  if (verbose)
    {  cout << "---------------------------------------------\n";
       cout << "LOCAL_SOL \n";
       cout << c[4] << " " << c[3] << " " << c[2] << " " << c[1] << " ";
       cout << c[0] << "      p=" << p << endl;
    }

  bigint r[2],t;
  long e,fl,i;
  bigint p2=sqr(p);
  Term term;
  Poly F;

  int zeromodp=1;
  for (i=0; (i<5) && zeromodp; i++) {zeromodp=div(p,c[i]);}
  if (zeromodp)
    { // Case II (Case I is below)
      if (verbose) cout << "f is 0 mod p: Case II" << endl;
      zeromodp=1;
      for (i=0; (i<5) && zeromodp; i++) { zeromodp=div(p2,c[i]); }
      bigint *dd=new bigint[5];
      if (zeromodp)
         { for (i=0; i<5; i++)   { dd[i]=c[i]/p2; }
	   if (verbose) cout << "f is 0 mod p^2, recursing" << endl;
           int res=local_sol(p,dd,verbose);
           delete[] dd;
           return res;
         }
      for (i=0; i<5; i++) { dd[i]=c[i]/p; }
      Factorization fact_f=fact_c(dd);
      // Is there a non-repeated root
      if (verbose) cout << "Factorization of f/p = "<<fact_f << endl;
      for (i=0; i<fact_f.length(); i++)
        { term=fact_f[i];
	  F=term.a;
	  e=term.b;
          if ((deg(F)==1) && (e==1))
	    { if (verbose)
	      cout << "Non-Repeated Root " << -ConstTerm(F)
		   <<", returning 1"<<endl;
              delete[] dd;
              return 1;
	    }
        }

      /* Go thru each repeated root and make the
         required transformation
      */
      bigint *d=new bigint[5];
      fl=0;
      for (i=0; i<fact_f.length() && (!fl); i++)
        { term=fact_f[i];
          F=term.a;
          e=term.b;
          if ((deg(F)==1) && (e!=1))
            { bigint r=-rep(ConstTerm(F));
              if (verbose)
		cout << "Repeated Root=" << r << ", recursing"<<endl;
 // Using f(pX+r)/p^2
              d[4]=dd[4]*p2*p;
              d[3]=p2*(dd[3]+4*dd[4]*r);
              d[2]=p*(dd[2]+6*dd[4]*r*r+3*dd[3]*r);
              d[1]=(2*dd[2]*r+4*dd[4]*r*r*r+3*dd[3]*r*r+dd[1]);
              d[0]=((((dd[4]*r+dd[3])*r+dd[2])*r+dd[1])*r+dd[0])/p;
              fl=local_sol(p,d,verbose);
            }
        }
      delete[] d;
      delete [] dd;
      return fl;
    }
// CASE I
  if (verbose) cout << "f is not 0 mod p: Case I" << endl;
  bigint unit;
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

  Factorization fact_f=fact_c(c);
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
  bigint te,g0,g1,g2;
  g2=rep(coeff(F,0)); // lift to ZZ
  g1=rep(coeff(F,1));
  g0=0;
  if (deg(F)==2) { g0=1; }
  if (verbose)
    cout << "g = " << g0 << " " << g1 << " " << g2 << endl;
  // Now determine h
  bigint h0,h1,h2,h3,h4;
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
  bigint *d=new bigint[5];
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

  delete[] d;
  return fl;
}


/* END OF FILE MLOCSOL.CC */
