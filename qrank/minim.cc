// minim.h: implementation of quartic minimization functions
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
 
#include "marith.h"
#include "unimod.h"
#include "points.h"
#include "mquartic.h"
#include "transform.h"
#include "msoluble.h"
#include "minim.h"

//#define DEBUG_MINIM

bigint root_p(const bigint& a, const bigint& b, const bigint& c, 
		 const bigint& d, const bigint& e, const bigint& p)
// assuming p|I, p|J, returns the unique alpha mod p 
// modulo which quartic has a root of multiplicity at least 3
// returns -1 if multiple root is at infinity (if a=b=0 mod p)
// (program does not actaully use this dubious feature)
{
  if(div(p,a)&&div(p,b)) return BIGINT(-1);
  if(div(p,e)&&div(p,d)) return BIGINT(0);
  
  // Now we have to find the multiple root, using invariant theory:
  if(p==2)
    {
      return BIGINT(1);  // the only other possibility
    }
  if(p==3)
    {
      if(div(p,a)) return mod(-b*e,p);
      else         return mod(-a*d,p);
    }
  bigint b2=sqr(b);
  bigint ac=a*c;
  bigint p_seminv = mod(3*b2-8*ac,p);
  if(is_zero(p_seminv)) // quadruple root
    {
      return mod(-b*invmod(4*a,p),p);
    }
  else // triple root only
    {
      if(div(p,a)) // fourth root is at infinity
	{
	  return mod(c*invmod(3*b,p),p);
	}
      bigint t=invmod(4*a*p_seminv,p);
      bigint r_seminv = b*b2+8*sqr(a)*d-4*ac*b;
      bigint alpha = mod((3*r_seminv-b*p_seminv)*t,p);
      return alpha;
    }
}

int minim_p(bigint& a, bigint& b, bigint& c, 
	      bigint& d, bigint& e, const bigint& p,
	      scaled_unimod& m)
// assuming p^4|I, p^6|J, (or stronger conditions when p=2 or p=3)
// returns an equivalent quartic with invariants divided by p^4, p^6;
// m holds the transformation matrix, must be initialized (say with identity)
// returns success, can be 0 only for p=2
{
  bigint a0,b0,c0,d0,e0,r;
  bigint p2=sqr(p), temp;

  int two   = (p==2);
  int three = (p==3);

  // First test for trivial case where p^2 divides all coeffs:

  int p2divall=::divides(a,p2,a0,r);
  if(p2divall) {p2divall=::divides(b,p2,b0,r);}
  if(p2divall) {p2divall=::divides(c,p2,c0,r);}
  if(p2divall) {p2divall=::divides(d,p2,d0,r);}
  if(p2divall) {p2divall=::divides(e,p2,e0,r);}

  if(p2divall) // trivial case, all coeffs divisible by p^2
    {
#ifdef DEBUG_MINIM
      cout<<"All coeffs divisible by p\n";
#endif
      a=a0; b=b0; c=c0; d=d0; e=e0; 
      m.u_scale(p);
      return 1;
    }

  // Next test for case where p divides all coeffs:

  int pdivall=::divides(a,p,a0,r);
  if(pdivall) { pdivall=::divides(b,p,b0,r);}
  if(pdivall) { pdivall=::divides(c,p,c0,r);}
  if(pdivall) { pdivall=::divides(d,p,d0,r);}
  if(pdivall) { pdivall=::divides(e,p,e0,r);}

  if(pdivall) // Case where all coeffs are divisible by p
                 // a0 etc hold coeffs/p
    {
#ifdef DEBUG_MINIM
      cout<<"All coeffs divisible by "<<p<<endl;
#endif
      if(div(p,a0)&&div(p,b0))  // mult root is at infty move to 0
	{
	  m_invert(a,b,c,d,e,m);
#ifdef DEBUG_MINIM
	  cout<<"Multiple root is at infinity, reversing coeffs\n";
#endif
	}
      else 
	{
	  if(ndiv(p,e0)||ndiv(p,d0)) // mult root is finite and not zero
        	                     // find it and shift it to 0
	    {
	      bigint alpha = root_p(a0,b0,c0,d0,e0,p);
#ifdef DEBUG_MINIM
	      cout<<"Multiple root is at "<<alpha<<" mod "<<p<<"\n";
	      cout<<"Shifting to 0...\n";
#endif
	      xshift(alpha,a,b,c,d,e,m);

	      if(ndiv(p2,e)||ndiv(p2,d)||ndiv(p2,c)) 
		{
		  cout<<"Error in c, d, e\n";
		  cout<<"root_p("<<a0<<","<<b0<<","<<c0<<","<<d0<<","<<e0<<","<<p<<") returns "<<alpha<<endl;
		}
	    }
	}
      // if ndiv(p2,b) must now shift second root to infty
      divide_exact(b,p,b0); // b0=b/p;
      if(ndiv(p2,a)&&ndiv(p,b0))
	{
	  bigint gamma = mod(-a0*invmod(b0,p),p);
#ifdef DEBUG_MINIM
	  cout<<"Triple root case, shifting fourth root to infinity\n";
	  cout<<"(gamma = "<<gamma<<")\n";
#endif
  // replace g(X,Z) by g(X, Z+gamma*X)
	  zshift(gamma,a,b,c,d,e,m);
	  divide_exact(b,p,b0); // b0=b/p;
	}
      if(two)
	{
	  if(ndiv(16,e)) // failure (cannot happen if 2^6|I, 2^7|J 
	                 //                        or 2^5|I & Q_2-soluble)
	    {
#ifdef DEBUG_MINIM
	      cout<<"Non-reducible\n";
#endif
	      return 0;
	    }
	}
      b=b0;
      divide_exact(c,p2,c);    // c=c/p2;
      divide_exact(d,p*p2,d);  // d=(d/p2)/p;
      divide_exact(e,p2*p2,e); // e=(e/p2)/p2;
      m.x_scale(p);
      m.u_scale(p2);
      return 1;
    }

  // Now the case where not all are divisible by p...
#ifdef DEBUG_MINIM
  cout<<"Not all coeffs divisible by "<<p<<endl;
#endif
  if(div(p,a)&&div(p,b)) // mult root is at infty,  move it to 0
    {
      m_invert(a,b,c,d,e,m);
#ifdef DEBUG_MINIM
      cout<<"Multiple root is at infinity, reversing coeffs\n";
#endif
    }
  else 
    {
      if(ndiv(p,e)||ndiv(p,d)) // mult root finite and not zero
	                       // find it and shift it to 0
	{
	  bigint alpha = root_p(a,b,c,d,e,p);
#ifdef DEBUG_MINIM
	  cout<<"Multiple root is at "<<alpha<<" mod "<<p<<"\n";
	  cout<<"Shifting to 0...\n";
#endif
	  xshift(alpha,a,b,c,d,e,m);

	  if(ndiv(p,c)||ndiv(p,d)||ndiv(p,e)) 
	    {
	      cout<<"Error in c, d or e\n";
	      cout<<"root_p("<<a<<","<<b<<","<<c<<","<<d<<","<<e<<","<<p<<") returns "<<alpha<<endl;
	    }
	}
    }
  // if ndiv(p,b) must now shift second root to infty
  if(ndiv(p,a)&&ndiv(p,b))
    {
      bigint gamma = mod(-a*invmod(b,p),p);
#ifdef DEBUG_MINIM
      cout<<"Triple root case, shifting fourth root to infinity\n";
      cout<<"(gamma = "<<gamma<<")\n";
#endif
  // replace g(X,Z) by g(X, Z+gamma*X)
      zshift(gamma,a,b,c,d,e,m);
    }
  if(div(p,a)) // triple root case
    {
      bigint beta; beta=0;
      if(three)
	{
	  long vpi = val(3,12*a*e-3*b*d+c*c);
	  if(vpi==4)
	    divide_exact(-e,BIGINT(27),temp);
	  else
	    divide_exact(-c,BIGINT(9),temp);
	  if(ndiv(3,temp)) beta = 3 * mod(temp * invmod(b,3) , 3);
	}
      else
	{
#ifdef DEBUG_MINIM
	  cout<<"fixing c to be divisible by p^2\n";
#endif
	  divide_exact(-c,p,temp);
	  if(ndiv(p,temp)) beta = p * mod(temp * invmod(3*b,p) , p);
	}
#ifdef DEBUG_MINIM
      cout<<"(beta = "<<beta<<")\n";
#endif
      xshift(beta,a,b,c,d,e,m);

      a=a*p2;
      divide_exact(c,p2,c);       // c=c/p2;
      divide_exact(d,p2*p2,d);    // d=(d/p2)/p2;
      divide_exact(e,p2*p2*p2,e); // e=((e/p2)/p2)/p2;
      m.x_scale(p2);
      m.u_scale(p*p2);
      return 1;
    }
  else  // quadruple root case
    {
      if(three)
	{
	  divide_exact(b,BIGINT(3),b0); // b0=b/3;
	  if(ndiv(p,b0))
	    {
	      bigint beta = 3 * mod(-b0*invmod(a,p),p);
#ifdef DEBUG_MINIM
	      cout<<"fixing b to be divisible by 9\n";
	      cout<<"(beta = "<<beta<<")\n";
#endif
	      xshift(beta,a,b,c,d,e,m);
	    }
	}
      if(two)
	{
	  if(ndiv(16,e)) // failure (cannot happen if 2^6|I, 2^7|J 
	                 //          unless "bad case quartic")
	    {
#ifdef DEBUG_MINIM
	      cout<<"Non-reducible, may be bad case\n";
#endif
	      return 0;
	    }
	}
      divide_exact(b,p,b);    // b=b/p;
      divide_exact(c,p2,c);    // c=c/p2;
      divide_exact(d,p*p2,d);  // d=(d/p2)/p;
      divide_exact(e,p2*p2,e); // e=(e/p2)/p2;
      m.x_scale(p);
      m.u_scale(p2);
      return 1;
    }
}

int is_nonmin(int smallp, long vpi, long vpj, long vpd, int assume_locsol)
// Given vpi = val(p,I) and vpj=val(p,J) returns 1 if non-minimal
// smallp = p if p=2,3 else =1.
// p=3: needs also vpd=val(p,disc)
// p=2: may or may not be minimizable, but worth a try
// (The commented out condition is sufficient but NOT necessary)
{
  if(!assume_locsol)   return (vpi>7)&&(vpj>11);
  //  if(smallp==2) return ((vpi>5)&&(vpj>8)&&(vpd>9));
  if(smallp==3) return  (((vpi>4)&&(vpj>8)) || ((vpi==4)&&(vpj==6)&&(vpd>14)));
  return (vpi>3)&&(vpj>5);
}

void minim_all(bigint& ga, bigint& gb, bigint& gc, bigint& gd, bigint& ge, 
	       bigint& I, bigint& J, const vector<bigint>& plist,
	       scaled_unimod& m,
	       int assume_locsol, int verb)
{
  unsigned long i; long j;
  for(i=0; i<plist.size(); i++)
    {
      bigint p=plist[i];
      long smallp=1;            // these save testing a (possibly big) p 
      if (p==2) smallp=2;       // all the time
      else if (p==3) smallp=3;
	  
      long vpi=1000; if(!is_zero(I)) vpi=val(p,I);
      long vpj=1000; if(!is_zero(J)) vpj=val(p,J);
      long vpd=0;
      int nonmin, success;
      if(smallp==3) vpd = val(p,4*I*sqr(I)-sqr(J));
      nonmin = is_nonmin(smallp,vpi,vpj,vpd, assume_locsol);
      if(!nonmin) 
	{
	  if(verb) cout<<"p="<<p<<": minimal already\n";
	  continue;
	}
      if(verb) 
	{
	  cout<<"p="<<p<<": ";
	  if(smallp==2) cout<<"(possibly) ";
	  cout<<"non-minimal (vp(I)="<<vpi<<", vp(J)="<<vpj<<")";
	}

      // Now nonminimal so do something
      if(verb) cout<<" minimalizing at "<<p<<"....\n";
      while(nonmin)
	{
	  success=minim_p(ga,gb,gc,gd,ge,p,m);
	  if(success) // can only fail for p=2
	    {
	      vpi-=4; vpj-=6; 
	      j=4; while(j--) divide_exact(I,p,I); // I/=p;
	      j=6; while(j--) divide_exact(J,p,J); // J/=p;
	      if(smallp==3) vpd-=12;
	      nonmin = is_nonmin(smallp,vpi,vpj,vpd, assume_locsol);
	    }
	  else nonmin=0;  // avoid looping when p=2!
	}
      if(verb) 
	{
	  cout<<"Finished minimalizing at "<<p<<", new coefficients: \n";
	  cout<<"("<<ga<<","<<gb<<","<<gc<<","<<gd<<","<<ge<<")"<<endl;
	  cout<<"transform = "<<m<<endl;
	}
      bigint newI = II(ga,gb,gc,gd,ge);
      bigint newJ = JJ(ga,gb,gc,gd,ge);
      if((I!=newI)||(J!=newJ))
	{
	  cout<<"Error in previous step: wrong I, J.  New quartic has\n";
	  cout<<"I = "<<newI<<"\nJ = "<<newJ<<endl;
	  cout<<"but should be\n";
	  cout<<"I = "<<I<<"\nJ = "<<J<<endl;
	}
      else
	if(verb) cout<<"I = "<<I<<"\nJ = "<<J<<endl;
    }
}
