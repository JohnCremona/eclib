// mquartic.h:   Declaration of class quartic and related functions
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
 
#ifndef _MQUARTIC_H
#define _MQUARTIC_H

#define NEW_EQUIV
#include "compproc.h"
#include "marith.h"
#include "points.h"
#include "realroots.h"

inline bigint II(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& d, const bigint& e)
{
  return 12*a*e - 3*b*d + c*c;
}

inline bigint JJ(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& d, const bigint& e)
{
  return (72*a*e + 9*b*d - 2*sqr(c)) * c - 27*(a*sqr(d) + sqr(b)*e);
}

inline bigint H_invariant(const bigint& a, const bigint& b, const bigint& c)
{
  return 8*a*c - 3*sqr(b);
}

inline bigint R_invariant(const bigint& a, const bigint& b, const bigint& c, 
		 const bigint& d)
{
  return pow(b,3) + 8*d*sqr(a) - 4*a*b*c;
}


class quartic {
public:
     // constructors (NOT inline as they all use "new"
        quartic();
        quartic(const bigint& qa, const bigint& qb, const bigint& qc, 
		const bigint& qd, const bigint& qe, 
                bigcomplex* qr,	int qt,
		const bigint& qi,const bigint& qj,const bigint& qdisc);
        quartic(const bigint& qa, const bigint& qb, const bigint& qc, 
		const bigint& qd, const bigint& qe);
  // latter calls set_roots_and_type()
	~quartic();
        quartic(const quartic& q);
  // member functions & operators
        void set_roots_and_type();
	void assign(const bigint& qa, const bigint& qb, const bigint& qc, 
		    const bigint& qd, const bigint& qe, 
		    bigcomplex *qr,    int qt,
		    const bigint& qi,const bigint& qj,
		    const bigint& qdisc);
	void assign(const bigint& qa, const bigint& qb, const bigint& qc, 
		    const bigint& qd, const bigint& qe);
  // latter calls set_roots_and_type()
        void operator=(const quartic& q);
        bigint geta() const {return a;}
        bigint getb() const {return b;}
        bigint getcc() const {return c;} //NB "getc" is a standard macro
        bigint getd() const {return d;}
        bigint gete() const {return e;}
	void getcoeffs(bigint& xa, bigint& xb, bigint& xc, bigint& xd, bigint& xe) const {xa=a; xb=b; xc=c; xd=d; xe=e;}
        int gettype() const {return type;}
        bigint getI() const {return ii;}
        bigint getJ() const {return jj;}
        bigint getH() const {return H_invariant(a,b,c);}
        bigint getdisc() const {return disc;}
        bigcomplex* getroots(void) const {return roots;}
        void doubleup()
          {
	    b*=2; c*=4; d*=8; e*=16; ii*=16; jj*=64; disc*=4096;
         }
        int trivial() const;     // Checks for a rational root
	long nrootsmod(long p) const;
        friend ostream& operator<<(ostream& s, const quartic& q);
	friend int new_equiv(quartic* q1, quartic* q2, int info);
	friend void qc(quartic& g,
		       const bigint& x0,  const bigint& y0,  const bigint& z0,
		       Curvedata * E, Curvedata* IJ_curve, 
		       const bigint& tr_u, const bigint& tr_r, 
		       const bigint& tr_s, const bigint& tr_t,  
		       Point& P, int verbose);
	void dump(ostream& s) const
	  {
	    s<<"Coeffs: ("<<a<<","<<b<<","<<c<<","<<d<<","<<e<<")\n"
	     <<"Roots("<<roots<<"): \n"<<roots[0]<<"\n"<<roots[1]<<"\n"
	     <<roots[2]<<"\n"<<roots[3]<<"\n"
	     <<"Type = "<<type<<", I="<<ii<<", J="<<jj<<endl;
	  }
// Implementation
private:
       bigint a,b,c,d,e; // coefficients
       bigcomplex* roots; // roots, array 0f 4 created in all constructors
       int type;       // 1, 2 or 3
       bigint ii,jj,disc;
// The following are used by new_equiv:  (NB p = -H)
       bigint p, r, psq, asq;
       void make_zpol();
       int have_zpol;
       unsigned long equiv_code;
      public:
       unsigned long set_equiv_code(const vector<long>& plist);
};

inline ostream& operator<<(ostream& s, const quartic& q)
{
 s<<"("<<q.a<<","<<q.b<<","<<q.c<<","<<q.d<<","<<q.e<<")"<<flush;
 return s;
}

quartic make_quartic(const bigint& a, const bigint& b, const bigint& c, const bigint& d, const bigint& e);

#endif
