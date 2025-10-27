// mquartic.h:   Declaration of class quartic and related functions
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
 
#ifndef _ECLIB_MQUARTIC_H
#define _ECLIB_MQUARTIC_H

#include "points.h"

inline ZZ II(const ZZ& a, const ZZ& b, const ZZ& c, 
	  const ZZ& d, const ZZ& e)
{
  return 12*a*e - 3*b*d + c*c;
}

inline ZZ JJ(const ZZ& a, const ZZ& b, const ZZ& c, 
	  const ZZ& d, const ZZ& e)
{
  return (72*a*e + 9*b*d - 2*sqr(c)) * c - 27*(a*sqr(d) + sqr(b)*e);
}

inline ZZ H_invariant(const ZZ& a, const ZZ& b, const ZZ& c)
{
  return 8*a*c - 3*sqr(b);
}

inline ZZ R_invariant(const ZZ& a, const ZZ& b, const ZZ& c, 
		 const ZZ& d)
{
  return pow(b,3) + 8*d*sqr(a) - 4*a*b*c;
}


class quartic {
public:
     // constructors
        quartic();
        quartic(const ZZ& qa, const ZZ& qb, const ZZ& qc, 
		const ZZ& qd, const ZZ& qe, 
                const vector<bigcomplex>& qr,	int qt,
		const ZZ& qi,const ZZ& qj,const ZZ& qdisc);
        quartic(const ZZ& qa, const ZZ& qb, const ZZ& qc, 
		const ZZ& qd, const ZZ& qe);
  // latter calls set_roots_and_type()
        quartic(const quartic& q);
  // member functions & operators
        void set_roots_and_type();
	void assign(const ZZ& qa, const ZZ& qb, const ZZ& qc, 
		    const ZZ& qd, const ZZ& qe, 
		    const vector<bigcomplex>& qr,    int qt,
		    const ZZ& qi,const ZZ& qj,
		    const ZZ& qdisc);
	void assign(const ZZ& qa, const ZZ& qb, const ZZ& qc, 
		    const ZZ& qd, const ZZ& qe);
  // latter calls set_roots_and_type()
        void operator=(const quartic& q);
        ZZ geta() const {return a;}
        ZZ getb() const {return b;}
        ZZ getcc() const {return c;} //NB "getc" is a standard macro
        ZZ getd() const {return d;}
        ZZ gete() const {return e;}
	void getcoeffs(ZZ& xa, ZZ& xb, ZZ& xc, ZZ& xd, ZZ& xe) const {xa=a; xb=b; xc=c; xd=d; xe=e;}
        int gettype() const {return type;}
        ZZ getI() const {return ii;}
        ZZ getJ() const {return jj;}
        ZZ getH() const {return H_invariant(a,b,c);}
        ZZ getdisc() const {return disc;}
        vector<bigcomplex> getroots(void) const {return roots;}
        void doubleup()
          {
	    b*=2; c*=4; d*=8; e*=16; ii*=16; jj*=64; disc*=4096;
         }
        vector<bigrational> rational_roots() const; // returns rational roots
        int trivial() const;     // Checks for a rational root
	long nrootsmod(long p) const;
        friend ostream& operator<<(ostream& s, const quartic& q);
	friend int new_equiv( quartic& q1, quartic& q2, int info);
	friend void qc(quartic& g,
		       const ZZ& x0,  const ZZ& y0,  const ZZ& z0,
		       Curvedata * E, Curvedata* IJ_curve, 
		       const ZZ& tr_u, const ZZ& tr_r, 
		       const ZZ& tr_s, const ZZ& tr_t,  
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
       ZZ a,b,c,d,e; // coefficients
       vector<bigcomplex> roots; // 4 roots, created in all constructors
       int type;       // 1, 2 or 3
       ZZ ii,jj,disc;
// The following are used by new_equiv:  (NB p = -H)
       ZZ p, r, psq, asq;
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

quartic make_quartic(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d, const ZZ& e);

#endif
