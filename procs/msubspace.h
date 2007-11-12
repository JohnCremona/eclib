// msubspace.h: declarations of multiprecision subspace class
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
 
#if     !defined(_MSUBSPACE_H)
#define _MSUBSPACE_H      1       //flags that this file has been included

#include "mmatrix.h"

class msubspace {

public:
     // constructors
        msubspace(int n=0) 
	  :pivots(iota(n)),basis(midmat(n)) {denom=1;}
        msubspace(const mat_m& b, const vec_i& p, const bigint& d)
	  :denom(d),pivots(p),basis(b) {}
        msubspace(const msubspace& s)
	  :denom(s.denom),pivots(s.pivots),basis(s.basis) {}
     // destructor
        ~msubspace() {;}
     // assignment
        void operator=(const msubspace& s) 
           {pivots=s.pivots; basis=s.basis; denom=s.denom;}

     // member functions & operators
        void clear() { pivots.init(); basis.init();}

     // non-member (friend) functions and operators
        friend int dim(const msubspace& s);      // the dimension
        friend bigint denom(const msubspace& s);   // the denominator
        friend vec_i pivots(const msubspace& s);// the pivot vector
        friend mat_m basis(const msubspace& s) ;// the basis matrix
        friend msubspace combine(const msubspace& s1, const msubspace& s2);
        friend mat_m restrict(const mat_m& m, const msubspace& s);
        friend msubspace pcombine(const msubspace& s1, const msubspace& s2, const bigint& pr);
        friend mat_m prestrict(const mat_m& m, const msubspace& s, const bigint& pr);
        friend msubspace lift(const msubspace& s, const bigint& pr, int trace);


// Implementation
private:
       bigint   denom;
       vec_i pivots;
       mat_m basis;
};


// Declarations of nonmember, nonfriend operators and functions:

msubspace kernel(const mat_m& mat, int method=0);
msubspace image(const mat_m& mat, int method=0);
msubspace eigenspace(const mat_m& mat, const bigint& lambda, int method=0);
msubspace subeigenspace(const mat_m& mat, const bigint& l, const msubspace& s, int method=0);


//The following work with msubspaces "mod p" using "echmodp" from
//mmatrix.h/cc to do gaussian elimination.  The "denom" of each is 1.

msubspace pkernel(const mat_m& mat, const bigint& pr);
msubspace pimage(const mat_m& mat, const bigint& pr);
msubspace peigenspace(const mat_m& mat, const bigint& lambda, const bigint& pr);
msubspace psubeigenspace(const mat_m& mat, const bigint& l, const msubspace& s, const bigint& pr);

inline int dim(const msubspace& s) {return ncols(s.basis);}  // the dimension
inline bigint denom(const msubspace& s) {return s.denom;}   // the denominator
inline vec_i pivots(const msubspace& s) {return s.pivots;} // the pivot vector
inline mat_m basis(const msubspace& s) {return s.basis;}  // the basis matrix

#endif
