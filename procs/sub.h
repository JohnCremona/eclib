// sub.h: declaration of class subspace
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
 
// SCALAR_OPTION must be set to 1 or 2 by including file

class subspace {

public:
     // constructors
        subspace(int n=0);
        subspace(const mat& b, const vec& p, scalar d);
	subspace(const subspace& s);
     // destructor
        ~subspace();
     // assignment
	void operator=(const subspace& s);

     // member functions & operators
        inline void clear() { pivots.init(); basis.init();}
        inline scalar den() const {return denom;}     // the denominator
        inline vec pivs() const {return pivots;} // the pivot vector
        inline mat bas() const {return basis;}   // the basis matrix

     // non-member (friend) functions and operators
        friend int dim(const subspace& s);      // the dimension
        friend scalar denom(const subspace& s);   // the denominator
        friend vec pivots(const subspace& s);// the pivot vector
        friend mat basis(const subspace& s) ;// the basis matrix
	friend subspace combine(const subspace& s1, const subspace& s2);
	friend mat restrict(const mat& m, const subspace& s, int cr);
	friend subspace pcombine(const subspace& s1, const subspace& s2, scalar pr);
	friend mat prestrict(const mat& m, const subspace& s, scalar pr, int cr);
	friend subspace lift(const subspace& s, scalar pr, int trace);


// Implementation
private:
       scalar   denom;
       vec pivots;
       mat basis;
};


// Declarations of nonmember, nonfriend operators and functions:

mat expressvectors(const mat& m, const subspace& s);
subspace kernel(const mat& m, int method=0);
subspace image(const mat& m, int method=0);
subspace eigenspace(const mat& m, scalar lambda, int method=0);
subspace subeigenspace(const mat& m, scalar l, const subspace& s, int method=0);


//The following work with subspaces "mod p" using "echmodp" from
//matrix.h/cc to do gaussian elimination.  The "denom" of each is 1.

subspace oldpkernel(const mat& m, scalar pr);
subspace pkernel(const mat& m, scalar pr);
subspace pimage(const mat& m, scalar pr);
subspace peigenspace(const mat& m, scalar lambda, scalar pr);
subspace psubeigenspace(const mat& m, scalar l, const subspace& s, scalar pr);


inline int dim(const subspace& s) {return ncols(s.basis);}  // the dimension
inline scalar denom(const subspace& s) {return s.denom;}    // the denominator
inline vec pivots(const subspace& s) {return s.pivots;}     // the pivot vector
inline mat basis(const subspace& s) {return s.basis;}       // the basis matrix
