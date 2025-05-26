// sub.h: declaration of class subspace
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
 
// Not to be included directly by user: use subspace.h which defines
// _ECLIB_SUBSPACE_H and includes this twice
//

class subspace {

public:
  // constructors
  subspace(int n=0) :denom(1),pivots(vec_i::iota(n)),basis(mat::identity_matrix(n)) {}
  subspace(const mat& b, const vec_i& p, const scalar& d) :denom(d),pivots(p),basis(b) {}
  subspace(const subspace& s) :denom(s.denom),pivots(s.pivots),basis(s.basis) {}

     // assignment
	void operator=(const subspace& s);

     // member functions & operators
        inline void clear() { pivots.init(); basis.init();}
        inline scalar den() const {return denom;}     // the denominator
        inline vec_i pivs() const {return pivots;} // the pivot vector
        inline mat bas() const {return basis;}   // the basis matrix

     // non-member (friend) functions and operators
        friend int dim(const subspace& s);      // the dimension
        friend scalar denom(const subspace& s);   // the denominator
        friend vec_i pivots(const subspace& s);// the pivot vector
        friend mat basis(const subspace& s) ;// the basis matrix
	friend subspace combine(const subspace& s1, const subspace& s2);
        friend mat restrict_mat(const mat& m, const subspace& s, int cr);
	friend subspace pcombine(const subspace& s1, const subspace& s2, const scalar& pr);
	friend mat prestrict(const mat& m, const subspace& s, const scalar& pr, int cr);
        friend int lift(const subspace& s, const scalar& pr, subspace& ans);

// Implementation
private:
       scalar   denom;
       vec_i pivots;
       mat basis;
};

// Declarations of nonmember, nonfriend operators and functions:

mat expressvectors(const mat& m, const subspace& s);
subspace kernel(const mat& m, int method=0);
subspace image(const mat& m, int method=0);
subspace eigenspace(const mat& m, const scalar& lambda, int method=0);
subspace subeigenspace(const mat& m, const scalar& l, const subspace& s, int method=0);


//The following work with subspaces "mod p" using "echmodp" from
//matrix.h/cc to do gaussian elimination.  The "denom" of each is 1.

subspace oldpkernel(const mat& m, const scalar& pr);
subspace pkernel(const mat& m, const scalar& pr);
subspace pimage(const mat& m, const scalar& pr);
subspace peigenspace(const mat& m, const scalar& lambda, const scalar& pr);
subspace psubeigenspace(const mat& m, const scalar& l, const subspace& s, const scalar& pr);


inline int dim(const subspace& s) {return s.basis.ncols();}  // the dimension
inline scalar denom(const subspace& s) {return s.denom;}    // the denominator
inline vec_i pivots(const subspace& s) {return s.pivots;}     // the pivot vector
inline mat basis(const subspace& s) {return s.basis;}       // the basis matrix
