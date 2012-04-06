// conic.h: declarations of functions for solving conics (see also legendre.h)
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
 
#ifndef _CONIC_H
#define _CONIC_H      1  //flags that this file has been included

#include "quadratic.h"

bigint cancel1(bigint& x, bigint& y, bigint& z);
     // cancels common factors only, return gcd

void cancel(bigint& x, bigint& y, bigint& z);
     // cancels common factors and leaves z>=0 or z=0 and x>=0

int solve_conic(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
		const vector<bigint>& factorbase,
                bigint& x, bigint& y, bigint& z, int method=4);
int solve_conic(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                bigint& x, bigint& y, bigint& z, int method=4);
     // Solves axx+bxz+czz=dyy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, c, d, bb-4ac non-zero
  
int solve_conic(const quadratic& q, const bigint& d,
		bigint& x, bigint& y, bigint& z, int method=4);

int solve_conic(const quadratic& q, const bigint& d,
		       const vector<bigint>& factorbase,
		bigint& x, bigint& y, bigint& z, int method=4);

int solve_conic_diag(const bigint& a, const vector<bigint>& aplist,
		     const bigint& b, const vector<bigint>& bplist,
                bigint& x, bigint& y, bigint& z, int method);
     // Solves xx-azz=byy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, b non-zero square-free, their prime divisors in aplist, bplist 

void conic_mordell_reduce(const bigint& a, const bigint& b, const bigint& c, bigint& x0, bigint& y0, bigint& z0, int verb=0);
     // Given a>0, b>0, c<0, abc square-free and ax^2+by^2+cz^2=0
     // reduces x, y, z in place using Mordell's method (page 48)
     // to achieve Holzer's bounds |z|<=sqrt(ab) etc.

void conic_diag_reduce(const bigint& a, const bigint& b, bigint& x, bigint& y, bigint& z, int verb=0);
     // As above but with a,b square-free only, calls conic_mordell_reduce

int solve_conic_param(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
		const vector<bigint>& factorbase,
                quadratic& qx, quadratic& qy, quadratic& qz, int method=4, int verb=0);

int solve_conic_param(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                quadratic& qx, quadratic& qy, quadratic& qz, int method=4, int verb=0);
     // Solves axx+bxz+czz=dyy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, c, d, bb-4ac non-zero
     // qx,qy,qz are arrays of coeffs of parametrizing quadratics
     //   with leading coeffs qx[0],qy[0],qz[0] one solution

int solve_conic_param(const quadratic& q, const bigint& d,
		      const vector<bigint>& factorbase,
		      quadratic& qx, quadratic& qy, quadratic& qz, 
		      int method=4, int verb=0);

int solve_conic_param(const quadratic& q, const bigint& d,
		      quadratic& qx, quadratic& qy, quadratic& qz, 
		      int method=4, int verb=0);

int testsol(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                const bigint& x, const bigint& y, const bigint& z, int verb=0);

// Tests to see if a given solution is a non-trivial solution
int testsol(const quadratic& q, const bigint& d,
	    const bigint& x, const bigint& y, const bigint& z, 
	    int verb=0);

int testlocsol(const bigint& a, 
	       const bigint& b, 
	       const bigint& c);
// tests if ax^2+by^2+cz^2=0 is soluble, where a, b, c are pairwise
// coprime and square-free

int testlocsol(const bigint& a, const vector<bigint>& alist, 
	       const bigint& b, const vector<bigint>& blist, 
	       const bigint& c, const vector<bigint>& clist);
// tests if ax^2+by^2+cz^2=0 is soluble, where a, b, c are pairwise
// coprime and square-free, their prime factors being in alist etc.

int testlocsol(const bigint& a, const vector<bigint>& alist, 
	       const bigint& b, const vector<bigint>& blist);
// tests if ax^2+by^2=z^2 is soluble, where a, b are
// square-free, their prime factors being in alist and blist.

int testparamsol(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                const quadratic& qx, const quadratic& qy, const quadratic& qz, int verb=0);

// Tests to see if a given parametrization is a solution
int testparamsol(const quadratic& q, const bigint& d,
		 const quadratic& qx, const quadratic& qy, const quadratic& qz, 
		 int verb=0);

//miscellaneous test functions:

void testmodsqrt();
void testsqf();
void testcancel();

// Output utilities:

void show_xyz(const bigint& x, const bigint& y, const bigint& z);
void show_cert(const bigint& p, const bigint& q, const bigint& r);
void show_eqn(const bigint& a, const bigint& b, const bigint& c);
void show_eqn_cert(const bigint& a, const bigint& b, const bigint& c,
		   const bigint& p, const bigint& q, const bigint& r);
void show_all(const bigint& a, const bigint& b, const bigint& c, 
	      const bigint& p, const bigint& q, const bigint& r, 
	      bigint& x, bigint& y, bigint& z);

#endif
