// conic.h: declarations of functions for solving conics (see also legendre.h)
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
 
#ifndef _ECLIB_CONIC_H
#define _ECLIB_CONIC_H      1  //flags that this file has been included

#include "quadratic.h"

ZZ cancel1(ZZ& x, ZZ& y, ZZ& z);
     // cancels common factors only, return gcd

void cancel(ZZ& x, ZZ& y, ZZ& z);
     // cancels common factors and leaves z>=0 or z=0 and x>=0

int solve_conic(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,
		const vector<ZZ>& factorbase,
                ZZ& x, ZZ& y, ZZ& z, int method=4);
int solve_conic(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,
                ZZ& x, ZZ& y, ZZ& z, int method=4);
     // Solves axx+bxz+czz=dyy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, c, d, bb-4ac non-zero
  
int solve_conic(const quadratic& q, const ZZ& d,
		ZZ& x, ZZ& y, ZZ& z, int method=4);

int solve_conic(const quadratic& q, const ZZ& d,
		       const vector<ZZ>& factorbase,
		ZZ& x, ZZ& y, ZZ& z, int method=4);

int solve_conic_diag(const ZZ& a, const vector<ZZ>& aplist,
		     const ZZ& b, const vector<ZZ>& bplist,
                ZZ& x, ZZ& y, ZZ& z, int method);
     // Solves xx-azz=byy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, b non-zero square-free, their prime divisors in aplist, bplist 

void conic_mordell_reduce(const ZZ& a, const ZZ& b, const ZZ& c, ZZ& x0, ZZ& y0, ZZ& z0, int verb=0);
     // Given a>0, b>0, c<0, abc square-free and ax^2+by^2+cz^2=0
     // reduces x, y, z in place using Mordell's method (page 48)
     // to achieve Holzer's bounds |z|<=sqrt(ab) etc.

void conic_diag_reduce(const ZZ& a, const ZZ& b, ZZ& x, ZZ& y, ZZ& z, int verb=0);
     // As above but with a,b square-free only, calls conic_mordell_reduce

int solve_conic_param(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,
		const vector<ZZ>& factorbase,
                quadratic& qx, quadratic& qy, quadratic& qz, int method=4, int verb=0);

int solve_conic_param(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,
                quadratic& qx, quadratic& qy, quadratic& qz, int method=4, int verb=0);
     // Solves axx+bxz+czz=dyy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, c, d, bb-4ac non-zero
     // qx,qy,qz are arrays of coeffs of parametrizing quadratics
     //   with leading coeffs qx[0],qy[0],qz[0] one solution

int solve_conic_param(const quadratic& q, const ZZ& d,
		      const vector<ZZ>& factorbase,
		      quadratic& qx, quadratic& qy, quadratic& qz, 
		      int method=4, int verb=0);

int solve_conic_param(const quadratic& q, const ZZ& d,
		      quadratic& qx, quadratic& qy, quadratic& qz, 
		      int method=4, int verb=0);

int testsol(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,
                const ZZ& x, const ZZ& y, const ZZ& z, int verb=0);

// Tests to see if a given solution is a non-trivial solution
int testsol(const quadratic& q, const ZZ& d,
	    const ZZ& x, const ZZ& y, const ZZ& z, 
	    int verb=0);

int testlocsol(const ZZ& a, 
	       const ZZ& b, 
	       const ZZ& c);
// tests if ax^2+by^2+cz^2=0 is soluble, where a, b, c are pairwise
// coprime and square-free

int testlocsol(const ZZ& a, const vector<ZZ>& alist, 
	       const ZZ& b, const vector<ZZ>& blist, 
	       const ZZ& c, const vector<ZZ>& clist);
// tests if ax^2+by^2+cz^2=0 is soluble, where a, b, c are pairwise
// coprime and square-free, their prime factors being in alist etc.

int testlocsol(const ZZ& a, const vector<ZZ>& alist, 
	       const ZZ& b, const vector<ZZ>& blist);
// tests if ax^2+by^2=z^2 is soluble, where a, b are
// square-free, their prime factors being in alist and blist.

int testparamsol(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d,
                const quadratic& qx, const quadratic& qy, const quadratic& qz, int verb=0);

// Tests to see if a given parametrization is a solution
int testparamsol(const quadratic& q, const ZZ& d,
		 const quadratic& qx, const quadratic& qy, const quadratic& qz, 
		 int verb=0);

//miscellaneous test functions:

void testmodsqrt();
void testsqf();
void testcancel();

// Output utilities:

void show_xyz(const ZZ& x, const ZZ& y, const ZZ& z);
void show_cert(const ZZ& p, const ZZ& q, const ZZ& r);
void show_eqn(const ZZ& a, const ZZ& b, const ZZ& c);
void show_eqn_cert(const ZZ& a, const ZZ& b, const ZZ& c,
		   const ZZ& p, const ZZ& q, const ZZ& r);
void show_all(const ZZ& a, const ZZ& b, const ZZ& c,
	      const ZZ& p, const ZZ& q, const ZZ& r,
	      const ZZ& x, const ZZ& y, const ZZ& z);

#endif
