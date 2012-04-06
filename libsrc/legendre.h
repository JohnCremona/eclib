// legendre.h: declarations of functions for solving legendre equations
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
 
#ifndef _LEGENDRE_H
#define _LEGENDRE_H      1 //flags that this file has been included

#include "quadratic.h"

void minv(const bigint& a1, const bigint& a2, 
	  const bigint& b1, const bigint& b2, 
	  const bigint& c1, const bigint& c2, 
	  bigint& xmin, const bigint& ymin);

// (e+f*i) is (a+b*i) mod (c+d*i) in Z[i]
void GIreduce(const bigint& a, const bigint& b, const bigint& c, 
	      const bigint& d, bigint& e, bigint& f);
// (e+f*i) = gcd ( (a+b*i) , (c+d*i) ) in Z[i]
void GIgcd(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& d, bigint& e, bigint& f);

// Solve Legendre's equation ax^2+by^2+cz^2=0 using Rusin's reduction, 
// without assuming a,b,c pairwise coprime
// returns 0 if not soluble
int legendre_solve(const bigint& a, const bigint& b, const bigint& c, 
		    bigint& x, bigint& y, bigint& z,
		    int use_lll=0);

int legendre_solve(const bigint& a, const bigint& b, const bigint& c, 
		   const vector<bigint>& factorbase,
		    bigint& x, bigint& y, bigint& z,
		    int use_lll=0);

// Solve Legendre's equation ax^2+by^2+cz^2=0 using Rusin's reduction, 
// given "certificate" (n,p,q) satisfying a|n^2+bc, b|p^2+ac, c|q^2+ab.
// assumes a,b,c pairwise coprime
void legendre_solve_cert(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& n, const bigint& p, const bigint& q, 
		    bigint& x, bigint& y, bigint& z);

//Ensure that input is valid
int checkin(const bigint& a,const bigint& b,const bigint& c,
	    const bigint& n,const bigint& p,const bigint& q);

// Check that purported solution is OK
int check_leg(const bigint& a, const bigint& b, const bigint& c,
	      bigint& x, bigint& y, bigint& z);

// Check that purported solution is OK & in correct lattice
int check_leg(const bigint& a, const bigint& b, const bigint& c,
	      const bigint& n, const bigint& p, const bigint& q, 
	      bigint& x, bigint& y, bigint& z);

//Peel off a known square factor  u  from coefficient  a
void lem2a(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   const bigint& u,
	   bigint& x, bigint& y, bigint& z);

//Peel off a known square factor  u  from coefficient  b
void lem2b(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   const bigint& u,
	   bigint& x, bigint& y, bigint& z);

//Peel off a known square factor  u  from coefficient  c
void lem2c(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   const bigint& u,
	   bigint& x, bigint& y, bigint& z);


void lem4(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& n, const bigint& p, const bigint& q, 
	  bigint& x, bigint& y, bigint& z);

// Finds a certificate or returns 0 if none exists:
int make_certificate(const bigint& a, const bigint& b, const bigint& c, 
		     bigint& n, bigint& p, bigint& q);
int make_certificate(const bigint& a, const vector<bigint>& apdivs, 
		     const bigint& b, const vector<bigint>& bpdivs, 
		     const bigint& c, const vector<bigint>& cpdivs, 
		     bigint& n, bigint& p, bigint& q);

// Check to see if  b is congruent to  +- c  mod a  (assumed positive!)
//  if not, how much to divide  a  by to ensure that congruence holds?
bigint should(const bigint& a, const bigint& b, const bigint& c);

void legendre_reduce(const bigint& a, const bigint& b, const bigint& c, 
		     bigint& x0, bigint& y0, bigint& z0, int verb=0);
     // Given a, b, c,  ax^2+by^2+cz^2=0
     // reduces x, y, z in place using Mordell's method (page 48)
     // to achieve Holzer's bounds |z|<=sqrt(ab) etc.
  // (just permutes & passes to conic_mordell_reduce())

void new_legendre_reduce(const bigint& a, const bigint& b, const bigint& c, 
		     bigint& x0, bigint& y0, bigint& z0, int verb=0);
     // Given a, b, c,  ax^2+by^2+cz^2=0
     // reduces x, y, z in place using quadratics

void legendre_via_lll(const bigint& a, const bigint& b, const bigint& c, 
		      const bigint& k1, const bigint& k2, const bigint& k3, 
		      bigint& x, bigint& y, bigint& z);

//
// Given one solution, returns quadratics parametrizing all solutions,
// with discriminants -4bc, -4ac, -4ab.  Here a, b, c are assumed
// pairwise coprime but not square-free, and a>0, b>0, c<0.

void legendre_param(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& x0, const bigint& y0, const bigint& z0, 
		    quadratic& qx, quadratic& qy, quadratic& qz);

//These versions EITHER return 0 with a solution in x,y,z in the lattice
//               OR return 1, 2, 3 and a square factor u of a, b, c (resp)
//               OR return -1 if something fails (should not happen)
int legendre_solve_cert_1(const bigint& a, const bigint& b, const bigint& c, 
			  const bigint& n, const bigint& p, const bigint& q, 
			  bigint& x, bigint& y, bigint& z,
			  bigint& u);
int lem4_1(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   bigint& x, bigint& y, bigint& z,
	   bigint& u);

bigfloat holzer_measure(const bigint& a, const bigint& b, const bigint& c, 
			const bigint& x, const bigint& y, const bigint& z);
// max{|a|x^2,|b|y^2,|c|z^2}/|abc|   ( < 1 for a Holzer-reduced solution)

#endif
