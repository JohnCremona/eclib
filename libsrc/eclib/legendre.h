// legendre.h: declarations of functions for solving legendre equations
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
 
#ifndef _ECLIB_LEGENDRE_H
#define _ECLIB_LEGENDRE_H      1 //flags that this file has been included

#include "quadratic.h"

void minv(const ZZ& a1, const ZZ& a2, 
	  const ZZ& b1, const ZZ& b2, 
	  const ZZ& c1, const ZZ& c2, 
	  ZZ& xmin, const ZZ& ymin);

// (e+f*i) is (a+b*i) mod (c+d*i) in Z[i]
void GIreduce(const ZZ& a, const ZZ& b, const ZZ& c, 
	      const ZZ& d, ZZ& e, ZZ& f);
// (e+f*i) = gcd ( (a+b*i) , (c+d*i) ) in Z[i]
void GIgcd(const ZZ& a, const ZZ& b, const ZZ& c, 
	   const ZZ& d, ZZ& e, ZZ& f);

// Solve Legendre's equation ax^2+by^2+cz^2=0 using Rusin's reduction, 
// without assuming a,b,c pairwise coprime
// returns 0 if not soluble
int legendre_solve(const ZZ& a, const ZZ& b, const ZZ& c, 
		    ZZ& x, ZZ& y, ZZ& z,
		    int use_lll=0);

int legendre_solve(const ZZ& a, const ZZ& b, const ZZ& c, 
		   const vector<ZZ>& factorbase,
		    ZZ& x, ZZ& y, ZZ& z,
		    int use_lll=0);

// Solve Legendre's equation ax^2+by^2+cz^2=0 using Rusin's reduction, 
// given "certificate" (n,p,q) satisfying a|n^2+bc, b|p^2+ac, c|q^2+ab.
// assumes a,b,c pairwise coprime
void legendre_solve_cert(const ZZ& a, const ZZ& b, const ZZ& c, 
		    const ZZ& n, const ZZ& p, const ZZ& q, 
		    ZZ& x, ZZ& y, ZZ& z);

//Ensure that input is valid
int checkin(const ZZ& a,const ZZ& b,const ZZ& c,
	    const ZZ& n,const ZZ& p,const ZZ& q);

// Check that purported solution is OK
int check_leg(const ZZ& a, const ZZ& b, const ZZ& c,
	      const ZZ& x, const ZZ& y, const ZZ& z);

// Check that purported solution is OK & in correct lattice
int check_leg(const ZZ& a, const ZZ& b, const ZZ& c,
	      const ZZ& n, const ZZ& p, const ZZ& q, 
	      const ZZ& x, const ZZ& y, const ZZ& z);

//Peel off a known square factor  u  from coefficient  a
void lem2a(const ZZ& a, const ZZ& b, const ZZ& c, 
	   const ZZ& n, const ZZ& p, const ZZ& q, 
	   const ZZ& u,
	   ZZ& x, ZZ& y, ZZ& z);

//Peel off a known square factor  u  from coefficient  b
void lem2b(const ZZ& a, const ZZ& b, const ZZ& c, 
	   const ZZ& n, const ZZ& p, const ZZ& q, 
	   const ZZ& u,
	   ZZ& x, ZZ& y, ZZ& z);

//Peel off a known square factor  u  from coefficient  c
void lem2c(const ZZ& a, const ZZ& b, const ZZ& c, 
	   const ZZ& n, const ZZ& p, const ZZ& q, 
	   const ZZ& u,
	   ZZ& x, ZZ& y, ZZ& z);


void lem4(const ZZ& a, const ZZ& b, const ZZ& c, 
	  const ZZ& n, const ZZ& p, const ZZ& q, 
	  ZZ& x, ZZ& y, ZZ& z);

// Finds a certificate or returns 0 if none exists:
int make_certificate(const ZZ& a, const ZZ& b, const ZZ& c, 
		     ZZ& n, ZZ& p, ZZ& q);
int make_certificate(const ZZ& a, const vector<ZZ>& apdivs, 
		     const ZZ& b, const vector<ZZ>& bpdivs, 
		     const ZZ& c, const vector<ZZ>& cpdivs, 
		     ZZ& n, ZZ& p, ZZ& q);

// Check to see if  b is congruent to  +- c  mod a  (assumed positive!)
//  if not, how much to divide  a  by to ensure that congruence holds?
ZZ should(const ZZ& a, const ZZ& b, const ZZ& c);

void legendre_reduce(const ZZ& a, const ZZ& b, const ZZ& c, 
		     ZZ& x0, ZZ& y0, ZZ& z0, int verb=0);
     // Given a, b, c,  ax^2+by^2+cz^2=0
     // reduces x, y, z in place using Mordell's method (page 48)
     // to achieve Holzer's bounds |z|<=sqrt(ab) etc.
  // (just permutes & passes to conic_mordell_reduce())

void new_legendre_reduce(const ZZ& a, const ZZ& b, const ZZ& c, 
		     ZZ& x0, ZZ& y0, ZZ& z0, int verb=0);
     // Given a, b, c,  ax^2+by^2+cz^2=0
     // reduces x, y, z in place using quadratics

void legendre_via_lll(const ZZ& a, const ZZ& b, const ZZ& c, 
		      const ZZ& k1, const ZZ& k2, const ZZ& k3, 
		      ZZ& x, ZZ& y, ZZ& z);

//
// Given one solution, returns quadratics parametrizing all solutions,
// with discriminants -4bc, -4ac, -4ab.  Here a, b, c are assumed
// pairwise coprime but not square-free, and a>0, b>0, c<0.

void legendre_param(const ZZ& a, const ZZ& b, const ZZ& c, 
		    const ZZ& x0, const ZZ& y0, const ZZ& z0, 
		    quadratic& qx, quadratic& qy, quadratic& qz);

//These versions EITHER return 0 with a solution in x,y,z in the lattice
//               OR return 1, 2, 3 and a square factor u of a, b, c (resp)
//               OR return -1 if something fails (should not happen)
int legendre_solve_cert_1(const ZZ& a, const ZZ& b, const ZZ& c, 
			  const ZZ& n, const ZZ& p, const ZZ& q, 
			  ZZ& x, ZZ& y, ZZ& z,
			  ZZ& u);
int lem4_1(const ZZ& a, const ZZ& b, const ZZ& c, 
	   const ZZ& n, const ZZ& p, const ZZ& q, 
	   ZZ& x, ZZ& y, ZZ& z,
	   ZZ& u);

bigfloat holzer_measure(const ZZ& a, const ZZ& b, const ZZ& c, 
			const ZZ& x, const ZZ& y, const ZZ& z);
// max{|a|x^2,|b|y^2,|c|z^2}/|abc|   ( < 1 for a Holzer-reduced solution)

#endif
