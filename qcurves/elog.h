// elog.h: declarations of elliptic logarithm functions
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
 
#ifndef _ELOG_H_
#define _ELOG_H_

#include "curve.h"
#include "points.h"
#include "cperiods.h"

//////////////////////////////////////////////////////////////////////////
//
//  Functions for passing between the complex torus C/L and a
//  Weierstrass model for an elliptic curve
//
//////////////////////////////////////////////////////////////////////////

//  1. Elliptic logarithm

// Given an elliptic curve and its (precomputed) periods, and a real
// or rational (but NOT complex) point P=(x,y), returns the unique
// complex number z such that

// (1) \wp(z)=x+b2/12, \wp'(z)=2y+a1*x+a3, 

// (2) either z is real and 0\le z\lt w1, or Delta>0, z-w2/2 is real
// and 0\le z-w2/2\le w1.  

// Here, [w1,w2] is the standard period lattice basis

// c.f. Cohen page 399

// First & second functions:  P=[x,y] with x,y real; maps to z mod lattice

bigcomplex ellpointtoz(const Curvedata& E, const Cperiods& per, 
		       const bigfloat& x, const bigfloat& y);

inline bigcomplex ellpointtoz(const Curvedata& E, const Cperiods& per, 
			      const vector<bigfloat> P) 
{return ellpointtoz(E,per,P[0],P[1]);}

// Third function:  P=rational point; maps to z mod lattice

inline bigcomplex elliptic_logarithm(const Curvedata& E, const Cperiods& per, 
				     const Point& P)
{
  if(P.iszero()) return bigcomplex(to_bigfloat(0));
  bigfloat xP, yP;
  realify_point(P,xP,yP);
  return ellpointtoz(E,per,xP,yP);
}

//  2. Weierstrass functions (interface to cperiods.h/cc)

// Cperiods is a class containing a basis for the period lattice L;
// it knows how to compute points from z mod L; so this function
// effectively does the same as PARI's ellztopoint()
//
// First function:  given z mod L, returns complex vector [x,y]

vector<bigcomplex> ellztopoint(Curvedata& E,  Cperiods& per, 
			       const bigcomplex& z);

// Second function, expects to return a rational point.
// User supplies a denominator for the point; if it doesn't work, the
// Point returned is 0 on the curve

Point ellztopoint(Curvedata& E,  Cperiods& per, const bigcomplex& z, 
		  const bigint& den);

// Returns a (possibly empty) vector of solutions to m*Q=P

// First version will compute the Cperiods itself, so best to use the
// second one if more than one call is to be made for the same curve

vector<Point> division_points(Curvedata& E,  const Point& P, int m);
vector<Point> division_points(Curvedata& E,  Cperiods& per, const Point& P, int m);

// Returns a vector of solutions to m*Q=0 (including Q=0)

// First version will compute the Cperiods itself, so best to use the
// second one if more than one call is to be made for the same curve

vector<Point> torsion_points(Curvedata& E,int m);
vector<Point> torsion_points(Curvedata& E,  Cperiods& per, int m);

#endif
