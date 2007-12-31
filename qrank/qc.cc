// qc.cc: implementation of function for mapping quartic point to curve
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
 
#include "marith.h"
#include "points.h"   // from qcurve library
#include "mquartic.h"
#include "qc.h"

//#define DEBUG

void qc(quartic& g,
        const bigint& x0,  const bigint& y0,  const bigint& z0,
        Curvedata * E, 
	Curvedata* IJ_curve, 
	const bigint& tr_u, const bigint& tr_r, 
	const bigint& tr_s, const bigint& tr_t,  
	Point& P, int verbose)
{
  bigint aa,bb,cc,dd,ee,p,q,r,t,xp,yp,zp;
  bigint a=g.a, b=g.b, c=g.c, d=g.d, e=g.e;

#ifdef DEBUG
  cout << "In qc(...) with:\n";
  cout << "IJ_curve = "<<(Curve)(*IJ_curve)<<"\n";
  cout << "min curve = "<<(Curve)(*E)<<"\n";
  cout << "[u,r,s,t] = ["<<tr_u<<","<<tr_r<<","<<tr_s<<","<<tr_t<<"]\n";
  cout << "(x0:y0:z0) = ("<<x0<<" : "<<y0<<" : "<<z0<<")\n";
  cout << "Quartic = (a,b,c,d,e) = ("<<a<<", "<<b<<", "<<c<<", "<<d<<", "<<e<<")\n";
#endif

  bigint z02=sqr(z0); 

  if(isqrt(a,q)) {z02=1;}  // else z0=0 which sets zp=0 below
  else if(isqrt(e,q)){t=a; a=e; e=t; t=b; b=d; d=t; }
    else 
      {
	const bigint& z03=z0*z02;  const bigint& z04=sqr(z02);
	const bigint& x02=sqr(x0); const bigint& x03=x0*x02;
	q=y0;
	e=z04*a;
	dd=z03*(4*a*x0 + b*z0);
	cc=z02*(6*a*x02 + 3*b*x0*z0 + c*z02);
	bb=z0*(4*a*x03 + 3*b*x02*z0 + 2*c*x0*z02 + d*z03);
	a=sqr(y0);
	b=bb; c=cc; d=dd;
      }

#ifdef DEBUG
  if(verbose) cout<<"Quartic transformed = (a,b,c,d,e) = ("<<a<<", "<<b<<", "<<c<<", "<<d<<", "<<e<<")\nq =  "<<q << "\n";
#endif

  p = 3*b*b-8*a*c;
  r = b*b*b+8*a*a*d-4*a*b*c;
  xp = 3*p, yp = 27*r, zp = 2*q*z02;  // The z0^2 factor since (I,J) have
                                      // been scaled up

  Point oldP(IJ_curve,xp*zp,yp,pow(zp,3));
  int valid;
#ifdef DEBUG
  valid = oldP.isvalid();
  if(verbose||!valid) cout << "Point "<<oldP<<" on IJ_curve "
         <<(Curve)(*IJ_curve);
  if(!valid) cout <<" --NOT OK\n";
  if(verbose) cout<<"\n";
#endif

  P = transform(oldP, E, tr_u, tr_r, tr_s, tr_t);

  valid = P.isvalid();
  if(verbose||!valid) cout<<"Point = "<<P;
  if(!valid) {cout << " -- warning: NOT on curve " << (Curve)(*E); abort();}
  if(verbose) cout << "\n\theight = " << height(P)<< "\n";
}


