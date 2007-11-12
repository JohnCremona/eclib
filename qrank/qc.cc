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

// New version:

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

  P = shift(oldP, E, tr_u, tr_r, tr_s, tr_t);

  valid = P.isvalid();
  if(verbose||!valid) cout<<"Point = "<<P;
  if(!valid) {cout << " -- warning: NOT on curve " << (Curve)(*E); abort();}
  if(verbose) cout << "\n\theight = " << height(P)<< "\n";
}

/* Old version

void qc(quartic& g,
        const bigint& x0,  const bigint& y0,  const bigint& z0,
        Curvedata * E, Point& P, int verbose, int have_curve)
{
#ifdef DEBUG
  cout  << "(x0:y0:z0) = ("<<x0<<" : "<<y0<<" : "<<z0<<")\n";
#endif
  bigint a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
#ifdef DEBUG
	if(verbose) cout<<"Quartic = (a,b,c,d,e) = ("<<a<<", "<<b<<", "<<c<<", "<<d<<", "<<e<<")\n";
	if(verbose) cout  << "(x0:y0:z0) = ("<<x0<<" : "<<y0<<" : "<<z0<<")\n";
#endif
  bigint q,f,aa,bb,cc,dd,ee;

  if(isqrt(a,q)) {ee=a; dd=b; cc=c; bb=d; aa=e; }
  else if(isqrt(e,q)){aa=a; bb=b; cc=c; dd=d; ee=e; }
    else 
      {
	q=y0;
	aa=pow(z0,4)*a;
	bb=pow(z0,3)*(4*a*x0 + b*z0);
	cc=pow(z0,2)*(6*a*x0*x0 + 3*b*x0*z0 + c*z0*z0);
	dd=z0*(4*a*pow(x0,3) + 3*b*x0*x0*z0 + 2*c*x0*z0*z0 + d*pow(z0,3));
	ee=y0*y0;
      }
#ifdef DEBUG
	if(verbose) cout<<"Quartic transformed = (a,b,c,d,e) = ("<<aa<<", "<<bb<<", "<<cc<<", "<<dd<<", "<<ee<<")\nq =  "<<q << "\n";
#endif
  bigint xp = -4*cc*q*q+dd*dd;

  bigint a1 = 2*dd, 
         a2 = 4*cc*ee-dd*dd, 
         a3 = 16*ee*ee*bb, 
         a4 = -64*ee*ee*ee*aa;
  bigint a6 = a2*a4;
  bigint b2 = a1*a1 + 4*a2; 
  bigint b4 = 2*a4 + a1*a3;
  bigint b6 = a3*a3 + 4*a6;
  bigint c4 = b2*b2 - 24*b4; 
  bigint c6 = -b2*b2*b2 + 36*b2*b4 - 216*b6;
  bigint newa1,newa2,newa3,newa4,newa6,newc4, newc6, newdiscr, u, u2;
#ifdef DEBUG
  cout<<"Original c4, c6 = " << c4 << ", " << c6 << "\n";
#endif
  minimise_c4c6(c4,c6,0,newc4,newc6,newdiscr,u);
#ifdef DEBUG
  cout<<"Minimal c4, c6 = " << newc4 << ", " << newc6 << "\n";
#endif
  c4c6_to_ai(newc4,newc6,newa1,newa2,newa3,newa4,newa6);
#ifdef DEBUG
  cout<<"Minimal ai = " << newa1 << ", " << newa2 << ", " 
      << newa3 << ", " << newa4 << ", " << newa6 << "\n";
#endif

  Curve newE(newa1,newa2,newa3,newa4,newa6);
  if(have_curve) // E already points to the curve point should be on
    {
      if(newE!=(*E))
	{
	  cout << "Warning: constructed curve " << newE <<"\n";
	  cout << " not the same as original curve " << (Curve)(*E) << endl;
	}
    }
  else
    {
      (*E) = newE;  // copy curve just constructed into given Curvedata*
    }

  bigint s = (u*newa1 - a1)/2;  mulx(u,u,u2);
  bigint r = (u2*newa2 - a2 + s*a1 + s*s)/3; 
  bigint t = (u2*u*newa3 - a3 - r*a1)/2;
#ifdef DEBUG
  if(verbose) cout<<"[u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]\n";
#endif

  bigint newx, newy, newz;
  raw_shift_point(xp, 0, 1, u, r, s, t, newx, newy, newz);

  P.init(E, newx*newz, newy, pow(newz,3));

  int valid = P.isvalid();
  if(verbose||!valid) cout<<"Point = "<<P << "\n";
  if(!valid) cout << " -- warning: NOT on curve \n";
}

*/ // end of old version


