// cperiods.h: declarations of class Cperiods & period lattice functions
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
 
#ifndef _CPERIODS_H_
#define _CPERIODS_H_

#include "compproc.h"
#include "curve.h"

#define TWOPI 2*Pi()
#define TWOPIEYE bigcomplex(to_bigfloat(0), TWOPI)


inline bigcomplex q(const bigcomplex& z) // q(z) = exp(2*pi*i * z)
{
  bigfloat twopix = TWOPI * real(z);
  return exp(-TWOPI * imag(z)) *
         bigcomplex( cos(twopix), sin(twopix) );
}

//
// the Cperiods class:
//

// w1, w2 are a lattice basis with tau=w1/w2 in fundamental region
// wR = least (positive) real period
// wI = least imaginary period
// wRI = either wI    (iff lattice_type==2) 
//         or (wR+wI)/2 (iff lattice_type==1) 

class Cperiods {
  bigcomplex w1, w2, tau;
  bigcomplex wR, wI, wRI;
  bigcomplex e1, e2, e3;  // 2-division values
  int lattice_type;           // 2 for rectangular, 1 for triangular 
  bigcomplex qtau, w1squared, w1cubed, sum3;
  void store_sums();  // sets quantities on previous line
public:
  Cperiods() 
    : w1(to_bigfloat(0)), w2(to_bigfloat(0)), tau(to_bigfloat(0)), 
      wR(to_bigfloat(0)), wI(to_bigfloat(0)), wRI(to_bigfloat(0)), 
      lattice_type(0) 
  {;}
  Cperiods(bigfloat x, bigfloat y, int type)
    : lattice_type(type)
  { 
    if (type==1) 
      {
	wR=2*x; 
	wI=bigcomplex(to_bigfloat(0),2*y); 
	wRI=bigcomplex(x,y); 
      }
    else 
      {
	wR=x; 
	wI=wRI=bigcomplex(to_bigfloat(0),y); 
      }
    w1=wR; w2=wRI;
    tau = normalize(w2,w1);  // NB reverse params;  from compproc.h
    store_sums();
  }
  Cperiods(const Curvedata& E); 

  // copying:
  Cperiods(const Cperiods& cp)
    : w1(cp.w1), w2(cp.w2), tau(cp.tau), 
      wR(cp.wR), wI(cp.wI), wRI(cp.wRI), 
      e1(cp.e1), e2(cp.e2), e3(cp.e3), 
      lattice_type(cp.lattice_type),
      qtau(cp.qtau), w1squared(cp.w1squared), w1cubed(cp.w1cubed), 
      sum3(cp.sum3)
    {;}
  void operator=(const Cperiods& cp)
    { w1=cp.w1; w2=cp.w2; tau=cp.tau; 
    wR=cp.wR; wI=cp.wI; wRI=cp.wRI; 
    e1=cp.e1; e2=cp.e2; e3=cp.e3; 
    lattice_type=cp.lattice_type; 
    qtau=cp.qtau; w1squared=cp.w1squared; w1cubed=cp.w1cubed; 
    sum3=cp.sum3;
    }

  // member access functions
  friend inline bigcomplex gettau(const Cperiods& cp) {return cp.tau; }
  friend inline int get_lattice_type(const Cperiods& cp) {return cp.lattice_type; }
  bigcomplex get_real_period() const {return wR;}
  int getwi(bigcomplex& ww1, bigcomplex& ww2) const
    { ww1=w1; ww2=w2; return lattice_type; }
  int getwRI(bigcomplex& wr, bigcomplex& wri) const
    { wr=wR; wri=wRI; return lattice_type; }

  friend inline ostream& operator<<(ostream& os, const Cperiods& cp)
    { 
      os<<"[w_1,w_2] = ["<<cp.w1<<","<<cp.w2<<"]\n";
      os<<"tau       = "<<cp.tau<<" (abs(tau)="<<abs(cp.tau)<<")\n";
      switch(cp.lattice_type) {
      case 1:
	{
	  os<<"w_R = "<<cp.wR<<"\tw_IR = "<<cp.wRI<<endl; break;
	}
      case 2:
	{
	  os<<"w_R = "<<cp.wR<<"\tw_I = "<<cp.wI<<endl;
	}
      }
      return os;
    }

 // Weierstrass functions:
  bigcomplex X_coord(const bigcomplex& qz); // qz=q(z), z modulo lattice,
  bigcomplex Y_coord(const bigcomplex& qz); // gives coords in Y^2=4X^3+... model
  void XY_coords(bigcomplex& X, bigcomplex& Y, const bigcomplex& z);

  void makec4c6(bigcomplex& cc4, bigcomplex& cc6) const
  {
    getc4c6(w2,w1,cc4,cc6);  // from compproc.h: note the order!
  }
  Curve trans_to_curve(void) const
  {
    bigcomplex cc4, cc6;
    //    cout<<"Calling getc4c6 with first w2="<<w2<<", second w1="<<w1<<endl;
    getc4c6(w2,w1, cc4, cc6);  // from compproc.h: note the order!
    return Curve(Iround( real(cc4) ), Iround( real(cc6) ));
  }

  bigfloat get_e3() const {return real(e3);} // largest or only real root
  vector<bigcomplex> ellztopoint(const bigcomplex& z, const bigcomplex& a1, const bigcomplex& a2, const bigcomplex& a3);
};  // end of class Cperiods def'n

bigcomplex* solve_nonsingular_cubic(const bigint& c1, const bigint& c2, const bigint& c3); //Returns an array
// Gets the 3 2-division points given the coefficients 
void getei(const Curvedata& E, bigcomplex& e1, bigcomplex& e2, bigcomplex& e3);
// Reorders 3 complex nos so real parts are decreasing
void reorder1(bigcomplex& a, bigcomplex& b, bigcomplex& c);
//reorders 3 complex nos so e1 is real if any (
void reorder2(bigcomplex& e1, bigcomplex& e2, bigcomplex& e3);

#endif

//ends file periods.h
