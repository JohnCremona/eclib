// sifter.h: declaration of class for sifting E(Q)/2E(Q)
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
 
// NB This is used for proving that points are independent; now
// largely obsolete, being superceded by general saturation algorithms


// allow for multiple includes
#ifndef _SIFTER_
#define _SIFTER_


class sifter {
private: 
  Curvedata *E;
  bigint I, J, disc;
  bigint r,s,t;    // tranforms E to E_{I,J} (u=6)
  int rank;
  int verbose;

  int num_aux, max_dim_im;
  int ** eps_mat;
  int * pivcols;
  long * auxs; long * all_p; int * nroots;  
  long ** thetamod;  int**squares;  
  void init();  // define  auxiliary moduli and squares
  void clear();     // free memory
public:
  sifter(Curvedata* EE, int na, int verb=0)
    :E(EE), rank(0), verbose(verb), num_aux(na)
    {
      I = getc4(*E);
      J = 2*getc6(*E);
      disc = getdiscr(*E);
      E->getai(s,r,t,r,r); // r is just a dummy here
      r = 3*getb2(*E);     // this is its real value
      s = 3*s;
      t = 108*t;
      init();
    }
  ~sifter()
    {
      clear();
    }
  int code(const bigint& x, const bigint& z2, int i);
  int * eps(const bigint& x, const bigint& z2);
  void process(const Point& P);
  void process(const vector<Point>& Plist);
  int getrank() {return rank;}
  void vecout(int* v);
};


#endif
