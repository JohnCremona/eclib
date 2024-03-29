// sifter.h: declaration of class for sifting E(Q)/2E(Q)
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
 
// NB This is used for proving that points are independent; now
// largely obsolete, being superceded by general saturation algorithms


// allow for multiple includes
#ifndef _ECLIB_SIFTER_H
#define _ECLIB_SIFTER_H

#include <eclib/curve.h>

class sifter {
private: 
  Curvedata *E;
  bigint I, J, disc;
  bigint r,s,t;    // tranforms E to E_{I,J} (u=6)
  int rank;
  int verbose;

  int num_aux, max_dim_im;
  vector<vector<int>> eps_mat, squares;
  vector<int> pivcols, nroots;
  vector<long> auxs, all_p;
  vector<vector<long>> thetamod;
public:
  sifter(Curvedata* EE, int na, int verb=0);
  int code(const bigint& x, const bigint& z2, int i);
  vector<int> eps(const bigint& x, const bigint& z2);
  void process(const Point& P);
  void process(const vector<Point>& Plist);
  int getrank() {return rank;}
  void vecout(const vector<int>& v);
};

#endif
