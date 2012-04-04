// xsplit.h: Declaration of class form_finder
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
 
#if     !defined(_XSPLIT_H)
#define _XSPLIT_H      1       //flags that this file has been included

#include "method.h"  // #defines form_finder=form_finder0/1/2/3/4
#include "splitbase.h"

// flags set on construction:

// dual:  use dual bases throughout (the usual case)
// plusflag: if 1 then targetdim=1 and conjmat not used
//           if 0 then targetdim=2 and conjmat is used

// bigmats: if 1 then uses opmat() to get operators on ambient space,
// does its own restriction to the current subspace; if 0 then uses
// opmat_restricted() to get operators restricted to current space
// directly -- ONLY possible when dual=1

class form_finder {
protected:
  splitter_base* h;
  int plusflag, dual, bigmats, verbose, targetdim;
  long maxdepth, mindepth, depth, subdim, dimen;
  SCALAR denom1;
  ssubspace** nest;     // array of pointers to subspaces
// "Current" subspace is *nest[depth] of dimension subdim
  vector<long> eiglist;
  vec bplus, bminus;

  int *havemat;
  char** opfilenames;  // temp filenames
  smat conjmat;  // only used if plus==0 and bigmats==1
  smat the_opmat;
  smat *submats;  // holds current restriction for i>0
  void make_opmat(long i);  // puts it in the_opmat
  void make_submat();
  void go_down(long eig, int last=0);
  void go_up();
  void make_basis();
  vec getbasis1(const ssubspace* s); //assuming dim(s)=1, get basis vector
public:
  form_finder(splitter_base* hh, int plus, int maxd, int mind=0, 
              int dualflag=1, int bigmatsflag=0, int v=0);
  ~form_finder(void); 
  void find();
  void recover(vector< vector<long> > eigs);
  void splitoff(const vector<long>& eigs);
  vec getbasis() const {return bplus;}
  vec getbasisplus() const {return bplus;}
  vec getbasisminus() const {return bminus;}
};

#endif
