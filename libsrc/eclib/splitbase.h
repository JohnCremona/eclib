// splitbase.h: Declaration of class splitter_base
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
 
#if     !defined(_ECLIB_SPLITBASE_H)
#define _ECLIB_SPLITBASE_H      1       //flags that this file has been included

#include "types.h"

// The following class must be a base class of any class wishing to 
// use the form_finder class; the using class must implement some or 
// all of these virtual functions before creating a form_finder with 
// "this" as first parameter in the constructor:
//
// In all cases implement opmat, matdim, matden, use;  implement eigrange
// iff you are going to use the recursive search function find()
// 
// opmat(i) returns the i'th operator for i>=0; opmat(-1) should return
// the conjugation matrix in the case where target dimension is 2 (i.e.
// plusflag=0) and basis vectors bplus, bminus are required.
//
// matdim() returns the size of the matrices, i.e. dimension of ambient space
//
// matden() returns the implicit denominator of all matrices, needed to scale 
//          eigenvalues
//
// eigrange(i) returns a list of possible eigenvalues for opmat(i)
//
// use(basis1, basis2, eiglist) provides whatever processing is done with the
//                              1-dimensional eigenspaces as found.
//                             (basis2 will not be used when plusflag=1)

#include "types.h"

template<class T>
class splitter_base {
public:
  virtual Zmat<T> opmat(int,int,int=0) = 0;
  virtual Zvec<T> opmat_col(int i, int j, int verb=0) = 0;
  virtual Zmat<T> opmat_cols(int i, const vec_i& jlist, int verb=0) = 0;
  virtual Zmat<T> opmat_restricted(int,const subZspace<T>& s, int,int=0) = 0;
  virtual sZmat<T> s_opmat(int,int,int=0) = 0;
  virtual sZvec<T> s_opmat_col(int i, int j, int verb=0) = 0;
  virtual sZmat<T> s_opmat_cols(int i, const vec_i& jlist, int verb=0) = 0;
  virtual sZmat<T> s_opmat_restricted(int,const ssubZspace<T>& s, int, int=0) = 0;
  virtual long matdim(void) = 0;
  virtual T matden(void) = 0;
  virtual vector<long> eigrange(int) = 0;
  virtual long dimoldpart(const vector<long>) = 0;
  virtual void use(const Zvec<T>&, const Zvec<T>&, const vector<long>) = 0;
  virtual ~splitter_base() {;}
};

#endif
