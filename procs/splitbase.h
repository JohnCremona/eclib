// splitbase.h: Declaration of class splitter_base
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
 
#if     !defined(_SPLITBASE_H)
#define _SPLITBASE_H      1       //flags that this file has been included

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

class splitter_base {
public: 
  virtual mat opmat(int,int,int=0) = 0;
  virtual mat opmat_restricted(int,const subspace& s, int,int=0) = 0;
  virtual smat s_opmat(int,int,int=0) = 0;
  virtual smat s_opmat_restricted(int,const ssubspace& s, int, int=0) = 0;
  virtual long matdim(void) = 0;
  virtual long matden(void) = 0;
  virtual vector<long> eigrange(int) = 0;
  virtual long dimoldpart(const vector<long>) = 0;
  virtual void use(const vec&, const vec&, const vector<long>) = 0;
  virtual ~splitter_base() {;}
};


#endif
