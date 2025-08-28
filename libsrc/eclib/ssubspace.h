// ssubspace.h: declaration of template class ssubZspace
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2025 John Cremona
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
 
#if     !defined(_ECLIB_SSUBSPACE_H)
#define _ECLIB_SSUBSPACE_H      1       //flags that this file has been included

// The ssubZspace classes differ from the subZspace classes in two ways:
// 1. the basis matrix is sparse (sZmat not Zmat);
// 2. the class has a prime modulus and represents a subspace of F_p^n

#include "smatrix.h"

template<class T>
class ssubZspace {

public:
  // constructors
  ssubZspace<T>(int n=0, T mod=default_modulus<T>());
  ssubZspace<T>(const sZmat<T>& b, const Zvec<int>& p, T mod);
  ssubZspace<T>(const ssubZspace<T>& s);
  // assignment
  void operator=(const ssubZspace<T>& s);

  // member functions & operators
  inline void clear() { pivots.init(); basis=sZmat<T>(0,0);}
  inline Zvec<int> pivs() const {return pivots;}  // the pivot vector
  inline sZmat<T> bas() const {return basis;}   // the basis matrix
  inline T mod() const {return modulus;}   // the (prime) modulus

  // non-member (friend) functions and operators
  friend ssubZspace<T> combine<>(const ssubZspace<T>& s1, const ssubZspace<T>& s2);
  friend sZmat<T> restrict_mat<>(const sZmat<T>& m, const ssubZspace<T>& s);

  // Implementation
private:
  T modulus;
  Zvec<int> pivots;
  sZmat<T> basis;
};

// Declarations of nonmember, nonfriend operators and functions:

template<class T> inline Zvec<int> pivots(const ssubZspace<T>& s)  {return s.pivs();}

template<class T>
ssubZspace<T> kernel(const sZmat<T>& sm, T m);
template<class T>
ssubZspace<T> eigenspace(const sZmat<T>& sm, T lambda, T m);
template<class T>
ssubZspace<T> subeigenspace(const sZmat<T>& sm, T l, const ssubZspace<T>& s, T m);
// construction of a 1-dimensional sparse subspace from a vector:
template<class T>
ssubZspace<T> make1d(const Zvec<T>& bas, T& piv, T m);

#endif
