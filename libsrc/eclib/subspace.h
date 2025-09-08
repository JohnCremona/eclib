// subspace.h: manage declarations of subspace classes
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
 
#if     !defined(_ECLIB_SUBSPACE_H)
#define _ECLIB_SUBSPACE_H      1       //flags that this file has been included

#include "matrix.h"

template<class T>
class subZspace {

public:
  // constructors
  subZspace(int n=0) :denom(1),pivots(Zvec<int>::iota(n)),basis(Zmat<T>::identity_matrix(n)) {}
  subZspace(const Zmat<T>& b, const Zvec<int>& p, const T& d) :denom(d),pivots(p),basis(b) {}
  subZspace(const subZspace<T>& s) :denom(s.denom),pivots(s.pivots),basis(s.basis) {}

  // assignment
  void operator=(const subZspace<T>& s);

  // member functions & operators
  inline void clear() { pivots.init(); basis.init();}
  inline T den() const {return denom;}     // the denominator
  inline Zvec<int> pivs() const {return pivots;} // the pivot vector
  inline Zmat<T> bas() const {return basis;}   // the basis matrix

  // non-member (friend) functions and operators
  friend int dim<>(const subZspace<T>& s);      // the dimension
  friend T denom<>(const subZspace<T>& s);   // the denominator
  friend Zvec<int> pivots<>(const subZspace<T>& s);// the pivot vector
  friend Zmat<T> basis<>(const subZspace<T>& s) ;// the basis matrix
  friend subZspace<T> combine<>(const subZspace<T>& s1, const subZspace<T>& s2);
  friend Zmat<T> restrict_mat<>(const Zmat<T>& m, const subZspace<T>& s, int cr);
  friend subZspace<T> pcombine<>(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr);
  friend Zmat<T> prestrict<>(const Zmat<T>& m, const subZspace<T>& s, const T& pr, int cr);
  friend int lift<>(const subZspace<T>& s, const T& pr, subZspace<T>& ans);

  // Implementation
private:
  T   denom;
  Zvec<int> pivots;
  Zmat<T> basis;
};

// Declarations of nonmember, nonfriend operators and functions:

template<class T>
Zmat<T> expressvectors(const Zmat<T>& m, const subZspace<T>& s);
template<class T>
subZspace<T> kernel(const Zmat<T>& m, int method=0);
template<class T>
subZspace<T> image(const Zmat<T>& m, int method=0);
template<class T>
subZspace<T> eigenspace(const Zmat<T>& m, const T& lambda, int method=0);
template<class T>
subZspace<T> subeigenspace(const Zmat<T>& m, const T& l, const subZspace<T>& s, int method=0);


//The following work with subspaces "mod p" using "echmodp" from
//matrix.h/cc to do gaussian elimination.  The "denom" of each is 1.

template<class T>
subZspace<T> oldpkernel(const Zmat<T>& m, const T& pr);
template<class T>
subZspace<T> pkernel(const Zmat<T>& m, const T& pr);
template<class T>
subZspace<T> pimage(const Zmat<T>& m, const T& pr);
template<class T>
subZspace<T> peigenspace(const Zmat<T>& m, const T& lambda, const T& pr);
template<class T>
subZspace<T> psubeigenspace(const Zmat<T>& m, const T& l, const subZspace<T>& s, const T& pr);


template<class T>
inline int dim(const subZspace<T>& s) {return s.basis.ncols();}  // the dimension
template<class T>
inline T denom(const subZspace<T>& s) {return s.denom;}    // the denominator
template<class T>
inline Zvec<int> pivots(const subZspace<T>& s) {return s.pivots;}     // the pivot vector
template<class T>
inline Zmat<T> basis(const subZspace<T>& s) {return s.basis;}       // the basis matrix

// return a basis for the orthogonal complement of a<2^r (viewed as a bit vector of length r)
vector<long> dotperp(long a, int r);
// return a basis for the orthogonal complement of the span of a in alist (viewed as bit vectors of length r)
vector<long> dotperp(const vector<long>& alist, int r);

#endif
