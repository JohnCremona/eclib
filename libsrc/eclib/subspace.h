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

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i
#define svec svec_i
#define smat smat_i
#define smat_elim smat_i_elim

#include "sub.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l
#define svec svec_l
#define smat smat_l
#define smat_elim smat_l_elim

#include "sub.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m

#include "sub.h"

#undef scalar
#undef vec
#undef mat
#undef subspace

///////////////////////////////////////////////////////////////////////////

template<class T> class vecT;
template<class T> class svecT;
template<class T> class smatT;
template<class T> class smatT_elim;
template<class T> class matT;
template<class T> class subspaceT;

template<class T> int dim(const subspaceT<T>& s);      // the dimension
template<class T> T denom(const subspaceT<T>& s);   // the denominator
template<class T> vecT<int> pivots(const subspaceT<T>& s);// the pivot vector
template<class T> matT<T> basis(const subspaceT<T>& s) ;// the basis matrix
template<class T> subspaceT<T> combine(const subspaceT<T>& s1, const subspaceT<T>& s2);
template<class T> matT<T> restrict_mat(const matT<T>& m, const subspaceT<T>& s, int cr);
template<class T> subspaceT<T> pcombine(const subspaceT<T>& s1, const subspaceT<T>& s2, const T& pr);
template<class T> matT<T> prestrict(const matT<T>& m, const subspaceT<T>& s, const T& pr, int cr);
template<class T> int lift(const subspaceT<T>& s, const T& pr, subspaceT<T>& ans);


template<class T>
class subspaceT {

public:
  // constructors
  subspaceT(int n=0) :denom(1),pivots(vecT<int>::iota(n)),basis(matT<T>::identity_matrix(n)) {}
  subspaceT(const matT<T>& b, const vecT<int>& p, const T& d) :denom(d),pivots(p),basis(b) {}
  subspaceT(const subspaceT<T>& s) :denom(s.denom),pivots(s.pivots),basis(s.basis) {}

  // assignment
  void operator=(const subspaceT<T>& s);

  // member functions & operators
  inline void clear() { pivots.init(); basis.init();}
  inline T den() const {return denom;}     // the denominator
  inline vecT<int> pivs() const {return pivots;} // the pivot vector
  inline matT<T> bas() const {return basis;}   // the basis matrix

  // non-member (friend) functions and operators
  friend int dim<>(const subspaceT<T>& s);      // the dimension
  friend T denom<>(const subspaceT<T>& s);   // the denominator
  friend vecT<int> pivots<>(const subspaceT<T>& s);// the pivot vector
  friend matT<T> basis<>(const subspaceT<T>& s) ;// the basis matrix
  friend subspaceT<T> combine<>(const subspaceT<T>& s1, const subspaceT<T>& s2);
  friend matT<T> restrict_mat<>(const matT<T>& m, const subspaceT<T>& s, int cr);
  friend subspaceT<T> pcombine<>(const subspaceT<T>& s1, const subspaceT<T>& s2, const T& pr);
  friend matT<T> prestrict<>(const matT<T>& m, const subspaceT<T>& s, const T& pr, int cr);
  friend int lift<>(const subspaceT<T>& s, const T& pr, subspaceT<T>& ans);

  // Implementation
private:
  T   denom;
  vecT<int> pivots;
  matT<T> basis;
};

// Declarations of nonmember, nonfriend operators and functions:

template<class T>
matT<T> expressvectors(const matT<T>& m, const subspaceT<T>& s);
template<class T>
subspaceT<T> kernel(const matT<T>& m, int method=0);
template<class T>
subspaceT<T> image(const matT<T>& m, int method=0);
template<class T>
subspaceT<T> eigenspace(const matT<T>& m, const T& lambda, int method=0);
template<class T>
subspaceT<T> subeigenspace(const matT<T>& m, const T& l, const subspaceT<T>& s, int method=0);


//The following work with subspaces "mod p" using "echmodp" from
//matrix.h/cc to do gaussian elimination.  The "denom" of each is 1.

template<class T>
subspaceT<T> oldpkernel(const matT<T>& m, const T& pr);
template<class T>
subspaceT<T> pkernel(const matT<T>& m, const T& pr);
template<class T>
subspaceT<T> pimage(const matT<T>& m, const T& pr);
template<class T>
subspaceT<T> peigenspace(const matT<T>& m, const T& lambda, const T& pr);
template<class T>
subspaceT<T> psubeigenspace(const matT<T>& m, const T& l, const subspaceT<T>& s, const T& pr);


template<class T>
inline int dim(const subspaceT<T>& s) {return s.basis.ncols();}  // the dimension
template<class T>
inline T denom(const subspaceT<T>& s) {return s.denom;}    // the denominator
template<class T>
inline vecT<int> pivots(const subspaceT<T>& s) {return s.pivots;}     // the pivot vector
template<class T>
inline matT<T> basis(const subspaceT<T>& s) {return s.basis;}       // the basis matrix


#endif
