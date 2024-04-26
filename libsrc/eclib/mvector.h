// mvector.h: declarations of multiprecision integer vector class
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
 
#if     !defined(_ECLIB_MVECTOR_H)
#define _ECLIB_MVECTOR_H      1       //flags that this file has been included

#include "vector.h"

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
#define svec svec_m
#define smat smat_m
#define smat_elim smat_m_elim

#include "vec.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

vec_m to_vec_m(const vec_i& v);
vec_m to_vec_m(const vec_l& v);
vec_i to_vec_i(const vec_m& v);
vec_l to_vec_l(const vec_m& v);

#if(0)

class vec_m {
friend class mat_m;
friend class msubspace;
public:
  // constructors
  explicit vec_m(long n=0);
  explicit vec_m(const vector<bigint>& array) :entries(array) {};
  vec_m(const vec_m& v) :entries(v.entries) {};           // copy constructor
  explicit vec_m(const vec_i&);                           // copy int vector
  explicit vec_m(const vec_l&);                           // copy long vector
  // member functions & operators
  void init(long n=0);                 // (re)-initializes
  vec_m& operator=(const vec_m&);         // assignment
  bigint& operator[](long i);             // the i'th component
  bigint operator[](long i) const;        // the i'th component
  vec_m& operator+=(const vec_m&);
  void addmodp(const vec_m&, const bigint&);
  vec_m& operator-=(const vec_m&);
  vec_m& operator*=(const bigint&);
  vec_m& operator/=(const bigint&);
  vec_m slice(long,long=-1) const;           // returns subvec_m
  vec_m operator[](const vec_i&) const;  // subscript composition
  vec_m operator[](const vec_l&) const; // subscript composition
  void set(long i, const bigint& x);                  // sets v[i]=x
  void add(long i, const bigint& x);                  // v[i]+=x
  bigint sub(long i) const;                    // same as v[i] (no ref)
  //converts to a vector of ints or longs if possible.
  //the parameter here is a dummy just to distinguish these
  vec_i shorten(int) const;
  vec_l shorten(long) const;

  // non-member (friend) functions and operators
  friend long dim(const vec_m&);                  // the dimension
  friend bigint operator*(const vec_m&, const vec_m&);   // dot product
  friend vec_m operator*(const mat_m& m, const vec_m& v);
  friend int operator==(const vec_m&, const vec_m&);
  friend int operator!=(const vec_m&, const vec_m&);
  friend int trivial(const vec_m&);                  // v==zerovec_m?
  friend ostream& operator<< (ostream&s, const vec_m&);
  friend istream& operator>> (istream&s, vec_m&);
  friend bigint content(const vec_m&);
  friend void swapvec(vec_m& v, vec_m& w);
  friend int member(const bigint& a, const vec_m& v);//tests if a=v[i] for some i
  friend mat_m restrict_mat(const mat_m& m, const msubspace& s);
  friend mat_m prestrict(const mat_m& m, const msubspace& s, const bigint& pr);
  friend bigint dotmodp(const vec_m& v1, const vec_m& v2, const bigint& pr);

// Implementation
private:
       vector<bigint> entries;
};

// Declaration of non-member, non-friend functions

vec_m operator+(const vec_m&);                   // unary
vec_m operator-(const vec_m&);                   // unary
vec_m operator+(const vec_m&, const vec_m&);
vec_m addmodp(const vec_m&, const vec_m&, const bigint&);
vec_m operator-(const vec_m&, const vec_m&);
vec_m operator*(const bigint&, const vec_m&);       // componentwise
vec_m operator*(long, const vec_m&);                 // componentwise
vec_m operator/(const vec_m&, const bigint&);       // componentwise
void makeprimitive(vec_m& v);
void elim(const vec_m& a, vec_m& b, long pos);
void elim1(const vec_m& a, vec_m& b, long pos);
void elim2(const vec_m& a, vec_m& b, long pos, const bigint& lastpivot);
vec_m express(const vec_m& v, const vec_m& v1, const vec_m& v2);
vec_m lift(const vec_m& v, const bigint& pr);
  //lifts a mod-p vec_m to a rational and scales to a primitive vec in Z.
int liftok(vec_m& v, const bigint& pr);
  //same, in place; returns success

// extra function definitions

long dim(const vec_m& v);
int operator!=(const vec_m& v, const vec_m& w);
vec_m operator+(const vec_m& v);
vec_m operator-(const vec_m& v);
vec_m operator+(const vec_m& v1, const vec_m& v2);
vec_m addmodp(const vec_m& v1, const vec_m& v2, const bigint& pr);
vec_m operator-(const vec_m& v1, const vec_m& v2);
vec_m operator*(const bigint& scal, const vec_m& v);
vec_m operator*(long scal, const vec_m& v);
vec_m operator/(const vec_m& v, const bigint& scal);
void makeprimitive(vec_m& v);
void elim(const vec_m& a, vec_m& b, long pos);
void elim1(const vec_m& a, vec_m& b, long pos);
void elim2(const vec_m& a, vec_m& b, long pos, const bigint& lastpivot);

#endif // end of #if(0)

#endif
