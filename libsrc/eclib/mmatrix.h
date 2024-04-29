// mmatrix.h: declarations of multiprecision integer matrix class
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
 
#if     !defined(_ECLIB_MMATRIX_H)
#define _ECLIB_MMATRIX_H      1       //flags that this file has been included

#include "mvector.h"
#include "matrix.h"

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

#include "mat.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

mat_m to_mat_m(const mat_i& v);
mat_m to_mat_m(const mat_l& v);
mat_i to_mat_i(const mat_m& v);
mat_l to_mat_l(const mat_m& v);

#if(0)

int liftmat(const mat_m& mm, const bigint& pr, mat_m& m, bigint& dd, int trace=0);
int lift(const msubspace& s, const bigint& pr, msubspace& ans, int trace=0);

class mat_m {
  friend class msubspace;
public:
  // constructors
  mat_m(long nr=0, long nc=0)   :nro(nr), nco(nc) {entries.resize(nr*nc, bigint(0));}
  mat_m(long nr, long nc, const vector<bigint> ent) : nro(nr), nco(nc), entries(ent) {}
  mat_m(const mat_m& m) :nro(m.nro), nco(m.nco){ entries = m.entries;}  // copy constructor
  explicit mat_m(const mat_i&);
  explicit mat_m(const mat_l&);

  // member functions & operators
  void init(long nr=0, long nc=0);
  mat_m& operator=(const mat_m&);            // assignment with copy
  bigint& operator()(long i, long j);        // returns ref to (i,j) entry
  bigint operator()(long i, long j) const;   // returns (i,j) entry
  mat_m slice(long,long,long=-1,long=-1) const;// returns submatrix
  bigint sub(long i, long j) const;             // returns the (i,j) entry
  vec_m row(long i) const;                // returns row i (as a vector)
  vec_m col(long j) const;                // returns col j (as a vector)
  void set(long i, long j, const bigint& x);     // sets the (i,j) entry to x
  void add(long i, long j, const bigint& x);  // adds x to the (i,j) entry
  void setrow(long i, const vec_m& v);
  void setcol(long i, const vec_m& v);
  void swaprows(long r1, long r2);
  void multrow(long r, const bigint& scal);
  void divrow(long r, const bigint& scal);
  bigint content() const;
  bigint row_content(long r) const;
  void clearrow(long r);
  void makeprimitive();
  void operator+=(const mat_m&);
  void operator-=(const mat_m&);
  void operator*=(const bigint&);
  void operator/=(const bigint&);
  // shortens to a matrix o ints or longs if possible
  //the parameter here is a dummy just to distinguish these
  mat_i  shorten(int) const;
  mat_l shorten(long) const;
  long nrows() const {return nro;}
  long ncols() const {return nco;}
  long rank() const;
  long nullity() const;
  bigint trace() const;
  vector<bigint> charpoly() const;
  bigint determinant() const;

  // non-member (friend) functions and operators
  friend mat_m operator*(const mat_m&, const mat_m&);
  friend vec_m operator*(const mat_m&, const vec_m&);
  friend int operator==(const mat_m&, const mat_m&);
  friend ostream& operator<< (ostream&s, const mat_m&);
  friend istream& operator>> (istream&s, mat_m&);
  friend mat_m colcat(const mat_m& a, const mat_m& b);
  friend mat_m rowcat(const mat_m& a, const mat_m& b);
  friend mat_m directsum(const mat_m& a, const mat_m& b);
  friend void elimrows(mat_m& m, long r1, long r2, long pos); //plain elimination, no clearing
  friend void elimrows1(mat_m& m, long r1, long r2, long pos); //elimination + clearing
  friend void elimrows2(mat_m& m, long r1, long r2, long pos, const bigint& last); //elimination + divide by last pivot
  friend mat_m echelon0(const mat_m& m, vec_i& pcols, vec_i& npcols,
                        long& rk, long& ny, bigint& d);
  friend void elimp(mat_m& m, long r1, long r2, long pos, const bigint& pr);
  friend mat_m echelonp(const mat_m& m1, vec_i& pcols, vec_i& npcols,
                        long& rk, long& ny, bigint& d, const bigint& pr);
  friend mat_m echmodp(const mat_m& m, vec_i& pcols, vec_i& npcols,
                       long& rk, long& ny, const bigint& pr);
  friend msubspace combine(const msubspace& s1, const msubspace& s2);
  friend mat_m restrict_mat(const mat_m& m, const msubspace& s);
  friend int liftmat(const mat_m& mm, const bigint& pr, mat_m& m, bigint& dd, int trace);
  friend int lift(const msubspace& s, const bigint& pr, msubspace& ans, int trace);
  friend msubspace pcombine(const msubspace& s1, const msubspace& s2, const bigint& pr);
  friend mat_m prestrict(const mat_m& m, const msubspace& s, const bigint& pr);
  friend mat_m matmulmodp(const mat_m&, const mat_m&, const bigint& pr);

  // Implementation
private:
  long nro,nco;
  vector<bigint> entries;  // stored in one array of size nro+nco, by rows
};

// Declaration of non-friend functions

mat_m operator+(const mat_m&);                   // unary
mat_m operator-(const mat_m&);                   // unary
mat_m operator+(const mat_m& m1, const mat_m& m2);
mat_m operator-(const mat_m& m1, const mat_m& m2);
mat_m operator*(const bigint& scal, const mat_m& m);
mat_m operator/(const mat_m& m, const bigint& scal);
int operator!=(const mat_m& m1, const mat_m& m2);
mat_m midmat(long n);  // = multi-idmat
mat_m transpose(const mat_m& m);
// submatrix of rows indexed by v, all columns
mat_m rowsubmat(const mat_m& m, const vec_i& v);
// submatrix of rows indexed by iv, columns indexed by jv
mat_m submatrix(const mat_m& m, const vec_i& iv, const vec_i& jv);
mat_m echelon(const mat_m& m, vec_i& pcols, vec_i& npcols,
                          long& rk, long& ny, bigint& d, int method=0);
mat_m echelon(const mat_m& m, vec_l& pcols, vec_l& npcols,
                          long& rk, long& ny, bigint& d, int method=0);
mat_m addscalar(const mat_m&, const bigint&);
vec_m apply(const mat_m&, const vec_m&);

#endif // end of #if(0)

#endif
