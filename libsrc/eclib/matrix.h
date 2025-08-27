// matrix.h:  manage declarations of integer matrix classes
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
 
#if     !defined(_ECLIB_MATRIX_H)
#define _ECLIB_MATRIX_H      1       //flags that this file has been included

#include "vector.h"

template <class T>
class Zmat {
friend class subZspace<T>;
friend class sZmat<T>;
friend class sZvec<T>;
friend class Zvec<T>;
friend class sZmat_elim<T>;
public:
  // constructors
  Zmat(long nr=0, long nc=0)   :nro(nr), nco(nc) {entries.resize(nr*nc, T(0));}
  Zmat(long nr, long nc, const vector<T>& ent) : nro(nr), nco(nc), entries(ent) {}
  Zmat(const Zmat<T>& m)   :nro(m.nro), nco(m.nco), entries(m.entries) {}

  // member functions & operators
  void init(long nr=0, long nc=0);
  Zmat& operator=(const Zmat<T>&);                // assignment with copy
  T& operator()(long i, long j);        // returns ref to (i,j) entry
  T operator()(long i, long j) const;   // returns (i,j) entry
  Zmat slice(long,long,long=-1,long=-1) const;// returns submatrix
  T sub(long i, long j) const;             // returns the (i,j) entry
  Zvec<T> row(long i) const;                // returns row i (as a vec)
  Zvec<T> col(long j) const;                // returns col j (as a vec)
  void set(long i, long j, const T& x);     // sets the (i,j) entry to x
  void add(long i, long j, const T& x);  // adds x to the (i,j) entry
  void setrow(long i, const Zvec<T>& v);
  void setcol(long i, const Zvec<T>& v);
  void swaprows(long r1, long r2);
  void multrow(long r, const T& scal);
  void divrow(long r, const T& scal);
  T content() const;
  T row_content(long r) const;
  void clearrow(long r);
  void make_primitive();
  void operator+=(const Zmat<T>&);
  void operator-=(const Zmat<T>&);
  void operator*=(const T&);
  void operator/=(const T&);
  const vector<T> get_entries()const{return entries;}
  long nrows() const {return nro;}
  long ncols() const {return nco;}
  long rank() const;
  long nullity() const;
  T trace() const;
  vector<T> charpoly() const;
  T determinant() const;
  void output(ostream&s=cout) const;
  void output_pari(ostream&s=cout) const;
  void output_pretty(ostream&s=cout) const;
  void reduce_mod_p(const T& p);

  static Zmat<T> scalar_matrix(long n, const T& a);
  static Zmat<T> identity_matrix(long n) {return scalar_matrix(n, T(1));}

  // non-member (friend) functions and operators
  friend void add_row_to_vec<>(Zvec<T>& v, const Zmat<T>& m, long i);
  friend void sub_row_to_vec<>(Zvec<T>& v, const Zmat<T>& m, long i);
  friend Zmat operator*<>(const Zmat<T>&, const Zmat<T>&);
  friend Zvec<T> operator*<>(const Zmat<T>&, const Zvec<T>&);
  friend int operator==<>(const Zmat<T>&, const Zmat<T>&);
  friend istream& operator>><> (istream&s, Zmat<T>&);
  friend Zmat colcat<>(const Zmat<T>& a, const Zmat<T>& b);
  friend Zmat rowcat<>(const Zmat<T>& a, const Zmat<T>& b);
  friend Zmat directsum<>(const Zmat<T>& a, const Zmat<T>& b);
  friend void elimrows<>(Zmat<T>& m, long r1, long r2, long pos); //plain elimination, no clearing
  friend void elimrows1<>(Zmat<T>& m, long r1, long r2, long pos); //elimination + clearing
  friend void elimrows2<>(Zmat<T>& m, long r1, long r2, long pos, const T& last); //elimination + divide by last pivot
  friend Zmat echelon0<>(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                          long& rk, long& ny, T& d);
  friend void elimp<>(Zmat<T>& m, long r1, long r2, long pos, const T& pr);
  friend void elimp1<>(Zmat<T>& m, long r1, long r2, long pos, const T& pr);
  friend Zmat echelonp<>(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                          long& rk, long& ny, T& d, const T& pr);
  friend Zmat echmodp<>(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const T& pr);
  friend Zmat echmodp_uptri<>(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const T& pr);
  friend Zmat ref_via_flint<>(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const T& pr);
  friend Zmat ref_via_ntl<>(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                             long& rk, long& ny, const T& pr);
  friend Zmat rref<>(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                     long& rk, long& ny, const T& pr);
  friend long rank_via_ntl<>(const Zmat<T>& M, const T& pr);
  friend T det_via_ntl<>(const Zmat<T>& M, const T& pr);
  friend subZspace<T> combine<>(const subZspace<T>& s1, const subZspace<T>& s2);
  friend Zmat restrict_mat<>(const Zmat<T>& m, const subZspace<T>& s, int cr);
  friend int liftmat<>(const Zmat<T>& mm, const T& pr, Zmat<T>& m, T& dd);
  friend int lift<>(const subZspace<T>& s, const T& pr, subZspace<T>& ans);
  friend subZspace<T> pcombine<>(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr);
  friend Zmat prestrict<>(const Zmat<T>& m, const subZspace<T>& s, const T& pr, int cr);
  friend Zmat matmulmodp<>(const Zmat<T>&, const Zmat<T>&, const T& pr);
  friend long population<>(const Zmat<T>& m); // #nonzero entries
  friend T maxabs<>(const Zmat<T>& m); // max entry
  friend double sparsity<>(const Zmat<T>& m); // #nonzero entries/#entries
  // Implementation
private:
  long nro,nco;
  vector<T> entries;  // stored in one array, by rows
};

// Declaration of non-friend functions

template<class T>
inline ostream& operator<< (ostream&s, const Zmat<T>&m)
{m.output(s); return s;}

template<class T>
Zmat<T> operator+(const Zmat<T>&);                   // unary
template<class T>
Zmat<T> operator-(const Zmat<T>&);                   // unary
template<class T>
Zmat<T> operator+(const Zmat<T>& m1, const Zmat<T>& m2);
template<class T>
Zmat<T> operator-(const Zmat<T>& m1, const Zmat<T>& m2);
template<class T>
Zmat<T> operator*(const T& scal, const Zmat<T>& m);
template<class T>
Zmat<T> operator/(const Zmat<T>& m, const T& scal);
template<class T>
int operator!=(const Zmat<T>& m1, const Zmat<T>& m2);
template<class T>
Zmat<T> transpose(const Zmat<T>& m);
// submatrix of only rows indexed by v, all columns
template<class T>
Zmat<T> rowsubmat(const Zmat<T>& m, const Zvec<int>& v);
template<class T>
Zmat<T> rowsubmat(const Zmat<T>& m, const Zvec<long>& v);
// submatrix of rows indexed by iv, columns indexed by jv
template<class T>
Zmat<T> submat(const Zmat<T>& m, const Zvec<int>& iv, const Zvec<int>& jv);
template<class T>
Zmat<T> submat(const Zmat<T>& m, const Zvec<long>& iv, const Zvec<long>& jv);

template<class T>
Zmat<T> echelon(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                          long& rk, long& ny, T& d, int method=0);
template<class T>
Zmat<T> addscalar(const Zmat<T>&, const T&);
template<class T>
Zvec<T> apply(const Zmat<T>&, const Zvec<T>&);

template<class T> Zmat<T> rref(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const T& pr)
{
#if FLINT
  return ref_via_flint(M, pcols, npcols, rk, ny, pr);
#else
  return ref_via_ntl(M, pcols, npcols, rk, ny, pr);
#endif
}

// Construct an NTL mat_lzz_p (matrix mod p) from a mat mod pr
template<class T>
mat_zz_p mat_zz_p_from_mat(const Zmat<T>& M, const T& pr);

// Construct a mat (T type same as pr) from an NTL mat_lzz_p

template<class T>
Zmat<T> mat_from_mat_zz_p(const mat_zz_p& A, const T& pr); // type of T fixes return type

// conversions between matrices of different scalar types

mat_m to_mat_m(const mat_i& m);
mat_m to_mat_m(const mat_l& m);
inline mat_m to_mat_m(const mat_m& m) {return m;}

mat_i to_mat_i(const mat_m& m);
mat_i to_mat_i(const mat_l& m);
inline mat_i to_mat_i(const mat_i& m) {return m;}

mat_l to_mat_l(const mat_m& m);
mat_l to_mat_l(const mat_m& m);
inline mat_l to_mat_l(const mat_l& m) {return m;}

#endif

