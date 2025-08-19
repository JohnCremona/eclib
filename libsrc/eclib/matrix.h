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

#include "marith.h"
#include "vector.h"
#include "limits.h" // MAX_INT gcc >= 4.3

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
#include "mat.h"
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
#include "mat.h"
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


///////////////////////////////////////////////////////////////////////////

template<class T> class vecT;
template<class T> class svecT;
template<class T> class smatT;
template<class T> class smatT_elim;
template<class T> class matT;
template<class T> class subspaceT;

template<class T> void add_row_to_vec(vecT<T>& v, const matT<T>& m, long i);
template<class T> void sub_row_to_vec(vecT<T>& v, const matT<T>& m, long i);
template<class T> matT<T> operator*(const matT<T>&, const matT<T>&);
template<class T> vecT<T> operator*(const matT<T>&, const vecT<T>&);
template<class T> int operator==(const matT<T>&, const matT<T>&);
template<class T> istream& operator>> (istream&s, matT<T>&);
template<class T> matT<T> colcat(const matT<T>& a, const matT<T>& b);
template<class T> matT<T> rowcat(const matT<T>& a, const matT<T>& b);
template<class T> matT<T> directsum(const matT<T>& a, const matT<T>& b);
template<class T> void elimrows(matT<T>& m, long r1, long r2, long pos); //plain elimination, no clearing
template<class T> void elimrows1(matT<T>& m, long r1, long r2, long pos); //elimination + clearing
template<class T> void elimrows2(matT<T>& m, long r1, long r2, long pos, const T& last); //elimination + divide by last pivot
template<class T> matT<T> echelon0(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                      long& rk, long& ny, T& d);
template<class T> void elimp(matT<T>& m, long r1, long r2, long pos, const T& pr);
template<class T> void elimp1(matT<T>& m, long r1, long r2, long pos, const T& pr);
template<class T> matT<T> echelonp(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                      long& rk, long& ny, T& d, const T& pr);
template<class T> matT<T> echmodp(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                     long& rk, long& ny, const T& pr);
template<class T> matT<T> echmodp_uptri(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                     long& rk, long& ny, const T& pr);
template<class T> matT<T> ref_via_flint(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                           long& rk, long& ny, const T& pr);
template<class T> matT<T> ref_via_ntl(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                         long& rk, long& ny, const T& pr);
template<class T> matT<T> rref(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                               long& rk, long& ny, const T& pr);
template<class T> long rank_via_ntl(const matT<T>& M, const T& pr);
template<class T> T det_via_ntl(const matT<T>& M, const T& pr);
template<class T> subspaceT<T> combine(const subspaceT<T>& s1, const subspaceT<T>& s2);
template<class T> matT<T> restrict_mat(const matT<T>& m, const subspaceT<T>& s, int cr);
template<class T> int liftmat(const matT<T>& mm, const T& pr, matT<T>& m, T& dd);
template<class T> int lift(const subspaceT<T>& s, const T& pr, subspaceT<T>& ans);
template<class T> subspaceT<T> pcombine(const subspaceT<T>& s1, const subspaceT<T>& s2, const T& pr);
template<class T> matT<T> prestrict(const matT<T>& m, const subspaceT<T>& s, const T& pr, int cr);
template<class T> matT<T> matmulmodp(const matT<T>&, const matT<T>&, const T& pr);
template<class T> matT<T> echmodp_d(const matT<T>& mat, vecT<int>& pcols, vecT<int>& npcols, long& rk, long& ny, double pr);
template<class T> long population(const matT<T>& m); // #nonzero entries
template<class T> T maxabs(const matT<T>& m); // max entry
template<class T> double sparsity(const matT<T>& m); // #nonzero entries/#entries

template <class T>
class matT {
friend class subspaceT<T>;
friend class smatT<T>;
friend class svecT<T>;
friend class vecT<T>;
friend class smatT_elim<T>;
public:
  // constructors
  matT(long nr=0, long nc=0)   :nro(nr), nco(nc) {entries.resize(nr*nc, T(0));}
  matT(long nr, long nc, const vector<T>& ent) : nro(nr), nco(nc), entries(ent) {}
  matT(const matT<T>& m)   :nro(m.nro), nco(m.nco), entries(m.entries) {}

  // member functions & operators
  void init(long nr=0, long nc=0);
  matT& operator=(const matT<T>&);                // assignment with copy
  T& operator()(long i, long j);        // returns ref to (i,j) entry
  T operator()(long i, long j) const;   // returns (i,j) entry
  matT slice(long,long,long=-1,long=-1) const;// returns submatrix
  T sub(long i, long j) const;             // returns the (i,j) entry
  vecT<T> row(long i) const;                // returns row i (as a vec)
  vecT<T> col(long j) const;                // returns col j (as a vec)
  void set(long i, long j, const T& x);     // sets the (i,j) entry to x
  void add(long i, long j, const T& x);  // adds x to the (i,j) entry
  void setrow(long i, const vecT<T>& v);
  void setcol(long i, const vecT<T>& v);
  void swaprows(long r1, long r2);
  void multrow(long r, const T& scal);
  void divrow(long r, const T& scal);
  T content() const;
  T row_content(long r) const;
  void clearrow(long r);
  void make_primitive();
  void operator+=(const matT<T>&);
  void operator-=(const matT<T>&);
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

  static matT<T> scalar_matrix(long n, const T& a);
  static matT<T> identity_matrix(long n) {return scalar_matrix(n, T(1));}

  // non-member (friend) functions and operators
  friend void add_row_to_vec<>(vecT<T>& v, const matT<T>& m, long i);
  friend void sub_row_to_vec<>(vecT<T>& v, const matT<T>& m, long i);
  friend matT operator*<>(const matT<T>&, const matT<T>&);
  friend vecT<T> operator*<>(const matT<T>&, const vecT<T>&);
  friend int operator==<>(const matT<T>&, const matT<T>&);
  friend istream& operator>><> (istream&s, matT<T>&);
  friend matT colcat<>(const matT<T>& a, const matT<T>& b);
  friend matT rowcat<>(const matT<T>& a, const matT<T>& b);
  friend matT directsum<>(const matT<T>& a, const matT<T>& b);
  friend void elimrows<>(matT<T>& m, long r1, long r2, long pos); //plain elimination, no clearing
  friend void elimrows1<>(matT<T>& m, long r1, long r2, long pos); //elimination + clearing
  friend void elimrows2<>(matT<T>& m, long r1, long r2, long pos, const T& last); //elimination + divide by last pivot
  friend matT echelon0<>(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                          long& rk, long& ny, T& d);
  friend void elimp<>(matT<T>& m, long r1, long r2, long pos, const T& pr);
  friend void elimp1<>(matT<T>& m, long r1, long r2, long pos, const T& pr);
  friend matT echelonp<>(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                          long& rk, long& ny, T& d, const T& pr);
  friend matT echmodp<>(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                         long& rk, long& ny, const T& pr);
  friend matT echmodp_uptri<>(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                               long& rk, long& ny, const T& pr);
  friend matT ref_via_flint<>(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                               long& rk, long& ny, const T& pr);
  friend matT ref_via_ntl<>(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                             long& rk, long& ny, const T& pr);
  friend matT rref<>(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                     long& rk, long& ny, const T& pr);
  friend long rank_via_ntl<>(const matT<T>& M, const T& pr);
  friend T det_via_ntl<>(const matT<T>& M, const T& pr);
  friend subspaceT<T> combine<>(const subspaceT<T>& s1, const subspaceT<T>& s2);
  friend matT restrict_mat<>(const matT<T>& m, const subspaceT<T>& s, int cr);
  friend int liftmat<>(const matT<T>& mm, const T& pr, matT<T>& m, T& dd);
  friend int lift<>(const subspaceT<T>& s, const T& pr, subspaceT<T>& ans);
  friend subspaceT<T> pcombine<>(const subspaceT<T>& s1, const subspaceT<T>& s2, const T& pr);
  friend matT prestrict<>(const matT<T>& m, const subspaceT<T>& s, const T& pr, int cr);
  friend matT matmulmodp<>(const matT<T>&, const matT<T>&, const T& pr);
  friend matT echmodp_d<>(const matT<T>& mat, vecT<int>& pcols, vecT<int>& npcols, long& rk, long& ny, double pr);
  friend long population<>(const matT<T>& m); // #nonzero entries
  friend T maxabs<>(const matT<T>& m); // max entry
  friend double sparsity<>(const matT<T>& m); // #nonzero entries/#entries
  // Implementation
private:
  long nro,nco;
  vector<T> entries;  // stored in one array, by rows
};

// Declaration of non-friend functions

template<class T>
inline ostream& operator<< (ostream&s, const matT<T>&m)
{m.output(s); return s;}

template<class T>
matT<T> operator+(const matT<T>&);                   // unary
template<class T>
matT<T> operator-(const matT<T>&);                   // unary
template<class T>
matT<T> operator+(const matT<T>& m1, const matT<T>& m2);
template<class T>
matT<T> operator-(const matT<T>& m1, const matT<T>& m2);
template<class T>
matT<T> operator*(const T& scal, const matT<T>& m);
template<class T>
matT<T> operator/(const matT<T>& m, const T& scal);
template<class T>
int operator!=(const matT<T>& m1, const matT<T>& m2);
template<class T>
matT<T> transpose(const matT<T>& m);
// submatrix of only rows indexed by v, all columns
template<class T>
matT<T> rowsubmat(const matT<T>& m, const vecT<int>& v);
template<class T>
matT<T> rowsubmat(const matT<T>& m, const vecT<long>& v);
// submatrix of rows indexed by iv, columns indexed by jv
template<class T>
matT<T> submat(const matT<T>& m, const vecT<int>& iv, const vecT<int>& jv);
template<class T>
matT<T> submat(const matT<T>& m, const vecT<long>& iv, const vecT<long>& jv);

template<class T>
matT<T> echelon(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                          long& rk, long& ny, T& d, int method=0);
template<class T>
matT<T> addscalar(const matT<T>&, const T&);
template<class T>
vecT<T> apply(const matT<T>&, const vecT<T>&);

template<class T> matT<T> rref(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
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
mat_zz_p mat_zz_p_from_mat(const matT<T>& M, const T& pr);

// Construct a mat (T type same as pr) from an NTL mat_lzz_p

template<class T>
matT<T> mat_from_mat_zz_p(const mat_zz_p& A, const T& pr); // type of T fixes return type

#endif

