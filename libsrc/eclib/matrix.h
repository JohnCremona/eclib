// matrix.h:  manage declarations of integer matrix classes
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

template<class T> Zmat<T> colcat(const Zmat<T>& a, const Zmat<T>& b);
template<class T> Zmat<T> rowcat(const Zmat<T>& a, const Zmat<T>& b);
template<class T> Zmat<T> directsum(const Zmat<T>& a, const Zmat<T>& b);

// Function to compute the REF of a matrix using FLINT:
template<class T> Zmat<T> ref(const Zmat<T>& M, T& d);
template<class T> Zmat<T> ref(const Zmat<T>& M, T& d,
                              Zvec<int>& pcols, Zvec<int>& npcols,
                              long& rk, long& ny);

// Function to compute the REF of a matrix modulo a prime, using FLINT:
template<class T> Zmat<T> ref_mod_p(const Zmat<T>& M, const T& pr);
template<class T> Zmat<T> ref_mod_p(const Zmat<T>& M, const T& pr,
                                    Zvec<int>& pcols, Zvec<int>& npcols,
                                    long& rk, long& ny);

// Divide out by common content, hence reducing A/d
template<class T> void cancel(Zmat<T>& A, T& d);
template<class T> T inverse(const Zmat<T>& A, Zmat<T>& Ainv);
template<class T> subZspace<T> combine(const subZspace<T>& s1, const subZspace<T>& s2);
template<class T> ssubZspace<T> combine(const ssubZspace<T>& s1, const ssubZspace<T>& s2);
template<class T> sZmat<T> restrict_mat(const sZmat<T>& m, const ssubZspace<T>& s);
template<class T> Zmat<T> restrict_mat(const Zmat<T>& m, const subZspace<T>& s, int cr=0);
template<class T> sZmat<T> restrict_mat(const sZmat<T>& m, const subZspace<T>& s);
template<class T> int liftmat(const Zmat<T>& mm, const T& pr, Zmat<T>& m, T& dd);
template<class T> int lift(const subZspace<T>& s, const T& pr, subZspace<T>& ans);
template<class T> subZspace<T> pcombine(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr);
template<class T> Zmat<T> matmulmodp(const Zmat<T>&, const Zmat<T>&, const T& pr);
template<class T> Zvec<T> matvecmulmodp(const Zmat<T>&, const Zvec<T>&, const T& pr);
template<class T> long population(const Zmat<T>& m);
template<class T> T maxabs(const Zmat<T>& m);
template<class T> double sparsity(const Zmat<T>& m);
template<class T> int is_permutation_matrix(const Zmat<T>& m);
template<class T> int is_signed_permutation_matrix(const Zmat<T>& m);
template<class T> int is_identity_matrix(const Zmat<T>& m);
template<class T> Zmat<T> addscalar(const Zmat<T>&, const T&);
template<class T> Zvec<T> apply(const Zmat<T>&, const Zvec<T>&);
template<class T> Zvec<T> operator*(const Zmat<T>& m, const Zvec<T>& v);
template<class T> Zmat<T> operator*(const Zmat<T>&, const Zmat<T>&);
template<class T> Zvec<T> operator*(const Zmat<T>&, const Zvec<T>&);
template<class T> int operator==(const Zmat<T>&, const Zmat<T>&);
template<class T> istream& operator>> (istream&s, Zmat<T>&);
// template<class T> Zmat<T> prestrict(const Zmat<T>& m, const subZspace<T>& s, const T& pr, int cr);

template <class T>
class Zmat {
friend class subZspace<T>;
friend class sZmat<T>;
friend class sZvec<T>;
friend class Zvec<T>;
friend class sZmat_elim<T>;
public:
  // constructors
  explicit Zmat(long nr=0, long nc=0)   :nro(nr), nco(nc) {entries.resize(nr*nc, T(0));}
  Zmat(long nr, long nc, const vector<T>& ent) : nro(nr), nco(nc), entries(ent) {}
  Zmat(const Zmat<T>& m)   :nro(m.nro), nco(m.nco), entries(m.entries) {}

  // member functions & operators
  void init(long nr=0, long nc=0);
  Zmat& operator=(const Zmat<T>&);                // assignment with copy
  T& operator()(long i, long j);        // returns reference to (i,j) entry
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
  ZZX charpoly() const; // returns an NTL poly
  T determinant() const;
  int is_zero() const;
  void output(ostream&s=cout) const;
  void output_raw(ostream&s=cout) const; // no "[", "]", commas or newlines
  void output_pari(ostream&s=cout) const;
  void output_pretty(ostream&s=cout) const;
  // same as output() except no newlines between rows
  void output_flat(ostream&s = cout) const;

  void reduce_mod_p(const T& p);

  static Zmat<T> scalar_matrix(long n, const T& a);
  static Zmat<T> identity_matrix(long n) {return scalar_matrix(n, T(1));}

  void append_rows(int n=1); // append n zero rows
  void append_row(const Zvec<T>& new_row); // append given new row
  void append_rows(const vector<Zvec<T>>& new_rows); // append given new rows
  void delete_rows(int n=1); // delete the last n rows

  // non-member (friend) functions and operators

  friend void cancel<>(Zmat<T>& A, T& d);
  friend T inverse<>(const Zmat<T>& A, Zmat<T>& Ainv);
  friend void add_row_to_vec<>(Zvec<T>& v, const Zmat<T>& m, long i);
  friend void sub_row_to_vec<>(Zvec<T>& v, const Zmat<T>& m, long i);
  friend Zmat operator*<>(const Zmat<T>&, const Zmat<T>&);
  friend Zvec<T> operator*<>(const Zmat<T>&, const Zvec<T>&);
  friend int operator==<>(const Zmat<T>&, const Zmat<T>&);
  friend istream& operator>><> (istream&s, Zmat<T>&);
  friend Zmat colcat<>(const Zmat<T>& a, const Zmat<T>& b);
  friend Zmat rowcat<>(const Zmat<T>& a, const Zmat<T>& b);
  friend Zmat directsum<>(const Zmat<T>& a, const Zmat<T>& b);
  friend subZspace<T> combine<>(const subZspace<T>& s1, const subZspace<T>& s2);
  friend Zmat<T> restrict_mat<>(const Zmat<T>& m, const subZspace<T>& s, int cr);
  friend int liftmat<>(const Zmat<T>& mm, const T& pr, Zmat<T>& m, T& dd);
  friend int lift<>(const subZspace<T>& s, const T& pr, subZspace<T>& ans);
  friend subZspace<T> pcombine<>(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr);
  // friend Zmat<T> prestrict<>(const Zmat<T>& m, const subZspace<T>& s, const T& pr, int cr);
  friend Zmat<T> matmulmodp<>(const Zmat<T>&, const Zmat<T>&, const T& pr);
  friend Zvec<T> matvecmulmodp<>(const Zmat<T>&, const Zvec<T>&, const T& pr);
  friend long population<>(const Zmat<T>& m); // #nonzero entries
  friend T maxabs<>(const Zmat<T>& m); // max entry
  friend double sparsity<>(const Zmat<T>& m); // #nonzero entries/#entries
  friend int is_permutation_matrix<>(const Zmat<T>& m); // test for being a permutation matrix
  friend int is_signed_permutation_matrix<>(const Zmat<T>& m);
  friend int is_identity_matrix<>(const Zmat<T>& m);

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

// Construct an NTL mat_lzz_p (matrix mod p) from a mat mod pr
template<class T>
mat_ZZ_p mat_ZZ_p_from_mat(const Zmat<T>& M, const T& pr);

// Construct a mat (T type same as pr) from an NTL mat_lzz_p
template<class T>
Zmat<T> mat_from_mat_ZZ_p(const mat_ZZ_p& A, const T& pr); // type of T fixes return type

// convert between matrices over different integer types
template<class T>
mat_m to_mat_m(const Zmat<T>& M);
template<class T>
mat_i to_mat_i(const Zmat<T>& M);
template<class T>
mat_l to_mat_l(const Zmat<T>& M);
template<class T>
mat_I to_mat_I(const Zmat<T>& M);

// convert an eclib Zmat to an NTL mat_ZZ:
template<class T>
mat_ZZ mat_to_mat_ZZ(Zmat<T> A);

// convert an eclib Zmat to an NTL mat_ZZ_p:
template<class T>
mat_ZZ_p mat_to_mat_ZZ_p(Zmat<T> A);

// compute char poly of A:
ZZX charpoly(const mat_ZZ& A);

// compute char poly of A:
template<class T>
ZZX charpoly(const Zmat<T>& A);

// compute char poly of A/den mod m:
ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den = to_ZZ(1), const ZZ& m = to_ZZ(0));

// return A mod m (or just A if m==0)
mat_ZZ reduce_mat(const mat_ZZ& A, const ZZ& m);

// evaluate f(A) (assumes f monic)
mat_ZZ evaluate(const ZZX& f, const mat_ZZ& A);
// evaluate f(A) (assumes f monic)
template<class T>
mat_m evaluate(const ZZX& f, const Zmat<T>& A);
// evaluate f(A) where the coeffs of f are in co
template<class T>
Zmat<T> evaluate(const Zvec<T>& co, const Zmat<T>& A);

// p should be monic:
mat_m CompanionMatrix(const ZZX& p);

// check that a matrix is a scaled involution (modulo m, if m!=0):
int check_involution(const mat_ZZ& A, const ZZ& den, const ZZ& m, int verbose=0);

// check that two matrices commute (modulo m, if m!=0):
int commute(const mat_ZZ& A, const mat_ZZ& B, const ZZ& m);

// check that a matrix commutes (modulo m, if m!=0) with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist, const ZZ& m);
// rank of an NTL matrix:
long rank(mat_ZZ A);

// nullity of an NTL matrix:
long nullity(mat_ZZ A);

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vec_m& co, const vector<mat_m>& mats);
// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vector<ZZ>& co, const vector<mat_m>& mats);

// create flint matrix (type fmpz_mat_t) copy of a Zmat<T>:
template<class T>
void flint_mat_from_mat(fmpz_mat_t& A, const Zmat<T>& M);
// convert a flint matrix (type fmpz_mat_t) to a Zmat<T>.  The dummy
// variable is needed to determine the return type
template<class T>
Zmat<T> mat_from_flint_mat(fmpz_mat_t& A, const T& dummy);

// Hermite Normal Form (via FLINT)
template<class T>
Zmat<T> HNF(const Zmat<T>& M);
template<class T>
Zmat<T> HNF(const Zmat<T>& M, Zmat<T>& U); // U*M=H with U unimodular

// Smith Normal Form (via FLINT)
template<class T>
Zmat<T> SNF(const Zmat<T>& M);
template<class T>
Zmat<T> SNF(const Zmat<T>& M, Zmat<T>& U, Zmat<T>& V); // U*M*V=S with U,V unimodular

// LLL row-reduction (via FLINT)
template<class T>
Zmat<T> LLL(const Zmat<T>& M);
template<class T>
Zmat<T> LLL(const Zmat<T>& M, Zmat<T>& U); // U*M=L with U unimodular

// From a matrix in REF, extract the pivotal and non-pivotal columns
template<class T>
void pnpcols(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols, long& rk, long& ny);

#endif
