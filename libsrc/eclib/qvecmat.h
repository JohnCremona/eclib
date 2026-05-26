// File QVECMAT.H: classes for working with rational vectors and matrices
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

// NB These classes do not provide all the functionality one might
// think of, just what is needed for the classes Field, FieldElement,
// FieldIso, Newspace and Newform.  Also, the base integer type is
// currently fixed to be ZZ (NTL integers) but could easily be
// templated to work also with int, long and INT (wrapping FLINT
// integers).  This would require a templatised Rational class to
// replace the classes rational (for long), bigrational (for ZZ) and
// RAT (wrapping FLINT fmpq).

#ifndef _QVECMAT_H
#define _QVECMAT_H      1

#include "linalg.h"
#include "bigrat.h"

class Qvec; // vec_m (=Zvec<ZZ>) with common denominator
class Qmat; // mat_m (=Zmat<ZZ>) with common denominator
class Field;
class FieldIso;
class FieldElement;

class Qvec {
  friend class Field;
  friend class FieldElement;
  friend class Order;
  friend class Qmat;
  vec_m numerator;
  ZZ denom;
  void cancel(); // divides through by gcd(content(numerator, denom))

public:
  Qvec() :numerator(vec_m()), denom(to_ZZ(1)) {;}
  explicit Qvec(int d) :numerator(vec_m(d)), denom(to_ZZ(1)) {;}
  explicit Qvec(const vec_m& c, const ZZ& d=to_ZZ(1)) :numerator(c), denom(d) { cancel();}
  explicit Qvec(const bigrational& r) :numerator(vec_m({r.num()})), denom(r.den()) {;}
  static Qvec unit_vector(long d, long i) {return Qvec(vec_m::unit_vector(d,i));}

  // Before calling this the size (dimension) must be set
  void read(istream& s) { s >> numerator >> denom;}

  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;
  void set_zero() {numerator.clear(); denom=ZZ(1);}
  void set_unit_vector(int i) {numerator.clear(); numerator[i]=denom=ZZ(1);}
  int is_zero() const {return trivial(numerator);}
  int is_integral() const {return denom==1;}
  vec_m get_numerator() const {return numerator;}
  ZZ get_denom() const {return denom;}
  int dim() const {return numerator.dim();}
  bigrational operator[](int i) const {return bigrational(numerator[i], denom);}
  int operator==(const Qvec& b) const {return (denom==b.denom) && (numerator==b.numerator);}
  int operator!=(const Qvec& b) const {return (denom!=b.denom) || (numerator!=b.numerator);}
  Qvec operator+(const Qvec& b) const; // add
  void operator+=(const Qvec& b); // add b to this
  Qvec operator-(const Qvec& b) const; // subtract
  void operator-=(const Qvec& b); // subtract b from this
  Qvec operator-() const {return Qvec(-numerator, denom);} // unary minus
  void operator *= (const ZZ& c) {numerator *= c; cancel();}
  void operator *= (int c) {numerator *= ZZ(c); cancel();}
  void operator *= (long c) {numerator *= ZZ(c); cancel();}
  inline friend Qvec reverse(const Qvec& v) {return Qvec(reverse(v.numerator), v.denom);}
  inline friend Qvec operator*(const ZZ& c, const Qvec& v) {return Qvec(c*v.numerator, v.denom);}
  inline friend ostream& operator<<(ostream& s, const Qvec& x) { s << x.str();  return s;}
  friend Qvec operator*(const Qmat&m, const Qvec& v);
  friend Qvec operator*(const mat_m&m, const Qvec& v);
};

inline istream& operator>>(istream& s, Qvec& x) {x.read(s); return s;}
inline Qvec operator*(int c, const Qvec& v) {return ZZ(c)*v;}
inline Qvec operator*(long c, const Qvec& v) {return ZZ(c)*v;}

class Qmat {
  friend class FieldElement;
  friend class Order;
  friend class Qvec;
  mat_m numerator;
  ZZ denom;
  void cancel(); // divides through by gcd(content(numerator, denom))

public:
  Qmat() :numerator(mat_m()), denom(to_ZZ(1)) {;}
  Qmat(long nr, long nc) :numerator(mat_m(nr,nc)), denom(to_ZZ(1)) {;}
  explicit Qmat(const mat_m& c, const ZZ& d=to_ZZ(1)) :numerator(c), denom(d) { cancel();}
  static Qmat identity(long d) {return Qmat(mat_m::identity_matrix(d));}
  static Qmat zero(long d) {return Qmat(mat_m(d,d));}

  // Before calling this the size (dimension) must be set
  void read (istream& s) { s >> numerator >> denom;}

  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;
  int is_zero() const {return numerator.is_zero();}
  int is_integral() const {return denom==1;}
  int is_identity() const {return denom==1 && numerator==mat_m::identity_matrix(numerator.nrows());}
  mat_m get_numerator() const {return numerator;}
  ZZ get_denom() const {return denom;}
  int nrows() const {return numerator.nrows();}
  int ncols() const {return numerator.ncols();}

  void setcol(int i, const Qvec& v);
  void setrow(int i, const Qvec& v);
  Qvec col(int i) const {return Qvec(numerator.col(i), denom);}
  Qvec row(int i) const {return Qvec(numerator.row(i), denom);}

  // append n 0 rows at bottom:
  void append_rows(int n = 1)
  {
    numerator.append_rows(n);
  }

  // append one 0 row at bottom:
  void append_row()
  {
    numerator.append_rows(1);
  }

  // append v as a new row at bottom:
  void append_row(const Qvec& v)
  {
    numerator.append_rows(1);
    setrow(numerator.nrows(), v);
  }

  // append v for v in vlist as new rows at bottom.
  void append_rows(const vector<Qvec>& vlist)
  {
    int old_nro = numerator.nrows();
    int extra_nro = vlist.size();
    append_rows(extra_nro);
    for (int i=0; i<extra_nro; i++)
      setrow(old_nro+i+1, vlist[i]);
  }

  // delete last row, recancel if requested (e.g. will not be
  // necessary if the last row is 0):
  void delete_row(int recancel=0)
  {
    numerator.delete_rows(1);
    if (recancel) cancel();
  }

  // delete last n rows, recancel if requested (e.g. will not be
  // necessary if the last rows are all 0):
  void delete_rows(int n, int recancel=0)
  {
    numerator.delete_rows(n);
    if (recancel) cancel();
  }

  // trace, det, charpoly, inverse only for square matrices (not checked)
  bigrational trace() const {return bigrational(numerator.trace(), denom);}
  bigrational det() const {return bigrational(numerator.determinant(), pow(denom, numerator.nrows()));}
  ZZX charpoly() const {return scaled_charpoly(mat_to_mat_ZZ(numerator), denom);}
  Qmat inverse() const;
  mat_m companion_transform(Qmat& B, Qmat& Binv) const;

  int operator==(const Qmat& b) const {return (denom==b.denom) && (numerator==b.numerator);}
  int operator!=(const Qmat& b) const {return (denom!=b.denom) || (numerator!=b.numerator);}
  Qmat operator+(const Qmat& b) const; // add
  void operator+=(const Qmat& b); // add b to this
  Qmat operator-(const Qmat& b) const; // subtract
  void operator-=(const Qmat& b); // subtract b from this
  Qmat operator-() const {return Qmat(-numerator, denom);} // unary minus
  void operator *= (const ZZ& c) {numerator *= c; cancel();}
  void operator *= (int c) {numerator *= ZZ(c); cancel();}
  void operator *= (long c) {numerator *= ZZ(c); cancel();}
  void operator *= (const Qmat& m) {numerator = numerator*m.numerator; denom *= m.denom; cancel();}
  inline friend Qmat operator*(const ZZ& c, const Qmat& x) {return Qmat(c*x.numerator, x.denom);}
  inline friend Qmat operator*(const Qmat&m1, const Qmat& m2) {return Qmat(m1.numerator*m2.numerator, m1.denom*m2.denom);}
  inline friend Qvec operator*(const Qmat&m, const Qvec& v) {return Qvec(m.numerator*v.numerator, m.denom*v.denom);}
  inline friend Qvec operator*(const mat_m&m, const Qvec& v) {return Qvec(m*v.numerator, v.denom);}
  inline friend ostream& operator<<(ostream& s, const Qmat& x) { s << x.str();  return s;}

  inline friend Qmat HNF(const Qmat& M) {return Qmat(HNF(M.numerator), M.denom);}
  inline friend Qmat HNF(const Qmat& M, mat_m& U) {return Qmat(HNF(M.numerator, U), M.denom);}

  inline friend Qmat SNF(const Qmat& M) {return Qmat(SNF(M.numerator), M.denom);}
  inline friend Qmat SNF(const Qmat& M, mat_m& U, mat_m& V) {return Qmat(SNF(M.numerator, U, V), M.denom);}

  inline friend Qmat LLL(const Qmat& M) {return Qmat(LLL(M.numerator), M.denom);}
  inline friend Qmat LLL(const Qmat& M, mat_m& U) {return Qmat(LLL(M.numerator, U), M.denom);}

  inline friend Qmat transpose(const Qmat& M) {return Qmat(transpose(M.numerator), M.denom);}
  inline friend int is_permutation_matrix(const Qmat& M)
  {return is_one(M.denom) && is_permutation_matrix(M.numerator);}
  inline friend int is_signed_permutation_matrix(const Qmat& M)
  {return is_one(M.denom) && is_signed_permutation_matrix(M.numerator);}
};

inline istream& operator>>(istream& s, Qmat& x) {x.read(s); return s;}
inline Qmat operator*(int c, const Qmat& v) {return ZZ(c)*v;}
inline Qmat operator*(long c, const Qmat& v) {return ZZ(c)*v;}
inline Qmat operator*(const mat_m& A, const Qmat& B) {return Qmat(A)*B;}
inline Qmat operator*(const Qmat& A, const mat_m& B) {return A*Qmat(B);}

#endif
