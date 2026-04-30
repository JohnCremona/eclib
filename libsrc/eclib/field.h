// File FIELD.H: class for working with number fields for Hecke eigenvalues
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

#ifndef _FIELD_H
#define _FIELD_H      1

#include "linalg.h"
#include "bigrat.h"

class Qvec; // vec_m (=Zvec<ZZ>) with common denominator
class Qmat; // mat_m (=Zmat<ZZ>) with common denominator
class Field;
class FieldIso;
class FieldElement;

class Qvec {
  friend class FieldElement;
  friend class Qmat;
  vec_m numerator;
  ZZ denom;
  void cancel(); // divides through by gcd(content(numerator, denom))

public:
  Qvec() :numerator(vec_m()), denom(to_ZZ(1)) {;}
  explicit Qvec(int d) :numerator(vec_m(d)), denom(to_ZZ(1)) {;}
  explicit Qvec(const vec_m& c, const ZZ& d=to_ZZ(1)) :numerator(c), denom(d) { cancel();}
  static Qvec unit_vector(long d, long i) {return Qvec(vec_m::unit_vector(d,i));}

  // Before calling this the sie (dimension) must be set
  void read (istream& s) { s >> numerator >> denom;}

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
  inline friend Qvec operator*(const ZZ& c, const Qvec& v) {return Qvec(c*v.numerator, v.denom);}
  inline friend ostream& operator<<(ostream& s, const Qvec& x) { s << x.str();  return s;}
};

inline istream& operator>>(istream& s, Qvec& x) {x.read(s); return s;}
inline Qvec operator*(int c, const Qvec& v) {return ZZ(c)*v;}
inline Qvec operator*(long c, const Qvec& v) {return ZZ(c)*v;}

class Qmat {
  friend class FieldElement;
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

  // Before calling this the sie (dimension) must be set
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
  Qvec col(int i) const {return Qvec(numerator.col(i), denom);}

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
  inline friend ostream& operator<<(ostream& s, const Qmat& x) { s << x.str();  return s;}
};

inline istream& operator>>(istream& s, Qmat& x) {x.read(s); return s;}
inline Qmat operator*(int c, const Qmat& v) {return ZZ(c)*v;}
inline Qmat operator*(long c, const Qmat& v) {return ZZ(c)*v;}

class Field {
  friend class FieldIso;
  friend class FieldElement;
  friend class Newform;
private:
  string var;   // name of generator
  int d;        // degree
  ZZX minpoly;  // irreducible poly of degree d
  ZZ denom; // minpoly is the (integral) min poly of A/denom
  vector<mat_m> Cpowers;  // C^i for i=0,1,...,d-1, where C =dxd companion matrix with min.poly. minpoly
public:
  ~Field() {minpoly.kill();}
  Field(); // defaults to Q
  explicit Field(const Qmat& A, Qmat& B, Qmat& Binv, string a="a", int verb=0);
  explicit Field(const Qmat& A, string a="a", int verb=0)
  {
    Qmat B, Binv;
    *this = Field(A, B, Binv, a, verb);
  }
  explicit Field(const ZZX& p, string a="a", int verb=0);

  FieldElement operator()(const bigrational& x) const;
  FieldElement operator()(const ZZ& x) const;
  FieldElement operator()(long x) const;
  FieldElement operator()(int x) const;
  FieldElement gen() const;
  FieldElement element(const vec_m& c, const ZZ& d=to_ZZ(1)) const;
  FieldElement element(const Qvec& v) const;
  int degree() const {return d;}
  int isQ() const {return d==1;}
  ZZX poly() const {return minpoly;}
  int operator==(const Field& F) const {return d==F.d && (d==1 || minpoly==F.minpoly);}
  int operator!=(const Field& F) const {return d!=F.d || (d!=1 && minpoly!=F.minpoly);}
  void display(ostream&s = cout) const;
  string get_var() const {return var;}
  void set_var(const string& v)  {var = v;}
  // String for pretty output, like "Q" or "Q(i) = Q[X]/(X^2+1)", or
  // (if raw) raw output, suitable for re-input, like "Q" or "i [1 0 1]":
  string str(int raw=0) const;
  void read(istream& s);
  friend istream& operator>>(istream& s, Field& F);

  // Apply polredabs (if canonical) or polredbest to the defining
  // polynomial, define a new field with that poly and return an
  // isomorphism from this to that.  If the poly was already
  // polredabsed, or if F is QQ, return the identity. Otherwise a new
  // Field is created with provided variable name.
  FieldIso reduction_isomorphism(string newvar, Field& Fred, int canonical=0) const;

  // Return an iso from this=Q(a) to Q(b) where b is in this field and generates
  FieldIso change_generator(const FieldElement& b, Field& Qb) const;

  // Return an embedding from K to K(sqrt(r)), optionally applying
  // polredabs (if reduce=2) or polredbest (if reduce=1) to the
  // codomain.  sqrt_r is set to sqrt(r) in the codomain, so sqrt_r^2
  // = image of r.  Caller supplies a reference to the target field
  // which will be overwritten (so caller already has the pointer to
  // it).
  FieldIso sqrt_embedding(const FieldElement& r, string newvar,
                          Field& F_sqrt_r, FieldElement& sqrt_r,
                          int reduce=1) const;

  // Return an embedding from K to K(sqrt(r1),sqrt(r2),...),
  // optionally applying polredabs (if reduce=2) or polredbest (if
  // reduce=1) to the codomain.  sqrts is set to [sqrt(r1), sqrt(r2),
  // ..], a list of elements of the codoimain, so (sqrt_ri)^2 = image
  // of ri.  Caller supplies a reference to the target field which
  // will be overwritten (so caller already has the pointer to it).
  FieldIso sqrt_embedding(const vector<FieldElement>& r_list, string newvar,
                          Field& F_sqrt_r, vector<FieldElement>& sqrt_r_list,
                          int reduce=1) const;
};

inline ostream& operator<<(ostream& s, const Field& F)
{ s << F.str();  return s; }

class FieldElement {
  friend class Field;
  friend class FieldIso;
  friend FieldElement evaluate(const ZZX& f, const FieldElement a);
private:
  Field const * F;
  Qvec v;             // only used for degree>1, i.e. not Q
  bigrational val;    // only used for degree=1, i.e. Q
  // In general the field element is the v-combination of power basis
  // of F, but when F is Q this class is just a wrapper round the
  // bigrational class and v is ignored.
public:
  FieldElement() {;}
  explicit FieldElement(const Field& HF)
    :F(&HF), v(vec_m(F->d)) {if (F->d==1) val = bigrational(0);}
  FieldElement(const Field& HF, const vec_m& c, const ZZ& d=to_ZZ(1))
    :F(&HF), v(c,d) {;}
  FieldElement(const Field& HF, const Qvec& c)
    :F(&HF), v(c) {;}

  // creation from a rational (general F)
  FieldElement(const Field& HF, const ZZ& a, const ZZ& d=to_ZZ(1))
    :F(&HF), v(a*vec_m::unit_vector(F->d, 1), d), val(bigrational(a,d)) {;}
  // creation from a rational
  FieldElement(const Field& HF, const bigrational& r)
    :F(&HF), v(r.num()*vec_m::unit_vector(F->d, 1), r.den()), val(r) {;}

  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;

  const Field* field_ptr() const {return F;}
  int field_degree() const {return F->d;}
  int field_is_Q() const {return F->d == 1;}

  mat_m matrix() const; // ignores denom
  // NB Since we do not have polynomials with rational coefficients,
  // both charpoly and minpoly are scaled to be primitive rather than
  // monic.
  ZZX charpoly() const;
  ZZX minpoly() const;
  int degree() const {return deg(minpoly());}
  bigrational norm() const;
  bigrational trace() const;
  int is_zero() const;
  int is_one() const;
  int is_minus_one() const;
  int is_generator() const {return degree()==F->d;}
  bigrational get_val() const {return val;}
  Qvec coords() const {return v;}
  Qmat power_matrix() const; // cols are coords of powers
  ZZ get_denom() const {return (field_is_Q()? val.den() : v.denom);}
  int in_same_field(const FieldElement& b) const {return (F==b.F) || (*F==*b.F);}
  int operator==(const FieldElement& b) const;
  int operator!=(const FieldElement& b) const;

  void set_zero();
  void set_one();

  FieldElement operator+(const FieldElement& b) const; // add
  FieldElement operator+(const ZZ& b) const {return operator+(F->operator()(b));} // add
  FieldElement operator+(const int& b) const {return operator+(F->operator()(b));} // add
  FieldElement operator+(const long& b) const {return operator+(F->operator()(b));} // add
  void operator+=(const FieldElement& b); // add b to this
  void operator+=(const ZZ& b) { operator+=(F->operator()(b));} // add b
  void operator+=(const int& b) { operator+=(F->operator()(b));} // add b
  void operator+=(const long& b) { operator+=(F->operator()(b));} // add b

  FieldElement operator-(const FieldElement& b) const; // subtract
  FieldElement operator-(const ZZ& b) const {return operator-(F->operator()(b));} // subtract
  FieldElement operator-(const int& b) const {return operator-(F->operator()(b));} // subtract
  FieldElement operator-(const long& b) const {return operator-(F->operator()(b));} // subtract
  void operator-=(const FieldElement& b); // subtract b from this
  void operator-=(const ZZ& b) { operator-=(F->operator()(b));} // subtract b
  void operator-=(const int& b) { operator-=(F->operator()(b));} // subtract b
  void operator-=(const long& b) { operator-=(F->operator()(b));} // subtract b
  FieldElement operator-() const;                           // unary minus

  FieldElement operator*(const FieldElement& b) const; // product
  FieldElement operator*(const ZZ& b) const {return operator*(F->operator()(b));} // product
  FieldElement operator*(const int& b) const {return operator*(F->operator()(b));} // product
  FieldElement operator*(const long& b) const {return operator*(F->operator()(b));} // product
  void operator*=(const FieldElement& b); // multiply by b
  void operator*=(const ZZ& b) { operator*=(F->operator()(b));} // multiply by b
  void operator*=(const int& b) { operator*=(F->operator()(b));} // multiply by b
  void operator*=(const long& b) { operator*=(F->operator()(b));} // multiply by b

  FieldElement inverse() const; // raise error if zero      // inverse
  FieldElement operator/(const FieldElement& b) const; // divide (raise error if b is zero)
  FieldElement operator/(const ZZ& b) const {return operator/(F->operator()(b));} // divide
  FieldElement operator/(const int& b) const {return operator/(F->operator()(b));} // divide
  FieldElement operator/(const long& b) const {return operator/(F->operator()(b));} // divide
  void operator/=(const FieldElement& b);                        // divide by b
  void operator/=(const ZZ& b) { operator/=(F->operator()(b));} // divide by b
  void operator/=(const int& b) { operator/=(F->operator()(b));} // divide by b
  void operator/=(const long& b) { operator/=(F->operator()(b));} // divide by b
  void negate(); // negate in place

  // NB for a in F, either [Q(sqrt(a))=Q(a)] or [Q(sqrt(a)):Q(a)]=2.
  // The first function only applies when a has maximal degree:
  // return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
  int is_absolute_square(FieldElement& r)  const;
  // Same as above if the denom is 1
  int is_absolute_integral_square(FieldElement& r)  const;
  // The second function applies in general: return 1 and r
  // s.t. r^2=this, with deg(r)=degree(), else 0. Here, ntries is the
  // number of squares this is multiplied by to get odd co-degree.
  int is_square(FieldElement& r, int ntries=100) const;

  // return 1 and set r to the rational value if the degree is 1
  int is_rational(bigrational& r) const;
  // Same with no value needed
  int is_rational() const {bigrational r;  return is_rational(r);}

  // return 1 iff this is an algebraic integer
  int is_integral() const;

  // x must be initialised with a Field before input to x
  void read (istream& s);
  friend istream& operator>>(istream& s, FieldElement& x);
};

inline ostream& operator<<(ostream& s, const FieldElement& x)
{ s << x.str();  return s;}

FieldElement evaluate(const ZZX& f, const FieldElement a);

class FieldIso {
  friend class Field;
  friend class FieldElement;
private:
  Field const * domain;   // pointer to const object
  Field const * codomain; // pointer to const object
  Qmat isomat;
  int id_flag; // flags that this is the identity (domain=codomain and isomat=id)
  void set_id_flag()
  {
    id_flag = (domain==codomain) && isomat.is_identity();
  }
public:
  // Dummy constructor
  FieldIso() {;}

  // Constructor from a known matrix
  FieldIso( const Field& F1, const Field& F2, const mat_m& M, const ZZ& d = ZZ(1), int id=-1)
    :domain(&F1), codomain(&F2), isomat(Qmat(M,d)), id_flag(id)
  {
    if (id_flag==-1)
      set_id_flag();
  }
  // Constructor from a known matrix
  FieldIso( const Field& F1, const Field& F2, const Qmat& M, int id=-1)
    :domain(&F1), codomain(&F2), isomat(M), id_flag(id)
  {
    if (id_flag==-1)
      set_id_flag();
  }
  // Partial constructor from two fields, used when inputting just the matrix and denominator.
  // This initializes isomat to the correct size as the 0 matrix
  FieldIso( const Field& F1, const Field& F2)
    :domain(&F1), codomain(&F2), isomat(Qmat(F2.d,F1.d)), id_flag(0) {;}
  // Identity
  explicit FieldIso( const Field& F1)
    :domain(&F1), codomain(&F1), isomat(mat_m::identity_matrix(F1.d)), id_flag(1) {;}

  // inverse isomorphism
  FieldIso inverse() const;
  // map x in domain to an element of the codomain
  FieldElement operator()(const FieldElement& x) const;
  // same to all in a list of elements of the domain
  vector<FieldElement> operator()(const vector<FieldElement>& x) const;

  // precompose this FieldIso with another (requires iso.codomain = domain)
  void precompose(const FieldIso& iso);
  // postcompose this FieldIso with another (requires iso.domain = codomain)
  void postcompose(const FieldIso& iso);
  // return postcomposion of this and iso (requires iso.domain = codomain)
  FieldIso operator*(const FieldIso& iso);
  // access data
  int is_identity() const {return id_flag;}
  int is_nontrivial() const {return !id_flag;}
  Field const * dom() const {return domain;}
  Field const * codom() const {return codomain;}
  Qmat matrix() const {return isomat;}

  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;

  // x must be initialised with domain and codomain, this just inputs
  // the matrix and denominator
  void read(istream& s);
  friend istream& operator>>(istream& s, FieldIso& x);
};

inline ostream& operator<<(ostream& s, const FieldIso& x)
{ s << x.str(); return s; }


#endif
