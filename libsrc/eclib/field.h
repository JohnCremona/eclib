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

class Field;
class FieldElement;
class FieldIso;
class Order;

////////////////////////////////////////////////////////////////////////
//
// FieldElement class
//
////////////////////////////////////////////////////////////////////////

class FieldElement {
  friend class Field;
  friend class Order;
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
  // creation of 0 in the given field
  explicit FieldElement(const Field& HF);
  // creation from an integer vector of coords with denominator
  FieldElement(const Field& HF, const vec_m& c, const ZZ& d=to_ZZ(1));
  // creation from a rational vector of coords
  FieldElement(const Field& HF, const Qvec& c);
  // creation from a rational (general F)
  FieldElement(const Field& HF, const ZZ& a, const ZZ& d=to_ZZ(1));
  // creation from a rational
  FieldElement(const Field& HF, const bigrational& r);

  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;

  const Field* field_ptr() const {return F;}
  int field_degree() const;
  int field_is_Q() const;

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
  int is_generator() const;
  bigrational get_val() const {return val;}
  Qvec coords() const {return v;}
  Qmat power_matrix() const; // cols are coords of powers
  ZZ get_denom() const {return (field_is_Q()? val.den() : v.denom);}
  int in_same_field(const FieldElement& b) const;
  int operator==(const FieldElement& b) const;
  int operator!=(const FieldElement& b) const;

  void set_zero();
  void set_one();

  FieldElement operator+(const FieldElement& b) const; // add
  FieldElement operator+(const ZZ& b) const;           // add
  FieldElement operator+(const int& b) const;          // add
  FieldElement operator+(const long& b) const;         // add
  void operator+=(const FieldElement& b); // add b to this
  void operator+=(const ZZ& b);           // add b to this
  void operator+=(const int& b);          // add b to this
  void operator+=(const long& b);         // add b to this

  FieldElement operator-(const FieldElement& b) const; // subtract
  FieldElement operator-(const ZZ& b) const;
  FieldElement operator-(const int& b) const;
  FieldElement operator-(const long& b) const;
  void operator-=(const FieldElement& b); // subtract b from this
  void operator-=(const ZZ& b);
  void operator-=(const int& b);
  void operator-=(const long& b);
  FieldElement operator-() const;                           // unary minus

  FieldElement operator*(const FieldElement& b) const; // product
  FieldElement operator*(const ZZ& b) const;
  FieldElement operator*(const int& b) const;
  FieldElement operator*(const long& b) const;
  void operator*=(const FieldElement& b); // multiply by b
  void operator*=(const ZZ& b);
  void operator*=(const int& b);
  void operator*=(const long& b);

  FieldElement inverse() const; // raise error if zero      // inverse
  FieldElement operator/(const FieldElement& b) const; // divide (raise error if b is zero)
  FieldElement operator/(const ZZ& b) const;
  FieldElement operator/(const int& b) const;
  FieldElement operator/(const long& b) const;
  void operator/=(const FieldElement& b);                        // divide by b
  void operator/=(const ZZ& b);
  void operator/=(const int& b);
  void operator/=(const long& b);
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

// List of x^i for i=0,1,...,deg(F)
vector<FieldElement> powers(const FieldElement& x);

// To coordinate matrix (by columns unless rows=1). If rev then
// rows/columns are reversed coordinate vectors.
Qmat coord_matrix(const vector<FieldElement>& alist, int rows=0, int rev=0);

// From coordinate matrix (by columns unless rows=1).  If rev then
// rows/columns are reversed coordinate vectors.
vector<FieldElement> from_coord_matrix(const Field& F, const Qmat& M, int rows=0, int rev=0);

inline ostream& operator<<(ostream& s, const FieldElement& x)
{ s << x.str();  return s;}

FieldElement evaluate(const ZZX& f, const FieldElement a);

////////////////////////////////////////////////////////////////////////
//
// Order class
//
////////////////////////////////////////////////////////////////////////

class Order{

public:
  // Constructors:
  Order() {rank=0;}
  explicit Order(const Field& F); // equation order
  Order(const vector<FieldElement>& v, int basis=1); // given a Z-basis or just a Z-spanning set
  Order(const vector<FieldElement>& v, const mat_m pcm); // same with known power_coords_matrix
  Order(const Field& F, const ZZ& i, const mat_m& M); // Order in F given pcm

  // Access data:
  int rk() const {return rank;}
  const vector<FieldElement> get_basis() const {return Zbasis;}
  mat_m get_pcm() const {return power_coords_matrix;}
  Qmat get_bm() const {return basis_matrix;}
  // index of equation order in this
  ZZ get_order_index() const {return index;}
  // discriminant of this order
  ZZ get_disc() const {return disc;}
  // discriminant of equation order
  ZZ get_poldisc() const {return poldisc;}
  string str(int raw=0) const;

  // coords w.r.t. Zbasis of an arbitrary element of F
  Qvec coords(const FieldElement& a) const  { return power_coords_matrix * a.v;}
  // denominator of Zbasis coords of an arbitrary element of F
  ZZ denom(const FieldElement& a) const { return coords(a).denom;}
  // coords w.r.t. Zbasis of an element of F in this order
  vec_m integral_coords(const FieldElement& a) const;
  // membership test
  int contains(const FieldElement& a) const;
  // membership test returning coords
  int contains(const FieldElement& a, vec_m& c) const;
  // containment test
  int contains(const Order& O2) const;
  // FieldElement from integer coords
  FieldElement operator()(const vec_m& coords) const;
  // FieldElement from rational coords
  FieldElement operator()(const Qvec& coords) const;

  // Functions to enlarge the order. If check==1, check that the
  // elements provided are algebraic integers.

  // Sum of two orders (in the same field!)
  Order operator+(const Order& O2);
  // Sum of this and power order of a
  Order operator+(const FieldElement& a);
  // Add another order to this:
  Order operator+=(const Order& O2);
  // Add the power order of a to this
  Order operator+=(const FieldElement& a);

  // The old versions add all products of the gens of the two orders
  // (slower):

  // Extend by a (which must be an algebraic integer), returning the
  // index of the extension
  ZZ extend_by_old(const FieldElement& a, int check=1);
  // Extend by all a in alist (which must be algebraic integers),
  // returning the index of the extension
  ZZ extend_by_old(const vector<FieldElement>& alist, int check=1);

  // These versions add just a and then check all products one by one,
  // adding any which are not yet contained  (faster):

  // Extend by a (which must be an algebraic integer), returning the
  // index of the extension: second version which uses the next two
  // internally.
  ZZ extend_by(const FieldElement& a, int check=1);
  // Extend by all a in alist (which must be algebraic integers),
  // returning the index of the extension
  ZZ extend_by(const vector<FieldElement>& alist, int check=1);

  // LLL-reduce basis (using the basis matrix to reduce)
  void LLL_reduce();
  // LLL-reduce basis (using the coord matrix of alist to reduce)
  void LLL_reduce(const vector<FieldElement>& alist);

private:
  // Add one algebraic integer to the Z-span: the result may not be an
  // order but this is only used internally.
  void add_one(const FieldElement& a, int check=0);
  // Check that the Z-span of the current Zbasis is closed under
  // multiplication.  If not, a will hold a missing product.
  int check_order(FieldElement& a) const;

private:
  vector<FieldElement> Zbasis; // Z-basis of the order
  mat_m power_coords_matrix; // columns are coords of power basis w.r.t. Zbasis
  Qmat basis_matrix; // columns are coords of Zbasis w.r.t. power basis
  ZZ index; // index of equation order in this
  ZZ disc;  // discriminant of this order
  ZZ poldisc;  // discriminant of ambient field's defining polynomial
  int rank; // of this order, i.e. degree of ambient field
};

// Compute Maximal Order (via lib)pari.  If bound>0 then the order may
// not be p-maximal for p>bound.
Order MaximalOrder(const Field* F, const ZZ& bound = ZZ(0));

// For x an algebraic integer, the order spanned by x^i for
// i=0..deg(F)-1.  If check then check that x is an algebraic integer.
Order PowerOrder(const FieldElement& x, int check=0);

////////////////////////////////////////////////////////////////////////
//
// Field class
//
////////////////////////////////////////////////////////////////////////

class Field {
  friend class FieldIso;
  friend class FieldElement;
  friend class Order;
private:
  string var;   // name of generator
  int d;        // degree
  ZZX minpoly;  // irreducible poly of degree d
  vector<mat_m> Cpowers;  // C^i for i=0,1,...,d-1, where C =dxd companion matrix of minpoly
  vec_m Cpower_traces; // traces of C^i

  int have_integral_basis; // 0 if not yet computed
  Order Integers; // ring of integers (if have_integral_basis==1)

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
  vector<FieldElement> power_basis() const;
  FieldElement element(const vec_m& c, const ZZ& d=to_ZZ(1)) const;
  FieldElement element(const Qvec& v) const;
  FieldElement operator()(const Qvec& v) const;

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

  // compute Integers order (via libpari).  See MaximalOrder() for bound parameter
  void make_integers(const ZZ& bound=ZZ(0));
  // recreate integral basis from index and base_change_matrix's inverse (after reading from file)
  void set_integers(const ZZ& i, const mat_m& M);
  // check whether we have an integral basis
  int has_integral_basis() const {return have_integral_basis;}
  // return maximal order
  const Order& maximal_order() const {return Integers;}
  const Order& integers() const {return Integers;}

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
  friend Order MaximalOrder(const Field* F, const ZZ& bound);
};

inline ostream& operator<<(ostream& s, const Field& F)
{ s << F.str();  return s; }

////////////////////////////////////////////////////////////////////////
//
// FieldIso class
//
////////////////////////////////////////////////////////////////////////

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
  // Partial constructor from two fields, used when inputting just the matrix.
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

inline ostream& operator<<(ostream& s, const Order& x)
{ s << x.str(); return s; }

#endif
