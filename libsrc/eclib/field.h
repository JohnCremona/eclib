// File FIELD.H: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#ifndef _FIELD_H
#define _FIELD_H      1

#include "linalg.h"
#include "bigrat.h"

class Field;
class FieldIso;
class FieldElement;

extern const Field* FieldQQ;

// Divide through by gcd of content(M) and d
void cancel_mat(mat_m& M, ZZ& d);

class Field {
  friend class FieldIso;
  friend class FieldElement;
  friend class Newform;
private:
  string var;   // name of generator
  int d;        // degree
  ZZX minpoly;  // irreducible poly of degree d
  ZZ denom; // minpoly is the (integral) min poly of A/denom
  mat_m A;        // dxd matrix with scaled min.poly. minpoly
  mat_m B, Binv;  // Binv*B = Bdet*I
  // B has first column Bdet1*e1, Binv has first column Bdet2*e1,
  // Bdet=Bdet1*Bdet2, Bdet3 = denom*Bdet2
  ZZ Bdet, Bdet1, Bdet2, Bdet3, Bcontent;
  // Binv*A*B = Bdet*denom*C
  // cols of Binv are Bdet * coeffs of basis w.r.t. a-powers
  mat_m C;        // dxd companion matrix with min.poly. minpoly
  vector<mat_m> Cpowers;  // C^i for i=0,1,...,d-1
public:
  //Field(const ZZX& p);
  Field(); // defaults to Q
  explicit Field(const mat_m& m, const ZZ& den = to_ZZ(1), string a="a", int verb=0);
  explicit Field(const ZZX& p, string a="a", int verb=0);
  FieldElement rational(const bigrational& x) const;
  FieldElement rational(const ZZ& x) const;
  FieldElement rational(long x) const;
  FieldElement rational(int x) const;
  FieldElement one() const;
  FieldElement minus_one() const;
  FieldElement two() const;
  FieldElement minus_two() const;
  FieldElement zero() const;
  FieldElement gen() const;
  FieldElement element(const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0) const;
  int degree() const {return d;}
  int isQ() const {return this==FieldQQ || d==1;}
  ZZX poly() const {return minpoly;}
  mat_m basis() const {return Binv;} // columns are Bfactor * coeffs of basis w.r.t. a-powers
  //  ZZ basis_factor() const {return Bdet;}
  mat_m inv_basis() const {return B;} // columns are coeffs of a-powers w.r.t. basis
  void display(ostream&s = cout, int raw=0) const; // if raw, also display raw basis
  void display_bases(ostream&s = cout) const; // display powers of A and C and bases in both embeddings
  string get_var() const {return var;}
  void set_var(const string& v)  {var = v;}
  // String for pretty output, like "Q" or "Q(i) = Q[X]/(X^2+1)", or
  // (if raw) raw output, suitable for re-input, like "Q" or "i [1 0 1]":
  string str(int raw=0) const;
  friend istream& operator>>(istream& s, Field* F);
  friend istream& operator>>(istream& s, Field** F);

  // Apply polredabs to the defining polynomial, define a new field
  // with that poly and return an isomorphism from this to that.  If
  // the poly was already polredabsed, or if F is QQ, return the
  // identity. Otherwise a new Field is created with provided variable
  // name.
  FieldIso reduction_isomorphism(string newvar) const;

  // Return an iso from this=Q(a) to Q(b) where B is in this field and generates
  FieldIso change_generator(const FieldElement& b) const;

  // Return an iso from this=Q(a) to Q(b) where b^2=r, optionally
  // applying polredabs to the codomain.  sqrt_r is set to sqrt(r) in
  // the codomain, so sqrt_r^2 = image of r
  FieldIso sqrt_embedding(const FieldElement& r, string newvar, FieldElement& sqrt_r, int reduce=1) const;
};

inline ostream& operator<<(ostream& s, const Field& F)
{ s << F.str();  return s; }

class FieldElement {
  friend class Field;
  friend class FieldIso;
  friend FieldElement evaluate(const ZZX& f, const FieldElement a);
private:
  const Field* F;

  // In general the field element is (1/denom)*coords-combination of power basis of F
  // NB On construction every element will be reduced using cancel()
  vec_m coords; // length F->d
  ZZ denom;     // >=1
  void cancel(); // divides through by gcd(content(coords, denom))

  // When F is Q this is just a wrapper round eclib's bigrational
  // class and coords and denom are ignored
  bigrational val;
public:
  FieldElement()
    :F(FieldQQ) {;}
  explicit FieldElement(const Field* HF)
    :F(HF), coords(vec_m(HF->d)), denom(to_ZZ(1))  {if (HF==FieldQQ) val = bigrational(0);}
  // raw means the given coords are w.r.t. the B-basis
  FieldElement(const Field* HF, const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0);
  // creation from a rational (general F)
  FieldElement(const Field* HF, const ZZ& a, const ZZ& d=to_ZZ(1))
    :F(HF), coords(a*vec_m::unit_vector(HF->d, 1)), denom(d), val(bigrational(a,d)) { cancel();}
  // creation from a rational (F=Q)
  explicit FieldElement(const bigrational& r)
    :F(FieldQQ), val(r) {;}
  // creation from a rational (general F)
  FieldElement(const Field* HF, const bigrational& r)
    :F(HF), coords(r.num()*vec_m::unit_vector(HF->d, 1)), denom(r.den()), val(r) {;}


  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;

  const Field* field() const {return F;}
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
  vec_m get_coords() const {return coords;}
  ZZ get_denom() const {return denom;}
  int operator==(const FieldElement& b) const;
  int operator!=(const FieldElement& b) const;

  FieldElement operator+(const FieldElement& b) const; // add
  FieldElement operator+(const ZZ& b) const {return operator+(FieldElement(F,b));} // add
  void operator+=(const FieldElement& b); // add b to this
  void operator+=(const ZZ& b) { operator+=(FieldElement(F,b));} // add b

  FieldElement operator-(const FieldElement& b) const; // subtract
  FieldElement operator-(const ZZ& b) const {return operator-(FieldElement(F,b));} // subtract
  void operator-=(const FieldElement& b); // subtract b from this
  void operator-=(const ZZ& b) { operator-=(FieldElement(F,b));} // subtract b
  FieldElement operator-() const;                           // unary minus

  FieldElement operator*(const FieldElement& b) const; // product
  FieldElement operator*(const ZZ& b) const {return operator*(FieldElement(F,b));} // product
  void operator*=(const FieldElement& b); // multiply by b
  void operator*=(const ZZ& b) { operator*=(FieldElement(F,b));} // multiply by b

  FieldElement inverse() const; // raise error if zero      // inverse
  FieldElement operator/(const FieldElement& b) const; // divide (raise error if b is zero)
  FieldElement operator/(const ZZ& b) const {return operator/(FieldElement(F,b));} // divide
  void operator/=(const FieldElement& b);                        // divide by b
  void operator/=(const ZZ& b) { operator/=(FieldElement(F,b));} // divide by b
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
  friend istream& operator>>(istream& s, FieldElement& x);
};

inline ostream& operator<<(ostream& s, const FieldElement& x)
{ s << x.str();  return s;}

FieldElement evaluate(const ZZX& f, const FieldElement a);

class FieldIso {
  friend class Field;
  friend class FieldElement;
private:
  const Field* domain;
  const Field* codomain;
  mat_m isomat;
  ZZ denom;
  int id_flag; // flags that this is the identity (domain=codomain and isomat=id)
  void set_id_flag()
  {
    id_flag = (domain==codomain) && is_one(denom) && (isomat==mat_m::identity_matrix(domain->d));
  }
public:
  // Dummy constructor
  FieldIso()
    :domain(FieldQQ), codomain(FieldQQ), isomat(mat_m::identity_matrix(1)), denom(ZZ(1)), id_flag(1) {;}
  // Constructor from a known matrix
  FieldIso(const Field* F1, const Field* F2, const mat_m& M, const ZZ& d = ZZ(1), int id=-1)
    :domain(F1), codomain(F2), isomat(M), denom(d), id_flag(id)
  {
    if (id_flag==-1)
      set_id_flag();
  }
  // Partial constructor from two fields, used when inputting just the matrix and denominator.
  // This initializes isomat to the correct size as the 0 matrix
  FieldIso(const Field* F1, const Field* F2)
    :domain(F1), codomain(F2), isomat(mat_m(F2->degree(),F1->degree())), denom(ZZ(1)), id_flag(0) {;}
  // Identity
  explicit FieldIso(const Field* F1)
    :domain(F1), codomain(F1), isomat(mat_m::identity_matrix(F1->d)), denom(ZZ(1)), id_flag(1) {;}
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
  const Field* dom() const {return domain;}
  const Field* codom() const {return codomain;}
  mat_m matrix() const {return isomat;}
  ZZ den() const {return denom;}

  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;

  // x must be initialised with domain and codomain, this just inputs
  // the matrix and denominator
  friend istream& operator>>(istream& s, FieldIso& x);
};

inline ostream& operator<<(ostream& s, const FieldIso& x)
{ s << x.str(); return s; }


#endif
