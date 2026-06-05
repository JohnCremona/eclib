// File FIELD.CC: classes for working with number fields for Hecke eigenvalues
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

#include "eclib/field.h"
#include "eclib/polred.h"


////////////////////////////////////////////////////////////////////////
//
// Implementation of Field methods
//
////////////////////////////////////////////////////////////////////////

//#define DEBUG_ARITH

//#define DEBUG_FIELD_CONSTRUCTOR

Field::Field(const ZZX& p, string a, int verb)
  :have_integral_basis(0)
{
#ifdef DEBUG_FIELD_CONSTRUCTOR
  cout << "In Field constructor with poly " << ::str(p) << ", name = " << a << endl;
  verb = 1;
#endif
  if (IsMonic(p) && IsIrreducible(p))
    {
      // Set m to be the companion matrix of p.  The function
      // CompanionMatrix(p) returns a mat_ZZ so we do this manually.
      d = deg(p);
      mat_m m(d,d);
      for(int i=1; i<d; i++)
        {
          m(i+1,i) = ZZ(1);
          m(i,d) = -coeff(p, i-1);
        }
      m(d,d) = -coeff(p, d-1);
      // Finally call the other constructor
#ifdef DEBUG_FIELD_CONSTRUCTOR
      cout << "Calling the other Field constructor " << endl;
#endif
      *this = Field(Qmat(m), a, verb);
    }
  else
    {
      cerr << "Error: poly should be monic and irreducible" << endl;
#ifdef DEBUG_FIELD_CONSTRUCTOR
      display_factors(p);
#endif
    }
#ifdef DEBUG_FIELD_CONSTRUCTOR
  cout << "Constructed field " << this << " --> " << *this << endl;
#endif
}

Field::Field() // defaults to Q
  :var(""), d(1), have_integral_basis(0)
{
  SetX(minpoly);
#ifdef DEBUG_FIELD_CONSTRUCTOR
  cout << "Constructed field " << this << " --> " << *this << endl;
#endif
}

Field::Field(const Qmat& A, Qmat& B, Qmat& Binv, string a, int verb)
  : var(a), d(A.nrows()), have_integral_basis(0)
{
  if (verb)
    {
      cout << "----------------------------"<<endl;
      cout << "In Field constructor (var = "<<a<<")"<< endl;
    }
  // NB We assume that A is integral in the sense that its monic char
  // poly is integral and irreducible.
  minpoly = A.charpoly();
  if (verb)
    cout << " - min poly = " << ::str(minpoly) << ", generator " << var << endl;
  //assert(IsMonic(minpoly));

  mat_m C = A.companion_transform(B, Binv);

  // NB The coords of a field element given a 'raw' coord Qvec v is Binv*v.

  Cpowers.resize(d);
  Cpower_traces.init(d);
  Cpowers[0] = mat_m::identity_matrix(d);
  Cpower_traces[1] = ZZ(d); // vec_m indices start at 1
  for (int i=1; i<d; i++)
    {
      Cpowers[i] = C*Cpowers[i-1];
      Cpower_traces[i+1] = Cpowers[i].trace();
    }
  if(verb)
    {
      cout<<"basis  matrix = " << Binv.str()<<endl;
      cout<<"inverse basis matrix = " << B.str() <<endl;
      if (verb>1)
        {
          cout<<"companion matrix  = ";
          C.output_flat(cout);
          cout<<endl;
        }
      cout << "Leaving Field constructor" << endl;
      cout << "----------------------------"<<endl;
    }
#ifdef DEBUG_FIELD_CONSTRUCTOR
  cout << "Constructed field " << this << " --> " << *this << endl;
#endif
}

// String for prettier output, like "Q" or "Q(i) = Q[X]/(X^2+1)" or
// raw output, suitable for re-input, like "Q" or "i [1 0 1]":
string Field::str(int raw) const
{
  // cout << "In Field::str() with field " << this << " ---> " << endl;
  ostringstream s;
  if (d==1)
    s << "Q";
  else
    {
      if (raw)
        s << var << " " << minpoly;
      else
        s << "Q("<<var<<") = Q[X]/(" << ::str(minpoly, "X")<<")";
    }
  return s.str();
}

//#define DEBUG_FIELD_INPUT
void Field::read(istream& s)
{
  s >> var;
  if (var=="Q")
    *this = Field();
  else
    {
      s >> minpoly;
      *this = Field(minpoly, var);
    }
}

istream& operator>>(istream& s, Field& F)
{
  F.read(s);
  return s;
}

void Field::display(ostream&s) const
{
  if (d==1)
    {
      s << "Q" << endl;
    }
  else
    {
      s << "Q(" << var << ") with defining polynomial "<< ::str(minpoly) <<" (of degree "<<d
        << " and discriminant " << discriminant(minpoly) << ")" << endl;
    }
}

FieldElement Field::element(const vec_m& c, const ZZ& d) const
{
  return FieldElement(*this, c, d);
}

FieldElement Field::element(const Qvec& v) const
{
  return FieldElement(*this, v);
}

FieldElement Field::operator()(const bigrational& x) const
{
  return FieldElement(*this, x);
}

FieldElement Field::operator()(const ZZ& x) const
{
  return operator()(bigrational(x));
}

FieldElement Field::operator()(long x) const
{
  return operator()(bigrational(x));
}

FieldElement Field::operator()(int x) const
{
  return operator()(bigrational(x));
}

FieldElement Field::operator()(const Qvec& v) const
{
  return element(v);
}

FieldElement Field::gen() const
{
  if (d==1)
    return operator()(1);
  else
    return FieldElement(*this, Qvec::unit_vector(d, 2));
}

vector<FieldElement> Field::power_basis() const
{
  return powers(gen());
}

// compute integral basis (via libpari), fill integral_basis and basis_change_matrix
void Field::make_integers(const ZZ& bound)
{
  if (!have_integral_basis)
    {
      Integers = MaximalOrder(this, bound);
      have_integral_basis = 1;
    }
}

// recreate integral basis from index and base_change_matrix's inverse (after reading from file)

void Field::set_integers(const ZZ& i, const mat_m& M)
{
  Integers = Order(*this, i, M);
  have_integral_basis = 1;
}

////////////////////////////////////////////////////////////////////////
//
// Implementation of FieldElement methods
//
////////////////////////////////////////////////////////////////////////

// creation of 0 in the given field
FieldElement::FieldElement(const Field& HF)
    :F(&HF), v(vec_m(F->d))
{
  if (F->d==1)
    val = bigrational(0);
}

// creation from an integer vector of coords with denominator
FieldElement::FieldElement(const Field& HF, const vec_m& c, const ZZ& d)
  :F(&HF), v(c,d)
{
  if (F->d==1)
    val = c[1];
}

// creation from a rational vector of coords
FieldElement::FieldElement(const Field& HF, const Qvec& c)
    :F(&HF), v(c)
{
  if (F->d==1)
    val = c[1];
}

// creation from a rational (general F)
FieldElement::FieldElement(const Field& HF, const ZZ& a, const ZZ& d)
  :F(&HF), v(a*vec_m::unit_vector(F->d, 1), d), val(bigrational(a,d))
{
  ;
}

// creation from a rational
FieldElement::FieldElement(const Field& HF, const bigrational& r)
  :F(&HF), v(r.num()*vec_m::unit_vector(F->d, 1), r.den()), val(r)
{
  ;
}

int FieldElement::field_degree() const
{
  return F->d;
}

int FieldElement::field_is_Q() const
{
  return F->d == 1;
}

int FieldElement::is_generator() const
{
  return degree()==F->d;
}

int FieldElement::in_same_field(const FieldElement& b) const
{
  return (F==b.F) || (*F==*b.F);
}

int FieldElement::is_zero() const
{
  if (field_is_Q())
    return val.is_zero();
  else
    return v.is_zero();
}

int FieldElement::is_one() const
{
  static const bigrational one(1);
  if (field_is_Q())
    return val == one;
  else
    return v == Qvec::unit_vector(F->d,1);
}

int FieldElement::is_minus_one() const
{
  static const bigrational minus_one(-1);
  if (field_is_Q())
    return val == minus_one;
  else
    return v == - Qvec::unit_vector(F->d,1);
}

void FieldElement::negate() // negate in place
{
  if (field_is_Q())
    val = -val;
  else
    v = -v;
}

// String for pretty printing, used in default <<
// or for raw output, suitable for re-input (with Field known):
string FieldElement::str(int raw) const
{
  ostringstream s;
  if (field_is_Q())
    {
      s << val; // if val is an integer this just outputs the integer with no "/" + denominator
      return s.str();
    }
  if (raw)
    s << v.str(1);
  else
    {
      string n = ::str(v.numerator, F->var);
      if (n[0]=='+')
        n.erase(0,1);
      if (v.denom==1)
        s << n;
      else
        s << "(" << n << ")/" << v.denom;
    }
  return s.str();
}

// If not Q, v must be initialised with the right size
void FieldElement::read(istream& s)
{
  if (field_is_Q())
    s >> val;
  else
    s >> v;
}

// x must be initialised with a Field before input to x
istream& operator>>(istream& s, FieldElement& x)
{
  x.read(s);
  return s;
}

int FieldElement::operator==(const FieldElement& b) const
{
  return in_same_field(b) && ( field_is_Q()? val==b.val : v==b.v);
}

int FieldElement::operator!=(const FieldElement& b) const
{
  return (!in_same_field(b)) || ( field_is_Q()? val!=b.val : v!=b.v);
}

mat_m FieldElement::matrix() const // ignores denom, not used for Q
{
  if (field_is_Q())
    return mat_m::scalar_matrix(1, num(val));
  return lin_comb_mats(v.numerator, F->Cpowers);
}

// NB Since we do not have polynomials with rational coefficients,
// both charpoly and minpoly are scaled to be primitive rather than
// monic.

ZZX FieldElement::charpoly() const
{
  ZZX cp;
  bigrational r;
  if (is_rational(r)) // then charpoly = minpoly^degree
    {
      ZZX mp = minpoly();
      cp = mp;
      int d = F->d - 1;
      while (d--)
        cp *= mp;
      return cp;
    }

  // Now this is not rational (and the field is not QQ)
  cp = ::charpoly(mat_to_mat_ZZ(matrix()));
  // now replace X by denom*X and make primitive
  ZZ dpow(v.denom);
  int d = F->d;
  for (int i=1; i<=d; i++)
    {
      SetCoeff(cp, i, dpow*coeff(cp, i));
      if (i!=d)
        dpow *= v.denom;
    }
  return PrimitivePart(cp);
}

// The charpoly is a power of the irreducible minpoly. NB This
// primitive integer polynomial is only the actual minpoly if it is
// monic.
ZZX FieldElement::minpoly() const
{
  ZZX cp;
  bigrational r;
  if (is_rational(r)) // then deg(minpoly)=1, whatever the field
    {
      SetX(cp);
      SetCoeff(cp, 0, -num(r));
      SetCoeff(cp, 1, den(r));
      return cp;
    }

  // Now this is not rational (and the field is not QQ)
  cp = charpoly();
  return factor(cp)[0].a;
}

bigrational FieldElement::norm() const
{
  int d = F->d;
  if (d==1) return val;
  const ZZ& pow_den = pow(v.denom, d);
  bigrational r;
  if (is_rational(r)) // then norm = r**degree
    return bigrational(pow(r.num(), d), pow_den);
  else
    return bigrational(matrix().determinant(), pow_den);
}

bigrational FieldElement::trace() const
{
  int d = F->d;
  if (d==1) return val;
  bigrational r;
  if (is_rational(r)) // then trace = r*degree
    return ZZ(d) * r;
  else
    return bigrational(v.numerator * F->Cpower_traces, v.denom);
}

void FieldElement::set_zero()
{
  val = bigrational(0);
  v.set_zero();
}

void FieldElement::set_one()
{
  val = bigrational(1);
  v.set_unit_vector(1);
}

FieldElement FieldElement::operator+(const ZZ& b) const {return operator+(F->operator()(b));} // add
FieldElement FieldElement::operator+(const int& b) const {return operator+(F->operator()(b));} // add
FieldElement FieldElement::operator+(const long& b) const {return operator+(F->operator()(b));} // add
void FieldElement::operator+=(const ZZ& b) { operator+=(F->operator()(b));} // add b
void FieldElement::operator+=(const int& b) { operator+=(F->operator()(b));} // add b
void FieldElement::operator+=(const long& b) { operator+=(F->operator()(b));} // add b

// add b to this
void FieldElement::operator+=(const FieldElement& b)
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      cerr << "In operator += with this = " << (*this) << " in field " << (*F)
           << " and that = " << b << " in field " << *(b.F) << endl;
      exit(1);
    }
  if (b.is_zero())
    return;
  if (field_is_Q())
    val += b.val;
  else
    v += b.v;
}

FieldElement FieldElement::operator+(const FieldElement& b) const
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  FieldElement a = *this;
  if (!b.is_zero())
    a += b;
  return a;
}

FieldElement FieldElement::operator-() const
{
  if (field_is_Q())
    return FieldElement(*F, -val);
  else
    return FieldElement(*F, -v);
}

FieldElement FieldElement::operator-(const ZZ& b) const {return operator-(F->operator()(b));} // subtract
FieldElement FieldElement::operator-(const int& b) const {return operator-(F->operator()(b));} // subtract
FieldElement FieldElement::operator-(const long& b) const {return operator-(F->operator()(b));} // subtract

// subtract b
void FieldElement::operator-=(const FieldElement& b)
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (b.is_zero())
    return;
  if (field_is_Q())
    val -= b.val;
  else
    v -=b.v;
}

void FieldElement::operator-=(const ZZ& b) { operator-=(F->operator()(b));} // subtract b
void FieldElement::operator-=(const int& b) { operator-=(F->operator()(b));} // subtract b
void FieldElement::operator-=(const long& b) { operator-=(F->operator()(b));} // subtract b

FieldElement FieldElement::operator-(const FieldElement& b) const
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  FieldElement a = *this;
  if (!b.is_zero())
    a -= b;
  return a;
}

void FieldElement::operator*=(const FieldElement& b) // multiply by b
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (is_zero()) return;
  if (b.is_zero()) {set_zero(); return;}
  if (field_is_Q())
    val *= b.val;
  else
    v = Qvec((matrix()*b.matrix()).col(1),  get_denom() * b.get_denom());
}

void FieldElement::operator*=(const ZZ& b) { operator*=(F->operator()(b));} // multiply by b
void FieldElement::operator*=(const int& b) { operator*=(F->operator()(b));} // multiply by b
void FieldElement::operator*=(const long& b) { operator*=(F->operator()(b));} // multiply by b

FieldElement FieldElement::operator*(const FieldElement& b) const
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (b.is_zero()) return b;
  FieldElement a = *this;
  if (!is_zero())
    a *= b;
  return a;
}

FieldElement FieldElement::operator*(const ZZ& b) const {return operator*(F->operator()(b));} // product
FieldElement FieldElement::operator*(const int& b) const {return operator*(F->operator()(b));} // product
FieldElement FieldElement::operator*(const long& b) const {return operator*(F->operator()(b));} // product

FieldElement FieldElement::inverse() const // raise error if zero
{
  if (is_zero())
    {
      cerr << "Attempt to invert zero!" << endl;
      exit(1);
    }
  if (is_one() || is_minus_one())
    return FieldElement(*this);
  if (field_is_Q())
    return FieldElement(*F, recip(val));

  mat_m M = matrix(), Minv;
  ZZ Mdet = ::inverse(M,Minv); // so M*Minv = Mdet*identity
  FieldElement ans = FieldElement(*F, v.denom*Minv.col(1), Mdet);
  // assert (operator*(ans).is_one());
  return ans;
}

void FieldElement::operator/=(const FieldElement& b)      // divide by b
{
  if (b.is_zero())
    {
      cerr << "Attempt to divide by zero!" << endl;
      exit(1);
    }
  if (!in_same_field(b))
    {
      cerr << "Attempt to divide elements of different fields!" << endl;
      cerr << "LHS: " << (*this) << ", " << flush;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS: " << b << ", " << flush;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (is_zero())
    return;
  if (field_is_Q())
    {
      val /= b.val;
      return;
    }
  operator*=(b.inverse());
}

void FieldElement::operator/=(const ZZ& b) { operator/=(F->operator()(b));} // divide by b
void FieldElement::operator/=(const int& b) { operator/=(F->operator()(b));} // divide by b
void FieldElement::operator/=(const long& b) { operator/=(F->operator()(b));} // divide by b

FieldElement FieldElement::operator/(const FieldElement& b) const // raise error if b is zero
{
  if (b.is_zero())
    {
      cerr << "Attempt to divide by zero!" << endl;
      exit(1);
    }
  if (!in_same_field(b))
    {
      cerr << "Attempt to divide elements of different fields!" << endl;
      cerr << "LHS: " << (*this) << ", " << flush;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS: " << b << ", " << flush;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  FieldElement a = *this;
  if (!is_zero())
    a /= b;
  return a;
}

FieldElement FieldElement::operator/(const ZZ& b) const {return operator/(F->operator()(b));} // divide
FieldElement FieldElement::operator/(const int& b) const {return operator/(F->operator()(b));} // divide
FieldElement FieldElement::operator/(const long& b) const {return operator/(F->operator()(b));} // divide

FieldElement evaluate(const ZZX& f, const FieldElement a)
{
  FieldElement fa = (*(a.F))(0);
  for(int i=deg(f); i>=0; i--)
    {
      fa += coeff(f,i);
      if(i)
        fa *= a;
    }
  return fa;

}

// return 1 and set r to the rational value if the degree is 1
int FieldElement::is_rational(bigrational& r) const
{
  r = val;
  if (F->d ==1)
    return 1;
  r = bigrational(v.numerator[1], v.denom);
  for (int i=2; i <= F->d; i++)
    {
      if (v.numerator[i]!=0)
        return 0;
    }
  return 1;
}

// return 1 iff this is an algebraic integer.  Since the field's
// defining polynomial is monic integral, it is sufficient, but not
// necessary that the coordinate vector be integral:
int FieldElement::is_integral() const
{
  return (get_denom()==ZZ(1)) || IsMonic(charpoly());
}

// NB for a in F, either [Q(sqrt(a))=Q(a)] or [Q(sqrt(a)):Q(a)]=2.
// The first function only applies when a has maximal degree:
// return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
int FieldElement::is_absolute_square(FieldElement& r) const
{
  if (field_is_Q())
    return (val.is_square(r.val));
  // field not Q, reduce to integral case if necessary
  const ZZ& vden = v.denom;
  if (::is_one(vden))
    return is_absolute_integral_square(r);
  FieldElement x = FieldElement(*F, vden*v.numerator);
  int res = x.is_absolute_integral_square(r);
  if (res)
    r /= vden;
  return res;
}

// Same as above if the min poly is known and denom=1
int FieldElement::is_absolute_integral_square(FieldElement& r)  const
{
  if (field_is_Q())
    return (val.is_square(r.val));
  // assert (::is_one(v.denom));
  ZZX g, g0, g1, f = minpoly();
  if (::is_square(f, g)) // Tests if f(x^2) = (+/-) g(x)*g(-x)
    {
      parity_split(g, g0, g1); // g(x) = g0(x^2) + x*g1(x^2)
      // Now 0 = g(r) = g0(a)+r*g1(a) and g1(a)!=0
      r = - evaluate(g0,*this)/evaluate(g1,*this);
      if (r*r== *this)
        return 1;
      else
        {
          cout << (*this) << " has scaled min poly " << ::str(f)
               << " whose double has factor " << ::str(g)
               << " with even part g0 = " << ::str(g0)
               << " and odd part g1 = " << ::str(g1) << endl;
          cout << "These evaluate to " << evaluate(g0,*this) << " and " << evaluate(g1,*this)
               << " with negative quotient/denom r = " << r
               << " but r*r = " << r*r << endl;
          return 0;
        }
    }
  return 0;
}

// The second function applies in general:
// return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
//#define DEBUG_IS_SQUARE
int FieldElement::is_square(FieldElement& r, int ntries) const
{
#ifdef DEBUG_IS_SQUARE
  cout << "Testing whether " << (*this) << " is a square" << endl;
#endif
  r = *this; // sets the field, and is a default
  if (is_zero() || is_one())
    return 1;

  if (field_is_Q())
    return (val.is_square(r.val));

  // If this has full degree, or if the relative degree is odd, just
  // call is_absolute_square()
  int d = F->d;
  ZZX m = minpoly();
  if ((d/deg(m))%2)
    {
#ifdef DEBUG_IS_SQUARE
      cout << "Is_square() works directly" << endl;
#endif
      int res = is_absolute_square(r);
#ifdef DEBUG_IS_SQUARE
      if (res)
        cout << (*this) << " is square with sqrt(" << (*this) << ") = " << r << endl;
      else
        cout << (*this) << " is not a square" << endl;
#endif
      return res;
    }
#ifdef DEBUG_IS_SQUARE
  cout << "Is_square() multipying by successive squares" << endl;
#endif

  // Otherwise we multiply this by a "random" square to have full degree first
  FieldElement b = F->gen();
  for (int i=0; i<ntries; i++, b+=to_ZZ(1))
    {
      FieldElement abb = (*this)*b*b, rb(*F);
      ZZX mb = abb.minpoly();
      if ((d/deg(mb))%2)
        {
#ifdef DEBUG_IS_SQUARE
          cout << "Is_square() succeeds after " << i+1 << " tries, using b = " << b << endl;
#endif
          int res = abb.is_absolute_square(rb);
          if (res)
            r = rb/b;
          return res;
        }
      // else keep trying
    } // end of loop over shifts i
  cout << "is_square() fails on " << (*this) << " after " << ntries << " tries" << endl;
  return 0;
}

Qmat FieldElement::power_matrix() const // cols are coords of powers
{
  return coord_matrix(powers(*this));
}

////////////////////////////////////////////////////////////////////////
//
// Implementation of FieldIso methods
//
////////////////////////////////////////////////////////////////////////

//#define DEBUG_COMPOSE

// precompose this FieldIso with another (requires iso.codomain = domain)
void FieldIso::precompose(const FieldIso& iso)
{
#ifdef DEBUG_COMPOSE
  cout << " precomposing " << *this << " with " << iso << endl;
#endif
  if (iso.is_identity())
    return;
  if (is_identity())
    {
      *this = iso;
      return;
    }
    if (domain==iso.codomain)
    {
      isomat = isomat * iso.isomat;
      domain = iso.domain;
#ifdef DEBUG_COMPOSE
      cout << " after precomposing, the iso is " << *this << endl;
#endif
    }
  else
    {
      cerr << "Cannot precompose " << *this << " with " << iso << endl;
      exit(1);
    }
}

// postcompose this FieldIso with another (requires iso.domain = codomain)
void FieldIso::postcompose(const FieldIso& iso)
{
#ifdef DEBUG_COMPOSE
  cout << " postcomposing " << *this << " with " << iso << endl;
#endif
  if (iso.is_identity())
    return;
  if (is_identity())
    {
      *this = iso;
      return;
    }
  if (codomain==iso.domain)
    {
      isomat = iso.isomat * isomat;
      codomain = iso.codomain;
#ifdef DEBUG_COMPOSE
      cout << " after postcomposing, the iso is " << *this << endl;
#endif
    }
  else
    {
      cerr << "Cannot postcompose " << *this << " with " << iso << endl;
      exit(1);
    }
}

// return postcomposion of this and iso (requires iso.domain = codomain)
FieldIso FieldIso::operator*(const FieldIso& iso)
{
  FieldIso comp(*this);
  comp.postcompose(iso);
  return comp;
}

string FieldIso::str(int raw) const
{
  if (raw)
    return isomat.str(raw);
  ostringstream s;
  //cout << "(domain = " << domain << ", codomain = " << codomain << "): " << endl;
  //s << "(domain = " << domain << ", codomain = " << codomain << "): " << endl;
  if (id_flag)
    s << "Identity automorphism of " << domain->str();
  else
    {
      if (domain==codomain)
        {
          s << "Automorphism of " << domain->str();
          //cout << "Automorphism of " << domain->str() << endl;
        }
      else
        {
          s << "Embedding of " << domain->str() << " into " << codomain->str();
          //cout << "Embedding of " << domain->str() << " into " << codomain->str();
          // s << " with matrix\n" << isomat;
          // if (!IsOne(denom)) s << "/ "<< denom;
        }
      s << " mapping " << domain->gen() << " to " << operator()(domain->gen());
      //cout << " mapping " << domain->gen() << " to " << operator()(domain->gen()) << endl;
    }
  return s.str();
}

// x must be initialised with domain and codomain, this just inputs
// the matrix and denominator
void FieldIso::read(istream& s)
{
  isomat = Qmat(codomain->degree(), domain->degree());
  s >> isomat;
  set_id_flag();
}

istream& operator>>(istream& s, FieldIso& x)
{
  x.read(s);
  return s;
}

// inverse isomorphism
FieldIso FieldIso::inverse() const
{
  if (id_flag)
    return *this;
  else
    return FieldIso(*codomain, *domain, isomat.inverse(), id_flag);
}

// map x in domain to an element of the codomain
FieldElement FieldIso::operator()(const FieldElement& x) const
{
  if (id_flag) return x;
  if ((x.field_ptr()==domain) || (*x.field_ptr()==*domain))
    {
      bigrational r;
      if (x.is_rational(r))
        return (*codomain)(r);

      FieldElement y(*codomain, isomat * x.v);

      // sanity check that the min poly has not changed
      if (! (x.minpoly()==y.minpoly()))
        {
          cout << "Error in applying field isomorphism\n"
               << *this
               << "\n to x = " << x << "\n --> y = " << y << endl;
          cout << "x has minpoly "<< ::str(x.minpoly()) << endl;
          cout << "y has minpoly "<< ::str(y.minpoly()) << endl;
          exit(1);
        }
      return y;
    }
  cerr << "Error in FieldIso::operator(): domain = " << domain << flush << " --> " << *domain
       << ", argument x = " << x << " in field " << x.field_ptr() << flush << " --> " << *x.field_ptr() <<endl;
  cerr << "Cannot apply FieldIso\n" << *this << "\n to " << x << " in " << *(x.field_ptr()) << endl;
  exit(1);
  return FieldElement(*codomain);
}

// same to all in a list of elements of the domain
vector<FieldElement> FieldIso::operator()(const vector<FieldElement>& x) const
{
  vector<FieldElement> y(x.size());
  std::transform(x.begin(), x.end(), y.begin(), [this](const FieldElement& a){return (*this)(a);});
  return y;
}

//#define DEBUG_REDUCE

// Apply polredabs (if canonical) or polredbest to the defining
// polynomial, define a new field with that poly and return an
// isomorphism from this to that.  If the poly was already
// polredabsed, or if F is QQ, return the identity. Otherwise a new
// Field is created with provided variable name.
FieldIso Field::reduction_isomorphism(string newvar, Field& Fred, int canonical) const
{
#ifdef DEBUG_REDUCE
  cout << "In Field::reduction_isomorphism(), minpoly = " << ::str(minpoly) << endl;
#endif
  Fred = *this;
  if (d==1)
    {
#ifdef DEBUG_REDUCE
      cout << " - Field is Q, so identity" << endl;
#endif
      return FieldIso(*this, Fred,  Qmat::identity(d));
    }
  ZZX h; ZZ denh;
  ZZX g = polred(minpoly, h, denh, canonical);

  if (minpoly==g)
    {
#ifdef DEBUG_REDUCE
      cout << " - " << ::str(minpoly) << " is already reduced, so identity" << endl;
#endif
      return FieldIso(*this, Fred,  Qmat::identity(d));
    }
#ifdef DEBUG_REDUCE
  cout << " - reduced minpoly = " << ::str(g) << endl;
#endif

  // construct the reduced field:
  Fred = Field(g, newvar);
#ifdef DEBUG_REDUCE
  cout << " - reduced field is\n" << Fred << endl;
#endif
  // Construct the isomorphism matrix from F to Fred.

  FieldElement a = evaluate(h,Fred.gen()) / denh;
#ifdef DEBUG_REDUCE
  cout << " - image of gen is (" << a << ") " << endl;
#endif
  Qmat M = a.power_matrix();
#ifdef DEBUG_REDUCE
  cout << " - iso matrix = \n" << M << endl;
#endif
  return FieldIso(*this, Fred, M, 0); // 0: not the identity
}

//#define DEBUG_CHANGE_GEN
// Return an iso from this=Q(a) to Q(b) where B is in this field and generates
FieldIso Field::change_generator(const FieldElement& b, Field& Qb) const
{
  // default (trivial extension)
  Qb = *this; // copy of this field
  FieldIso iso(*this);

  if (b.field_ptr() != this)
    {
      cerr << "Cannot change generator of " << *this << " to " << b
           << " which is in a different field " << *b.field_ptr() << endl;
      return iso;
    }
  ZZX b_pol(b.minpoly());
#ifdef DEBUG_CHANGE_GEN
  cout << "Constructing isomorphism from " << *this << " to Q(b) with b = "
       << b << ", minpoly(b) = " << b_pol << endl;
#endif
  if (deg(b_pol)!=d)
    {
      cerr << "Cannot change generator of " << *this << " to " << b
           << " which only has degree " << deg(b_pol) << endl;
      return iso;
    }
  if (!IsMonic(b_pol))
    {
      cerr << "Cannot change generator of " << *this << " to " << b
           << " which is not integral " << endl;
      return iso;
    }
  if ((d==1) || (b_pol==minpoly)) // then the identity works but the
                                  // codomain must be a pointer to the
                                  // supplied field
    return FieldIso(*this, Qb,  Qmat::identity(d));

  // Create the new field:
  Qb = Field(b_pol, var+string("1"));

  // Let B be the matrix whose columns are the powers of b expressed
  // in the original a-power basis; then its inverse is the matrix of
  // the isomorphism.

  Qmat B = b.power_matrix();
  iso = FieldIso(*this, Qb, B.inverse(), 0);

  // check:
  FieldElement isob = iso(b);
  if (isob != Qb.gen())
    {
      cerr << "Error in Field::change_generator(b) with b = " << b << "\n";
      cerr << "b has minpoly " << ::str(b.minpoly()) << "\n";
      cerr << "iso(b) = " << isob << " with minpoly " << ::str(isob.minpoly()) << "\n";
      exit(1);
    }
#ifdef DEBUG_CHANGE_GEN
  cout << "Field::change_generator() with this = " << this
       << " returns an iso with domain " << iso.domain << " and codomain " << iso.codomain << endl;
#endif
  return iso;
}

// Return an embedding from K to K(sqrt(r)), optionally applying
// polredabs (if reduce=2) or polredbest (if reduce=1) to the
// codomain.  sqrt_r is set to sqrt(r) in the codomain, so sqrt_r^2 =
// image of r.  Caller supplies a reference to the target field which
// will be overwritten (so caller already has the pointer to it).

//#define DEBUG_SQRT_EMBEDDING
FieldIso Field::sqrt_embedding(const FieldElement& r, string newvar, Field& F_sqrt_r, FieldElement& sqrt_r, int reduce) const
{
#ifdef DEBUG_SQRT_EMBEDDING
  cout << "In simple sqrt_embedding() with base field " << this << " --> " << *this << " and r = " << r
       << ", target field pointer = " << &F_sqrt_r << endl;
#endif
  if (r.field_ptr() != this)
    {
      cerr << "Cannot adjoin sqrt(" << r << ") to " << *this
           << " as it is in a different field " << *r.field_ptr() << endl;
      return FieldIso(*this);
    }
  if (r.is_square(sqrt_r))
    {
      cout << "Adjoining sqrt(" << r << ") to " << *this << " is trivial since "
           << r << " is already a square, with root " << sqrt_r << endl;
      F_sqrt_r = *this;
      return FieldIso(*this, F_sqrt_r, Qmat::identity(d));
    }

  // If r has degree < d we replace it with an equivalent element
  // (modulo squares) of maximal degree by multiplying by a square.
  // Then the sqrt field is generated by f(X^2) where f is the char
  // poly of r.
  FieldElement one(operator()(1));
  FieldElement rss(r), s(one);
  if (rss.degree()<d)
    {
      s = gen();
      rss = r*s*s;
      while (rss.degree()<d)
        {
          s += one;
          rss = r*s*s;
        }
    }
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " s = " << s << ", rss = " << rss << endl;
#endif

  ZZX sqrt_rss_pol = XtoX2(rss.charpoly());
  // assert (IsIrreducible(sqrt_rss_pol)); // must be else r is a square

  Field F_sqrt_r_orig = Field(sqrt_rss_pol, newvar); // the extension field before reduction
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " extended field (before any reduction) is " << &F_sqrt_r_orig << " --> " << F_sqrt_r_orig << endl;
#endif

  // Now we embed the original field into the new field in three steps:
  // Q(a) -~-> Q(rss) --> Q(sqrt(rss)) -~-> Q(b)
  // The first and last are isomorphisms, the last (optional) is polredabs reduction.

  Field Qrss; // this field expressed as Q(r*s*s), for s such that r*s*s has full degree
  FieldIso iso(change_generator(rss, Qrss));  // Q(a) -> Q(rss)
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " After changing generator to rss, Qrss = " << &Qrss << " --> " << Qrss << endl;
  cout << " first map (isomorphism) is " << iso << endl;
#endif
  s = iso(s); // image of s in Q(rss)
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " iso(s) = " << s << endl;
#endif

  // the images of the powers of rss are the even powers of sqrt_rss:
  mat_m isomat(2*d, d);
  for (int j=0; j<d; j++)
    isomat.setrow(2*j+1, vec_m::unit_vector(d,j+1));
  FieldIso emb(Qrss, F_sqrt_r_orig, isomat); // Q(rss) -> Q(sqrt(rss))
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " second map (embedding) is " << emb << endl;
#endif
  s = emb(s); // image of s in Q(sqrt(rss))
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " emb(s) = " << s << endl;
#endif
  sqrt_r = F_sqrt_r_orig.gen()/s; // sqrt(r) = sqrt(rss)/s
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " sqrt(r) = " << sqrt_r << endl;
#endif
  emb.precompose(iso);
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " embedding (without reduction) is " << emb << endl;
#endif
  if (reduce)
    {
#ifdef DEBUG_SQRT_EMBEDDING
      cout << " reducing via "
           << (reduce>1? "polredabs" : "polredbest")
           << "..." << endl;
#endif
      FieldIso red = F_sqrt_r_orig.reduction_isomorphism(newvar, F_sqrt_r, reduce>1);
      emb.postcompose(red);
      sqrt_r = red(sqrt_r);
#ifdef DEBUG_SQRT_EMBEDDING
      cout << " third map (reduction isomorphism) is " << red << endl;
#endif
    }
  else // No reduction
    {
      FieldIso red(F_sqrt_r_orig, F_sqrt_r, Qmat::identity(F_sqrt_r_orig.degree()));
      F_sqrt_r = F_sqrt_r_orig;
      emb.postcompose(red);
      sqrt_r = red(sqrt_r);
#ifdef DEBUG_SQRT_EMBEDDING
      cout << " third map (identity isomorphism) is " << red << endl;
#endif
    }
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " final embedding is " << emb << endl;
  cout << " In the extension, sqrt(r) = " << sqrt_r << endl;
#endif
  return emb;
}

// Return an embedding from K to K(sqrt(r1),sqrt(r2),...), optionally
// applying polredabs (if reduce=2) or polredbest (if reduce=1) to the
// codomain.  sqrts is set to [sqrt(r1), sqrt(r2), ..], a list of
// elements of the codoimain, so (sqrt_ri)^2 = image of ri.  Caller
// supplies a reference to the target field which will be overwritten
// (so caller already has the pointer to it).
FieldIso Field::sqrt_embedding(const vector<FieldElement>& r_list, string newvar,
                               Field& F_sqrt_r, vector<FieldElement>& sqrt_r_list,
                               int reduce) const
{
#ifdef DEBUG_SQRT_EMBEDDING
  cout << "In multiple (" << r_list.size() << ") sqrt_embedding() with base field "
       << this << " --> " << *this << " and r_list = " << r_list
       << ", target field pointer = " << &F_sqrt_r << endl;
  for (const auto& r: r_list)
    {
      cout << r.field_ptr() << " --> " << r << endl;
      assert (this == r.field_ptr());
    }
#endif

  // Trivial case where r_list is empty:
  if (r_list.empty())
    {
#ifdef DEBUG_SQRT_EMBEDDING
      cout << " sqrt_embedding() returning trivial embedding into self " << endl;
#endif
      sqrt_r_list.clear();
      F_sqrt_r = *this;
      return FieldIso(*this, F_sqrt_r, Qmat::identity(d));
    }

  // First call the special case to construct K(sqrt(r1)), then recurse

  Field F_sqrt_r1;
  FieldElement sqrt_r1(F_sqrt_r1);
  FieldIso emb1 = sqrt_embedding(r_list.front(), newvar, F_sqrt_r1, sqrt_r1, 0); // no reduction here
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " first embedding is " << emb1 <<endl;
#endif

  // Embed the remaining ri (if any) into F_sqrt_r1:
  vector<FieldElement> r_list_1(r_list.size()-1);
  std::transform(r_list.begin()+1, r_list.end(), r_list_1.begin(),
                 [emb1](const FieldElement& r){ return emb1(r);});

  // Use recursion to extend F_sqrt_r1:
  FieldIso emb2 = F_sqrt_r1.sqrt_embedding(r_list_1, newvar, F_sqrt_r, sqrt_r_list, reduce);

#ifdef DEBUG_SQRT_EMBEDDING
  cout << " second embedding is " << emb1 <<endl;
#endif

  // Insert emb2(sqrt_r1) into sqrt_r_list:
  sqrt_r_list.insert(sqrt_r_list.begin(), emb2(sqrt_r1));

  // Compose the two embeddings:
  return emb1*emb2;
}


////////////////////////////////////////////////////////////////////////
//
// Implementation of Order methods
//
////////////////////////////////////////////////////////////////////////

Order::Order(const Field& K) // equation order
  : F(&K),
    rank(K.d),
    power_coords_matrix(mat_m::identity_matrix(K.d)),
    basis_matrix(Qmat(power_coords_matrix)),
    index(1),
    poldisc(discriminant(K.minpoly)),
    disc(poldisc)
{
  ;
}

Order::Order(const Field& K, const vector<FieldElement>& v, int basis)
  : F(&K),
    rank(F->d),
    poldisc(discriminant(F->minpoly))
{
  if (basis)
    {
      basis_matrix = coord_matrix(v);
      power_coords_matrix = basis_matrix.inverse().get_numerator();
      index = abs(power_coords_matrix.determinant());
      disc = poldisc/(index*index);
    }
  else
    {
      Qmat M = HNF(coord_matrix(v, 1, 1)); // 1: rows, 1: reverse coords
      M.delete_rows(v.size()-rank);   // rows are coords of new Z-basis
      const auto& bas = from_coord_matrix(*F, M, 1, 1); // 1: rows, 1: reverse coords
      basis_matrix = coord_matrix(bas); // = M with cols reversed
      power_coords_matrix = basis_matrix.inverse().get_numerator();
      index = abs(power_coords_matrix.determinant());
      disc = poldisc/(index*index);
    }
}

Order::Order(const Field& K, const vector<FieldElement>& v, const mat_m pcm)
  : F(&K),
    rank(F->d),
    power_coords_matrix(pcm),
    basis_matrix(coord_matrix(v)),
    index(abs(pcm.determinant())),
    poldisc(discriminant(F->minpoly)),
    disc(poldisc/(index*index))
{
  ;
}

Order::Order(const Field& K, const ZZ& i, const mat_m& M) // Order in F given pcm
  : F(&K),
    rank(F->d),
    power_coords_matrix(M),
    basis_matrix(Qmat(M).inverse()),
    index(abs(M.determinant())),
    poldisc(discriminant(F->minpoly)),
    disc(poldisc/(index*index))
{
  ;
}

vector<FieldElement> Order::get_basis() const
{
  return from_coord_matrix(*F, basis_matrix);
}

FieldElement Order::basis_elt(int i) const
{
  return FieldElement(*F, basis_matrix.col(i));
}

// coords w.r.t. Z-basis of an arbitrary element of F
Qvec Order::coords(const FieldElement& a) const
{
  return power_coords_matrix * a.coords(); // not a.v in case field is Q
}

// coords w.r.t. Z-basis of an element of F in this order
vec_m Order::integral_coords(const FieldElement& a) const
{
  Qvec c = coords(a);
  if (!is_one(c.denom))
    cout << "Error in computing integral_coords(a) for a = " << a
         << " in order " << (*this)
         << ": it has denominator " << c.denom
         << " so is not in the order" << endl;
  return c.numerator;
}

// membership test
int Order::contains(const FieldElement& a) const
{
  return (F==a.field_ptr()) && is_one(denom(a));
}

// membership test returning coords
int Order::contains(const FieldElement& a, vec_m& c) const
{
  if (F!=a.field_ptr())
    return 0;
  Qvec qc = coords(a);
  c = qc.numerator;
  return is_one(qc.denom);
}

// containment test
int Order::contains(const Order& O2) const
{
  auto M = basis_matrix * O2.power_coords_matrix;
  return M.is_integral();
}

// FieldElement from integer coords
FieldElement Order::operator()(const vec_m& coords) const
{
  return (*F)(basis_matrix*Qvec(coords));
}

// FieldElement from rational coords
FieldElement Order::operator()(const Qvec& coords) const
{
  return operator()(coords.numerator) / coords.denom;
}

string Order::str(int raw) const
{
  if (rank==0) // only possible after a default construction woth no field
    return string("Empty order: no field");

  ostringstream s;

  auto bas = get_basis();
  if (raw) // only output the integral basis coords
    {
      for (auto bi: bas)
        s << bi.str(1) << "\n";
    }
  else
    {
      s << "Order in " << *F << " with Z-basis " << bas;
    }
  return s.str();
}

// Sum of two orders (in the same field!)
Order Order::operator+(const Order& O2)
{
  if (F != O2.F)
    {
      cerr << "Cannot add orders in defferent fields" << endl;
      exit(1);
    }
  vector<FieldElement> gens;
  gens.reserve(rank*rank);
  for (auto a: O2.get_basis())
    for (auto b: get_basis())
      gens.push_back(a*b);
  return Order(*F, gens, 0); // 0 means not a basis, just Z-gens
}

// Add another order to this:
Order Order::operator+=(const Order& O2)
{
  *this = *this + O2;
  return *this;
}

// Sum of this and the power order of a
Order Order::operator+(const FieldElement& a)
{
  if (contains(a))
    return (*this);
  else
    return (*this) + PowerOrder(a);
}

// Add the power order of a to this
Order Order::operator+=(const FieldElement& a)
{
  if (!contains(a))
    *this = *this + a;
  return *this;
}

// Check if this is just the equation order (with standard basis in standard order)
int Order::is_equation_order() const
{
  return is_identity_matrix(power_coords_matrix);
}

// Functions to enlarge the order

// Add one a (algebraic integer) to the Z-span: the result may not be
// an order but this is only used internally.  NB This does not
// recompute index or disc: save these until we have an actual order.
void Order::add_one(const FieldElement& a, int check)
{
  if (contains(a))
    return;
  if (check && !a.is_integral())
    return;
  add_one(a.coords());
}

// Same given just the field coords of the new element
void Order::add_one(const Qvec& v)
{
  if (is_one(v.denom))
    return;
  // assert ((*F)(v).is_integral());

  Qmat M(transpose(basis_matrix)); // rows are field coords of current gens
  M.append_row(v);
  M = HNF(M);           // rows are coords of new Z-basis, last row is 0
  // assert(trivial(M.get_numerator().row(rank+1)));
  M.delete_row();
  M = LLL(M);

  basis_matrix = transpose(M);
  Qmat A = basis_matrix.inverse();
  // assert (is_one(A.get_denom()));
  power_coords_matrix = A.get_numerator();
  index = abs(power_coords_matrix.determinant());
}

// Check that the Z-span of the current Z-basis is closed under
// multiplication.  If not, a will hold a missing product.
int Order::check_order(FieldElement& a) const
{
  for (int i=1; i<=rank; i++)
    {
      const FieldElement& x = basis_elt(i);
      for (int j=i; j<=rank; j++)
        {
          a = x * basis_elt(j);
          if (!contains(a))
            return 0;
        }
    }
  // if we reach here, all products are in the Z-span
  return 1;
}

// Same, v will hold field coords of a missing product
int Order::check_order(Qvec& v) const
{
  for (int i=1; i<=rank; i++)
    {
      const FieldElement& x = basis_elt(i);
      // assert (x.is_integral());
      for (int j=i; j<=rank; j++)
        {
          FieldElement y = basis_elt(j);
          // assert (y.is_integral());
          FieldElement z = x*y;
          v = z.coords(); // field coords
          if (!contains(z)) // not in equation order
            return 0;
        }
    }
  // if we reach here, all products are in the Z-span
  return 1;
}

// Extend by a (an algebraic integer), returning the index of the
// extension: incremental version using the previous two internal
// functions.
ZZ Order::extend_by_one(const FieldElement& a, int check)
{
  if (contains(a))
    {
      return ZZ(1);
    }
  if (check && !a.is_integral())
    {
      cerr << "extend_by(): a is not integral" << endl;
      cerr << "char poly is " << ::str(a.charpoly()) << endl;
      return ZZ(1);
    }
  return extend_by_one(a.coords());
}

// Extend by v, field coords of an algebraic integer, returning the
// index of the extension.
ZZ Order::extend_by_one(const Qvec& v)
{
  if (is_one(v.denom))
    return ZZ(1);
  auto a = (*F)(v);
  // assert (a.is_integral());
  if (contains(a))
    return ZZ(1);

  ZZ old_index = index;
  int nsteps=0;
  Qvec w = v;
  while (nsteps < rank*rank)
    {
      nsteps += 1;
      // cout << "nsteps = " << nsteps << ", index = " << index << endl;
      add_one(w);
      if (check_order(w))
        {
          disc = poldisc/(index*index);
          // cout << "extend_by() returning after " << nsteps << " steps " << endl;
          return index/old_index;
        }
    }
  // code should never reach here
  cout << "extend_by() failed after " << nsteps << " steps, current index is " << index << endl;
  disc = poldisc/(index*index);
  return index/old_index;
}

// Extend by all a in alist (which must be algebraic integers),
// returning the index of the extension
ZZ Order::extend_by(const vector<FieldElement>& alist, int check)
{
  if (check &&
      !all_of(alist.begin(), alist.end(), [](const FieldElement&a){return a.is_integral();}))
    {
      cout << "Cannot extend order " << *this << " by " << alist
           << " which are not all algebraic integers" << endl;
      return ZZ(1);
    }
  ZZ index_gain(1);
  for (auto a: alist)
    index_gain *= extend_by_one(a, 0);
  return index_gain;
}

// Extend by all v in vlist (which must be field coords of algebraic
// integers), returning the index of the extension.
ZZ Order::extend_by(const vector<Qvec>& vlist)
{
  ZZ index_gain(1);
  for (auto v: vlist)
    index_gain *= extend_by_one(v);
  return index_gain;
}

// Extend by all columns v of M (as in previous method), returning
// the index of the extension
ZZ Order::extend_by(const Qmat& M)
{
  ZZ index_gain(1);
  for (int i=1; i<=M.ncols(); i++)
    index_gain *= extend_by_one(M.col(i));
  return index_gain;
}

// LLL-reduce basis (using the basis matrix to reduce) If the basis
// matrix is a signed permutation matrix then (1) the order is the
// equation order, and (2) the basis is a permutation of the power
// basis (up to sign), in which we do nothing.
void Order::LLL_reduce()
{
  mat_m U;
  LLL_reduce(U);
}

// Same, returning the unimodular basis change matrix
void Order::LLL_reduce(mat_m& U)
{
  Qmat M = LLL(transpose(basis_matrix), U);
  if (is_signed_permutation_matrix(M))
    {
      U = transpose(M.get_numerator()) * U;  //transpose=inverse for a signed perm. matrix
      M = Qmat::identity(rank);
    }
  basis_matrix = transpose(M);
  power_coords_matrix = basis_matrix.inverse().get_numerator();
  // index and disc are unchanged
}

// LLL-reduce basis (using the coord matrix of alist to reduce)
void Order::LLL_reduce(const vector<FieldElement>& alist)
{
  vector<FieldElement> alist1(alist);
  alist1.push_back((*F)(1));
  const mat_m& C = (power_coords_matrix * coord_matrix(alist1)).numerator;
  LLL_reduce(C);
}

// Same, returning the unimodular basis change matrix
void Order::LLL_reduce(const vector<FieldElement>& alist, mat_m& U)
{
  vector<FieldElement> alist1(alist);
  alist1.push_back((*F)(1));
  const mat_m& C = (power_coords_matrix * coord_matrix(alist1)).numerator;
  LLL_reduce(C, U);
}

// LLL-reduce basis (using the coord matrix of alist to reduce)
void Order::LLL_reduce(const mat_m& C)
{
  mat_m U;
  LLL_reduce(C, U);
}

// Same, returning the unimodular basis change matrix
void Order::LLL_reduce(const mat_m& C, mat_m& U)
{
  mat_m V;    // to hold inverse unimodular transform
  const mat_m& L = LLL(C, U);  // L = U*C = LLL-reduced coords of alist
  ZZ d = inverse(U, V);
  if (is_one(-d))
    V = -V;

  // Now V*U = id.  We post-multiply basis_matrix by V and premultiply
  // power_coords_matrix by U, preserving the relation basis_matrix *
  // power_coords_matrix = id.  If the new matrices are signed
  // permutation matrices, we replace both with the identity.

  const auto& M = basis_matrix * V;     // potential new basis matrix
  if (is_signed_permutation_matrix(M))
    {
      power_coords_matrix = mat_m::identity_matrix(rank);
      basis_matrix = Qmat(power_coords_matrix);
      d = inverse(power_coords_matrix, U);
    }
  else
    {
      basis_matrix = M;
      power_coords_matrix = U * power_coords_matrix;
    }
  // cout << "After  LLL, coords of alist are:\n" << transpose(L) << endl;
  // index and disc are unchanged
}

// Compute Maximal Order (via lib)pari.  If bound>0 then the order may
// not be p-maximal for p>bound.
Order MaximalOrder(const Field* F, const ZZ& bound)
{
  // cout << "In MaximalOrder(F) with F = " << F << " -> " << *F << endl;
  ZZ ind;
  vector<Qvec> zbc;
  mat_m pcm;
  nfinit(F->minpoly, ind, zbc, pcm, bound);
  // cout << "nfinit() returns ind = " << ind << "\nzbc = " <<zbc << "\npcm = " << pcm
  //      << "\nwith det = " << pcm.determinant() << endl;
  vector<FieldElement> bas(F->d);
  std::transform(zbc.begin(), zbc.end(), bas.begin(),
                 [F](const Qvec& v) {return (*F)(v);});
  // cout << "Now calling Order constructor with basis = " << bas << endl;
  Order O(*F, bas, pcm);
  // cout << "... O = " << O << endl;
  // cout << "End of MaximalOrder(F) with F = " << F << " -> " << *F << endl;
  return O;
}

// List of x^i for i=0,1,...,deg(F)
vector<FieldElement> powers(const FieldElement& x)
{
  const Field& F = *x.field_ptr();
  int d = F.degree();
  vector<FieldElement> bas(d);
  bas[0] = F(1);
  for (int i=1; i<d; i++)
    bas[i] = x*bas[i-1];
  return bas;
}

// To coordinate matrix (by columns unless rows=1).  If rev then
// rows/columns are reversed coordinate vectors.
Qmat coord_matrix(const vector<FieldElement>& alist, int rows, int rev)
{
  Qmat M(alist[0].field_degree(), alist.size());
  int i = 1;
  for (auto& a: alist)
    {
      Qvec v = a.coords();
      if (rev)
        v = reverse(v);
      M.setcol(i++, v);
    }
  return (rows? transpose(M) : M);
}

// From coordinate matrix (by columns unless rows=1).  If rev then
// rows/columns are reversed coordinate vectors.
vector<FieldElement> from_coord_matrix(const Field& F, const Qmat& M, int rows, int rev)
{
  int n = (rows? M.nrows(): M.ncols());
  vector<FieldElement> alist(n);
  for (int i=1; i<=n; i++)
    {
      Qvec v =  (rows? M.row(i): M.col(i));
      if (rev)
        v = reverse(v);
      alist[i-1] = FieldElement(F, v);
    }
  return alist;
}

// For x an algebraic integer, the order spanned by x^i for
// i=0..deg(F)-1.  If check then check that x is an algebraic integer.
Order PowerOrder(const FieldElement& x, int check)
{
  if (check && !x.is_integral())
    {
      cout << "Cannot creat power order from " << x
           << " which is not an algenbraic integer!" << endl;
      return Order();
    }
  else
    return Order(*x.field_ptr(), powers(x));
}

// The old versions add all products of the gens of the two orders (slower):

// ZZ Order::extend_by_old(const FieldElement& a, int check)
// {
//   if (check && !a.is_integral())
//     {
//       cout << "Cannot extend order " << *this << " by " << a
//            << " which is not an algenbraic integer!" << endl;
//       return ZZ(1);
//     }
//   ZZ old_index = index;
//   if (!contains(a))
//     operator+=(a);
//   return index/old_index;
// }

// Extend by all a in alist (which must be algebraic integers),
// returning the index of the extension
// ZZ Order::extend_by_old(const vector<FieldElement>& alist, int check)
// {
//   ZZ index_gain(1);
//   for (auto a: alist)
//     {
//       if (check && !a.is_integral())
//         {
//           cout << "Cannot extend order " << *this << " by " << a
//                << " which is not an algenbraic integer!" << endl;
//           return ZZ(1);
//         }
//       index_gain *= extend_by_old(a, 0); // 0: no need to check a again
//     }
//   return index_gain;
// }
