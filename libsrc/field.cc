// File FIELD.CC: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#include "eclib/field.h"
#include "eclib/polred.h"

//#define DEBUG_ARITH

const Field* FieldQQ = new Field();

//#define DEBUG_FIELD_CONSTRUCTOR
Field::Field(const ZZX& p, string a, int verb)
{
#ifdef DEBUG_FIELD_CONSTRUCTOR
  cout << "In Field constructor with poly " << ::str(p) << ", name = " << a << endl;
  verb = 1;
#endif
  if (IsMonic(p) && IsIrreducible(p))
    {
      // Set A to be the companion matrix of p.  The function
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
      *this = Field(m, ZZ(1), a, verb);
    }
  else
    {
      cerr << "Error: poly should be monic and irreducible" << endl;
#ifdef DEBUG_FIELD_CONSTRUCTOR
      display_factors(p);
#endif
      *this = Field();
    }
}

Field::Field() // defaults to Q
{
  d=1;
  A.init(1,1); // zero matrix
  B.init(1,1); // zero matrix
  C.init(1,1); // zero matrix
  denom = B(1,1) = Bdet = Bdet1 = Bdet2 = Bdet3 = ZZ(1);
  Binv = B;
  var = ""; // will not be used anyway
  SetX(minpoly);
}

Field::Field(const mat_m& m, const ZZ& den, string a, int verb)
  : var(a), d(m.nrows()), denom(den), A(m)
{
  if (verb)
    {
      cout << "----------------------------"<<endl;
      cout << "In Field constructor (var = "<<a<<")"<< endl;
    }
  // NB the assumption is that A/denom is integral, i.e. its (monic)
  // char poly is integral and irreducible.
  minpoly = scaled_charpoly(mat_to_mat_ZZ(A), denom);
  if (verb)
    cout << " - min poly = " << ::str(minpoly) << ", generator " << var << endl;

  // Compute change of basis matrix B, with column j equal to
  // denom^(n-j)*A^(j-1)v for j from 1 to d
  vec_m v(d);
  v[1] = pow(denom,d-1); // so v=[1,0,...,0]*denom^(d-1)
  // cout<<"v = "<<v<<endl;
  B.init(d,d);
  B.setcol(1,v);
  for(int i=2; i<=d; i++)
    {
      v = A*v / denom;
      B.setcol(i,v);
    }
  Bcontent = B.content();
  // cout << "Content of original B is " << Bcontent << endl;
  B /= Bcontent;
  Binv.init(d,d);
  Bdet = inverse(B,Binv); // so B*Binv = Bdet*identity
  mat_m I = mat_m::identity_matrix(d);
  assert (B*Binv == Bdet*I);
  Bdet1 = B(1,1);
  Bdet2 = Binv(1,1);
  Bdet3 = Bdet2*denom;
  // Now we have Binv*A*B = denom*Bdet * C, where C = companion
  // matrix of minpoly. i.e. the B-conjugate of A can be divided by
  // denom to give the integral matrix C.
  C = (Binv*A*B) / (denom*Bdet);
  assert (CompanionMatrix(minpoly) == mat_to_mat_ZZ(C));

  Cpowers.resize(d);
  Cpowers[0] = I;
  for (int i=1; i<d; i++)
    Cpowers[i] = C*Cpowers[i-1];
  if(verb)
    {
      cout<<"basis  matrix = ";
      output_flat_matrix(Binv);
      cout<<endl;
      cout<<"inverse basis matrix = ";
      output_flat_matrix(B);
      cout<<endl;
      cout << "basis factor  = " << Bdet1 << "*" << Bdet2 << "=" << Bdet <<endl;
      if (verb>1)
        {
          cout<<"companion matrix  = ";
          output_flat_matrix(C);
          cout<<endl;
        }
      cout << "Leaving Field constructor" << endl;
      cout << "----------------------------"<<endl;
    }
}

// String for prettier output, like "Q" or "Q(i) = Q[X]/(X^2+1)" or
// raw output, suitable for re-input, like "Q" or "i [1 0 1]":
string Field::str(int raw) const
{
  ostringstream s;
  if (isQ())
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
istream& operator>>(istream& s, Field* F)
{
#ifdef DEBUG_FIELD_INPUT
  cout << "In operator>>(istream& s, Field* F)..." << endl;
#endif
  string var;
  s >> var;
#ifdef DEBUG_FIELD_INPUT
  cout << "- input var = " <<var << endl;
#endif
  if (var=="Q")
    {
#ifdef DEBUG_FIELD_INPUT
      cout << "- setting F to QQ" << endl;
#endif
      *F = Field();
    }
  else
    {
      ZZX f;
      s >> f;
#ifdef DEBUG_FIELD_INPUT
      cout << "- not QQ" <<endl;
      cout << "- input f = " << ::str(f) << endl;
#endif
      *F = Field(f, var);
#ifdef DEBUG_FIELD_INPUT
      cout << *F << endl;
      F->display();
#endif
    }
  return s;
}

istream& operator>>(istream& s, Field** F)
{
  string var;
  s >> var;
  if (var=="Q")
    *F = (Field*)FieldQQ;
  else
    {
      ZZX f;
      s >> f;
      *F = new Field(f, var);
    }
  return s;
}

void Field::display_bases(ostream&s) const
{
  if (isQ())
    return;
  mat_m I(mat_m::identity_matrix(d));

  s << "Powers of A (i.e. powers of " << var << " in A-embedding):\n";
  vector<mat_m> Apowers(d);
  Apowers[0] = I;
  for (int i=1; i<d; i++)
    Apowers[i] = A*Apowers[i-1];
  for (auto Apow: Apowers)
    {
      Apow.output(s);
      s<<endl;
      s<<endl;
    }
  ZZ fac = pow(denom, d-1) * Bdet;
  s << "Basis in A-embedding, scaled by "<< denom << "^" << d-1 << " * " << Bdet << " = " << fac <<":\n";
  fac = Bdet*Bcontent;
  s << "(first columns should be standard basis vectors * "<< fac <<")\n";
  for(int i=1; i<=d; i++)
    {
      vec_m coli = Binv.col(i);
      for(int j=1; j<=d; j++)
        coli[j] *= pow(denom, d-j);
      mat_m M = lin_comb_mats(coli, Apowers);
      M.output(s);
      s<<endl;
      s<<endl;
      if (!(M.col(1)==fac*I.col(i)))
        {
          cout<<"Scaling problem:  column(1) = "<<M.col(1)
              <<" but fac*e_"<<i<<" = "<<fac*I.col(i)<<endl;
        }
    }
  s << "Powers of C (i.e. powers of " << var << " in C-embedding):\n";
  for (auto Cpow: Cpowers)
    {
      Cpow.output(s);
      s<<endl;
      s<<endl;
    }
  s << "Basis in C-embedding, scaled by "<<Bdet<<":\n";
  s << "(first columns should be columns of basis())\n";
  for(int i=1; i<=d; i++)
    {
      mat_m M = lin_comb_mats(Binv.col(i), Cpowers);
      M.output(s);
      s<<endl;
      s<<endl;
      assert (M.col(1) == Binv.col(i));
    }
}

void Field::display(ostream&s, int raw) const
{
  string fpol = ::str(minpoly);
  if (isQ() || (d==1))
    {
      s << "Q" << endl;
      return;
    }
  s << "Q(" << var << ") with defining polynomial "<< ::str(minpoly) <<" (of degree "<<d
    << " and discriminant " << discriminant(minpoly) << ")" << endl;
  if(raw && !isQ())
    {
      s << "   Raw basis with respect to alpha-power basis:\n";
      for(int i=1; i<=d; i++)
        {
          FieldElement bi = element(vec_m::unit_vector(d,i), ZZ(1), 1);
          s << "   #"<<i<<": "<< bi << endl;
          s << "   (with char poly = " << ::str(bi.charpoly())
            << " and minpoly = " << ::str(bi.minpoly()) << "\n";
        }
    }
}

////////////////////////////////////////////////////////////////////////

// raw means the given coords are w.r.t. the B-basis
FieldElement::FieldElement(const Field* HF, const vec_m& c, const ZZ& d, int raw)
    :F(HF), coords(c), denom(d)
{
  // cout<<"Constructing a FieldElement in Q("<<F->var<<")\n";
  if (raw)
    {
      coords = F->Binv *coords;
      denom *= F->Bdet3;
    }
  cancel();
  // cout<<" - finished constructing " << (*this) << "\n";
}

FieldElement Field::element(const vec_m& c, const ZZ& d, int raw) const
{
  return FieldElement(this, c, d, raw);
}

FieldElement Field::rational(const bigrational& x) const
{
  return FieldElement(this, x.num(), x.den());
}

FieldElement Field::rational(const ZZ& x) const
{
  return rational(bigrational(x));
}

FieldElement Field::rational(long x) const
{
  return rational(bigrational(x));
}

FieldElement Field::rational(int x) const
{
  return rational(bigrational(x));
}

FieldElement Field::zero() const
{
  return rational(0);
}

FieldElement Field::one() const
{
  return rational(1);
}

FieldElement Field::minus_one() const
{
  return rational(-1);
}

FieldElement Field::two() const
{
  return rational(2);
}

FieldElement Field::minus_two() const
{
  return rational(-2);
}

FieldElement Field::gen() const
{
  if (isQ())
    return rational(1);
  else
    return FieldElement(this, vec_m::unit_vector(d, 2));
}

void FieldElement::cancel() // divides through by gcd(content(coords, denom))
{
  if (IsOne(denom))
    return;
  if (denom<0)
    {
      denom = -denom;
      coords = -coords;
    }
  ZZ g = gcd(content(coords), denom);
  if (IsOne(g))
    return;
  denom /=g;
  coords /= g;
}

int FieldElement::is_zero() const
{
  if (F->isQ())
    return val.is_zero();
  return trivial(coords);
}

int FieldElement::is_one() const
{
  static const bigrational one(1);
  if (F->isQ())
    return val==one;
  return IsOne(denom) && coords == vec_m::unit_vector(F->d,1);
}

int FieldElement::is_minus_one() const
{
  static const bigrational minus_one(-1);
  if (F->isQ())
    return val==minus_one;
  return IsOne(denom) && coords == -vec_m::unit_vector(F->d,1);
}

void FieldElement::negate() // negate in place
{
  if (F->isQ())
    val = -val;
  else
    coords = -coords;
}

// String for pretty printing, used in default <<
// or for raw output, suitable for re-input (with Field known):
string FieldElement::str(int raw) const
{
  // cout << "In FieldElement::str() with var = " << F->var << endl;
  // cout << "coords = " << coords << ", denom = " << denom << endl;
  ostringstream s;
  if (F->isQ())
    {
      s << val; // if val is an integer this just outputs the integer with no "/" + denominator
      return s.str();
    }
  if (raw)
    s << coords << " " << denom;
  else
    {
      string n = ::str(coords, F->var);
      if (n[0]=='+')
        n.erase(0,1);
      if (denom==1)
        s << n;
      else
        s << "(" << n << ")/" << denom;
    }
  return s.str();
}

// x must be initialised with a Field before input to x
istream& operator>>(istream& s, FieldElement& x)
{
  if (x.field()->isQ())
    s >> x.val;
  else
    s >> x.coords >> x.denom;
  return s;
}

int FieldElement::operator==(const FieldElement& b) const
{
  return in_same_field(b) && (F->isQ()? val==b.val : (denom==b.denom) && (coords==b.coords));
}

int FieldElement::operator!=(const FieldElement& b) const
{
  return (!in_same_field(b)) || (F->isQ()? val!=b.val : (denom!=b.denom) || (coords!=b.coords));
}

mat_m FieldElement::matrix() const // ignores denom, not used for Q
{
  if (F->isQ())
    return mat_m::scalar_matrix(1, num(val));
  return lin_comb_mats(coords, F->Cpowers);
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
      int d = field()->degree() - 1;
      while (d--)
        cp *= mp;
      return cp;
    }

  // Now this is not rational (and the field is not QQ)
  cp = ::charpoly(mat_to_mat_ZZ(matrix()));
  // now replace X by denom*X and make primitive
  ZZ dpow(denom);
  int d = F->d;
  for (int i=1; i<=d; i++)
    {
      SetCoeff(cp, i, dpow*coeff(cp, i));
      if (i!=d)
        dpow *= denom;
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
  int d = field()->degree();
  if (d==1) return val;
  bigrational r;
  if (is_rational(r)) // then norm = r**degree
    return bigrational(pow(r.num(), d), pow(r.den(), d));
  else
    return bigrational(matrix().determinant(), pow(denom, field()->degree()));
}

bigrational FieldElement::trace() const
{
  int d = field()->degree();
  if (d==1) return val;
  bigrational r;
  if (is_rational(r)) // then norm = r**degree
    return ZZ(d) * r;
  else
    return bigrational(matrix().trace(), denom);
}

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
  if (F->isQ())
    {
      val += b.val;
      return;
    }
  coords = b.denom*coords + denom*b.coords;
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator+(const FieldElement& b) const
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      exit(1);
    }
  FieldElement a = *this;
  if (!b.is_zero())
    a += b;
  return a;
}

FieldElement FieldElement::operator-() const
{
  if (F->isQ())
    return FieldElement(-val);
  return FieldElement(F, -coords, denom);
}

// subtract b
void FieldElement::operator-=(const FieldElement& b)
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      exit(1);
    }
  if (b.is_zero())
    return;
  if (F->isQ())
    {
      val -= b.val;
      return;
    }
  coords = b.denom*coords - denom*b.coords;
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator-(const FieldElement& b) const
{
  if (!in_same_field(b))
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
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
  if (b.is_zero()) {*this = b; return;}
  if (F->isQ())
    {
      val *= b.val;
      return;
    }
  coords = (matrix()*b.matrix()).col(1);
  denom *= b.denom;
  cancel();
}

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

FieldElement FieldElement::inverse() const // raise error if zero
{
  if (is_zero())
    {
      cerr << "Attempt to invert zero!" << endl;
      exit(1);
    }
  if (is_one() || is_minus_one())
    return FieldElement(*this);
  if (F->isQ())
    return FieldElement(recip(val));

  mat_m M = matrix(), Minv;
  ZZ Mdet = ::inverse(M,Minv); // so M*Minv = Mdet*identity
  FieldElement ans = FieldElement(F, denom*Minv.col(1), Mdet);
  assert (operator*(ans).is_one());
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
      exit(1);
    }
  if (is_zero())
    return;
  if (F->isQ())
    {
      val /= b.val;
      return;
    }
  operator*=(b.inverse());
}

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
      exit(1);
    }
  FieldElement a = *this;
  if (!is_zero())
    a /= b;
  return a;
}

FieldElement evaluate(const ZZX& f, const FieldElement a)
{
  FieldElement fa(a.F, to_ZZ(0));
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
  if (F->isQ())
    return 1;
  r = bigrational(coords[1],denom);
  for (int i=2; i <= F->d; i++)
    {
      if (coords[i]!=0)
        return 0;
    }
  return 1;
}

// return 1 iff this is an algebraic integer
int FieldElement::is_integral() const
{
  return IsMonic(charpoly());
}

// NB for a in F, either [Q(sqrt(a))=Q(a)] or [Q(sqrt(a)):Q(a)]=2.
// The first function only applies when a has maximal degree:
// return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
int FieldElement::is_absolute_square(FieldElement& r) const
{
  if (F->isQ())
    return (val.is_square(r.val));
  // field not Q, reduce to integral case if necessary
  if (::is_one(denom))
    return is_absolute_integral_square(r);
  FieldElement x = FieldElement(F, denom*coords);
  int res = x.is_absolute_integral_square(r);
  if (res)
    r /= denom;
  return res;
}

// Same as above if the min poly is known and denom=1
int FieldElement::is_absolute_integral_square(FieldElement& r)  const
{
  if (F->isQ())
    return (val.is_square(r.val));
  assert (::is_one(denom));
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

  if (F->isQ())
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
      FieldElement abb = (*this)*b*b, rb(F);
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

void cancel_mat(mat_m& M, ZZ& d)
{
  if (IsOne(d))
    return;
  ZZ g = gcd(M.content(), d);
  if (IsOne(g))
    return;
  M /= g;
  d /= g;
}

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
      denom *= iso.denom;
      cancel_mat(isomat,denom);
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
      denom *= iso.denom;
      cancel_mat(isomat,denom);
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
  ostringstream s;
  if (raw)
    {
      for (auto mij: isomat.get_entries())
        s << mij << " ";
      s << denom;
    }
  else
    {
      if (id_flag)
        s << "Identity automorphism of " << domain->str();
      else
        {
          if (domain==codomain)
            s << "Automorphism of " << domain->str();
          else
            s << "Isomorphism from " << domain->str() << " to " << codomain->str();
          // s << " with matrix\n" << isomat;
          // if (!IsOne(denom)) s << "/ "<< denom;
          s << " mapping " << domain->gen() << " to " << operator()(domain->gen());
        }
    }
  return s.str();
}

// x must be initialised with domain and codomain, this just inputs
// the matrix and denominator
istream& operator>>(istream& s, FieldIso& x)
{
  // cout << "Reading a FieldIso into " << x << endl;
  // cout << "Domain has degree " << x.domain->degree() << endl;
  // cout << "Codomain has degree " << x.codomain->degree() << endl;
  s >> x.isomat >> x.denom;
  x.set_id_flag();
  // cout << "After reading, the FieldIso is " << x << endl;
  // cout << "Matrix = " << x.isomat << ", denominator = " << x.denom << ", id_flag = " << x.id_flag << endl;
  return s;
}

// inverse isomorphism
FieldIso FieldIso::inverse() const
{
  if (id_flag) return *this;
  mat_m inversemat;
  ZZ d = ::inverse(isomat, inversemat);
  inversemat *= denom;
  cancel_mat(inversemat, d);
  return FieldIso(codomain, domain, inversemat, d, 0); // 0: not the identity
}

// map x in domain to an element of the codomain
FieldElement FieldIso::operator()(const FieldElement& x) const
{
  if (id_flag) return x;
  if ((x.field()==domain) || (*x.field()==*domain))
    {
      if (x.field()->isQ())
        return codomain->rational(x.val);

      FieldElement y(codomain, isomat*x.coords, denom*x.denom);
      // sanity check that the min poly has not changed
      if (! (x.minpoly()==y.minpoly()))
        {
          cout << "Error in applying field isomorphism\n" << *this << "\n to x = " << x << "\n --> y = " << y << endl;
          cout << "x has minpoly "<< ::str(x.minpoly()) << endl;
          cout << "y has minpoly "<< ::str(y.minpoly()) << endl;
        }
      return y;
    }
  cerr << "Cannot apply FieldIso\n" << *this << "\n to " << x << " in " << *(x.field()) << endl;
  exit(1);
  return FieldElement(codomain);
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
FieldIso Field::reduction_isomorphism(string newvar, int canonical) const
{
#ifdef DEBUG_REDUCE
  cout << "In Field::reduction_isomorphism(), minpoly = " << ::str(minpoly) << endl;
#endif
  if (isQ())
    {
#ifdef DEBUG_REDUCE
      cout << " - Field is Q, so identity" << endl;
#endif
      return FieldIso(this); // identity
    }
  ZZX h; ZZ denh;
  ZZX g = polred(minpoly, h, denh, canonical);

  if (minpoly==g)
    {
#ifdef DEBUG_REDUCE
      cout << " - " << ::str(minpoly) << " is already reduced, so identity" << endl;
#endif
      return FieldIso(this); // identity
    }
#ifdef DEBUG_REDUCE
  cout << " - reduced minpoly = " << ::str(g) << endl;
#endif
  // construct the reduced field:
  const Field* Fred = new Field(g, newvar);
#ifdef DEBUG_REDUCE
  cout << " - reduced field is\n" << *Fred << endl;
#endif
  // construct the isomorphism matrix from F to Fred:
  mat_m M(d,d);
  // denh * image of F's gen in Fred:
  FieldElement a = evaluate(h,Fred->gen()); // / Fred->rational(denh);
#ifdef DEBUG_REDUCE
  cout << " - image of gen is (" << a << ") / " << denh << endl;
#endif
  FieldElement apow = a; // power of a
  ZZ denhpowmax = pow(denh, d-1);
  ZZ denhpow = denhpowmax;
  M.setcol(1, denhpow * vec_m::unit_vector(d, 1));
  denhpow /= denh;
  M.setcol(2, denhpow * a.get_coords());
  for (int i=3; i<=d; i++)
    {
      denhpow /= denh;
      apow *= a;
      M.setcol(i, denhpow * apow.get_coords());
    }
  cancel_mat(M, denhpowmax);
#ifdef DEBUG_REDUCE
  cout << " - iso matrix = \n" << M;
  if (denh>1)
    cout << " / " << denhpowmax;
  cout << endl;
#endif
  FieldIso iso(this, Fred, M, denhpowmax, 0); // 0: not the identity
  return iso;
}

//#define DEBUG_CHANGE_GEN
// Return an iso from this=Q(a) to Q(b) where B is in this field and generates
FieldIso Field::change_generator(const FieldElement& b) const
{
  FieldIso iso(this); // default

  if (b.field() != this)
    {
      cerr << "Cannot change generator of " << *this << " to " << b
           << " which is in a different field " << *b.field() << endl;
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
  if ((d==1) || (b_pol==minpoly)) // then the identity will do
    {
      return iso;
    }

  // Create the new field:
  const Field* b_field = new Field(b_pol, var+string("1"));

  // To define the map to the new field we need to express a as a
  // polynomial in b.  The coordinates of b^j w.r.t. a are the columns
  // of the matrix M with first column m_1=e_1 and j'th column m_j =
  // (B/dB)*m_{j-1}.  Instead of multiplying M on the right by
  // diag(1,dB,dB^2,...)^-1 (scaling its columns down) we will
  // multiply Minv on the left by diag(1,dB,dB^2,...), (scaling its
  // rows up).

  mat_m bmat(b.matrix()), M(d,d),  Minv(d,d);
  ZZ dbmat(b.denom);
  vec_m v(vec_m::unit_vector(d,1));
  // Set the volumns of M in turn, multiplying by bmat
  M.setcol(1,v);
  for(int j=2; j<=d; j++)
    {
      v = bmat*v;
      M.setcol(j,v);
    }
#ifdef DEBUG_CHANGE_GEN
  cout << "M = " << M << endl;
#endif
  ZZ da = inverse(M,Minv); // so M*Minv = da*identity
#ifdef DEBUG_CHANGE_GEN
  cout << "Before scaling by dbmat = " << dbmat << ", Minv = " << Minv
       << " and denom(Minv) = " << da << endl;
#endif
  // Multiply the rows of Minv by successive posers of dbmat
  if (!IsOne(dbmat))
    {
      ZZ dbmatpow(dbmat);
      for(int i=2; i<=d; i++)
        {
          Minv.multrow(i, dbmatpow);
          if (i<d)
            dbmatpow *= dbmat;
        }
    }
#ifdef DEBUG_CHANGE_GEN
  cout << "After scaling, Minv = " << Minv
       << " and denom(Minv) = " << da << endl;
#endif
  cancel_mat(Minv, da);
#ifdef DEBUG_CHANGE_GEN
  cout << "After cancelling, Minv = " << Minv
       << " and denom(Minv) = " << da << endl;
#endif
  // The coeffs of 1,a,a^2,... as polynomials in b are the columns
  // 1,2,3,... of Minv/da.
  iso =  FieldIso(this, b_field, Minv, da, 0); // 0: not the identity
  // check:
  FieldElement isob = iso(b);
  if (isob!=b_field->gen())
    {
      cerr << "Error in Field::change_generator(b) with b = " << b << "\n";
      cerr << "b has minpoly " << ::str(b.minpoly()) << "\n";
      cerr << "iso(b) = " << isob << " with minpoly " << ::str(isob.minpoly()) << "\n";
      exit(1);
    }
  delete b_field;
  return iso;
}

// Return an iso from this=Q(a) to Q(b) where b^2=r, optionally
// applying polredabs (if reduce=2) or polredbest (if reduce=1) to the
// codomain.  sqrt_r is set to sqrt(r) in the codomain, so sqrt_r^2 =
// image of r.

//#define DEBUG_SQRT_EMBEDDING
FieldIso Field::sqrt_embedding(const FieldElement& r, string newvar, FieldElement& sqrt_r, int reduce) const
{
#ifdef DEBUG_SQRT_EMBEDDING
  cout << "In sqrt_embedding() with base field " << *this << " and r = " << r << endl;
#endif
  if (r.field() != this)
    {
      cerr << "Cannot adjoin sqrt(" << r << ") to " << *this
           << " as it is in a different field " << *r.field() << endl;
      return FieldIso(this);
    }
  if (r.is_square(sqrt_r))
    {
      cout << "Adjoining sqrt(" << r << ") to " << *this << " is trivial since "
           << r << " is already a square, with root " << sqrt_r << endl;
      return FieldIso(this);
    }

  // If r has degree < d we replace it with an equivalent element of
  // maximal degree by multiplying by a square.  Then the sqrt field
  // is generated by f(X^2) where f is the char poly.
  FieldElement s(rational(1));
  FieldElement rss = r;
  if (rss.degree()<d)
    {
      s = gen();
      rss = r*s*s;
      while (rss.degree()<d)
        {
          s += ZZ(1);
          rss = r*s*s;
        }
    }
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " s = " << s << ", rss = " << rss << endl;
#endif
  // Now we adjoin sqrt(rss) instead
  ZZX sqrt_rss_pol = XtoX2(rss.charpoly());
  assert (IsIrreducible(sqrt_rss_pol)); // must be else r is a square
  Field* F_sqrt_rss = new Field(sqrt_rss_pol, newvar);
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " extended field is " << *F_sqrt_rss << endl;
#endif
  // Now we embed this into the new field in three steps:
  // Q(a) -~-> Q(rss) c-> Q(sqrt(rss)) -~-> Q(b)
  // The first and last are isomorphisms, the last (optional) is polredabs reduction.
  FieldIso iso(change_generator(rss));  // Q(a) -> Q(rss)
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " first map (isomorphism) is " << iso << endl;
#endif
  s = iso(s); // image of s in Q(rss)
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " iso(s) = " << s << endl;
#endif
  const Field* Qrss = iso.codomain; // Q(rss)
  // the images of the powers of rss are the even powers of sqrt_rss:
  mat_m isomat(2*d, d);
  for (int j=0; j<d; j++)
    isomat.setrow(2*j+1, vec_m::unit_vector(d,j+1));
  FieldIso emb(Qrss, F_sqrt_rss, isomat, ZZ(1)); // Q(rss) -> Q(sqrt(rss))
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " second map (embedding) is " << emb << endl;
#endif
  s = emb(s); // image of s in Q(sqrt(rss))
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " emb(s) = " << s << endl;
#endif
  sqrt_r = F_sqrt_rss->gen()/s; // sqrt(r) = sqrt(rss)/s
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " sqrt_r = " << sqrt_r << endl;
#endif
  emb.precompose(iso);
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " embedding (before reduction) is " << emb << endl;
#endif
  if (reduce)
    {
#ifdef DEBUG_SQRT_EMBEDDING
      cout << " reducing via "
           << (reduce>1? "polredabs" : "polredbest")
           << "..." << endl;
#endif
      FieldIso red = F_sqrt_rss->reduction_isomorphism(newvar, reduce>1);
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " third map (reduction isomorphism) is" << red << endl;
#endif
  emb.postcompose(red);
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " final embedding is " << emb << endl;
#endif
  sqrt_r = red(sqrt_r);
#ifdef DEBUG_SQRT_EMBEDDING
  cout << " In the extension, sqrt(r) = " << sqrt_r << endl;
#endif
    }
  return emb;
}
