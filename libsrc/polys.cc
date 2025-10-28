// polys.cc : implements interface to NTL polynomials
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
 
#include "eclib/polys.h"

ZZ_pX reduce(const ZZX& f, const galois_field& Fq)
{
  ZZ_pX fmodq;
  for(int i=0; i<=deg(f); i++)
    SetCoeff(fmodq,i,ZtoGF(Fq,coeff(f,i)));
  return fmodq;
}

vector<gf_element> roots(const ZZ_pX& f)
{
  // make f monic:
  ZZ_pX f1=f;
  MakeMonic(f1);
  // reduce to distinct roots case:
  ZZ_pX X; SetX(X); 
  ZZ_pX g = PowerXMod(ZZ_p::modulus(),f1)-X;
  vec_ZZ_p r; FindRoots(r,GCD(f1,g)); 
  vector<gf_element>ans;
  for(int i=0; i<r.length(); i++) ans.push_back(r[i]);
  return ans;
}

vector<ZZ> rootsmod(const vector<ZZ>& coeffs, ZZ q)
{
  galois_field Fq(q);
  ZZ_pX f;
  unsigned long i, deg = coeffs.size()-1;
  for (i=0; i<=deg; i++) SetCoeff(f,i,ZtoGF(Fq,coeffs[i]));

  vector<gf_element> r = roots(f);
  vector<ZZ>ans;
  for(i=0; i<r.size(); i++) ans.push_back(LiftGF(r[i]));

  sort(ans.begin(),ans.end());
  return ans;
}

//#define TRACE_ROOTS
vector<bigrational> roots(const ZZX& f)
{
#ifdef TRACE_ROOTS
  cout<<"Finding rational roots of polynomial f =  "<<f<<endl;
#endif
  vector<bigrational> ans;
  int i;
  ZZX g;
  bigrational root;
  ZZ c;
  vec_pair_ZZX_long factors;
  factor(c,factors,f);

#ifdef TRACE_ROOTS
  cout<<"f has " << factors.length() << " factors" << endl;
#endif

  for(i=0; i<factors.length(); i++)
    {
      g = factors[i].a;
#ifdef TRACE_ROOTS
      cout<<"factor "<<g<<" has degree "<<deg(g)<<endl;
#endif
      if(deg(g)==1)
        {
          root = bigrational(-coeff(g,0),coeff(g,1));
#ifdef TRACE_ROOTS
          cout<<"root "<<root<<endl;
#endif
          ans.push_back(root);
        }
    }
  sort(ans.begin(), ans.end());
  return ans;
}

// Return the list of *integral* roots of an integral polynomial.
// Intended for monic polys, but that is not a requirement
vector<ZZ> introots(const ZZX& f)
{
#ifdef TRACE_ROOTS
  cout<<"Finding integer roots of polynomial f =  "<<f<<endl;
#endif
  vector<bigrational> ratroots = roots(f);
  vector<ZZ> ans;
  if (ratroots.size()==0)
    return ans;
  for( const auto& r : ratroots)
    if (den(r)==1)
      ans.push_back(num(r));
  sort(ans.begin(), ans.end());
  return ans;
}

vector<bigrational> roots(const vector<ZZ>& coeffs)
{
#ifdef TRACE_ROOTS
  cout<<"Finding rational roots of polynomial f with coefficients "<<coeffs<<endl;
#endif
  ZZX f;
  vector<bigrational> ans;
  int i, d = coeffs.size()-1;  // degree
  if(d<1)
    return ans;
  for(i=0; i<=d; i++)
    SetCoeff(f,d-i,coeffs[i]);
#ifdef TRACE_ROOTS
  cout<<"f = "<<f<<endl;
#endif
  ans = roots(f);
#ifdef TRACE_ROOTS
  cout<<"roots of f: "<< ans << endl;
#endif
  return ans;
}



// root-finding functions for monic integer cubics and quartics
//
// With NTL we factor the polynomial in Z[X] and pick out degree 1 factors

vector<ZZ> Introotscubic(const ZZ& a, const ZZ& b, const ZZ& c)
{
  ZZX f;
  SetCoeff(f,3);   // sets it to 1
  SetCoeff(f,2,a);
  SetCoeff(f,1,b);
  SetCoeff(f,0,c);
  return introots(f);
}

vector<ZZ> Introotsquartic(const ZZ& a, const ZZ& b, const ZZ& c, const ZZ& d)
{
  ZZX f; vec_pair_ZZX_long factors; ZZ cont;
  SetCoeff(f,4);   // sets it to 1
  SetCoeff(f,3,a);
  SetCoeff(f,2,b);
  SetCoeff(f,1,c);
  SetCoeff(f,0,d);
  return introots(f);
}

// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)
int nrootscubic(const ZZ& b, const ZZ& c, const ZZ& d, const ZZ& p)
{
  vector<ZZ> coeffs;
  coeffs.push_back(d);
  coeffs.push_back(c);
  coeffs.push_back(b);
  coeffs.push_back(ZZ(1));
  return rootsmod(coeffs,p).size();
}

// factor a primitive (e.g. monic) polynomial
vec_pair_ZZX_long factor(const ZZX& f)
{
  vec_pair_ZZX_long factors;
  ZZ cont;
  factor(cont,factors,f);
  ::sort(factors.begin(), factors.end(), fact_cmp);
  return factors;
}

factor_comparison fact_cmp;
poly_comparison poly_cmp;
factor_modp_comparison fact_modp_cmp;
poly_modp_comparison poly_modp_cmp;

// pretty output for integer polynomials
string monomial_string(int i, const string& var)
{
  ostringstream s;
  if (i>0) s << var;
  if (i>1) s << "^" << i;
  return s.str();
}

string polynomial_string(const vector<ZZ>& coeffs, const string& var)
{
  //  cout<<"\npolynomial_string("<<coeffs<<")\n";
  if (std::all_of(coeffs.begin(), coeffs.end(), [](const ZZ&c){return IsZero(c);}))
    return "0";
  int d = coeffs.size()-1;
  ZZ c;
  ostringstream s;
  if (d==0)
    {
      s << coeffs[0];
      return s.str();
    }
  // All non-constant terms:
  for (int i=d; i>0; i--)
    {
      c = coeffs[i];
      if (c==0)
        continue;
      if (c>1)
        {
          if (i<d) // no + needed on leading term
            s << "+";
          s << c << "*";
        }
      else if (c==1 && i<d) s << "+";
      else if (c==-1) s << "-";
      else if (c<-1) s << "-" << abs(c) << "*";
      s << monomial_string(i, var);
      //cout<<" - after i="<<i<<": "<<s.str()<<endl;
    }
  // Constant term:
  c = coeffs[0];
  if (c>0) s << "+" << c;
  else if (c<0) s << "-" <<abs(c);
  //cout<<" - after i=0: "<<s.str()<<endl;

  return s.str();
}

string polynomial_string(const Zvec<ZZ> coeffs, const string& var)
{
  if (trivial(coeffs))
    return "0";
  int d = dim(coeffs); // one less than 'degree'
  vector<ZZ> co(d);
  for(int i=0; i<d; i++)
    co[i] = coeffs[i+1];
  return polynomial_string(co, var);
}

vector<ZZ> coeffs(const ZZX& p)
{
  int d = deg(p);
  vector<ZZ> v(d+1);
  for(int i=0; i<=d; i++)
    v[i] = coeff(p, i);
  return v;
}

vector<ZZ> coeffs(const ZZ_pX& p)
{
  int d = deg(p);
  vector<ZZ> v(d+1);
  for(int i=0; i<=d; i++)
    v[i] = rep(coeff(p, i));
  return v;
}

string polynomial_string(const ZZX& p, const string& var)
{
  return polynomial_string(coeffs(p), var);
}

string polynomial_string(const ZZ_pX& p, const string& var)
{
  return polynomial_string(coeffs(p), var);
}

// display factors of a polynomaial:
void display_factor(const pair_ZZX_long& f)
{
  ZZX p = f.a;
  string pol = polynomial_string(p);
  int d = deg(p), e = f.b;
  cout << "(degree " << d << ")\t"
       << pol
       << "\t to power " << e;
  //cout << " (coefficients " << p << ")";
}

void display_factors(const ZZX& f)
{
  ZZ cont; vec_pair_ZZX_long factors;
  factor(cont, factors, f);
  ::sort(factors.begin(), factors.end(), fact_cmp);
  long nf = factors.length();
  for(int i=0; i<nf; i++)
    {
      cout << (i+1) << ":\t";
      display_factor(factors[i]);
      cout<<endl;
    }
}

void display_factor(const pair_ZZ_pX_long& f)
{
  ZZ_pX p = f.a;
  string pol = polynomial_string(p);
  int d = deg(p), e = f.b;
  cout << "(degree " << d << ")\t"
       << pol
       << "\t to power " << e;
  //cout << " (coefficients " << p << ")";
}

void display_factors(const ZZ_pX& f)
{
  vec_pair_ZZ_pX_long factors = berlekamp(f);
  ::sort(factors.begin(), factors.end(), fact_modp_cmp);
  long nf = factors.length();
  for(int i=0; i<nf; i++)
    {
      cout << (i+1) << ":\t";
      display_factor(factors[i]);
      cout<<endl;
    }
}

// return f(X/c)*c^d: multiply coeff(f,i) by c^(d-i)
ZZX scale_poly_up(const ZZX& f, const ZZ& c)
{
  if (c==1) return f;
  ZZX g = f;
  ZZ cpow(1);
  long d = deg(f);
  for (int i=0; i<=d; i++)
    {
      SetCoeff(g, d-i, cpow*coeff(g, d-i));
      if (i<d)
        cpow *= c;
    }
  return g;
}

// return f(c*X)/c^d: divide coeff(f,i) by c^(d-i)
// NB only use when divisions are exact
ZZX scale_poly_down(const ZZX& f, const ZZ& c)
{
  if (c==1) return f;
  ZZX g = f;
  ZZ cpow(1);
  long d = deg(f);
  for (int i=0; i<=d; i++)
    {
      SetCoeff(g, d-i, coeff(g, d-i)/cpow);
      if (i<d)
        cpow *= c;
    }
  return g;
}

// return f(X) mod m
ZZX reduce_poly(const ZZX& f, const ZZ& m)
{
  if (m==0) return f;
  ZZX g = f;
  long d = deg(f);
  for (int i=0; i<=d; i++)
    SetCoeff(g, d-i, mod(coeff(g, d-i), m));
  return g;
}

// Coprime test
int AreCoprime(const ZZX& f, const ZZX& g)
{
  return deg(GCD(f, g)) == 0;
}

// Squarefree test
int IsSquareFree(const ZZX& f)
{
  return AreCoprime(f, diff(f));
}

// Irreducibility test (ignoring content)
int IsIrreducible(const ZZX& f)
{
  ZZ cont; vec_pair_ZZX_long factors;
  factor(cont, factors, f);
  return factors.length()==1;
}

// return f(X^2)
ZZX XtoX2(const ZZX& f)
{
  int n = deg(f);
  ZZX f2(2*n, LeadCoeff(f));
  for (int i=0; i<n; i++)
    SetCoeff(f2, 2*i, coeff(f,i));
  return f2;
}

// write f(x) = f0(X^2)+X*f1(X^2)
void parity_split(const ZZX& f, ZZX& f0, ZZX& f1)
{
  int n = deg(f);
  // if n is even, deg(f0)=n/2 and deg(f1)<=(n-2)/2
  // if n is odd, deg(f0)<=(n-1)/2 and deg(f1)=(n-1)/2
  // so both have degree <=n/2
  f0.SetLength(1 + n/2); // reserve enough space
  f1.SetLength(1 + n/2); // reserve enough space
  for (int i=0; i<=n; i++)
    {
      if (i%2) // i odd
        SetCoeff(f1, (i-1)/2, coeff(f,i));
      else // i even
        SetCoeff(f0, i/2, coeff(f,i));
    }
  f0.normalize(); // strip any leading 0s
  f1.normalize(); // strip any leading 0s
  if (! (f==XtoX2(f0)+LeftShift(XtoX2(f1),1)))
    cout << "parity_split of " << polynomial_string(f) << " gives \n"
         << "f0 = "<< polynomial_string(f0) << "\n"
         << "f1 = "<< polynomial_string(f1) << "\n"
         << " --> " << polynomial_string(XtoX2(f0)+LeftShift(XtoX2(f1),1)) << endl;
}

// assuming f irreducible:
// return 0 if f(x^2) is irreducible; else
// return 1 and set g where f(x^2)=g(x)g(-x) (*-1 if degree odd)
int is_square(const ZZX& f, ZZX& g)
{
  vec_pair_ZZX_long factors = factor(XtoX2(f));
  if (factors.length()==1)
    return 0;
  g = factors[0].a;
  return 1;
}
