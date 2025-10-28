// polys.h : defines interface to NTL polynomials
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

// allow for multiple includes
#ifndef _ECLIB_POLYS_H
#define _ECLIB_POLYS_H

// NB --warning!  In NTL there is a universal modulus set by a single
// call to ZZ_p::init(q).

#include "gf.h"
#include "bigrat.h"
#include "vector.h"

vector<ZZ> rootsmod(const vector<ZZ>& coeffs, ZZ p);
vector<ZZ_p> roots(const ZZ_pX& f);
vector<bigrational> roots(const vector<ZZ>& coeffs);
vector<bigrational> roots(const ZZX& f);
vector<ZZ> introots(const ZZX& f);

vector<ZZ> Introotscubic(const ZZ& a, const ZZ& b, const ZZ& c);
vector<ZZ> Introotsquartic(const ZZ& a, const ZZ& b, const ZZ& c,
                            const ZZ& d);

// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)
int nrootscubic(const ZZ& bb, const ZZ& cc, const ZZ& dd, const ZZ& p);

ZZ_pX reduce(const ZZX& f, const galois_field& Fq);

inline ZZ_pX reduce(const ZZX& f, const ZZ& q)
{return reduce(f,galois_field(q));}

// function to sort factorizations (lists of (factor,exponent) pairs),
// first by degree of factor then exponent of factor then
// lexicographically

struct factor_comparison {
  bool operator()(pair_ZZX_long& fac1, pair_ZZX_long& fac2)
  {
    // first sort by degree of the factor
    int s = deg(fac1.a) - deg(fac2.a);
    if(s) return (s<0); // true if fac1 has smaller degree

    // then sort by exponent of the factor
    s = fac1.b - fac2.b;
    if(s) return (s<0); // true if fac1 is to a lower exponent

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(fac1.a.rep.begin(), fac1.a.rep.end(),
                                        fac2.a.rep.begin(), fac2.a.rep.end());
  }
};

struct factor_modp_comparison {
  bool operator()(pair_ZZ_pX_long& fac1, pair_ZZ_pX_long& fac2)
  {
    auto cmp = [](const ZZ_p& a, const ZZ_p& b) {return rep(a)<rep(b);};
    // first sort by degree of the factor
    int s = deg(fac1.a) - deg(fac2.a);
    if(s) return (s<0); // true if fac1 has smaller degree

    // then sort by exponent of the factor
    s = fac1.b - fac2.b;
    if(s) return (s<0); // true if fac1 is to a lower exponent

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(fac1.a.rep.begin(), fac1.a.rep.end(),
                                        fac2.a.rep.begin(), fac2.a.rep.end(),
                                        cmp);
  }
};

// function to sort of polynomials, first by degree of factor
// then lexicographically

struct poly_comparison {
  bool operator()(ZZX& pol1, ZZX& pol2)
  {
    // first sort by degree of the factor
    int s = deg(pol1) - deg(pol2);
    if(s) return (s<0); // true if pol1 has smaller degree

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(pol1.rep.begin(), pol1.rep.end(),
                                        pol2.rep.begin(), pol2.rep.end());
  }
};

struct poly_modp_comparison {
  bool operator()(ZZ_pX& pol1, ZZ_pX& pol2)
  {
    auto cmp = [](const ZZ_p& a, const ZZ_p& b) {return rep(a)<rep(b);};
    // first sort by degree of the factor
    int s = deg(pol1) - deg(pol2);
    if(s) return (s<0); // true if pol1 has smaller degree

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(pol1.rep.begin(), pol1.rep.end(),
                                        pol2.rep.begin(), pol2.rep.end(),
                                        cmp);
  }
};

extern factor_comparison fact_cmp;
extern poly_comparison poly_cmp;
extern factor_modp_comparison fact_modp_cmp;
extern poly_modp_comparison poly_modp_cmp;

vector<ZZ> coeffs(const ZZX& p);
vector<ZZ> coeffs(const ZZ_pX& p);
string polynomial_string(const vector<ZZ>& coeffs, const string& var="X");
string polynomial_string(const Zvec<ZZ>& coeffs, const string& var="X");
string polynomial_string(const ZZX& p, const string& var="X");
string polynomial_string(const ZZ_pX& p, const string& var="X");

// factor a primitive (e.g. monic) polynomial
NTL::vec_pair_ZZX_long factor(const ZZX& f);

// display factors of a polynomial:
void display_factors(const ZZX& f);

// display factors of a polynomial mod p:
void display_factors(const ZZ_pX& f);

// return f(X/c)*c^d
ZZX scale_poly_up(const ZZX& f, const ZZ& c);
// return f(c*X)/c^d
ZZX scale_poly_down(const ZZX& f, const ZZ& c);

// return f(X) mod m (or just f if m==0)
ZZX reduce_poly(const ZZX& f, const ZZ& m);

// Coprime test
int AreCoprime(const ZZX& f, const ZZX& g);

// Monic test
inline int IsMonic(const ZZX& f)
{ return IsOne(LeadCoeff(f));}

// Irreducibility test (ignoring content)
int IsIrreducible(const ZZX& f);

// Squarefree test
int IsSquareFree(const ZZX& f);

// return f(X^2)
ZZX XtoX2(const ZZX& f);

// write f(x) = f0(X^2)+X*f1(X^2)
void parity_split(const ZZX& f, ZZX& f0, ZZX& f1);

// assuming f irreducible:
// return 0 if f(x^2) is irreducible; else
// return 1 and set g where f(x^2)=g(x)g(-x) (*-1 if degree odd)
int is_square(const ZZX& f, ZZX& g);

#endif // #define _POLYS_
