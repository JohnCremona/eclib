// polys.cc : implements uniform interface to LiDIA/NTL polynomials
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include "marith.h"
#include "polys.h"

FqPoly reduce(const ZPoly& f, const galois_field& Fq)
{
  NewFqPoly(Fq,fmodq);
  SetDegree(fmodq,Degree(f));
  for(int i=0; i<=Degree(f); i++)
    SetCoeff(fmodq,i,ZtoGF(Fq,PolyCoeff(f,i)));
  return fmodq;
}

#ifdef LiDIA_INTS

vector<gf_element> roots(const FqPoly& f)
{
  // make f monic:
  FqPoly f1;
  divide(f1,f,f.lead_coeff());
  // reduce to distinct roots case:
  galois_field F=f.get_field();
  FqPoly g,h; g.assign_x(F);         //X
  power(h,g,F.number_of_elements()); //X^q
  h=h-g;                             //X^q-X
  gcd(f1,f1,h);
  base_vector<gf_element> r = find_roots(f1);
  int i, nr = r.size();  
  vector<gf_element>ans;
  ans.reserve(nr);
  for(i=0; i<nr; i++) ans.push_back(r[i]);
  return ans;
}

#endif // LiDIA_INTS

#ifdef NTL_INTS

vector<gf_element> roots(const FqPoly& f)
{
  // make f monic:
  FqPoly f1=f;
  MakeMonic(f1);
  // reduce to distinct roots case:
  ZZ_pX X; SetX(X); 
  ZZ_pX g = PowerXMod(ZZ_p::modulus(),f1)-X;
  vec_ZZ_p r; FindRoots(r,GCD(f1,g)); 
  vector<gf_element>ans;
  for(int i=0; i<r.length(); i++) ans.push_back(r[i]);
  return ans;
}

/*
vector<bigint> rootsmod(const vector<bigint>& coeffs, bigint q)
{
  ZZ_p::init(q);
  ZZ_pX f; ZZ_p ci;
  for (long i=0; i<coeffs.size(); i++) 
    { ci=to_ZZ_p(coeffs[i]); SetCoeff(f,i,ci); }

  vector<gf_element> r = roots(f);
  int i, nr = r.length();

  vector<bigint> ans;
  ans.reserve(nr);
  for(int i=0; i<nr; i++) ans.push_back(rep(r[i]));
  sort(ans.begin(),ans.end());
  return ans;
}
*/

#endif // NTL_INTS

vector<bigint> rootsmod(const vector<bigint>& coeffs, bigint q)
{
  galois_field Fq(q);
  NewFqPoly(Fq,f);
  unsigned long i, deg = coeffs.size()-1;
  SetDegree(f,deg);
  for (i=0; i<=deg; i++) SetCoeff(f,i,ZtoGF(Fq,coeffs[i]));

  vector<gf_element> r = roots(f);
  vector<bigint>ans;
  for(i=0; i<r.size(); i++) ans.push_back(LiftGF(r[i]));

  sort(ans.begin(),ans.end());
  return ans;
}

// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)
int nrootscubic(const bigint& b, const bigint& c, const bigint& d, const bigint& p)
{
  vector<bigint> coeffs;
  coeffs.push_back(d);
  coeffs.push_back(c);
  coeffs.push_back(b);
  coeffs.push_back(BIGINT(1));
  return rootsmod(coeffs,p).size();
}
