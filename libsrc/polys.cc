// polys.cc : implements interface to NTL polynomials
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
 
#include <eclib/polys.h>

FqPoly reduce(const ZPoly& f, const galois_field& Fq)
{
  NewFqPoly(Fq,fmodq);
  SetDegree(fmodq,Degree(f));
  for(int i=0; i<=Degree(f); i++)
    SetCoeff(fmodq,i,ZtoGF(Fq,PolyCoeff(f,i)));
  return fmodq;
}

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

//#define TRACE_ROOTS
vector<bigrational> roots(const ZPoly& f)
{
#ifdef TRACE_ROOTS
  cout<<"Finding rational roots of polynomial f =  "<<f<<endl;
#endif

  vector<bigrational> ans;
  int i;
  ZPoly g;
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

vector<bigrational> roots(const vector<bigint>& coeffs)
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
