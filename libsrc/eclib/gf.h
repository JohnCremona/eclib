// gf.h:  interface for NTL's ZZ_p
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
#ifndef _ECLIB_GF_H_
#define _ECLIB_GF_H_

#include "marith.h"

extern map<ZZ,ZZ_pContext> ZZ_pContextCache;

class galois_field {
  ZZ q;
 public:
  galois_field(void);
  explicit galois_field(const ZZ& qq);
  ZZ characteristic() const {return q;}
};

// NB Here caller must ensure that a is a square;  q odd
inline ZZ_p sqrt(const galois_field& F, const ZZ_p& a)
{
  ZZ rd;
  sqrt_mod_p(rd,rep(a),F.characteristic());
  return to_ZZ_p(rd);
}

// Returns 1 if a is a square (root in r), else 0
inline int sqrt(const galois_field& F, const ZZ_p& a, ZZ_p& r)
{
  ZZ rd = to_ZZ(0),  repa=rep(a),  q = F.characteristic();  
  switch(legendre(repa,q))
    {
    case -1: return 0;
    case 1:
      sqrt_mod_p(rd,repa,q); // & carry through to next lines
    case 0: 
      r= to_ZZ_p(rd);
    }
  return 1;
}

inline ZZ_p root_of_unity(const galois_field& F, int n)
{
  ZZ qm1 = F.characteristic()-1;
  if(qm1%n !=0) return to_ZZ_p(0);
  qm1/=n;
  while(1)
    {
      ZZ_p mu = NTL::random_ZZ_p(); 
      if(mu==to_ZZ_p(0)) continue;
      power(mu,mu,qm1);
      if(mu!=to_ZZ_p(1))  return mu;
    }
}

inline ZZ order(const ZZ_p& z)
{
  ZZ_p one=to_ZZ_p(1);
  ZZ_p zn=z;  ZZ n(1);
  while (zn!=one) {zn*=z; n+=1;}
  return n;
}

inline vector<ZZ_p> roots_of_unity(const galois_field& Fq, int p)
{
  ZZ_p mu = root_of_unity(Fq,p); // =0 if p ndiv q-1
  vector<ZZ_p>mu_p;
  mu_p.resize(p);
  mu_p[0] = to_ZZ_p(1);
  for(int i=1; i<p; i++) mu_p[i]=mu_p[i-1]*mu;
  return mu_p;
}

#endif // #define _GF_H_
