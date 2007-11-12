// gf.h:  common interface for LiDIA's galois_field/gf_element & NTL's ZZ_p
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
 
// allow for multiple includes
#ifndef _GF_H_
#define _GF_H_

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)

#include "LiDIA/gf_element.h"

#define NewGF(field,name) gf_element name(field)
#define GFinit(field,name) name.assign_zero(field)
#define GFSetZ(name,Zvalue) name.assign(Zvalue)
#define LiftGF(name) name.lift_to_Z()
#define GFrandomize(name) name.randomize()
#define IsZero(name) name.is_zero()

inline gf_element ZtoGF(const galois_field& F, const bigint& a)
{
  gf_element b(F); b.assign(a); return b;
}

inline gf_element ItoGF(const galois_field& F, int a)
{
  bigint aa; aa.assign(a);
  gf_element b(F); b.assign(aa); return b;
}

inline gf_element root_of_unity(const galois_field& F, int n)
{
  gf_element mu(F); 
  bigint qm1 = F.number_of_elements()-1;
  if(!((qm1%n).is_zero())) return mu; // =0
  qm1/=n;
  mu.assign_primitive_element(F);
  power(mu,mu,qm1);
  return mu;
}

inline gf_element sqrt(const galois_field& F, const gf_element& a)
{
  return sqrt(a);
}

// Returns 1 if a is a square (root in r), else 0
inline int sqrt(const galois_field& F, const gf_element& a, gf_element& r)
{
  if ( a.is_square() ) {r=sqrt(a); return 1;} 
  return 0;
}

inline bigint order(const gf_element& z)
{
  return z.order();
}

#else // NTL

#include <NTL/ZZ_p.h>

class galois_field {
  ZZ q;
 public:
  galois_field(void) :q(to_ZZ(2))  {ZZ_p::init(q);} //dummy
  galois_field(const ZZ& qq) :q(qq) {ZZ_p::init(q);}
  bigint characteristic() const {return q;}
};

#define gf_element ZZ_p

#define NewGF(field,name) gf_element name
#define GFinit(field,name) name=to_ZZ_p(0)
#define GFSetZ(name,Zvalue) name=to_ZZ_p(Zvalue)
#define LiftGF(name) rep(name)
#define GFrandomize(name) random(name)

inline gf_element ZtoGF(const galois_field& F, const bigint& a)
{
  return to_ZZ_p(a);
}

inline gf_element ItoGF(const galois_field& F, int a)
{
  return to_ZZ_p(a);
}

// NB Here caller must ensure that a is a square;  q odd
inline gf_element sqrt(const galois_field& F, const gf_element& a)
{
  bigint rd;  
  ressol(rd,rep(a),F.characteristic());
  return ZtoGF(F,rd);
}

// Returns 1 if a is a square (root in r), else 0
inline int sqrt(const galois_field& F, const gf_element& a, gf_element& r)
{
  bigint rd = to_ZZ(0),  repa=rep(a),  q = F.characteristic();  
  switch(legendre(repa,q))
    {
    case -1: return 0;
    case 1:
      ressol(rd,repa,q); // & carry through to next lines
    case 0: 
      r= ZtoGF(F,rd);
    }
  return 1;
}

inline gf_element root_of_unity(const galois_field& F, int n)
{
  ZZ qm1 = F.characteristic()-1;
  if(qm1%n !=0) return to_ZZ_p(0);
  qm1/=n;
  while(1)
    {
      ZZ_p mu = random_ZZ_p(); 
      if(mu==to_ZZ_p(0)) continue;
      power(mu,mu,qm1);
      if(mu!=to_ZZ_p(1))  return mu;
    }
}

inline bigint order(const gf_element& z)
{
  gf_element one=z/z;
  gf_element zn=z;  bigint n=BIGINT(1);
  while (zn!=one) {zn*=z; n+=1;}
  return n;
}

#endif

inline vector<gf_element> roots_of_unity(const galois_field& Fq, int p)
{
  gf_element mu = root_of_unity(Fq,p); // =0 if p ndiv q-1
  vector<gf_element>mu_p;
  mu_p.resize(p);
  GFinit(Fq,mu_p[0]);
  GFSetZ(mu_p[0],1);
  for(int i=1; i<p; i++) mu_p[i]=mu_p[i-1]*mu;
  return mu_p;
}

#endif
