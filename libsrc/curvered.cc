// curvered.cc -- implementation of CurveRed class etc.
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

#include <eclib/curve.h>
#include <eclib/polys.h>    // for nrootscubic
#include <eclib/ffmod.h>
#include <eclib/parifact.h> // for ellap

ostream& operator<<(ostream& os, const Kodaira_code& c)
{
  int code=c.code;
  switch (code%10){
  case 0:
    os<<"I"<<((code)/10); break;
  case 1:
    os<<"I*"<<((code - 1)/10);  break;
  case 2:
    os<<"II   "; break;
  case 3:
    os<<"III  "; break;
  case 4:
    os<<"IV   "; break;
  case 5:
    os<<"IV*  "; break;
  case 6:
    os<<"III* "; break;
  case 7:
    os<<"II*  "; break;
  default:
    os<<"???? "; break;
  };
  return os;
}

//
// Tate's algorithm -- constructor function for a CurveRed
//
// subroutines -- not general purpose (.cc only)

bigint root(const bigint& aa, int e, const bigint & p)
// the e'th root of aa, mod p
{
  bigint a = aa%p, r, b, xe(e);
  if (is_zero(a)) return a;
  if (e==2)
    {
      sqrt_mod_p(r, a, p);
      return r;
    }
  for (r = 1; r<p; ++r)
    {
      power_mod(b, r, xe, p);
      if (b==a) return r;
    }
  return bigint(0);
}

// test if quadratic aX^2 + bX + c = 0 (mod p) has roots
int rootsexist(const bigint& aa, const bigint& bb, const bigint& cc, const bigint& p)
{
  if (even(p)) return is_zero((aa*bb*cc)%p);
  const bigint& a = aa % p;
  if (is_zero(a)) return 1;
  const bigint& b = bb % p;
  const bigint& c = cc % p;
  if (is_zero(c)) return 1;
  return legendre(b*b - 4*a*c, p) >= 0;      //ie true if legendre(d,p)=0,1
}

//monic version
int rootsexist(const bigint& bb, const bigint& cc, const bigint& p)
{
  static const bigint one(1);
  return rootsexist(one,bb,cc,p);
}

CurveRed::~CurveRed()
{
}

CurveRed::CurveRed(const CurveRed& E)    : Curvedata(E), N(E.N)
{
  factor_discr(); // will only do anything if not already factored
  reduct_array = E.reduct_array;
}

void CurveRed::operator=(const CurveRed& E)
{
  if(this==&E) return;
  a1=E.a1;a2=E.a2;a3=E.a3;a4=E.a4;a6=E.a6;b2=E.b2;b4=E.b4;b6=E.b6;b8=E.b8;
  c4=E.c4;c6=E.c6;discr=E.discr;minimal_flag=1;
  discr_factored=E.discr_factored;
  the_bad_primes=E.the_bad_primes;
  conncomp=E.conncomp;
  ntorsion=E.ntorsion;
  reduct_array = E.reduct_array;
  N=E.N;
}

CurveRed::CurveRed(const Curvedata& E)
  : Curvedata(E, 1), //minimalize in constructor
    N(1)
{
  static const bigint one(1), three(3), twelve(12);
  // constructor stuff
  if (discr==0) {N = 0; return; }
  factor_discr(); // will only do anything if not already factored

  // local variables
  Curvedata C(*this);
  bigint temp, r, s, t, b, c, bb, cc, bc, d, w, x, mx, my,
          a2t, a3t, a4t, a6t, zero;
  int c_p=1, sw, loop, ix, iy;
  zero=0;
  // main loop - for each of the prime divisors of the discriminant.
  // Because the curve is minimal, Tate's algorithm reduce-loop is not needed

//cout<<"Running Tate's algorithm"<<endl;

  for (const auto& p : the_bad_primes)
    {
    int ord_p_discr = val(p,discr);
    int ord_p_j = ord_p_discr - 3*val(p,c4);
    if (ord_p_j < 0) ord_p_j = 0;
    bigint halfmodp = (p+1) >>1;
    int pdiv2 = even(p);
    int pdiv3 = (p==3);

    //change coords so that p|C.a3,C.a4,C.a6
    if ( pdiv2 )
      { if ( div(p,C.b2) )
          { r = root(C.a4,2,p);
            t = root(((r+C.a2)*r+C.a4)*r+C.a6,2,p); }
        else { temp=invmod(C.a1,p);
               r = temp*C.a3;
               t = temp*(C.a4 + r*r); }
      }
    else if ( pdiv3 )
      { if ( div(p,C.b2) ) r = root(-C.b6,3,p); else r = -invmod(C.b2,p)*C.b4;
        t  =  C.a1*r + C.a3; }
    else
      { if ( div(p,c4) ) r = -invmod(twelve,p)*C.b2;
        else r = -invmod(12*c4,p)*(c6+C.b2*c4);
        t =  -halfmodp*(C.a1*r+C.a3); }
    r = mod(r,p);
    t = mod(t,p);
    C.transform(r,zero,t);

    // test for Types In, II, III, IV
    if ( ndiv(p,c4) )
      {temp=-C.a2;
       if (rootsexist(C.a1,temp,p) ) c_p = ord_p_discr;
       else { if ( odd(ord_p_discr) ) c_p = 1; else c_p = 2; }
       reduct_array[p] = Reduction_type
         (ord_p_discr, 1, ord_p_j, 10*ord_p_discr, c_p);
       continue; }  // Type In (n=ord_p_discr)
    else if ( val(p,C.a6) < 2 )
      {reduct_array[p] = Reduction_type
         (ord_p_discr, ord_p_discr, ord_p_j, 2, 1);
       continue; } // Type II
    else if ( val(p,C.b8) < 3 )
      {reduct_array[p] = Reduction_type
         (ord_p_discr, ord_p_discr - 1, ord_p_j, 3, 2);
       continue; } // Type III
    else if ( val(p,C.b6) < 3 )
      {temp = -(C.a6/p)/p;
      bigint temp2 = C.a3/p;
      if ( rootsexist(temp2,temp,p) ) c_p = 3;
      else c_p = 1;
      reduct_array[p] = Reduction_type
	(ord_p_discr, ord_p_discr - 2, ord_p_j, 4, c_p);
      continue; } // Type IV

    // else change coords so that p|C.a1,C.a2, p^2|C.a3,C.a4, p^3|C.a6
    if ( pdiv2 )
      { s = root(C.a2,2,p);
        t = p*root((C.a6/p)/p,2,p); }
    else if ( pdiv3 )
      { s = C.a1;
        t = C.a3; }
    else
      { s = -C.a1*halfmodp;
        t = -C.a3*halfmodp; }
    C.transform(zero,s,t);

    //                             3     2
    // Analyse roots of the cubic T  + bT  + cT + d = 0, where
    // b=C.a2/p, c=(C.a4/p)/p, d=((C.a6/p)/p)/p
    b=C.a2/p;
    c=(C.a4/p)/p;
    d=((C.a6/p)/p)/p;
    bb=b*b; cc=c*c; bc=b*c;
    w = 27*d*d - bb*cc + 4*b*bb*d - 18*bc*d + 4*c*cc;
    x = 3*c - bb;

    if ( div(p,w) )
      {if ( div(p,x) ) sw = 3; else sw = 2; }
    else sw = 1;

//cout << "Analysing roots of cubic; case " << sw << endl;

    switch ( sw ) {
    case 1:
      //Three distinct roots - Type I*0
      reduct_array[p] = Reduction_type
        (ord_p_discr, ord_p_discr - 4, ord_p_j, 1, 1+nrootscubic(b,c,d,p) );
      break;

    case 2:
      // One double root - Type I*m for some m
      // Change coords so that the double root is T=0 mod p
      if ( pdiv2 ) r = root(c,2,p);
      else if ( pdiv3 ) r = c*invmod(b,p);
      else r = (bc - 9*d)*invmod(2*x,p);
      r = p * mod(r,p);
      C.transform(r,zero,zero);

      ix = 3; iy = 3; mx = p*p; my = p*p;
      loop = 1;
      while (loop)
      {
        a2t = C.a2/p;
        a3t = C.a3/my;
        a4t = (C.a4/p)/mx;
        a6t = (C.a6/mx)/my;
	temp = a3t*a3t + 4*a6t;
        if ( div(p,temp ) )
          {if ( pdiv2 ) t = my*root(a6t,2,p);
           else t = my*mod(-a3t*halfmodp, p);
           C.transform(zero,zero,t);
           my = my*p;
           iy++;
           a2t = C.a2/p;
           a3t = C.a3/my;
           a4t = (C.a4/p)/mx;
           a6t = (C.a6/mx)/my;
	   temp = a4t*a4t - 4*a6t*a2t;
           if ( div(p,temp) )
             {if ( pdiv2 ) r = mx*root( a6t*invmod(a2t,p), 2, p);
              else r = mx*mod( -a4t*invmod(2*a2t,p), p);
              C.transform(r,zero,zero);
              mx = mx*p;
              ix++;       // and stay in loop
             }
           else
             {if ( rootsexist(a2t,a4t,a6t,p) ) c_p = 4;
              else c_p = 2;
              loop = 0; }  // and exit loop
          }
        else
          { temp = -a6t;
	    if ( rootsexist(a3t,temp,p) ) c_p = 4;
           else c_p = 2;
           loop = 0; }
      }
      reduct_array[p] = Reduction_type
        (ord_p_discr, ord_p_discr - ix - iy + 1, ord_p_j,
         10 * (ix + iy) - 49, c_p );
      break;  // Type I*m

    case 3:
      // Triple root
      // change coords so that T=0 mod p
      if ( pdiv2 ) r = b;
      else if ( pdiv3 ) r = root(-d,3,p);
      else r = -b*invmod(three,p);
      r = p*mod(r,p);
      C.transform(r,zero,zero);

      a3t = (C.a3/p)/p;
      a6t = (((C.a6/p)/p)/p)/p;

      // test for Type IV*
      temp = a3t*a3t + 4*a6t;
      if ( ndiv(p,temp ) )
        {
	  temp = -a6t;
	  if ( rootsexist(a3t,temp,p) ) c_p = 3;
	  else c_p = 1;
	  reduct_array[p] = Reduction_type
	    (ord_p_discr, ord_p_discr - 6, ord_p_j, 5, c_p);
	  break; }

      // change coords so that p^3|C.a3, p^5|C.a6
      if ( pdiv2 ) t = -p*p*root(a6t,2,p);
      else t = p*p*mod(-a3t*halfmodp, p);
      C.transform(zero,zero,t);

      // test for types III*, II*
      if ( val(p,C.a4) < 4 )
        {reduct_array[p] = Reduction_type
           (ord_p_discr, ord_p_discr - 7, ord_p_j, 6, 2);
         break; }  // Type III*
      else if ( val(p,C.a6) < 6 )
        {reduct_array[p] = Reduction_type
           (ord_p_discr, ord_p_discr - 8, ord_p_j, 7, 1);
         break; } // Type II*

      else
        cout<<" ## Tate's algorithm reached end of loop !!!"<<endl;

      // at this point (only if the input curve were not minimal)
      // one would divide each ai by p^i and start again
    };  // end switch
  }     // end primes for-loop
  N = one;
  for( const auto& ri : reduct_array)
    N *= pow((ri.first), (ri.second).ord_p_N);
  for (auto p: the_bad_primes)
    setLocalRootNumber(p);
  return;
}     // end of Tate's algorithm

// Sort key for sorting lists of curves (LMFDB ordering):
// (1) conductor
// (2) list of ap for p < NP_SORT
// (3) a1,a2,a3,4,a6
vector<bigint> CurveRed::sort_key(const int NP_SORT) const
{
  vector<bigint> key;
  key.push_back(N);
  for(primevar pr(NP_SORT); pr.ok(); ++pr)
    key.push_back(bigint(ap(pr)));
  key.push_back(a1);
  key.push_back(a2);
  key.push_back(a3);
  key.push_back(a4);
  key.push_back(a6);
  return key;
}


// The local Tamagawa exponent -- same as Tamagawa number unless the
// component group is (2,2).  Use p=0 for reals
bigint local_Tamagawa_exponent(const CurveRed& c, const bigint& p)
{
  static const bigint one(1), two(2), four(4);
  if (is_zero(p)) return bigint(c.conncomp);
  auto ri = c.reduct_array.find(p);
  if (ri == c.reduct_array.end())
    return one;
  Reduction_type info = ri->second;
  int cp = info.c_p;
  if (cp!=4)
    return bigint(cp);
  // see if we have C4 or C2xC2
  int code = info.Kcode.code;
  return (code%20==1? two: four); // Type I*m, m even: [2,2], else 4
}

// The global Tamagawa exponent, i.e. the lcm of the exponents of
// the component groups at all bad primes (including infinity if
// real_too is 1), which is the lcm of the local Tamagawa exponents.
// So (with no further knowledge of the MW group) we know that m*P
// is in the good-reduction subgroup for all P, with this m.

bigint global_Tamagawa_exponent(const CurveRed& c, int real_too)
{
  static const bigint one(1);
  static const bigint two(2);
  bigint ans = ((real_too && (getconncomp(c)==2))? two: one);

  for( const auto&  ri : c.reduct_array)
    {
      Reduction_type info = ri.second;
      int code = info.Kcode.code;
      int ep = (code%20==1? 2: info.c_p); // Type I*m, m even: [2,2]
      ans = lcm(ans,bigint(ep));
    }
  return ans;
}

// Tamagawa primes: primes dividing any Tamagawa number
vector<long> tamagawa_primes(const CurveRed& C, int real_too)
{
  vector<bigint> T = pdivs(global_Tamagawa_exponent(C, real_too));
  vector<long> t(T.size());
  std::transform(T.begin(), T.end(), t.begin(), [](const bigint& x) {return I2long(x);});
  return t;
}

// CurveRed member access friend functions:

// NB If p is not a bad prime this will not cause an error, but will
// return the appropriate result for a good prime.  We do have to
// check that p is a valid key in the reduct_array map first,
// otherwise we would cause inclusion of a new entry in that map,
// which is impossible without removing the "const" qualifier from the
// CurveRed argument!

int CurveRed::ord_p_discr(const bigint& p) const
{
  auto ri = reduct_array.find(p);
  return (ri==reduct_array.end()? 0 : (ri->second).ord_p_discr);
}

int CurveRed::ord_p_N(const bigint& p) const
{
  auto ri = reduct_array.find(p);
  return (ri==reduct_array.end()? 0 : (ri->second).ord_p_N);
}

int CurveRed::ord_p_j_denom(const bigint& p) const
{
  auto ri = reduct_array.find(p);
  return (ri==reduct_array.end()? 0 : (ri->second).ord_p_j_denom);
}

int CurveRed::c_p(const bigint& p) const
{
  auto ri = reduct_array.find(p);
  return (ri==reduct_array.end()? 1 : (ri->second).c_p);
}

vector<bigint> CurveRed::all_cp() const
{
  vector<bigint> allcp(reduct_array.size());
  auto cp = [] (const pair<bigint,Reduction_type>& x) {return bigint(x.second.c_p);};
  std::transform(reduct_array.begin(), reduct_array.end(), allcp.begin(), cp);
  return allcp;
}

bigint CurveRed::prodcp() const
{
  static const bigint one(1);
  vector<bigint> allcp = all_cp();
  return std::accumulate(allcp.begin(), allcp.end(), one, std::multiplies<>());
}

// The local Tamagawa number.  Use p=0 for reals
bigint local_Tamagawa_number(const CurveRed& c, const bigint& p)
{
  return bigint(is_zero(p)? getconncomp(c): getc_p(c,p));
}

// The global Tamagawa number, = product of local ones.
bigint global_Tamagawa_number(const CurveRed& c, int real_too)
{
  return bigint(prodcp(c) * (real_too ? getconncomp(c) : 1));
}

Kodaira_code getKodaira_code(const CurveRed& c, const bigint& p)
{
  auto ri = c.reduct_array.find(p);
  return (ri==c.reduct_array.end()? Kodaira_code() : (ri->second).Kcode);
}

void CurveRed::output(ostream& os) const
{
  Curvedata::output(os);
  if(isnull()) return;
  os << "Conductor = " << N << endl;
}

ostream& operator<<(ostream& os, const Reduction_type& R)
{
  os << R.ord_p_discr << "\t";
  os << R.ord_p_N << "\t";
  os << R.ord_p_j_denom <<"\t";
  os << R.Kcode << "\t";
  os << R.c_p << "\t";
  os << R.local_root_number;
  return os;
}

void CurveRed::display(ostream& os)
{
  CurveRed::output(os);
  if(isnull()) return;
  os << "Global Root Number = " << GlobalRootNumber() << endl;
  os << "Reduction type at bad primes:\n";
  os <<"p\tord(d)\tord(N)\tord(j)\tKodaira\tc_p\troot_number\n";
  for( const auto& ri : reduct_array)
    {
      if((ri.second).local_root_number==0)
	setLocalRootNumber(ri.first);
      os << (ri.first) << "\t" << (ri.second) << endl;
    }
}

// Sign of functional equation for an elliptic curve over Q
//
// Taken partly from Odile Lecacheux's GP code, partly translated by
// JEC from the pari C code in elliptic.c from pari-2.1.3
//
// For the theory, see
//
// Halberstadt, Emmanuel. Signes locaux des courbes elliptiques en 2 et 3.
// (French) [Local root numbers of elliptic curves for $p=2$ or $3$]
// C. R. Acad. Sci. Paris S�r. I Math. 326 (1998), no. 9, 1047--1052.
//

// The following methods return local and global root numbers, just
// looking up the local numbers from the
// Reduction_type::local_root_number field, computing them if not
// already set (i.e. field contains 0)

int CurveRed::LocalRootNumber(const bigint& p) const
{
  if(is_zero(p)) return -1;  // the infinite prime
  auto ri = reduct_array.find(p);
  if(ri==reduct_array.end()) return 1; // good reduction case
  return (ri->second).local_root_number;
}

int CurveRed::GlobalRootNumber() const
{
  auto rn = [](const std::pair<bigint,Reduction_type>& ri){return (ri.second).local_root_number;};
  return std::transform_reduce(reduct_array.cbegin(), reduct_array.cend(),
                               -1, std::multiplies<>(), rn);
}

int kro(const bigint& d, const bigint& n);
int kro(const bigint& d, long n);
int kro(long d, long n);
int kro_m1(long x); // kronecker(-1,x) with x>0 odd
int kro_p2(long x); // kronecker(2,x) with x>0 odd
int kro_m2(long x); // kronecker(-2,x) with x>0 odd
int kro_3(long x); // kronecker(x,3)

// sets the "local root number" or local factor at 2 to +/-1

void CurveRed::setLocalRootNumber2()
{
  static const bigint two(2);
  auto ri = reduct_array.find(two);
  if(ri==reduct_array.end()) return;
  Reduction_type& info = ri->second;
  int kod = PariKodairaCode(info.Kcode);
  int n2  = neron(2,kod);

  bigint mu,mv; long u,v; int v4,v6;

  if (is_zero(c4)) {v4=12; u=0;}
  else {mu=c4; v4=divide_out(mu,two); u = posmod(mu,64);}

  if (is_zero(c6)) {v6=12; v=0;}
  else {mv=c6; v6=divide_out(mv,two); v = posmod(mv,64);}

  if (kod > 4)
    {
      info.local_root_number = div(2,a2+a3)? -1: 1;
      return;
    }

  if (kod < -9)
    {
      info.local_root_number = (n2==2)? -kro_m1(v) : -1;
      return;
    }

  bigint tmp = discr;
  divide_out(tmp,two);
  long d1=posmod(tmp,64);

  long x1=u+v+v, y1;

  switch(kod)
  {
  case 1: info.local_root_number = 1;
    return;
  case 2:
    switch(n2)
      {
      case 1:
	switch(v4)
	  {
	  case 4: info.local_root_number = kro_m1(u);
	    return;
	  case 5: info.local_root_number =  1;
	    return;
	  default: info.local_root_number =  -1;
	    return;
	  }
      case 2: info.local_root_number = (v6==7) ? 1 : -1;
	return;
      case 3: info.local_root_number = (v%8==5 || (u*v)%8==5) ? 1 : -1;
	return;
      case 4: if (v4>5)
	{
	  info.local_root_number = kro_m1(v);
	  return;
	}
      else
	{
	  info.local_root_number = (v4==5) ? -kro_m1(u) : -1;
	  return;
	}
      }
  case 3:
    switch(n2)
      {
      case 1: info.local_root_number = -kro_p2(u*v);
	return;
      case 2: info.local_root_number = -kro_p2(v);
	return;
      case 3: y1=posmod((u-(c6 >> 5)) , 16);
	info.local_root_number = (y1==7 || y1==11) ? 1 : -1;
	return;
      case 4: info.local_root_number = (v%8==3 || (2*u+v)%8==7) ? 1 : -1;
	return;
      case 5: info.local_root_number = v6==8 ? kro_p2(x1) : kro_m2(u);
	return;
      }
  case -1:
    switch(n2)
      {
      case 1: info.local_root_number = -kro_p2(x1);
	return;
      case 2: info.local_root_number = (v%8==7) || (x1%32==11) ? 1 : -1;
	return;
      case 3: info.local_root_number = v4==6 ? 1 : -1;
	return;
      case 4: if (v4>6)
	{
	  info.local_root_number = kro_m1(v);
	  return;
	}
      else
	{
	  info.local_root_number = v4==6 ? -kro_m1(u*v) : -1;
	  return;
	}
      }
  case -2: info.local_root_number = n2==1 ? kro_m2(v) : kro_m1(v);
    return;
  case -3:
    switch(n2)
      {
      case 1: y1=posmod((u-2*v),64);
	info.local_root_number = (y1==3) || (y1==19) ? 1 : -1;
	return;
      case 2:
	if(kro_m1(u)==1)
	  {
	    info.local_root_number = kro_p2(v);
	    return;
	  }
	else
	  {
	    info.local_root_number = kro_m2(v);
	    return;
	  }
      case 3:
	if(kro_m1(u)==1)
	  {
	    info.local_root_number = -kro_m2(u*v);
	    return;
	  }
	else
	  {
	    info.local_root_number = kro_p2(u*v);
	    return;
	  }
      case 4: info.local_root_number = v6==11 ? kro_m2(x1) : -kro_m2(u);
	return;
      }
  case -5:
    if (n2==1)
      {
	info.local_root_number = x1%32==23 ? 1 : -1;
	return;
      }
    else
      {
	info.local_root_number = -kro_p2(2*u+v);
	return;
      }
  case -6:
    switch(n2)
      {
      case 1: info.local_root_number = 1;
	return;
      case 2: info.local_root_number = v6==10 ? 1 : -1;
	return;
      case 3: info.local_root_number = (u%16==11) || ((u+4*v)%16==3) ? 1 : -1;
	return;
      }
  case -7:
    if (n2==1)
      {
	info.local_root_number = 1;
	return;
      }
    else
      {
        y1= posmod((u+(c6 >> 8)) , 16);
	if (v6==10)
	  {
	    info.local_root_number = (y1==9) || (y1==13) ? 1 : -1;
	    return;
	  }
	else
	  {
	    info.local_root_number = (y1==9) || (y1==5) ? 1 : -1;
	    return;
	  }
      }
  case -8: info.local_root_number = n2==2 ? kro_m1(v*d1) : -1;
    return;
  case -9: info.local_root_number = n2==2 ? -kro_m1(d1) : -1;
    return;
  default: info.local_root_number = -1;
    return;
  }
}

// sets the "local root number" or local factor at 3 to +/-1

void CurveRed::setLocalRootNumber3()
{
  static const bigint three(3);
  auto ri = reduct_array.find(three);
  if(ri==reduct_array.end()) return;
  Reduction_type& info = ri->second;
  int kod = PariKodairaCode(info.Kcode);
  int n2  = neron(3,kod);

  bigint mu,mv; long u,v; int v4;

  if (is_zero(c4)) {
    v4=12;
    u=0;
  }
  else {
    mu=c4;
    v4=divide_out(mu,three);
    u = posmod(mu,81);
  }

  if (is_zero(c6)) {
    v=0;
  }
  else {
    mv=c6;
    divide_out(mv,three);
    v = posmod(mv,81);
  }

  bigint tmp = discr;
  divide_out(tmp,three);
  long d1=posmod(tmp,81);
  long r6 = posmod(v,9);
  long K4=kro_3(u), K6=kro_3(v);

  if (kod > 4)
    {
      info.local_root_number = K6;
      return;
    }

  switch(kod)
    {
    case 1: case 3: case -3: info.local_root_number = 1;
      return;
    case 2:
      switch(n2)
	{
	case 1: info.local_root_number = (r6==4 || r6>6) ? 1 : -1;
	  return;
	case 2: info.local_root_number = -K4*K6;
	  return;
	case 3: info.local_root_number = 1;
	  return;
	case 4: info.local_root_number = -K6;
	  return;
	}
    case 4:
      switch(n2)
	{
	case 1: info.local_root_number = K6*kro_3(d1);
	  return;
	case 2: info.local_root_number = -K4;
	  return;
	case 3: info.local_root_number = -K6;
	  return;
	}
    case -2: info.local_root_number = n2==2 ? 1 : K6;
      return;
    case -4:
      switch(n2)
	{
	case 1:
	  if (v4==4)
	    {
	      info.local_root_number = (r6==4 || r6==8) ? 1 : -1;
	      return;
	    }
	  else
	    {
	      info.local_root_number = (r6==1 || r6==2) ? 1 : -1;
	      return;
	    }
	case 2: info.local_root_number = -K6;
	  return;
	case 3: info.local_root_number = (r6==2 || r6==7) ? 1 : -1;
	  return;
	case 4: info.local_root_number = K6;
	  return;
	}
    default: info.local_root_number = -1;
      return;
    }
}

// Given a prime p not 2 or 3, sets to +1 or -1, the "local root
// number" or local factor in the sign of the functional equation of
// L(E,s).

void CurveRed::setLocalRootNumber_not_2_or_3(const bigint& p)
{
  auto ri = reduct_array.find(p);
  if(ri==reduct_array.end()) return;
  Reduction_type& info = ri->second;
  if (info.ord_p_N == 1)
    {
      info.local_root_number = -kro(-c6,p);
      return;
    }

  long sp=posmod(p,24);
  if (info.ord_p_j_denom >0)
    {
      info.local_root_number = kro_m1(sp);
      return;
    }

  long ep=12 / gcd(12,info.ord_p_discr);
  if(ep==4)
    {
      info.local_root_number =  kro_m2(sp);
      return;
    }
  if(odd(ep))
    {
      info.local_root_number =  kro_3(sp);
      return;
    }
  info.local_root_number =  kro_m1(sp);
}


// Given a prime p, sets to +1 or -1, the "local root number" or local
// factor in the sign of the functional equation of L(E,s).
//
// This function just delegates to subsidiary ones for the cases
// p=2, p=3, and p>=5.
//

void CurveRed::setLocalRootNumber(const bigint& p)
{
  if (is_zero(p)) return;
  if (p==2) setLocalRootNumber2();
  else if (p==3) setLocalRootNumber3();
  else  setLocalRootNumber_not_2_or_3(p);
}


int kro(const bigint& d, const bigint& n)
{
  return kronecker(d,n);
}

int kro(const bigint& d, long n)
{
  return kronecker(d,n);
}

int kro(long d, long n)
{
  return kronecker(d,n);
}

int kro_m1(long x) // kronecker(-1,x) with x>0 odd
{
  static const int kro_m1_tab[4] = {0,1,0,-1};
  return kro_m1_tab[x&3];
}

int kro_p2(long x) // kronecker(2,x) with x>0 odd
{
  static const int kro_p2_tab[8] = {0,1,0,-1,0,-1,0,1};
  return kro_p2_tab[x&7];
}

int kro_m2(long x) // kronecker(-2,x) with x>0 odd
{
  static const int kro_m2_tab[8] = {0,1,0,1,0,-1,0,-1};
  return kro_m2_tab[x&7];
}

int kro_3(long x) // kronecker(x,3)
{
  static const int kro_3_tab[3] = {0,1,-1};
  return kro_3_tab[x%3];
}


int PariKodairaCode(Kodaira_code Kod)
{
  int ans, code=Kod.code;
  if (code==0) ans = 1;
  else if (code==1) ans = -1;
  else if (code%10 ==0) ans = (code/10)+4;
  else if (code%10 ==1) ans = ((1-code)/10)-4;
  else if (code>4) ans = code-9;
  else ans=code;
#ifdef DEBUG_ESIGN
  cout<<"PariKodairaCode("<<Kod<<") returns "<<ans<<endl;
#endif
  return ans;
}

//  p = 2 or 3 for the neron function

int CurveRed::neron(long p, int kod)
{
  bigint d=discr;
  int v4=val(p,c4);
  int v6=val(p,c6);
  int vd=val(p,d);
#ifdef DEBUG_ESIGN
  cout<<"In neron with p="<<p<<", v4="<<v4<<", v6="<<v6<<", vd="<<vd<<", kod="<<kod<<endl;
#endif
  if (p==3)
  {
    if (abs(kod)>4) return 1;
    else
      {
	switch(kod)
	  {
	  case -1: case 1: return v4&1 ? 2 : 1;
	  case -3: case 3: return (2*v6>vd+3) ? 2 : 1;
	  case -4: case 2:
	    switch (vd%6)
	      {
	      case 4: return 3;
	      case 5: return 4;
	      default: return v6%3==1 ? 2 : 1;
	      }
	  default: /* kod = -2 et 4 */
	    switch (vd%6)
	      {
	      case 0: return 2;
	      case 1: return 3;
	      default: return 1;
	      }
	  }
      }
  }
  if(p==2)
    {
      if (kod>4) return 1;
      else
	{
	  switch(kod)
	    {
	    case 1: return (v6>0) ? 2 : 1;
	    case 2:
	      if (vd==4) return 1;
	      else
		{
		  if (vd==7) return 3;
		  else return v4==4 ? 2 : 4;
		}
	    case 3:
	      switch(vd)
		{
		case 6: return 3;
		case 8: return 4;
		case 9: return 5;
		default: return v4==5 ? 2 : 1;
		}
	    case 4: return v4>4 ? 2 : 1;
	    case -1:
	      switch(vd)
		{
		case 9: return 2;
		case 10: return 4;
		default: return v4>4 ? 3 : 1;
		}
	    case -2:
	      switch(vd)
		{
		case 12: return 2;
		case 14: return 3;
		default: return 1;
		}
	    case -3:
	      switch(vd)
		{
		case 12: return 2;
		case 14: return 3;
		case 15: return 4;
		default: return 1;
		}
	    case -4: return v6==7 ? 2 : 1;
	    case -5: return (v6==7 || v4==6) ? 2 : 1;
	    case -6:
	      switch(vd)
		{
		case 12: return 2;
		case 13: return 3;
		default: return v4==6 ? 2 : 1;
		}
	    case -7: return (vd==12 || v4==6) ? 2 : 1;
	    default:
	      return v4==6 ? 2 : 1;
	    }
	}
    }
  cout<<"neron() returns 0 -- should not happen!"<<endl;
  return 0; /* should not occur */
}

// Trace of Frobenius (via pari) if p good
// (or 0 for additive reduction, +1 for split multiplicative, -1 for nonsplit)
long CurveRed::ap(long p) const
{
  bigint P(p);
  int f = min(2, ord_p_N(P));

  switch (f) {
    // Bad primes: return the p'th coefficient of the L-series
  case 2: return 0;
  case 1: return -LocalRootNumber(P);
  case 0: default: // good primes
    // NB p must be good here else the ai mod p may define a singular curve
    return ellap(posmod(a1,p), posmod(a2,p), posmod(a3,p), posmod(a4,p), posmod(a6,p), p);
  }
}

// Trace of Frobenius (via pari) if p good
// (or 0 for additive reduction, +1 for split multiplicative, -1 for nonsplit)
bigint CurveRed::ap(const bigint& p) const
{
  if (is_long(p))
    return bigint(ap(I2long(p)));

  // now p does not fit in a long. until we implement conversion from
  // NTL ZZ to pari integers we use our own point-counting

  static const bigint zero(0);
  static const bigint one(1);
  int f = min(2, ord_p_N(p));
  switch (f) {
    // Bad primes: return the p'th coefficient of the L-series
  case 2: return zero;
  case 1: return bigint(-LocalRootNumber(p));
  case 0: default: // good primes
    return one + p - reduce_curve(*this,p).group_order();
  }
}

// Quadratic twist of an elliptic curve (returns minimal model)
CurveRed QuadraticTwist(const CurveRed& E, const bigint& D)
{
  static bigint zero(0);
  bigint c4, c6, D2=D*D;
  E.getci(c4, c6);
  Curvedata ED(zero, zero, zero,-27*D2*c4,-54*D*D2*c6, 1); // 1 means minimise
  return CurveRed(ED);
}

// Given a list of elliptic curves E, and one discriminant D, return the
// list of twists of the curves by D
vector<CurveRed> QuadraticTwists(const vector<CurveRed>& EE, const bigint& D)
{
  vector<CurveRed> ans;
  std::transform(EE.begin(), EE.end(), std::back_inserter(ans),
                 [D] (const CurveRed& E) {return QuadraticTwist(E,D);});
  return ans;
}

// Given a list of elliptic curves E, and one prime p, return the
// list of twists of the curves by:
// +p if p=1 (mod 4)
// -p if p=3 (mod 4)
// -4, 8 and -8 if p=2

vector<CurveRed> PrimeTwists(const vector<CurveRed>& EE, const bigint& p)
{
  long p4 = posmod(p,4);
  if (p4==1)
    return QuadraticTwists(EE,p);
  if (p4==3)
    return QuadraticTwists(EE,-p);
  vector<CurveRed> ans;
  static const vector<bigint> D2 = {bigint(-4), bigint(-8), bigint(8)};
  for (auto D: D2)
    {
      vector<CurveRed> ans1 = QuadraticTwists(EE,D);
      ans.insert(ans.end(), ans1.begin(), ans1.end());
    }
  return ans;
}


// Given a list of elliptic curves, and a list of primes, return a
// list of all quadratic twists of the curves by discriminants supported on
// those primes (including the original curves)

vector<CurveRed> AllTwists(const vector<CurveRed>& EE, const vector<bigint>& PP)
{
  vector<CurveRed> ans = EE;
  for (auto p: PP)
    {
      vector<CurveRed> p_twists = PrimeTwists(ans, p);
      for (auto E: p_twists)
        {
          if (std::find(ans.begin(), ans.end(), E) == ans.end())
            ans.push_back(E);
        }
    }
  return ans;
}
