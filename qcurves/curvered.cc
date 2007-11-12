// curvered.cc -- implementation of CurveRed class etc.
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
 
#include "curve.h"
#include "polys.h" // for nrootscubic
#include "curvemod.h" // for point counting to get ap
#include "ffmod.h"

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
  bigint i, ans, temp; int found=0;
  const bigint& a = mod(aa, p);
  for (i = 1; (! found); ++i)
    {ans = i;
     if (e==2) temp = ans*ans - a;
     else temp = ans*ans*ans - a;
     found = div(p, temp); }
//   bigint ans;
//   int found = 0;
//   bigint a = aa % p;
//   for ( ans = 1; (ans<=p) && (! found); ++ans)
//     found = div(p, pow(ans,e) - a);     //N.B. pow(0,p) causes problems.
  return ans;
}

// test if quadratic aX^2 + bX + c = 0 (mod p) has roots
int rootsexist(const bigint& aa, const bigint& bb, const bigint& cc, const bigint& p)
{
  const bigint& a = aa % p;
  const bigint& b = bb % p;
  const bigint& c = cc % p;
  const bigint& temp = (a*b*c)%p;
  if (even(p)) return sign(temp)==0;
  if (sign(a)==0) return 1;
  const bigint& d = b*b - 4*a*c;
  return legendre(d, p) >= 0;      //ie true if legendre(d,p)=0,1
}

//monic version
int rootsexist(const bigint& bb, const bigint& cc, const bigint& p)
{
  static bigint one; one=1;
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
: Curvedata(E, 1) //minimalize in constructor
{
  // constructor stuff
  N=1;
  if (discr==0) {N = 0; return; }
  factor_discr(); // will only do anything if not already factored

  // local variables
  Curvedata C(*this);
  bigint p, halfmodp, temp, r, s, t, b, c, bb, cc, bc, d, w, x, mx, my,
          a2t, a3t, a4t, a6t, zero;
  int ord_p_discr, ord_p_j, c_p=1, pdiv2, pdiv3, sw, loop, ix, iy;
  zero=0;
  // main loop - for each of the prime divisors of the discriminant.
  // Because the curve is minimal, Tate's algorithm reduce-loop is not needed

//cout<<"Running Tate's algorithm"<<endl;
  vector<bigint>::const_iterator pi= the_bad_primes.begin();
  while(pi!=the_bad_primes.end())
    {
    p = *pi++;
    ord_p_discr = val(p,discr);
    ord_p_j = ord_p_discr - 3*val(p,c4);
    if (ord_p_j < 0) ord_p_j = 0;
    halfmodp = (p+1) >>1;
    pdiv2 = even(p);
    pdiv3 = (p==3);

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
      { if ( div(p,c4) ) r = -invmod(BIGINT(12),p)*C.b2;
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
      else r = -b*invmod(BIGINT(3),p);
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
        cerr<<" ## Tate's algorithm reached end of loop !!!"<<endl;

      // at this point (only if the input curve were not minimal)
      // one would divide each ai by p^i and start again
    };  // end switch
  }     // end primes for-loop
  N = BIGINT(1);
  map<bigint,Reduction_type>::const_iterator ri; 
  for(ri = reduct_array.begin(); ri!=reduct_array.end(); ri++)
    {
      N *= pow((ri->first), (ri->second).ord_p_N);
    }
  return;
}     // end of Tate's algorithm


// CurveRed member access friend functions:

// NB If p is not a bad prime this will not cause an error, but will
// return the appropriate result for a good prime.  We do have to
// check that p is a valid key in the reduct_array map first,
// otherwise we would cause inclusion of a new entry in that map,
// which is impossible without removing the "const" qualifier from the
// CurveRed argument!

int getord_p_discr(const CurveRed& c, const bigint& p)
{
  map<bigint,Reduction_type>::const_iterator ri = c.reduct_array.find(p);
  if(ri==c.reduct_array.end()) return 0;
  return (ri->second).ord_p_discr;
}

int getord_p_N(const CurveRed& c, const bigint& p)
{
  map<bigint,Reduction_type>::const_iterator ri = c.reduct_array.find(p);
  if(ri==c.reduct_array.end()) return 0;
  return (ri->second).ord_p_N;
}

int getord_p_j_denom(const CurveRed& c, const bigint& p)
{
  map<bigint,Reduction_type>::const_iterator ri = c.reduct_array.find(p);
  if(ri==c.reduct_array.end()) return 0;
  return (ri->second).ord_p_j_denom;
}

int getc_p(const CurveRed& c, const bigint& p)
{
  map<bigint,Reduction_type>::const_iterator ri = c.reduct_array.find(p);
  if(ri==c.reduct_array.end()) return 1;
  return (ri->second).c_p;
}

int prodcp(const CurveRed& c)
{
  int ans=1;
  map<bigint,Reduction_type>::const_iterator ri; 
  for(ri = c.reduct_array.begin(); ri!=c.reduct_array.end(); ri++)
    {
      ans *= (ri->second).c_p;
    }
  return ans;
}

Kodaira_code getKodaira_code(const CurveRed& c, const bigint& p)
{
  map<bigint,Reduction_type>::const_iterator ri = c.reduct_array.find(p);
  if(ri==c.reduct_array.end()) return Kodaira_code();
  return (ri->second).Kcode;
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
  os << "Global Root Number = " << GlobalRootNumber(*this) << endl;
  os << "Reduction type at bad primes:\n";
  os <<"p\tord(d)\tord(N)\tord(j)\tKodaira\tc_p\troot_number\n";
  map<bigint,Reduction_type>::const_iterator ri;  
  for(ri = reduct_array.begin(); ri!=reduct_array.end(); ri++)
    {
      if((ri->second).local_root_number==0) 
	setLocalRootNumber(ri->first);     
      os << (ri->first) << "\t" << (ri->second) << endl;
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
// C. R. Acad. Sci. Paris Sér. I Math. 326 (1998), no. 9, 1047--1052.
//

// The following functions return local and global root numbers, just
// looking up the local numbers from the
// Reduction_type::local_root_number field, computing them if not
// already set (i.e. field contains 0)

int LocalRootNumber(CurveRed& c, const bigint& p)
{
  if(is_zero(p)) return -1;  // the infinite prime
  map<bigint,Reduction_type>::const_iterator ri = c.reduct_array.find(p);
  if(ri==c.reduct_array.end()) return 1; // good reduction case
  if((ri->second).local_root_number==0) 
    c.setLocalRootNumber(p);
  return (ri->second).local_root_number;
}

int GlobalRootNumber(CurveRed& c)
{
  int ans=-1;
  map<bigint,Reduction_type>::const_iterator ri; 
  for(ri = c.reduct_array.begin(); ri!=c.reduct_array.end(); ri++)
    {
      if((ri->second).local_root_number==0)  
	c.setLocalRootNumber(ri->first);
      ans *= (ri->second).local_root_number;
    }
  return ans;
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
  static const bigint two = BIGINT(2);
  map<bigint,Reduction_type>::iterator ri = reduct_array.find(two);
  if(ri==reduct_array.end()) return;
  Reduction_type& info = ri->second;
  int kod = PariKodairaCode(info.Kcode);
  int n2  = neron(2,kod); 

#ifdef DEBUG_ESIGN
  cout<<"\nIn LocalRootNumber2(), n2 = "<<n2<<"..."<<endl;
#endif

  bigint mu,mv; long u,v; int v4,v6;

  if (is_zero(c4)) {v4=12; u=0;}
  else {mu=c4; v4=divide_out(mu,two); u = posmod(mu,64);}
#ifdef DEBUG_ESIGN
  cout<<"c4="<<c4<<", v4="<<v4<<", u="<<u<<endl;
#endif

  if (is_zero(c6)) {v6=12; v=0;}
  else {mv=c6; v6=divide_out(mv,two); v = posmod(mv,64);}
#ifdef DEBUG_ESIGN
  cout<<"c6="<<c6<<", v6="<<v6<<", v="<<v<<endl;
#endif

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
  static const bigint three = BIGINT(3);
  map<bigint,Reduction_type>::iterator ri = reduct_array.find(three);
  if(ri==reduct_array.end()) return;
  Reduction_type& info = ri->second;
  int kod = PariKodairaCode(info.Kcode);
  int n2  = neron(3,kod); 

#ifdef DEBUG_ESIGN
  cout<<"\nIn LocalRootNumber3()..."<<endl;
#endif
  bigint mu,mv; long u,v; int v4,v6;

  if (is_zero(c4)) { v4=12; u=0;}
  else {mu=c4; v4=divide_out(mu,three); u = posmod(mu,81);}

#ifdef DEBUG_ESIGN
  cout<<"c4="<<c4<<", v4="<<v4<<", u="<<u<<endl;
#endif

  if (is_zero(c6)) { v=0;}
  else {mv=c6; v6=divide_out(mv,three); v = posmod(mv,81);}

#ifdef DEBUG_ESIGN
  cout<<"c6="<<c6<<", v6="<<v6<<", v="<<v<<endl;
#endif

  bigint tmp = discr;
  divide_out(tmp,three);
  long d1=posmod(tmp,81);

#ifdef DEBUG_ESIGN
  cout<<"d1="<<d1<<endl;
#endif

  long r6 = posmod(v,9);
  long K4=kro_3(u), K6=kro_3(v);
#ifdef DEBUG_ESIGN
  cout<<"r6="<<r6<<endl;
  cout<<"K4="<<K4<<endl;
  cout<<"K6="<<K6<<endl;
#endif

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
  map<bigint,Reduction_type>::iterator ri = reduct_array.find(p);
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
#ifdef DEBUG_ESIGN
  cout<<"kro("<<d<<","<<n<<") returns "<<flush;
#endif
  int ans = kronecker(d,n);
#ifdef DEBUG_ESIGN
  cout<<ans<<endl;
#endif
  return ans;
}

int kro(const bigint& d, long n)
{
#ifdef DEBUG_ESIGN
  cout<<"kro("<<d<<","<<n<<") returns "<<flush;
#endif
  int ans = kronecker(d,n);
#ifdef DEBUG_ESIGN
  cout<<ans<<endl;
#endif
  return ans;
}

int kro(long d, long n)
{
#ifdef DEBUG_ESIGN
  cout<<"kro("<<d<<","<<n<<") returns "<<flush;
#endif
  int ans = kronecker(d,n);
#ifdef DEBUG_ESIGN
  cout<<ans<<endl;
#endif
  return ans;
}

int kro_m1(long x) // kronecker(-1,x) with x>0 odd
{
  static int kro_m1_tab[4] = {0,1,0,-1};
#ifdef DEBUG_ESIGN
  if (!((x>0)&&(odd(x))))  
    { 
      cout<<"kro_m1() called with x="<<x<<endl;
      return 0;
    }
#endif
  return kro_m1_tab[x&3];
}

int kro_p2(long x) // kronecker(2,x) with x>0 odd
{
  static int kro_p2_tab[8] = {0,1,0,-1,0,-1,0,1};
#ifdef DEBUG_ESIGN
  if (!((x>0)&&(odd(x))))  
    {
      cout<<"kro_p2() called with x="<<x<<endl;
      return 0;
    }
#endif
  return kro_p2_tab[x&7];
}

int kro_m2(long x) // kronecker(-2,x) with x>0 odd
{
  static int kro_m2_tab[8] = {0,1,0,1,0,-1,0,-1};
#ifdef DEBUG_ESIGN
  if (!((x>0)&&(odd(x))))  
    {
      cout<<"kro_m2() called with x="<<x<<endl;
      return 0;
    }
#endif
  return kro_m2_tab[x&7];
}

int kro_3(long x) // kronecker(x,3)
{
  static int kro_3_tab[3] = {0,1,-1};
#ifdef DEBUG_ESIGN
  if (!(x>0))  
    {
      cout<<"kro_3() called with x="<<x<<endl;
      return 0;
    }
#endif
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
  cerr<<"neron() returns 0 -- should not happen!"<<endl;
  return 0; /* should not occur */
}

// Here the CurveRed parameter is not const since the call to
// LocalRootNumber may have to compute and store it

bigint Trace_Frob(CurveRed& c, const bigint& p)
{
  const bigint zero=BIGINT(0);
  const bigint one=BIGINT(1);
  const bigint two=BIGINT(2);
  const bigint three=BIGINT(3);
  //  cout<<"Trace_Frob at "<<p<<endl;

  int f = getord_p_N(c,p);
  // Bad primes: for convenience returns the p'th coefficient of the L-series
  if(f>=2)  return zero;
  if(f==1)  return BIGINT(-LocalRootNumber(c,p));

  int x,a,b,d;
  bigint n=zero; 
  if(p==two) // curvemodq class only in characteristic > 3 
    {
      // Count points naively
      // y^2+(a1*x+a3)*y-(x^3+a2*x^2+a4*x+a6) = y^2+ay+b
      int a1=bigint_mod_long(c.a1,2), a2=bigint_mod_long(c.a2,2), 
	a3=bigint_mod_long(c.a3,2),  a4=bigint_mod_long(c.a4,2), 
	a6=bigint_mod_long(c.a6,2);   
      // x=0:
      a = odd(a3);        // 1 if odd else 0
      b = odd(a6);               
      n += (a?(b?0:2):1);
      // x=1:
      a = odd(a1+a3); 
      b = odd(1+a2+a4+a6); 
      n += (a?(b?0:2):1);
      return two-n;
    }
  if(p==three) // curvemodq class only in characteristic > 3 
    {
      // Count points naively
      // y^2+(a1*x+a3)*y-(x^3+a2*x^2+a4*x+a6) = y^2+ay+b
      int a1=bigint_mod_long(c.a1,3), a2=bigint_mod_long(c.a2,3), 
	a3=bigint_mod_long(c.a3,3),  a4=bigint_mod_long(c.a4,3), 
	a6=bigint_mod_long(c.a6,3);   
      for(x=-1; x<2; x++)
	{
	  a = (((x+a2)*x+a4)*x+a6)%3;
	  b = (a1*x+a3)%3;
	  d = (b*b+a)%3;  
	  if(d==2)d=-1;
	  if(d==-2)d=1;
	  n += (d+1);
	}
      return three-n;
    }
#if defined(LiDIA_INTS) || defined(LiDIA_ALL)
      // Count points naively
      // y^2+(a1*x+a3)*y-(x^3+a2*x^2+a4*x+a6) = y^2+ay+b
      int a1=bigint_mod_long(c.a1,3), a2=bigint_mod_long(c.a2,3), 
	a3=bigint_mod_long(c.a3,3),  a4=bigint_mod_long(c.a4,3), 
	a6=bigint_mod_long(c.a6,3);   
      for(x=0; x<p; x++)
	{
	  bigint a = (a1*x+a3)%p;
	  bigint b = (-((x+a2)*x+a4)*x+a6)%p;
	  bigint d = (b*b-a)%p;  
	  n += (1+legendre(d,p));
	}
#else // NTL
  curvemodq Cq = reduce_curve(c,p); 
  n = Cq.group_order();
#endif
  bigint ans = one+p-n;
  return ans;
}

