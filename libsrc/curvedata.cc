// curvedata.cc -- implementation of Curvedata class
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
#include <eclib/cubic.h>

Curvedata::Curvedata(const bigint& aa1, const bigint& aa2, 
		     const bigint& aa3, const bigint& aa4, 
		     const bigint& aa6, int min_on_init)
  : Curve(aa1,aa2,aa3,aa4,aa6),
    b2(a1*a1 + 4*a2), b4(2*a4 + a1*a3), b6(a3*a3 + 4*a6),
    minimal_flag(0), ntorsion(0)
{
  b8 =  (b2*b6 - b4*b4) / 4;
  c4 = b2*b2 - 24*b4;
  c6 = -b2*b2*b2 + 36*b2*b4 - 216*b6;
  discr = (c4*c4*c4 - c6*c6) / 1728;
  discr_factored=0;
  if(sign(discr)==0)  // singular curve, replace by null
    {
      a1=0;a2=0;a3=0;a4=0;a6=0;
      b2=0;b4=0;b6=0;b8=0;
      c4=0;c6=0;
      conncomp=0;
  }
  else
    {
      conncomp = sign(discr)>0 ? 2 : 1;
      if (min_on_init) minimalize(); // which sets discr
    }
}

//#define DEBUG_Q_INPUT
Curvedata::Curvedata(const vector<bigrational>& qai, bigint& scale)
: minimal_flag(0), ntorsion(0)
{
  bigrational qa1(qai[0]), qa2(qai[1]), qa3(qai[2]), qa4(qai[3]), qa6(qai[4]);
  scale=BIGINT(1);
  a1=num(qa1);  a2=num(qa2);  a3=num(qa3);  a4=num(qa4);  a6=num(qa6);
#ifdef DEBUG_Q_INPUT
  cout<<"In Curvedata constructor with ["<<qa1<<","<<qa2<<","<<qa3<<","<<qa4<<","<<qa6<<"]"<<endl;
#endif
  vector<bigint> plist=pdivs(den(qa1));
  plist=vector_union(plist,pdivs(den(qa2)));
  plist=vector_union(plist,pdivs(den(qa3)));
  plist=vector_union(plist,pdivs(den(qa4)));
  plist=vector_union(plist,pdivs(den(qa6)));
#ifdef DEBUG_Q_INPUT
  cout<<"Denominator primes: "<<plist<<endl;
#endif
  for( const auto& p : plist)
    {
#ifdef DEBUG_Q_INPUT
      cout<<"p =  "<<p<<endl;
#endif
      long e=val(p,den(qa1));
      e=max(e,ceil(rational(val(p,den(qa2)),2)));
      e=max(e,ceil(rational(val(p,den(qa3)),3)));
      e=max(e,ceil(rational(val(p,den(qa4)),4)));
      e=max(e,ceil(rational(val(p,den(qa6)),6)));
#ifdef DEBUG_Q_INPUT
      cout<<"e =  "<<e<<endl;
#endif
      if(e>0)
	{
	  bigint pe=pow(p,e);
	  scale *= pe;
	  bigint pei=pe;
	  a1*=pei; pei*=pe;
	  a2*=pei; pei*=pe;
	  a3*=pei; pei*=pe;
	  a4*=pei; pei*=pe; pei*=pe;
	  a6*=pei;
	}
    }
  a1/=den(qa1);  a2/=den(qa2);  a3/=den(qa3);  a4/=den(qa4);  a6/=den(qa6);
#ifdef DEBUG_Q_INPUT
  cout<<"After scaling, coeffs are ["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"]"<<endl;
#endif
  b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3;
  b6 = a3*a3 + 4*a6; b8 =  (b2*b6 - b4*b4) / 4;
  c4 = b2*b2 - 24*b4; c6 = -b2*b2*b2 + 36*b2*b4 - 216*b6;
  discr = (c4*c4*c4 - c6*c6) / 1728;
  discr_factored=0;
  if(sign(discr)==0)  // singular curve, replace by null
    {
      a1=0;a2=0;a3=0;a4=0;a6=0;
      b2=0;b4=0;b6=0;b8=0;
      c4=0;c6=0;
      conncomp=0;
    }
  else
    {
      conncomp = sign(discr)>0 ? 2 : 1;
    }
}

Curvedata::Curvedata(const Curve& c, int min_on_init)
: Curve(c),
  b2(a1*a1 + 4*a2), b4(2*a4 + a1*a3), b6(a3*a3 + 4*a6),
  minimal_flag(0), ntorsion(0)
{
  b8 =  (b2*b6 - b4*b4) / 4;
  c4 = b2*b2 - 24*b4;
  c6 = -b2*b2*b2 + 36*b2*b4 - 216*b6;
  discr = (c4*c4*c4 - c6*c6) / 1728;
  discr_factored=0;
  if(sign(discr)==0)  // singular curve, replace by null
    {
      a1=0;a2=0;a3=0;a4=0;a6=0;
      b2=0;b4=0;b6=0;b8=0;
      c4=0;c6=0;
      conncomp=0;
    }
  else
    {
      conncomp = sign(discr)>0 ? 2 : 1;
      if (min_on_init) minimalize(); // which sets discr
    }
}

Curvedata::Curvedata(const Curvedata& c)
: Curve(c), b2(c.b2), b4(c.b4), b6(c.b6), b8(c.b8), c4(c.c4),
  c6(c.c6), discr(c.discr), minimal_flag(c.minimal_flag),
  discr_factored(c.discr_factored), conncomp(c.conncomp), 
  ntorsion(c.ntorsion)
{
  if(discr_factored) the_bad_primes=c.the_bad_primes;
}

Curvedata::Curvedata(const Curvedata& c, int min_on_init)
: Curve(c), b2(c.b2), b4(c.b4), b6(c.b6), b8(c.b8), c4(c.c4),
  c6(c.c6), discr(c.discr), minimal_flag(c.minimal_flag),
  discr_factored(c.discr_factored), conncomp(c.conncomp), 
  ntorsion(c.ntorsion)
{
  if(discr_factored) the_bad_primes=c.the_bad_primes;
  if (min_on_init) minimalize(); 
}

Curvedata::Curvedata(const bigint& cc4, const bigint& cc6, int min_on_init)
  :minimal_flag(0), discr_factored(0), ntorsion(0)
{
  if (valid_invariants(cc4, cc6))
    {
      c4=cc4; c6=cc6;
      c4c6_to_ai(cc4, cc6, a1, a2, a3, a4, a6, b2, b4, b6, b8);
      /*
      cout<<"a1="<<a1<<"\t";
      cout<<"a2="<<a2<<"\t";
      cout<<"a3="<<a3<<"\t";
      cout<<"a4="<<a4<<"\t";
      cout<<"a6="<<a6<<"\n";
      */
      if (min_on_init) minimalize(); 
      else 
	{
	  discr = (c4*c4*c4 - c6*c6) / 1728;
	}
      conncomp = sign(discr)>0 ? 2 : 1;
    }
  else 
    {
      cout << " ## attempt to call Curve constructor\n"
	   << "    with invalid invariants c4 = "<<cc4<<", cc6 = "<<c6
	   << ": reading as null curve\n";
      a1=0; a2=0; a3=0; a4=0; a6=0;
      b2=0; b4=0; b6=0; b8=0; 
      c4=0; c6=0; discr=0;
    }
}

void Curvedata::operator=(const Curvedata& c)
{
  a1 = c.a1, a2 = c.a2, a3 = c.a3, a4 = c.a4, a6 = c.a6;
  b2 = c.b2, b4 = c.b4, b6 = c.b6, b8 = c.b8;
  c4 = c.c4, c6 = c.c6, discr = c.discr;
  minimal_flag = c.minimal_flag;
  discr_factored = c.discr_factored;
  if(discr_factored) the_bad_primes = c.the_bad_primes;
  conncomp = c.conncomp;
  ntorsion = c.ntorsion;
}

// Minimalizing, reducing, and standardizing a general curve
// based on Laska--Kraus--Connell algorithm

void Curvedata::minimalize()
{
  if (minimal_flag) return;
  if ( isnull() ) {minimal_flag = 1; return; }
  
// else we are ready for Laska-Kraus-Connell reduction

  bigint newc4, newc6, newdiscr, u;
  //cout<<"minimising c4, c6 = "<<c4<<", "<<c6<<"\n";
  minimise_c4c6(c4,c6,discr,newc4,newc6,newdiscr,u);
  //cout<<"minimal c4, c6 = "<<newc4<<", "<<newc6<<"\t"<<"( u = "<<u<<")\n";
  
  // now compute minimal equation
  if ( u > 1)     {       c4 = newc4; c6 = newc6;     }
  discr = newdiscr;

  if(discr_factored)
    { 
      if(u>1) // trim list of bad primes
	{
	  vector<bigint> new_bad_primes;
          for (const auto& p : the_bad_primes)
	    {
	      if(div(p,discr))
                new_bad_primes.push_back(p);
	    }
	  the_bad_primes=new_bad_primes;
	}
    }
  else
    the_bad_primes = pdivs(discr);

  //  cout<<"After Curvedata::minimalize(): discr = "<<discr<<", ";
  //  cout<<"with bad primes "<<the_bad_primes<<endl;
  c4c6_to_ai(c4,c6,a1,a2,a3,a4,a6,b2,b4,b6,b8);
  minimal_flag = 1;
}

Curvedata Curvedata::minimalize(bigint& u, bigint& r, bigint& s, bigint& t) const
{
  if (minimal_flag) 
    {
      Curvedata newc(*this); 
      r=0; s=0; t=0; u=1; 
      return newc;
    }
  // else we use Laska-Kraus-Connell reduction

  bigint newc4, newc6, newdiscr, u2;
  //cout<<"minimising c4, c6 = "<<c4<<", "<<c6<<"\n";
  minimise_c4c6(c4,c6,discr,newc4,newc6,newdiscr,u);
  //cout<<"minimal c4, c6 = "<<newc4<<", "<<newc6<<"\t"<<"( u = "<<u<<")\n";

  Curvedata newc(newc4, newc6, 0);   // no need to re-minimalise
  //cout<<"minimal curve = "<<newc<<endl;
  
  s = (u*newc.a1 - a1)/2;  u2=u*u;
  r = (u2*newc.a2 - a2 + s*a1 + s*s)/3; 
  t = (u2*u*newc.a3 - a3 - r*a1)/2;
  //  cout<<"r,s,t="<<r<<","<<s<<","<<t<<endl;
  return newc;
}

void Curvedata::transform(const bigint& r, const bigint& s, const bigint& t)   //NB u=1;
{
  a6 += r*(a4 + r*(a2 + r)) - t*(a3 + r*a1 + t);
  a4 += -s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r*r - 2*s*t;
  a3 += r*a1 +t+t;
  a2 += -s*a1 + 3*r - s*s;
  a1 += s+s;
  b2 = a1*a1 + 4*a2;
  b4 = a4+a4 + a1*a3;
  b6 = a3*a3 + 4*a6;
  b8 = (b2*b6 - b4*b4) / 4;
}
       
void Curvedata::input(istream& is)
{
  Curve::input(is);
  b2=a1*a1 + 4*a2; b4=2*a4 + a1*a3;
  b6=a3*a3 + 4*a6; b8= (b2*b6 - b4*b4) / 4;
  c4=b2*b2 - 24*b4; c6=-b2*b2*b2 + 36*b2*b4 - 216*b6;
  discr= (c4*c4*c4 - c6*c6) / 1728;
  minimal_flag=0;
  discr_factored=0;
  conncomp= sign(discr)>0 ? 2 : 1;
  ntorsion=0;
}

void Curvedata::output(ostream& os) const
{
  Curve::output(os);
  if (isnull()) {os << " --singular\n"; return; }
  if (minimal_flag) os << " (reduced minimal model)";
  os << endl;
  os << "b2 = " << b2 << "\t " 
     << "b4 = " << b4 << "\t " 
     << "b6 = " << b6 << "\t " 
     << "b8 = " << b8 << endl;
  os << "c4 = " << c4 << "\t\t" 
     << "c6 = " << c6 << endl;
  os << "disc = " << discr << "\t(";
  if (minimal_flag&&discr_factored) 
    os << "bad primes: " << the_bad_primes << ";\t";
  os << "# real components = " << conncomp << ")" << endl;
  if (ntorsion) os<<"#torsion = "<<ntorsion<<endl;
  else os<<"#torsion not yet computed"<<endl;
}

Curvedata opt_x_shift(const Curvedata& C, bigint& k)
{
  bigint b2,b4,b6,b8,four,zero; four=BIGINT(4); zero=BIGINT(0);
  C.getbi(b2,b4,b6,b8);
  cubic b_cubic(four,b2,2*b4,b6);
  k = b_cubic.shift_reduce();
  Curvedata CD(C);
  CD.transform(k,zero,zero);
  return CD;
}


