// curvedata.cc -- implementation of Curvedata class
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
#include "cubic.h"

Curvedata::Curvedata(const bigint& aa1, const bigint& aa2, 
		     const bigint& aa3, const bigint& aa4, 
		     const bigint& aa6, int min_on_init)
  : Curve(aa1,aa2,aa3,aa4,aa6),  minimal_flag(0), ntorsion(0)
{
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
      if (min_on_init) minimalize(); // which sets discr
    }
}

Curvedata::Curvedata(const Curve& c, int min_on_init)
: Curve(c), minimal_flag(0), ntorsion(0)
{
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
      if (min_on_init) minimalize(); // which sets discr
    }
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
      cerr << " ## attempt to call Curve constructor\n"
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
	  bigint p;
	  vector<bigint> new_bad_primes;
	  vector<bigint>::iterator pr = the_bad_primes.begin();
	  while(pr!=the_bad_primes.end())
	    {
	      bigint p=*pr++;
	      if(div(p,discr)) new_bad_primes.push_back(p);
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
  if (minimal_flag) os << " (minimal reduced form)";
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
  Curvedata CC(C);
  CC.transform(k,zero,zero);
  return CC;
}

#if(0)

// next the invariants package

bigfloat agm( bigfloat a, bigfloat b)
{
  bigfloat a0, b0 ;
  do{
    a0 = (a + b) / 2.0 ;
    b0 = sqrt( a * b ) ;
    a = a0 ; b = b0 ;
  } while ( !is_approx_zero(a-b) ) ;
  return ( a ) ;
}
 
bigfloat cuberoot(bigfloat x ) /* computes x^(1.0/3.0) for all x */
{ 
  bigfloat third(1); third/=3;
  if (x >= 0.0) 
    return (pow(x, third) ) ;
  else 
    return( - pow(-x, third) );
}

void  sort(bigfloat roots[3] )
{
        bigfloat x0, x1, x2, temp ;
        x0 = roots[0] ; x1 = roots[1] ; x2 = roots[2] ;
        if (x0 > x1) { temp = x0; x0 = x1 ; x1 = temp ; }  /* now x0 < x1 */
        if (x1 > x2) { temp = x1; x1 = x2 ; x2 = temp ; }  /* now x1 < x2, x0 < x2 */
        if (x0 > x1) { temp = x0; x0 = x1 ; x1 = temp ; }  /* now x0 < x1 < x2 */
        roots[0] = x0 ; roots[1] = x1; roots[2] = x2 ;
}

CurvedataExtra::CurvedataExtra(const Curvedata& c)
: Curvedata(c)
{
  if ( isnull() ) {roots[0]=roots[1]=roots[2]=period=0; return; }
  bigfloat fb2 = I2bigfloat(b2);
  bigfloat fc6 = I2bigfloat(c6);
  bigfloat fdiscr = I2bigfloat(discr);
  if (conncomp == 1) {
    bigfloat sqdisc  = 24.0 * sqrt(-3.0 * fdiscr) ;
    roots[2] = ( -fb2 + cuberoot(fc6 - sqdisc) +
                cuberoot(fc6 + sqdisc) ) / 12.0 ;
  } else {
    bigfloat theta = atan2( 24.0 * sqrt(3.0*fdiscr),  fc6) ;
    bigfloat fc4 = I2bigfloat(c4);
    bigfloat sqc4 = sqrt( fc4 ) ;
    roots[0] = (-fb2 + 2.0 * sqc4 * cos( theta/3.0 )) / 12.0 ;
    roots[1] = (-fb2 + 2.0 * sqc4 * cos( (theta+2.0*Pi())/3.0)) / 12.0;
    roots[2] = (-fb2 + 2.0 * sqc4 * cos( (theta+4.0*Pi())/3.0)) / 12.0;
    sort(roots) ;
  }
  // note that roots[2] (the THIRD one!) is always the largest root 
  // now compute the period
  bigfloat a, b ; // parameters to feed to agm()
  if (conncomp == 2){
    a = sqrt(roots[2] - roots[0]) ;
    b = sqrt(roots[2] - roots[1]) ;
  } else {
    bigfloat fb4 = I2bigfloat(b4); // bigfloat version of b4
    bigfloat fp = sqrt( fb4/2.0 + roots[2] * (fb2/2.0 + 3.0 * roots[2]) );
                                  // f', where f is the cubic 
    a = 2.0 * sqrt( fp ) ;
    b = sqrt( 3.0 * roots[2] + fb2/4.0 + 2.0 * fp ) ;
  }
  period =  (2.0 * Pi() / agm(a, b) );
}

void CurvedataExtra::output(ostream& os) const
{
  Curvedata::output(os);
  if(conncomp == 1) 
     os << "The real two-division point is " << roots[2] << endl ;
  else 
  if(conncomp == 2) 
     os << "The real two-division points are " 
        << roots[0] << ", " << roots[1] << ", "<< roots[2] << endl ; 
  os << "The real period is " << period << endl;
}
#endif

