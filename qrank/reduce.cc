// reduce.cc:  implementation of quartic reduction functions
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
#include "unimod.h"
#include "points.h"
#include "mquartic.h"
#include "transform.h"
#include "msoluble.h"
#include "minim.h"
#include "reduce.h"

#define REDUCE_B
//#define DEBUG_REDUCE
//#define DEBUG_ONESTEP
//#define DEBUG_FIRSTSTEP
//#define USE_OLD_BSD  // to use bsd-reduction for type 3 (worse!)

void reduce_b(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	      unimod& m)
{
  bigint a4=a<<2;
  bigint bmod4a = mod(b,a4);
  bigint alpha;
  divide_exact((bmod4a-b),a4,alpha);
  if(is_zero(alpha)) return;
  xshift(alpha,a,b,c,d,e,m);
  return;
}

bigfloat show(bigfloat x) {cout<<x<<endl; return x;}

bigfloat* types12_covar(const bigint& a, const bigint& b, const bigint& c, const bigint& d, 
			const bigfloat& xH, const bigfloat& phi)
{
  // The following covariant quadratic is the unique real quadratic factor 
  // of g6 which is positive definite
  bigfloat xa=I2bigfloat(a), xb=I2bigfloat(b), xc=I2bigfloat(c);
  bigfloat* ans = new bigfloat[3];
  ans[0] = 3*(4*xa*phi-xH);
  ans[1] = 6*(xb*phi+I2bigfloat(b*c-6*a*d));
  ans[2] = 2*phi*(xc-phi)+I2bigfloat(4*c*c-9*b*d);
#ifdef DEBUG_REDUCE
  cout<<"Before scaling, quadratic has coefficients: "
      <<ans[0]<<", "<<ans[1]<<", "<<ans[2]<<endl;
#endif
  ans[1]/=ans[0];
  ans[2]/=ans[0];
  ans[0]=1;
  return ans;
}

bigfloat* type3_covar(const bigfloat& xa, const bigfloat& xb, const bigfloat& xH, 
		      const bigfloat& rphi, const bigcomplex& cphi, int Risneg)
{
  bigfloat* ans = new bigfloat[3];
  bigfloat a4=4*xa;
  static bigfloat three(to_bigfloat(3));
  bigcomplex r1 = sqrt((a4*cphi-xH)/three);
  bigfloat r3   = sqrt(abs(a4*rphi-xH)/three);
  if(Risneg)  r3 = -r3;
#ifdef DEBUG_REDUCE
  cout<<"r1, r3 = "<<r1<<", "<<r3<<endl;
#endif
  bigfloat rr1 = abs(real(r1));
  bigfloat ir1 = abs(imag(r1));
  bigfloat rbeta   = ( r3 - xb )        / a4;
  bigfloat ibeta   = ( 2*ir1 )          / a4;
  bigfloat alpha_1 = ( 2*rr1 - r3 - xb) / a4;
  bigfloat alpha_2 = (-2*rr1 - r3 - xb) / a4;
  // roots are alpha_1, alpha_2 (real), rbeta+/-i*ibeta (non-real)
  bigfloat p      = -2*rbeta;
  bigfloat lambda =  2*ibeta;
  bigfloat q      = (p*p+lambda*lambda)/4;
#ifdef DEBUG_REDUCE
  cout<<"phi = "<<rphi<<endl;
  cout<<"Complex phi = "<<cphi<<endl;
  cout<<"Roots are "<<alpha_1<<", "<<alpha_2<<", "<<rbeta<<" +/-i* "<<ibeta<<endl;
  cout<<"p = "<<p<<", q = "<<q<<endl;
  cout << "lambda = " << lambda << endl;
#endif
#ifdef USE_OLD_BSD
  // B&SD's covariant:
  ans[0]=1; ans[1] = p; ans[2] = q;
#else
  //
  // Julia's covariant:
  bigfloat ar = abs(r1+r3);
  bigfloat br = abs(r1-r3);
  bigfloat t1sq = ir1*ar*ar;
  bigfloat t2sq = ir1*br*br;
  bigfloat usq  = rr1*ar*br;
  
  ans[0] = t1sq + t2sq + 2*usq;
  ans[1] = -2*alpha_1*t1sq -2*alpha_2*t2sq + 2*p*usq;
  ans[2] = alpha_1*alpha_1*t1sq + alpha_2*alpha_2*t2sq + 2*q*usq;

#ifdef DEBUG_REDUCE
  cout << "quadratic has coefficients "<<ans[0]<<", "<<ans[1]<<", "<<ans[2]<<endl;
  cout << "-disc = 4*q0*q2-q1^2 = " << (4*ans[0]*ans[2]-ans[1]*ans[1]) 
       << " (should be positive!)\n";
#endif
#endif
      
  ans[1] /= ans[0];
  ans[2] /= ans[0];
  ans[0] = 1;
  
  return ans;
}

// Compute the quadratic covariant of a real quartic:
bigfloat* quadratic_covariant(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e)
{
  bigint ii = II(a,b,c,d,e);
  bigint jj = JJ(a,b,c,d,e);
  bigint disc = 4*pow(ii,3)-sqr(jj);
  bigint  H = H_invariant(a,b,c), R = R_invariant(a,b,c,d);
  bigfloat xH = I2bigfloat(H);
  bigfloat xii = I2bigfloat(ii), xjj=I2bigfloat(jj);
  bigcomplex c1(to_bigfloat(0)), c2(-3*xii), c3(xjj);
  vector<bigcomplex> cphi = solvecubic( c1, c2, c3);
#ifdef DEBUG_REDUCE
      cout<<"Three roots phi are initially "<<cphi<<"\n";
#endif
  bigfloat * hcoeffs; // will hold coeffs of covariant quadratic
  bigfloat realphi, phi;       // will hold specific real roots
  long type;
  
  if(is_positive(disc))
    {
      // all the phi are real;  order them so that a*phi[i] decreases
      bigfloat phi1 = real(cphi[0]);
      bigfloat phi2 = real(cphi[1]);
      bigfloat phi3 = real(cphi[2]);
      if(is_positive(a))      orderreal(phi1,phi2,phi3); 
      else                    orderreal(phi3,phi2,phi1); 
#ifdef DEBUG_REDUCE
      cout<<"Three real phi are "<<phi1<<", "<<phi2<<", "<<phi3<<"\n";
#endif
      // So now a*phi1>a*phi2>a*phi3
      
      if((is_negative(H))&&(H*H>16*a*a*ii)) 
	{ 
	  type = 2; 
	  realphi = phi2; phi = phi2;
#ifdef DEBUG_REDUCE
	  cout<<"Type = 2, phi = "<<phi<<"\n";
#endif
	}
      else
	{ 
	  type = 1; 
	  realphi = phi1; phi = phi3;
#ifdef DEBUG_REDUCE
	  cout<<"Type = 1, phi = "<<phi<<"\n";
#endif
	}
      hcoeffs = types12_covar(a,b,c,d,xH,phi);
    }
  else // disc < 0
    {
      type=3;
#ifdef DEBUG_REDUCE
      cout<<"Type = 3, ";
#endif
      if (is_real(cphi[1])) 
	{
	  realphi=real(cphi[1]);
	  cphi[1]=cphi[2];
	  cphi[2]=realphi;
	}
      else 
	if (is_real(cphi[2])) 
	  {
	    realphi=real(cphi[2]);
	  }
	else 
	  {
	    realphi=real(cphi[0]);
	    cphi[0]=cphi[2];
	    cphi[2]=realphi;
	  }
#ifdef DEBUG_REDUCE
      cout<<"real phi = "<<realphi<<"\n";
      cout<<"complex phi = "<<cphi[0]<<"\n";
#endif
      bigfloat xa=I2bigfloat(a),xb=I2bigfloat(b);
      hcoeffs = type3_covar(xa,xb,xH,realphi,cphi[0],(sign(R)<0));
    }
  return hcoeffs;
}

void reduce(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	    unimod& m)
     // Construct a covariant quadratic, and reduce this
{
  bigfloat* hcoeffs = quadratic_covariant(a,b,c,d,e);
  unimod m1 = reduce_quad(hcoeffs[1],hcoeffs[2]);
  delete [] hcoeffs;
  // m1 contains the transform to reduce the quadratic;
  // now update the input transform and quartic:
  m*=m1;
  apply_transform(a,b,c,d,e,m1);

#ifdef DEBUG_REDUCE
  quartic newg(a,b,c,d,e);
  cout << "reduced quartic has coefficients " << newg << endl;
#endif

#ifdef REDUCE_B
  // Now reduce b so -2a < b <= 2a:
  bigint newa4=a<<2;
  bigint bmod4a = mod(b,newa4);
  bigint alpha;
  divide_exact((bmod4a-b),newa4,alpha);
  if(!is_zero(alpha))
    {
      xshift(alpha,a,b,c,d,e,m);
#ifdef DEBUG_REDUCE
      newg.assign(a,b,c,d,e);
      cout << "after reducing b (shifting by "<<alpha<<"), reduced quartic has coefficients " 
	   << newg << endl;
#endif
    }
#endif
}

// Finds a good unimodular matrix (a,b;c,d) 
// which raises z=x0+i*y0 when y0 is small
int first_step(const bigfloat& x0, const bigfloat& y0,
	       bigint& a, bigint& b, bigint& c, bigint& d);

// Finds a unimodular matrix (a,b;c,d) 
// which raises z=x0+i*y0 so im(z)>h
int one_step(const bigfloat& x0, const bigfloat& y0, const bigfloat& h,
	     bigint& a, bigint& b, bigint& c, bigint& d);

unimod reduce_quad_1(const bigfloat& b, const bigfloat& c);
unimod reduce_quad_2(const bigfloat& b, const bigfloat& c);

unimod reduce_quad(const bigfloat& b, const bigfloat& c)
{
  return reduce_quad_1(b,c);
}

// Given a pos. def. quadratic x^2+b*x+c, returns a unimod which
// reduces it (whose inverse takes its root into the fundamental
// region).
unimod reduce_quad_1(const bigfloat& bb, const bigfloat& cc)
{
  bigfloat b=bb, c=cc;
  bigfloat dz = 4*c-b*b;  // should be positive!
  bigfloat xz = -b/2;
  bigfloat yz = sqrt(abs(dz))/2;
#ifdef DEBUG_REDUCE_QUAD
  bigfloat az = c;
  cout << "Before any reduction:\n"; 
  cout << "Covariant quadratic has coefficients 1, "<<b<<", "<<c<<endl;
  cout << "-disc = " << dz << " (should be positive!)\n";
  cout << "Root has real part = " << xz << endl;
  cout << "imaginary part     = " << yz << endl;
  cout << "modulus^2          = " << az << endl;
#endif

  // Now do reduction.

  // Special first step: should raise so Im(z)>height
  bigfloat h(to_bigfloat(0.1));
  bigint ma, mb, mc, md;
  ma=1; mb=0; mc=0; md=1;
  one_step(xz,yz,h,ma,mb,mc,md);
#ifdef DEBUG_REDUCE_QUAD
  cout<<"one_step returns [a,b;c,d] = ["<<ma<<","<<mb<<";"<<mc<<","<<md<<"]\n";
#endif
  bigfloat xd=I2bigfloat(ma), xb=-I2bigfloat(mb), 
           xc=-I2bigfloat(mc), xa=I2bigfloat(md);
  bigfloat a1=(xa*xa+(b*xa+c*xc)*xc);
  bigfloat b1=2*(xa*xb+c*xc*xd)+b*(xa*xd+xb*xc);
  bigfloat c1=(xb*xb+(b*xb+c*xd)*xd);
  b=b1/a1;
  c=c1/a1;
  unimod m1(md,-mb,-mc,ma);

#ifdef DEBUG_REDUCE_QUAD
  dz = 4*a1*c1-b1*b1;  // should be positive!
  xz = -b1/(2*a1);
  yz = sqrt(abs(dz))/(2*a1);
  az = c1/a1;
  cout << "After first step of reduction:\n"; 
  cout << "Covariant quadratic has coefficients "<<a1<<", "<<b1<<", "<<c1<<endl;
  cout << "-disc = " << dz << " (should be positive!)\n";
  cout << "Root has real part = " << xz << endl;
  cout << "imaginary part     = " << yz << endl;
  cout << "modulus^2          = " << az << endl;
#endif

  bigint s = Iround(xz);
  bigfloat xk=I2bigfloat(s);
  m1.x_shift(s);
 
  c+=xk*(xk+b);  b+=2*xk;  xz-=xk;

  int reduced = (c>0.99999);
#ifdef DEBUG_REDUCE_QUAD
  dz = 4*c-b*b;  // should be positive!
  cout << "After preliminary shift by "<<s<<":\n"; 
  cout << "Covariant quadratic has coefficients 1, "<<b<<", "<<c<<endl;
  cout << "-disc = = " << dz << " (should be positive!)\n";
  cout << "Root has real part = " << xz << endl;
  cout << "imaginary part     = " << yz << endl;
  cout << "modulus^2          = " << c << endl;
#endif
  while(!reduced) 
    {
      // invert:      
      c=1/c; b*=(-c);
      m1.invert();

      // shift:
      xz = -b/2;
      s = Iround(xz);
      xk=I2bigfloat(s);
      c+=xk*(xk+b);  b+=2*xk;   xz -= xk;
      m1.x_shift(s);
      reduced = (c>0.99999);
#ifdef DEBUG_REDUCE_QUAD
      //recompute root:
      dz = 4*c-b*b;  // should be positive!
      yz = sqrt(abs(dz))/2;
      cout << "After inversion and shift by "<<s<<":\n"; 
      cout << "Covariant quadratic has coefficients "<<a<<", "<<b<<", "<<c<<endl;
      cout << "-disc = = " << dz << " (should be positive!)\n";
      cout << "Root has real part = " << xz << endl;
      cout << "imaginary part     = " << yz << endl;
      cout << "modulus^2          = " << c << endl;
#endif
    }
  return m1;
}

unimod reduce_quad_2(const bigfloat& b, const bigfloat& c)
{
  bigfloat xz = -b/2;
  bigfloat dz = 4*c-b*b;  // should be positive!
  bigfloat az = c;
  bigfloat yz = sqrt(abs(dz))/2;
  bigcomplex z(xz,yz);
#ifdef DEBUG_REDUCE_QUAD
  bigcomplex z0=z;
  cout << "Before any reduction:\n"; 
  cout << "Covariant quadratic has coefficients 1, "<<b<<", "<<c<<endl;
  cout << "-disc = " << dz << " (should be positive!)\n";
  cout << "Root = " << z << " with modulus^2 = " << az << endl;
#endif

  // Now do reduction.

  // Preliminary step shifts real part of root to [-1/2,1/2]
  //    gives better stability in low precision

  bigint s = Iround(xz);
  bigfloat xk=I2bigfloat(s);
  unimod m1; // default constructor initializes to identity
  m1.x_shift(s);
  xz-=xk;
  // yz is unchanged literally but we recompute it for better precision
  bigfloat b1=b+2*xk;
  bigfloat c1=xk*(xk+b)+c;
  dz = 4*c1-b1*b1;  // should be positive!
  az = c1;
  yz = sqrt(abs(dz))/2;
  z = bigcomplex(xz,yz);
  int reduced = (az>0.999);
#ifdef DEBUG_REDUCE_QUAD
  cout << "After preliminary shift:\n"; 
  cout << "Covariant quadratic has coefficients 1, "<<b1<<", "<<c1<<endl;
  cout << "-disc = = " << dz << " (should be positive!)\n";
  cout << "Root = " << z << " with modulus^2 = " << az << endl;
#endif
  if(reduced) return m1;

  // First stage: crude but guaranteed to double im(z):

  bigfloat a11, a12, a21, a22;
  bigint m11, m12, m21, m22;
  static bigfloat one(to_bigfloat(1));
  int changed=1;
  first_step(xz, yz, m11, m12, m21, m22);
  while(changed)
    {
      a11=I2bigfloat(m11); a12=I2bigfloat(m12); 
      a21=I2bigfloat(m21); a22=I2bigfloat(m22);
      z = (a11*z+a12) / (a21*z+a22);
      xz=real(z); yz=imag(z);
      az  = xz*xz + yz*yz;
      
      unimod n(Iround(a22),Iround(-a12),
	       Iround(-a21),Iround(a11));  // rounding the inverse matrix
      m1 *= n;
#ifdef DEBUG_REDUCE_QUAD
      cout << "z = "<<z<<endl;
      cout << "one matrix  = "<<n<<"\n"; 
      cout << "cumulative transform matrix = "<<m1<<"\n"; 
#endif
      changed=first_step(xz, yz, m11, m12, m21, m22);
    }

  // Second stage: standard algorithm

#ifdef DEBUG_REDUCE_QUAD
  cout << "After first stage , "<<endl;
  cout << "z     = "<<z<<endl;
#endif
  reduced=0;
  while(!reduced)
    {
      s = Iround(xz);
      xk=I2bigfloat(s);
      m1.x_shift(s);
      xz-=xk;
      z -=xk;
      az  = xz*xz + yz*yz;
#ifdef DEBUG_REDUCE_QUAD
      cout << "After shift by "<<s<<", "<<endl;
      cout << "z     = "<<z<<endl;
      cout << "|z|^2 = "<<az<<endl;
#endif
      reduced = (az>0.999);
      if(!reduced) 
	{
	  z=-one/z;      
	  xz=real(z); yz=imag(z);
#ifdef DEBUG_REDUCE_QUAD
	  cout << "After inverting, "<<endl;
	  cout << "z     = "<<z<<endl;
	  cout << "|z|^2 = "<<(1/az)<<endl;
#endif
	  m1.invert();
	}
    }

  // final shift of x-coord
  s = Iround(xz);
  m1.x_shift(s);

#ifdef DEBUG_REDUCE_QUAD
  bigfloat xs = I2bigfloat(s);
  xz -= xs;
  az  = xz*xz + yz*yz;
  z = bigcomplex(xz,yz);
  cout << "After reduction, root = " << z << " with modulus^2 = " 
       << az << endl;
  cout<<"Transform matrix = "<<m1<<endl;
  cout<<"This (inverted) applied to original root z0 = "<<z0<<endl;
  a11=I2bigfloat(m1(1,1)); a12=I2bigfloat(m1(1,2));  
  a21=I2bigfloat(m1(2,1)); a22=I2bigfloat(m1(2,2));
  bigcomplex z1=(a22*z0-a12) / (-a21*z0+a11);
  cout<<" gives z = "<<z1<<endl;
#endif
  return m1;
}

int first_step(const bigfloat& x0, const bigfloat& y0,
	       bigint& a, bigint& b, bigint& c, bigint& d)
// Finds a good unimodular matrix (a,b;c,d) 
// which raises z=x0+i*y0 when y0 is small
{
#ifdef DEBUG_FIRSTSTEP
  cout<<"In first_step with x0="<<x0<<", y0="<<y0<<endl;
#endif

  a=1; b=0; c=0; d=1;
  bigfloat xc = to_bigfloat(1)/(2*y0);
  c=Ifloor(xc);
  if(c<10) // return without doing anything
    {
      c=0; 
#ifdef DEBUG_FIRSTSTEP
      cout<<"first_step returns with no change"<<endl;
#endif
      return 0;
    } 
  d=-Iround(xc*x0);
  bigint g = bezout(-c,d,b,a); // a*d-b*c=g, may be >1
  if(g>1) {c/=g; d/=g;}
#ifdef DEBUG_FIRSTSTEP
  cout<<"first_step returns (a,b;c,d) = ("<<a<<","<<b<<";"<<c<<","<<d<<")"<<endl;
#endif
  return 1;
} 

int one_step(const bigfloat& x0, const bigfloat& y0, const bigfloat& h,
	     bigint& a, bigint& b, bigint& c, bigint& d)
// Finds a good unimodular matrix (a,b;c,d) 
// which raises z=x0+i*y0
{
#ifdef DEBUG_ONESTEP
  cout<<"In one_step with x0="<<x0<<", y0="<<y0<<endl;
#endif
  bigint k, newc, newd;
  bigfloat x(x0), s0(y0/h);
  bigfloat xk, x2, s1, s2, s=to_bigfloat(1), news;
  int i, ans=0;

  a=1; b=0; c=0; d=1;

  for(i=0; ; i++)
    {
      k = Iround(x); xk=I2bigfloat(k); x2=xk-x;
#ifdef DEBUG_ONESTEP
      cout<<i<<":\tk = "<<k<<", xk-x = "<<x2<<endl;
#endif
      newc=k*c-a; newd=k*d-b;
#ifdef DEBUG_ONESTEP
      cout<<"c'= "<<newc<<", d' = "<<newd<<endl;
#endif
      s1 = x0*I2bigfloat(newc)+I2bigfloat(newd);
      s2 = y0*I2bigfloat(newc);
      news = s1*s1+s2*s2;  // =Im(z)/Im(z')
#ifdef DEBUG_ONESTEP
      cout<<"s = "<<news<<endl;
      cout<<"Im(z) = "<<y0/news<<endl;
#endif
      if(news>s)      return ans; // new z lower than previous
      s=news;
      a= c;  c= newc;  b= d;  d= newd;
      if(s0>s)               return ans; // Im(new z)>0.1
      if(is_approx_zero(x2)) return ans;
      x = 1/x2;
#ifdef DEBUG_ONESTEP
      cout<<"New x = "<<x<<endl;
#endif
      ans=1;
    }
  return ans;
} 

