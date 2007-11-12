// cubic.cc:  implementation of integer cubic class
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
 //

#include "marith.h"
#include "unimod.h"
#include "cubic.h"
#include "realroots.h"

void cubic::init()
{
  coeffs = new bigint[4];
}

cubic::~cubic()
{
  delete [] coeffs;
}

void cubic::transform(const unimod& m)
{
  bigint m112=sqr(m(1,1)); bigint m113=m112*m(1,1);
  bigint m212=sqr(m(2,1)); bigint m213=m212*m(2,1);
  bigint m222=sqr(m(2,2)); bigint m223=m222*m(2,2);
  bigint m122=sqr(m(1,2)); bigint m123=m122*m(1,2);

  coeffs[0] = m113*a() + m(2,1)*m112*b() + m212*m(1,1)*c() + m213*d();
  coeffs[1] = m123*a() + m(2,2)*m122*b() + m222*m(1,2)*c() + m223*d();

  coeffs[3] = 3*m(1,2)*m112*a() + (m(2,2)*m112 + 2*m(2,1)*m(1,2)*m(1,1))*b() 
    + (2*m(2,2)*m(2,1)*m(1,1) + m212*m(1,2))*c() + 3*m(2,2)*m212*d();

  coeffs[2] = 3*m122*m(1,1)*a() + (2*m(2,2)*m(1,2)*m(1,1) + m(2,1)*m122)*b() 
    + (m222*m(1,1) + 2*m(2,2)*m(2,1)*m(1,2))*c() + 3*m222*m(2,1)*d();

}

void cubic::x_shift(const bigint& e, unimod& m)
{
  coeffs[3] += e*(c()+e*(  b()+  e*a()));
  coeffs[2] +=        e*(2*b()+3*e*a());
  coeffs[1] +=                 3*e*a();
  m.x_shift(e);
}

void cubic::y_shift(const bigint& e, unimod& m)
{
  coeffs[0] += e*(b()+e*(  c()+  e*d()));
  coeffs[1] +=        e*(2*c()+3*e*d());
  coeffs[2] +=                 3*e*d();
  m.y_shift(e);
}

void cubic::invert(unimod& m)
{
  swap(coeffs[0],coeffs[3]); ::negate(coeffs[0]);
  swap(coeffs[1],coeffs[2]); ::negate(coeffs[2]);
  m.invert();
}

bigint cubic::j_c1() const
{
  bigint b2=sqr(b());
  bigint b3=b()*b2;
  bigint b4=b()*b3;
  bigint b5=b()*b4;
  bigint b6=b()*b5;
  bigint a2=sqr(a());
  bigint a3=a()*a2;
  bigint a4=a()*a3;
  bigint c2=sqr(c());
  bigint c3=c()*c2;
  bigint c4=c()*c3;
  bigint c5=c()*c4;
  bigint c6=c()*c5;
  bigint d2=sqr(d());
  bigint d3=d()*d2;
  bigint d4=d()*d3;
  bigint ac=a()*c(), bd=b()*d();
  return - 108*b3*a2*d() - 3*b4*c2 + 54*a2*c4 + 18*b5*d() + 243*a2*d2*b2 -
    54*b3*ac*d() - 162*bd*c2*a2 - 54*a3*c3 + 486*a3*bd*c() + 3*c4*b2 -
      18*c5*a() + 54*c3*a()*bd - 243*d2*a2*c2 + 162*d2*ac*b2 + 2*c6 -
	729*a4*d2 - 2*b6 + 18*b4*ac - 27*a2*b2*c2 + 729*d4*a2 + 54*b3*d3 +
	  108*c3*d2*a() - 18*c4*bd + 27*d2*c2*b2 - 486*d3*ac*b() - 54*d2*b4;

}

bigint cubic::j_c2() const
{
  bigint b2=sqr(b());
  bigint b3=b()*b2;
  bigint b4=b()*b3;
  bigint b5=b()*b4;
  bigint b6=b()*b5;
  bigint a2=sqr(a());
  bigint a3=a()*a2;
  bigint a4=a()*a3;
  bigint c2=sqr(c());
  bigint c3=c()*c2;
  bigint c4=c()*c3;
  bigint c5=c()*c4;
  bigint c6=c()*c5;
  bigint d2=sqr(d());
  bigint d3=d()*d2;
  bigint d4=d()*d3;
  bigint ac=a()*c(), bd=b()*d();

  return - 108*b3*a2*d() + 12*b4*c2 - 216*a2*c4 - 72*b5*d() - 486*a3*c2*d()
    + 270*a2*c3*b() - 90*b3*c2*a() - 972*a2*d2*b2 + 216*b3*ac*d() +
      648*bd*c2*a2 - 54*a3*c3 + 486*a3*bd*c() - 16*c3*b3 + 216*d2*b3*a() +
	72*d()*b4*c() + 72*c4*b()*a() + 216*d()*c3*a2 - 432*d()*b2*a()*c2 
    - 729*a4*d2 - 2*b6
	  + 18*b4*ac - 27*a2*b2*c2 + 6*b5*c() - 648*b2*c()*a2*d() 
    + 162*a()*d()*b4 +  1458*a3*d2*b();
}

bigint cubic::j_c3() const
{ 
  bigint a = coeffs[0], b=coeffs[1], c=coeffs[2], d=coeffs[3];
  bigint b2=b*b;
  bigint b3=b*b2;
  bigint b4=b*b3;
  bigint b5=b*b4;
  bigint b6=b*b5;
  bigint a2=a*a;
  bigint a3=a*a2;
  bigint a4=a*a3;
  bigint c2=c*c;
  bigint c3=c*c2;
  bigint c4=c*c3;
  bigint c5=c*c4;
  bigint c6=c*c5;
  bigint d2=d*d;
  bigint d3=d*d2;
  bigint d4=d*d3;

  return 108*b3*a2*d - 12*b4*c2 + 216*a2*c4 + 72*b5*d - 486*a3*c2*d +
    270*a2*c3*b - 90*b3*c2*a + 972*a2*d2*b2 - 216*b3*c*a*d - 648*b*c2*a2*d
      + 54*a3*c3 - 486*a3*d*c*b - 16*c3*b3 + 216*d2*b3*a + 72*d*b4*c +
	72*c4*b*a + 216*d*c3*a2 - 432*d*b2*a*c2 + 729*a4*d2 + 2*b6 - 18*b4*a*c
	  + 27*a2*b2*c2 + 6*b5*c - 648*b2*c*a2*d + 162*a*d*b4 + 1458*a3*d2*b;
}

//#define DEBUG

bigcomplex cubic::hess_root() const
{
  bigfloat discr = I2bigfloat(disc());
  if(!is_positive(disc())) 
    {
      cout<<"Error: hess_root called with negative dicriminant!\n";
      return to_bigfloat(0);
    }
  bigfloat P = I2bigfloat(p_semi());
  bigfloat Q = I2bigfloat(q_semi());
  bigfloat delta = sqrt(3*discr);
  bigcomplex gamma(-Q,delta); gamma/=(2*P);
  return gamma;  
}

void cubic::hess_reduce(unimod& m)
{
  int s=1;  bigint k;
  m.reset();

  while(s)
    {
      s=0;
      k = roundover(-q_semi(),2*p_semi());
      if(!is_zero(k))
	{
	  s=1;	  x_shift(k,m);
#ifdef DEBUG
	  cout << "Shift by " << k << ": " << (*this) << endl;
#endif
	}
      if(p_semi()>r_semi())
	{
	  s=1;	  invert(m);
#ifdef DEBUG
	  cout << "invert: " << (*this) << endl;
#endif
	}
    }
  if(a()<0) {::negate(coeffs[0]); ::negate(coeffs[2]);}
}

void cubic::mathews_reduce(unimod& m)
{
  int s=1;  bigint k; bigfloat alpha;
  m.reset();
  while(s)
    {
      s=0;
      if(mat_c1()<0) 
	{
	  s=1; invert(m);
#ifdef DEBUG
	  cout << "invert: " << (*this) << endl;
#endif
	}
      alpha = real_root();
      k = Iround(-alpha/2 - I2bigfloat(b())/I2bigfloat(2*a()));
      x_shift(k,m);
#ifdef DEBUG
      cout << "Shift by "<<k<<": "<<(*this)<<endl;
#endif
      bigint plus1, minus1;  plus1=1; minus1=-1;
      while(mat_c2()>0) 
	{
	  s=1; x_shift(plus1,m);
#ifdef DEBUG
	  cout << "Shift by +1: "<<(*this)<<endl;
#endif
	}
      while(mat_c3()<0) 
	{
	  s=1; x_shift(minus1,m);
#ifdef DEBUG
	  cout << "Shift by -1: "<<(*this)<<endl;
#endif
	}
    }
  if(a()<0) {::negate(coeffs[0]); ::negate(coeffs[2]);}
}

void cubic::jc_reduce(unimod& m)
{
  int s=1;   bigint k, jc2, jc3;
  bigint plus1, minus1;  plus1=1; minus1=-1;
  m.reset();

  while(s)
    {
      s=0;
      if(j_c1()<0) 
	{
	  s=1;  invert(m);
#ifdef DEBUG
	  cout << "invert: " << (*this) << endl;
#endif
	}
      bigfloat alpha = real_root();
      bigfloat ra = I2bigfloat(a());
      bigfloat rb = I2bigfloat(b());
      bigfloat rc = I2bigfloat(c());
      bigfloat h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
      bigfloat h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc; 
      k = Iround(-h1/(2*h0));
      x_shift(k,m);
#ifdef DEBUG
      cout << "Shift by "<<k<<": "<<(*this)<<endl;
#endif
      while(j_c2()>0) 
	{
	  s=1; x_shift(minus1,m);
#ifdef DEBUG
	  cout << "Shift by -1: "<<(*this)<<endl;
#endif
	}
      while(j_c3()<0) 
	{
	  s=1;   x_shift(plus1,m);
#ifdef DEBUG
	  cout << "Shift by +1: "<<(*this)<<endl;
#endif
	}
    }
  if(a()<0) {::negate(coeffs[0]); ::negate(coeffs[1]);}
}

  // Just shifts x:
bigint cubic::shift_reduce()
{
  unimod m; bigint k;
  if(is_positive(disc()))
    {
      k = roundover(-q_semi(),2*p_semi());
    }
  else
    {
      bigfloat alpha = real_root();
      bigfloat ra = I2bigfloat(a());
      bigfloat rb = I2bigfloat(b());
      bigfloat rc = I2bigfloat(c());
      bigfloat h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
      bigfloat h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc; 
      k = Iround(-h1/(2*h0));
    }
  x_shift(k,m);
  return k;
}

bigfloat cubic::real_root() const
{
  bigfloat discr = I2bigfloat(disc());
  if(discr>=0) 
    {
      cout<<"Error: real_root called with positive dicriminant!\n";
      return to_bigfloat(0);
    }
  bigfloat P = I2bigfloat(p_semi());
  bigfloat Q = I2bigfloat(q_semi());
  bigfloat A = I2bigfloat(a());

  if(is_zero(P)) 
    {
      bigfloat Q = I2bigfloat(q_semi());
      bigfloat R = I2bigfloat(r_semi())/Q;
      bigfloat eta3  = I2bigfloat(d())/A - (I2bigfloat(c())*R)/(3*A);
      bigfloat eta   = cube_root(eta3);
      bigfloat alpha = -eta - R;
      return alpha;
    }

  bigfloat U = I2bigfloat(u_semi());
  bigfloat delta = sqrt(-3*discr);
  bigfloat gamma1 = (-Q+delta)/(2*P);  // roots of Hessian
  bigfloat gamma2 = (-Q-delta)/(2*P);  //
  bigfloat eta3  = (U-3*A*delta)/(U+3*A*delta);
  bigfloat eta   = cube_root(eta3);
  bigfloat alpha = (eta*gamma1-gamma2)/(eta-1);
  return alpha;
}
