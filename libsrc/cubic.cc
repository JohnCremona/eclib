// cubic.cc:  implementation of integer cubic class
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
 //

#include <eclib/realroots.h>
#include <eclib/cubic.h>
#include <eclib/marith.h>
#include <eclib/polys.h>
#include <cassert>

bigint zero(0), one(1);

// comparison operator for sorting
int operator<(const cubic& F1, const cubic& F2)
{
  return std::lexicographical_compare(F1.coeffs.begin(), F1.coeffs.end(), F2.coeffs.begin(), F2.coeffs.end());
}

vector<bigint> transform_helper(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                                const unimod& m)
{
  bigint m112=sqr(m(1,1)); bigint m113=m112*m(1,1);
  bigint m212=sqr(m(2,1)); bigint m213=m212*m(2,1);
  bigint m222=sqr(m(2,2)); bigint m223=m222*m(2,2);
  bigint m122=sqr(m(1,2)); bigint m123=m122*m(1,2);

  bigint A = m113*a + m(2,1)*m112*b + m212*m(1,1)*c + m213*d;
  bigint B = 3*m(1,2)*m112*a + (m(2,2)*m112 + 2*m(2,1)*m(1,2)*m(1,1))*b
    + (2*m(2,2)*m(2,1)*m(1,1) + m212*m(1,2))*c + 3*m(2,2)*m212*d;
  bigint C = 3*m122*m(1,1)*a + (2*m(2,2)*m(1,2)*m(1,1) + m(2,1)*m122)*b
    + (m222*m(1,1) + 2*m(2,2)*m(2,1)*m(1,2))*c + 3*m222*m(2,1)*d;
  bigint D = m123*a + m(2,2)*m122*b + m222*m(1,2)*c + m223*d;

  return {A,B,C,D};
}

vector<bigint> transform_helper(const vector<bigint>& abcd, const unimod& m)
{
  return transform_helper(abcd[0],abcd[1],abcd[2],abcd[3],m);
}

void cubic::transform(const unimod& m)
{
  coeffs = transform_helper(coeffs, m);
}

cubic transform(const cubic& F, const unimod& m)
{
  return cubic(transform_helper(F.coeffs, m));
}

void cubic::sl2_reduce(unimod& m)
{
  if (disc()<0)
    jc_reduce(m);
  else
    hess_reduce(m);
}

//#define DEBUG_NORMALISE

// - for an sl2-reduced cubic, normalise w.r.t. <-I> (default) or <S> or <TS>
// - updates m by multiplying by normalising transformation
void cubic::normalise(unimod& m)
{
  int nautos = 2;

#ifdef DEBUG_NORMALISE
  cout<<"Normalising "<<(*this)<<endl;
#endif

  if (disc()<zero)
    {
      // Coeffs of covariant are [h0,h1,h2].  We have extra autos if
      // h0=h2, i.e. C1=0.  The covariant is const*(X^2+Y^2) if also
      // h1=0, i.e. C4=0, or const*(X^2+XY+Y^2) if also h1=h0,
      // i.e. C2=0.
      bigint C1=j_c1(), C2=j_c2(), C4=j_c4();
      nautos = (C1==0? (C4==0? 4: (C2==0? 6: 2)): 2);
#ifdef DEBUG_NORMALISE
      cout<<"Covariant quantities are C1="<<C1<<", C2="<<C2<<", C3="<<j_c3()<<", C4="<<C4<<endl;
#endif
    }
  else
    {
      // Coeffs of covariant are [P,Q,R].  We have extra autos if P=R.
      // The covariant is const*(X^2+Y^2) if also Q=0, or
      // const*(X^2+XY+Y^2) if also Q=R.
      bigint P=p_semi(), Q=q_semi(), R=r_semi();
      nautos = (P==R? (Q==0? 4: (Q==R? 6: 2)): 2);
#ifdef DEBUG_NORMALISE
      cout<<"Hessian coefficients are P="<<P<<", Q="<<Q<<", R="<<R<<endl;
      cout<<"Hessian root is "<<hess_root()<<endl;
#endif
    }

#ifdef DEBUG_NORMALISE
      cout<<"Number of automorphisms = "<<nautos<<endl;
#endif

  vector<cubic> Flist; // will hold the candidates
  vector<unimod> transforms; // will hold the transformations taking original cubic to candidates
  cubic F1 = *this;
  Flist.push_back(F1);
  unimod autpower; // initialised to identity
  transforms.push_back(autpower);
  unimod auto2(-one,zero,zero,-one), auto4(zero,-one,one,zero), auto6(one,one,-one,zero);
  unimod aut = (nautos==4? auto4: (nautos==6? auto6: auto2));
  for(int i=1; i<nautos; i++)
    {
      F1.transform(aut);
      Flist.push_back(F1);
      autpower*=aut;
      transforms.push_back(autpower);
    }

#ifdef DEBUG_NORMALISE
  cout<<"Comparing "<<Flist.size()<<" candidates for the reduced representative of "<<(*this)<<": "<<Flist<<endl;
#endif
  auto biggest = std::max_element(Flist.begin(),Flist.end());
  coeffs = biggest->coeffs;
  autpower = *(transforms.begin() + (biggest-Flist.begin()));
#ifdef DEBUG_NORMALISE
  cout<<"Largest one after sorting is "<<(*this)<<", the transform by "<<autpower<<endl;
#endif
  m *= autpower;
  return;
}

// Tests for sl2/gl2-equivalence:
int cubic::sl2_equivalent(const cubic& G) const
{
  unimod m; // not used but the reduction functions need one
  cubic F1(*this), G1(G);
  if (F1==G1)
    return 1;
  F1.negate(m);
  if (F1==G1)
    return 1;
  if (F1.disc()!=G1.disc())
    return 0;
  F1.sl2_reduce(m);
  G1.sl2_reduce(m);
  return F1==G1;
}

int cubic::gl2_equivalent(const cubic& G) const
{
  unimod m(-one,zero,zero,one);
  return (sl2_equivalent(G) || sl2_equivalent(::transform(G,m)));
}

// Test for sl2/gl2-equivalence to one in a list:
int cubic::sl2_equivalent_in_list(const vector<cubic>& Glist) const
{
  for (auto Gi=Glist.begin(); Gi!=Glist.end(); ++Gi)
    if (sl2_equivalent(*Gi))
      return 1;
  return 0;
}

int cubic::gl2_equivalent_in_list(const vector<cubic>& Glist) const
{
  for (auto Gi=Glist.begin(); Gi!=Glist.end(); ++Gi)
    if (gl2_equivalent(*Gi))
      return 1;
  return 0;
}

// affine roots of F mod q.
// NB rootsmod requires a non-constant polynomial

// affine roots of F mod q, assuming leading coefficient a() is nonzero:
vector<bigint> cubic::roots_mod(const bigint& q) const
{
  bigint aq(a()%q), bq(b()%q), cq(c()%q), dq(d()%q);
  if (is_zero(aq) && is_zero(bq) && is_zero(cq))
    return {};
  return rootsmod({dq,cq,bq,aq}, q);
}

// Return 1 iff F has a projective root mod q:
int cubic::has_roots_mod(const bigint& q) const
{
  return div(q,a()) || (roots_mod(q).size() > 0);
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

void cubic::negate(unimod& m)
{
  for (int i=0; i<4; i++)
    ::negate(coeffs[i]);
  m.negate();
}

void cubic::seminegate(unimod& m)
{
  for (int i=0; i<2; i++)
    ::negate(coeffs[2*i+1]);
  m.seminegate();
}

// The quantity called C_1 in the paper, = Norm(h2-h0) and should be
// NON-NEGATIVE for a reduced form:

bigint cubic::j_c1() const
{
  bigint a = coeffs[0], b=coeffs[1], c=coeffs[2], d=coeffs[3];
  bigint b2=sqr(b);
  bigint b3=b*b2;
  bigint b4=b*b3;
  bigint b5=b*b4;
  bigint b6=b*b5;
  bigint a2=sqr(a);
  bigint a3=a*a2;
  bigint a4=a*a3;
  bigint c2=sqr(c);
  bigint c3=c*c2;
  bigint c4=c*c3;
  bigint c5=c*c4;
  bigint c6=c*c5;
  bigint d2=sqr(d);
  bigint d3=d*d2;
  bigint d4=d*d3;
  bigint ac=a*c, bd=b*d;
  return - 108*b3*a2*d - 3*b4*c2 + 54*a2*c4 + 18*b5*d + 243*a2*d2*b2 -
    54*b3*ac*d - 162*bd*c2*a2 - 54*a3*c3 + 486*a3*bd*c + 3*c4*b2 -
      18*c5*a + 54*c3*a*bd - 243*d2*a2*c2 + 162*d2*ac*b2 + 2*c6 -
	729*a4*d2 - 2*b6 + 18*b4*ac - 27*a2*b2*c2 + 729*d4*a2 + 54*b3*d3 +
	  108*c3*d2*a - 18*c4*bd + 27*d2*c2*b2 - 486*d3*ac*b - 54*d2*b4;

}

// The quantity called C_2 in the paper, = Norm(h0-h1) and should be
// NON-NEGATIVE for a reduced form:

bigint cubic::j_c2() const
{
  bigint a = coeffs[0], b=coeffs[1], c=coeffs[2], d=coeffs[3];
  bigint b2=sqr(b);
  bigint b3=b*b2;
  bigint b4=b*b3;
  bigint b5=b*b4;
  bigint b6=b*b5;
  bigint a2=sqr(a);
  bigint a3=a*a2;
  bigint a4=a*a3;
  bigint c2=sqr(c);
  bigint c3=c*c2;
  bigint c4=c*c3;
  bigint c5=c*c4;
  bigint c6=c*c5;
  bigint d2=sqr(d);
  bigint d3=d*d2;
  bigint d4=d*d3;
  bigint ac=a*c, bd=b*d;

  return 108*b3*a2*d - 12*b4*c2 + 216*a2*c4 + 72*b5*d + 486*a3*c2*d
    - 270*a2*c3*b + 90*b3*c2*a + 972*a2*d2*b2 - 216*b3*ac*d -
      648*bd*c2*a2 + 54*a3*c3 - 486*a3*bd*c + 16*c3*b3 - 216*d2*b3*a -
	72*d*b4*c - 72*c4*b*a - 216*d*c3*a2 + 432*d*b2*a*c2 
    + 729*a4*d2 + 2*b6
	  - 18*b4*ac + 27*a2*b2*c2 - 6*b5*c + 648*b2*c*a2*d 
    - 162*a*d*b4 -  1458*a3*d2*b;
}

// The quantity called C_3 in the paper, = Norm(h0+h1) and should be
// POSITIVE for a reduced form:

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

// The quantity C_4 (not in the paper), = Norm(h1)/8 and should be
// NON-NEGATIVE for a reduced form with C1=0 (i.e. when h0=h2 we want h1>=0).

bigint cubic::j_c4() const
{
  bigint a = coeffs[0], b=coeffs[1], c=coeffs[2], d=coeffs[3];
  bigint b2=b*b;
  bigint b3=b*b2;
  bigint b4=b*b3;
  bigint a2=a*a;
  bigint c2=c*c;
  bigint c3=c*c2;
  bigint c4=c2*c2;
  bigint d2=d*d;

  return 27*d*c3*a2 + (27*d2*b3 - 54*d*c2*b2 + 9*c4*b)*a + 9*d*c*b4 - 2*c3*b3;
}

//#define DEBUG_REDUCE

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

int cubic::is_hessian_reduced()
// for positive discriminant only
// The condition is -P < Q <= P < R or 0 <= Q <= P=R.
{
  bigint P = p_semi();
  bigint R = r_semi();
  if (P>R) return 0;
  // now P<=R
  bigint Q = q_semi();
  if (Q>P) return 0;
  // now Q<=P<=R
  if (P==R) return (Q>=0);
  return (Q>-P);
}

void cubic::hess_reduce(unimod& m)
{
  int s=1;  bigint k;
  m.reset();
#ifdef DEBUG_REDUCE
  cout<<"Using hess_reduce() on "<<(*this)<<endl;
#endif
  while(s)
    {
      s=0;
      // NB roundover(a,b) returns c such that a/b=c+x and -1/2 < x <= 1/2,
      // so after the shift (when P>0) we have -P <= Q < P.
      k = roundover(-q_semi(),2*p_semi());
      if(!is_zero(k))
	{
	  s=1;	  x_shift(k,m);
#ifdef DEBUG_REDUCE
	  cout << "Shift by " << k << ": " << (*this) << endl;
#endif
	}
      if(p_semi()>r_semi())
	{
	  s=1;	  invert(m);
#ifdef DEBUG_REDUCE
	  cout << "invert: " << (*this) << endl;
#endif
	}
    }
  // Now we have -P <= Q < P <= R and test for boundary condition
  if ((p_semi()==r_semi()) && (q_semi()<0))
    {
      invert(m);
#ifdef DEBUG_REDUCE
      cout << "Final inversion: " << (*this) << endl;
#endif
    }
  normalise(m);
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
#ifdef DEBUG_REDUCE
	  cout << "invert: " << (*this) << endl;
#endif
	}
      alpha = real_root();
      k = Iround(-alpha/2 - I2bigfloat(b())/I2bigfloat(2*a()));
      if (k!=0)
        {
          s=1;
          x_shift(k,m);
#ifdef DEBUG_REDUCE
      cout << "Shift by "<<k<<": "<<(*this)<<endl;
#endif
        }
      bigint plus1, minus1;  plus1=1; minus1=-1;
      while(mat_c2()>0)
	{
	  s=1; x_shift(plus1,m);
#ifdef DEBUG_REDUCE
	  cout << "Shift by +1: "<<(*this)<<endl;
#endif
	}
      while(mat_c3()<0)
	{
	  s=1; x_shift(minus1,m);
#ifdef DEBUG_REDUCE
	  cout << "Shift by -1: "<<(*this)<<endl;
#endif
	}
    }
  if(a()<0) negate(m);
}

int cubic::is_jc_reduced() // for negative discriminant only
{
  if (is_zero(a())) // we want the quadratic form (b,c,d) to be reduced
    {
      bigint b(coeffs[1]), c(coeffs[2]), d(coeffs[3]);
      if (b==d)
        return ((0<=c) && (c<=b));
      else
        return ((-b<c) && (c<=b) && (b<d));
    }
  bigint C1 =  j_c1();
  if (C1<0) return 0;
  // now C1>=0, i.e. h0<=h2
  bigint C2 = j_c2();
  if (C2<0) return 0;
  // now C1, C2 >=0, i.e. h1<=h0<=h2
  if (is_zero(C1)) // i.e. h0=h2
    {
      bigint C4 =  j_c4(); // = N(h1)/8, not in JCM paper
      return (C4>=0); // i.e. h1 >= 0
    }
  bigint C3 =  j_c3();
  return (C3>0); // i.e. h1 > -h0
}

void cubic::jc_reduce(unimod& m)
{
  int s=1;   bigint k, jc2, jc3;
  bigint plus1(one), minus1(-one);

  bigfloat alpha, ra, rb, rc, rd, h0, h1, h2;

  m.reset();
#ifdef DEBUG_REDUCE
      cout << "\nJC-reducing " << (*this) << "...\n";
      cout<<"C1="<<j_c1()<<", C2="<<j_c2()<<", C3="<<j_c3()<<", C4="<<j_c4()<<endl;
      alpha = real_root();
      cout<<"alpha = "<<alpha<<endl;
      ra = I2bigfloat(a());
      rb = I2bigfloat(b());
      rc = I2bigfloat(c());
      rd = I2bigfloat(d());
      h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
      h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc;
      h2 = 3*(ra*rc*alpha + rb*rc-3*ra*rd)*alpha + 2*rc*rc - 3*rb*rd;
      cout << "(h0,h1,h2) = ("<<h0<<", " << h1 << ", "<<h2<<")"<<endl;
#endif

      if (is_zero(a()))
    {
      bigint bb=b(), cc=c(), q,r;
      if (bb<0)
        {
          bb = -bb;
          cc = -cc;
          m.negate();
        }
      ::divides(-cc,2*bb,q,r);
      if (r>=bb)
        q+=1;
#ifdef DEBUG_REDUCE
      cout << "[a=0] shift "<< (*this);
#endif
      x_shift(q,m);
#ifdef DEBUG_REDUCE
      cout << " by "<<q<< " to get " << (*this) << endl;
      cout<<"b="<<b()<<", c="<<c()<<", d="<<d()<<endl;
#endif
      // assert (is_jc_reduced());
      normalise(m);
      return;
    }

  while(s)
    {
      s=0;
      if (is_zero(a()))
        {
          jc_reduce(m);
          return;
        }
      if(j_c1()<0) // then h0 <= h2 fails, so invert:
	{
	  s=1;  invert(m);
#ifdef DEBUG_REDUCE
	  cout << "inverting --> " << (*this) << endl;
          cout<<"C1="<<j_c1()<<", C2="<<j_c2()<<", C3="<<j_c3()<<", C4="<<j_c4()<<endl;
          alpha = real_root();
          cout<<"alpha = "<<alpha<<endl;
          ra = I2bigfloat(a());
          rb = I2bigfloat(b());
          rc = I2bigfloat(c());
          rd = I2bigfloat(d());
          h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
          h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc;
          h2 = 3*(ra*rc*alpha + rb*rc-3*ra*rd)*alpha + 2*rc*rc - 3*rb*rd;
          cout << "(h0,h1,h2) = ("<<h0<<", " << h1 << ", "<<h2<<")"<<endl;
#endif
	}
      if ((j_c2()<0) || (j_c3()<=0)) // then -h0 < h1 <= h0 fails, so shift:
        {
          s=1;
          alpha = real_root();
#ifdef DEBUG_REDUCE
          cout<<"alpha = "<<alpha<<endl;
#endif
          ra = I2bigfloat(a());
          rb = I2bigfloat(b());
          rc = I2bigfloat(c());
          h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
          h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc;
          h2 = 3*(ra*rc*alpha + rb*rc-3*ra*rd)*alpha + 2*rc*rc - 3*rb*rd;
          k = Iround(-h1/(2*h0)); // this is the amount to shift by
          if (k!=0)
            {
              x_shift(k,m);
#ifdef DEBUG_REDUCE
              cout << "Shift by "<<k<<"--> "<<(*this)<<endl;
              cout<<"C1="<<j_c1()<<", C2="<<j_c2()<<", C3="<<j_c3()<<", C4="<<j_c4()<<endl;
              alpha = real_root();
              cout<<"alpha = "<<alpha<<endl;
              ra = I2bigfloat(a());
              rb = I2bigfloat(b());
              rc = I2bigfloat(c());
              rd = I2bigfloat(d());
              h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
              h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc;
              h2 = 3*(ra*rc*alpha + rb*rc-3*ra*rd)*alpha + 2*rc*rc - 3*rb*rd;
              cout << "(h0,h1,h2) = ("<<h0<<", " << h1 << ", "<<h2<<")"<<endl;
#endif
            }
          // Two loops to guard against rounding error in computing k:
          while(j_c2()<0) // h1>h0 so shift by -1
            {
              x_shift(minus1,m);
#ifdef DEBUG_REDUCE
              cout << "Shift by -1 --> "<<(*this)<<endl;
#endif
            }
          while(j_c3()<=0) // h1<=-h0 so shift by +1
            {
              x_shift(plus1,m);
#ifdef DEBUG_REDUCE
              cout << "Shift by +1--> "<<(*this)<<endl;
#endif
            }
        }
      if (is_zero(j_c1()) && (j_c4()<0)) // h0=h2 and h1<0, so invert
        {
          s=1;  invert(m);
#ifdef DEBUG_REDUCE
          cout << "final inversion--> " << (*this) << endl;
          cout<<"C1="<<j_c1()<<", C2="<<j_c2()<<", C3="<<j_c3()<<", C4="<<j_c4()<<endl;
          alpha = real_root();
          cout<<"alpha = "<<alpha<<endl;
          ra = I2bigfloat(a());
          rb = I2bigfloat(b());
          rc = I2bigfloat(c());
          rd = I2bigfloat(d());
          h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
          h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc;
          h2 = 3*(ra*rc*alpha + rb*rc-3*ra*rd)*alpha + 2*rc*rc - 3*rb*rd;
          cout << "(h0,h1,h2) = ("<<h0<<", " << h1 << ", "<<h2<<")"<<endl;
#endif
        }
#ifdef DEBUG_REDUCE
      cout<<"C1="<<j_c1()<<", C2="<<j_c2()<<", C3="<<j_c3()<<", C4="<<j_c4()<<endl;
      alpha = real_root();
      cout<<"alpha = "<<alpha<<endl;
      ra = I2bigfloat(a());
      rb = I2bigfloat(b());
      rc = I2bigfloat(c());
      rd = I2bigfloat(d());
      h0 = (9*ra*ra*alpha + 6*ra*rb)*alpha  + 6*ra*rc-rb*rb;
      h1 = 6*(ra*rb*alpha + (rb*rb-ra*rc))*alpha + 2*rb*rc;
      h2 = 3*(ra*rc*alpha + rb*rc-3*ra*rd)*alpha + 2*rc*rc - 3*rb*rd;
      cout << "(h0,h1,h2) = ("<<h0<<", " << h1 << ", "<<h2<<")"<<endl;
#endif
    }
  normalise(m);
  assert (is_jc_reduced());
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

  if(is_zero(A)) 
    {
      return A;
    }

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

vector<bigrational> cubic::rational_roots() const
{
  return roots(coeffs);
}


vector<cubic> reduced_cubics(const bigint& disc, int include_reducibles, int gl2, int verbose)
{
  bigint a, b, c, d;
  bigint amax, bmin, bmax, a2, b2, b3, b4, cmin, cmax, r;
  bigint P, U, absU, U2, Ud, Db2;
  int sU;
  bigfloat i3a, a23, ra2, rb2, D, D2, D3, D32, D4, Pmax, Pmin;

  unimod m;

  bigfloat third = 1/to_bigfloat(3);
  bigfloat const1 = 2 / sqrt(to_bigfloat(27));
  bigfloat const2 = 3 / pow(to_bigfloat(4), third);
  bigfloat const3 = sqrt(to_bigfloat(8)/to_bigfloat(27));
  bigfloat const4 = to_bigfloat(27)/to_bigfloat(8);
  bigfloat const5 = sqrt(pow(to_bigfloat(2), third));

  int neg=(disc<0);
  vector<cubic> glist; // will hold all cubics found
  vector<cubic> reduced_glist; // will hold unique reduced cubics found
  if(verbose>1) cout << "Discriminant = " << disc << endl;
  D = I2bigfloat(disc);
  D2 = sqrt(abs(D));
  D3 = pow(abs(D),third);
  if (neg) D3=-D3;
  D32 = D2*D2*D2;
  D4 = sqrt(D2);
  // Upper bound on a:
  amax = Ifloor((neg? const3: const1) * D4);
  if(verbose>1) cout<<"Upper bound on a: " << amax << endl;
  if (include_reducibles)
    {
      //
      // Code for a=0:
      //
      bmax = (neg? Ifloor(const5*D4): Ifloor(D4));
      if(verbose>1) cout<<"a=0, b<="<<bmax<<endl;
      for (b=1; b<=bmax; b++)
        {
          b2=b*b;
          if (::divides(disc,b2,Db2,r))
            {
              b4 = 4*b;
              for (c=1-b; c<=b; c++)
                {
                  if (::divides(c*c-Db2,b4,d,r))
                    {
                      if ((verbose>1) && glist.size()==0) cout<<disc<<" :\n";
                      cubic g(a,b,c,d);
                      assert(g.disc()==disc);
                      g.sl2_reduce(m);
                      glist.push_back(g);
                      if(verbose>1)
                        {
                          cout<<"found "<<g;
                          cout<<" (reducible: leading coefficient 0)\n";
                        }
                    }
                }
            }
        }
    } // end of a=0 code
  //
  // Code for a>0:
  //
  for(a=1; a<=amax; a++)
    {
      a2=a*a;
      ra2 = I2bigfloat(a2);
      a23 = pow(ra2, third);
      i3a = 1/I2bigfloat(3*a);
      Pmax = (neg? pow(2*D32+const4*D*ra2, third): D2);
      Pmin = const2*D3*a23;
      if(verbose>1) cout<<"bounds on P: ["<<Pmin<<","<<Pmax<<"]"<<endl;
      bmax=(3*a)/2;
      bmin=-bmax;
      if (2*bmax==3*a) bmin++;
      if(verbose>1) cout<<"a="<<a<<"; bounds on b: ["<<bmin<<","<<bmax<<"]"<<endl;
      for(b=bmin; b<=bmax; b++)
        {
          b2=b*b; b3=b*b2;
          rb2 = I2bigfloat(b2);
          cmin =  Ifloor((rb2-Pmax)*i3a); // round down for safety
          cmax = Iceil((rb2-Pmin)*i3a);   // round up for safety
          if(verbose>1) cout<<"(a,b)=("<<a<<","<<b<<"); bounds on c: ["<<cmin<<","<<cmax<<"]"<<endl;
          for(c=cmin; c<=cmax; c++)
            {
              P = b2-3*a*c;
              U2 = 4*P*P*P-27*disc*a2;
              Ud = 2*b3-9*a*b*c;
              if(verbose>1) cout<<"(a,b,c)=("<<a<<","<<b<<","<<c<<"): P="<<P<<", U^2="<<U2<<endl;
              if(isqrt(U2,absU))
                {
                  for (sU=0; sU<2; sU++)
                    {
                      U = (sU? -absU: absU);
                      if(::divides(U-Ud,27*a2,d,r))
                        {
                          if (verbose>1 && glist.size()==0) cout<<disc<<" :\n";
                          cubic g(a,b,c,d);
                          assert(g.disc()==disc);
                          g.sl2_reduce(m);
                          if(verbose>1) cout<<"found "<<g;
                          int irred = g.is_irreducible();
                          if (verbose>1)
                            {
                              if (irred)
                                cout<<" (irrreducible)\n";
                              else
                                cout<<" (reducible)\n";
                            }
                          if (irred or include_reducibles)
                              glist.push_back(g);
                        } // d integral test
                    }     // sign(U) loop
                }         // U square test
            } // c loop
        }     // b loop
    }         // a loop

  if (verbose)
    {
      cout << glist.size() << " cubics found with discriminant " << disc << ".";
      if (glist.size()>0) cout << glist << "\n Now reducing and eliminating repeats...";
      cout<<endl;
    }

  for (auto gi=glist.begin(); gi!=glist.end(); gi++)
    {
      cubic g = *gi;
      if(verbose) cout<<g;
      if (neg)
        {
          if (verbose>1) cout<<": testing JC-reduction..."<<endl;
          if (g.is_jc_reduced())
            {
              if(verbose) cout<<"\t (JC-reduced)";
            }
          else
            {
              if(verbose) cout<<"\t---(reduces to)--->\t";
              g.jc_reduce(m);
              if(verbose) cout<<g;
            }
        }
      else
        {
          if (g.is_hessian_reduced())
            {
              if(verbose) cout<<"\t (Hessian-reduced)";
            }
          else
            {
              if(verbose) cout<<"\t---(reduces to)--->\t";
              g.hess_reduce(m);
              if(verbose) cout<<g;
            }
        }

      // Check to see if this reduced cubic is already in the list:
      int equiv;
      if (gl2)
        equiv = g.gl2_equivalent_in_list(reduced_glist);
      else
        equiv = g.sl2_equivalent_in_list(reduced_glist);
      if (equiv)
        {
          if(verbose)
            {
              cout<<" -REPEAT";
              if(gl2)
                cout<<" (GL(2,Z)-equivalent)";
              else
                cout<<" (SL(2,Z)-equivalent)";
            }
        }
      else
        {
          reduced_glist.push_back(g);
          if(verbose)
            {
              cout<<" -NEW";
              if(gl2)
                cout<<" (not GL(2,Z)-equivalent)";
              else
                cout<<" (not SL(2,Z)-equivalent)";
            }
        }
      if(verbose) cout<<endl;
    }
  if (verbose && glist.size()>0)
    {
      cout << reduced_glist.size();
      cout << (gl2? " GL":" SL");
      cout << "(2,Z)-inequivalent cubics found with discriminant " << disc << "." << endl;
    }
  return reduced_glist;
}
