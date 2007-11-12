// pointsmod.cc: implementation of classes pointmodq and curvemodqbasis
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
 // and functions for point counting and elliptic curve discrete log


// Under LiDIA pointmodq is a wrapper for point<gf_element>; under NTL
// it is self-contained.  We provided a common interface.

// curvemodqbasis is derived from curvemodq (see file curvemod.h) and
// contains a Z-basis for the group of points

// The baby-step-giant step algorithm in my_bg_algorithm is adapted
// from LiDIA's bg_algorithm() with few changes.

// The point-counting and group structure algorithm in
// my_isomorphism_type() provide the same functionality as LiDIA's
// isomorphism_type() but has been rewritten from scratch by JEC; a
// main difference from the LiDIA version is the use of Weil pairing
// when the group is not cyclic.  This is only intended for use when q
// is small-medium sized (NOT cryptographic!) -- as is also true of
// LiDIA's isomorphism_type().  The current implementation is only for
// prime fields, but the same strategy would work over arbitrary
// finite fields.

#include "curve.h"
#include "points.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "ffmod.h"

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)

void set_order_point(pointmodq& P, const bigint& n)
{;}

#else // NTL

void set_order_point(pointmodq& P, const bigint& n)
{P.set_order(n);}

pointmodq::pointmodq(const gf_element&x, const curvemodq& EE)  // a point with X=x or oo if none
  : order(BIGINT(0)), E(EE)
{
  set_x_coordinate(x);
}

// make a point with given x & return true, or return false if none
int pointmodq::set_x_coordinate(const gf_element& x)
{
  is0flag=1; order=0;
  gf_element two=to_ZZ_p(2);
  gf_element four=to_ZZ_p(4);
  gf_element a1,a2,a3,a4,a6; E.get_ai(a1,a2,a3,a4,a6);
  //  cout<<"E = "<<E<<endl;
  //  cout<<"Trying x = "<<x<<endl;
  const bigint q=E.q;
  gf_element b2 = a1*a1 + four*a2; 
  gf_element b4 = two*a4 + a1*a3;
  gf_element b6 = a3*a3 + four*a6; 
  gf_element d = ((four*x+b2)*x+(two*b4))*x+b6;
  switch(legendre(rep(d),q)) // NTL has no modular sqrt!?
    {
    case -1: return 0;
    case 0:     case 1:
      is0flag=0;
      X=x;
      Y=(sqrt(galois_field(q),d)-(a1*x+a3))/two;
      if(!(on_curve()))
	{
	  cerr<<"Error in pointmodq::set_x_coordinate("<<x<<"): result "
	      <<(*this)<<" is not a valid point on "<<E<<endl;
	  cerr<<"b2,b4,b6 = "<<b2<<","<<b4<<","<<b6<<" mod "<<q<<endl;
	  cerr<<"d = "<<d<<" mod "<<q<<endl;
	  abort();
	}
    }
  return 1;
}

bigint pointmodq::get_order()
{
  if(order==BIGINT(0)) order=my_order_point(*this);
  return order;
}

bigint pointmodq::get_order(const bigint& mult)
{
  if(order==BIGINT(0)) order=my_order_point(*this,mult);
  return order;
}

bigint pointmodq::get_order(const bigint& lower, const bigint& upper)
{
  if(order==BIGINT(0)) order=my_order_point(*this,lower,upper);
  return order;
}

void pointmodq::output(ostream& os) const
{ 
  if(is0flag) 
    os<<"oo mod "<<(E.q); 
  else 
    os<<"("<<X<<","<<Y<<") mod "<<(E.q);
}

// addition of points, etc

pointmodq pointmodq::operator+(const pointmodq& Q) const // add Q to this
{
  pointmodq ans(Q.get_curve()); // initialized to oo
  if(is0flag) return Q;
  if(Q.is0flag) return *this;
  gf_element XQ=Q.X, YQ=Q.Y;
  if(X==XQ) 
    {
      if(Y==YQ)	return this->twice();
      else      return ans;             // =oo
    }
  gf_element lambda = (Y-YQ)/(X-XQ);
  gf_element mu     =  Y-lambda*X;
  ans.X = lambda*(lambda+(E.a1))-(E.a2)-X-XQ;
  ans.Y = lambda*(ans.X)+mu;
  ans.is0flag=0;
  ans.order=0;
  if(!(ans.on_curve())) cerr<<"error in pointmodq::operator+() adding "<<(*this)<<" to "<<Q<<endl;
  return ans.negate();
}

pointmodq pointmodq::operator-(const pointmodq& Q)  const // sub Q from this
{
  return *this + Q.negate();
}

pointmodq pointmodq::negate(void) const  // negates P
{
  if(is0flag) return pointmodq(E);
  return pointmodq(X,-Y-(E.a1)*X-(E.a3),E);
}

pointmodq pointmodq::operator-(void) const  // -P
{
  return this->negate();
}

pointmodq pointmodq::twice(void) const // doubles P
{
  pointmodq ans(E);
  if(is0flag) return ans;
  // Do NOT make these static as the modulus might change!
  gf_element two=to_ZZ_p(2);
  gf_element three=to_ZZ_p(3);
  gf_element a1,a2,a3,a4,a6; E.get_ai(a1,a2,a3,a4,a6);
  gf_element den = two*Y+a1*X+a3;
  if(den==0) return ans;
  gf_element lambda=(three*X*X+two*a2*X+a4-a1*Y)/den;
  gf_element mu    = Y-lambda*X;
  ans.X = lambda*(lambda+a1)-a2-two*X;
  ans.Y = lambda*(ans.X)+mu;
  ans.is0flag=0;
  ans.order=0;
  if(!(ans.on_curve())) 
    cerr<<"\nerror in pointmodq::twice() with P = "<<(*this)<<": "<<(ans)<<" not on "<<E<<endl;
  return ans.negate();
}

// calculates nP for long n
pointmodq operator*(long n, const pointmodq& P)  // n*P
{
  pointmodq ans(P.get_curve());
  if(P.is0flag || n == 0) return ans;
  int negative = (n < 0) ;
  if(negative) n = - n ;
  if(n == 1) {
    return (negative? -P : P);
  }
  // now n >= 2
  if(n == 2){
    ans = P.twice() ;
    return (negative? -ans : ans);
  }
  // now n >= 3
  if(n&1) ans = P ;    // (else ans is still 0 from initialization)
  pointmodq  Q = P ; 
  while(n > 1){
    Q = Q.twice() ; // 2^k P
    n /= 2 ;
    if(n&1) ans = ans + Q ;
  }
  return (negative? -ans : ans);
}

// calculates nP for bigint n
pointmodq operator*(const bigint& n, const pointmodq& P)  // n*P
{
  static bigint one = BIGINT(1);
  static bigint two = BIGINT(2);
  pointmodq ans(P.get_curve());
  if(P.is0flag || is_zero(n)) return ans;
  int negative = (is_negative(n)) ;
  bigint nn = n;
  if(negative) nn = - n ;
  if(nn == one) {
    return (negative? -P : P);
  }
  // now nn >= 2
  if(nn == two){
    ans = P.twice() ;
    return (negative? -ans : ans);
  }
  // now nn >= 3
  if(odd(nn)) ans = P ;    // (else ans is still 0 from initialization)
  pointmodq  Q = P ; 
  while(nn > one){
    Q = Q.twice() ; // 2^k P
    nn >>= 1 ;
    if(odd(nn)) ans = ans + Q ;
  }
  return (negative? -ans : ans);
}

pointmodq curvemodq::random_point()
{
  gf_element x;
  pointmodq ans(*this);
  while(ans.is_zero())
    {
      random(x);
      ans=pointmodq(x,*this);
    }
  return ans;
}

#endif //LiDIA/NTL split

pointmodq reduce_point(const Point& P,  const curvemodq& Emodq)
{
  //  cout<<"Reducing "<<P<<" mod q -> "<<flush;
  galois_field Fq = get_field(Emodq);
  NewGF(Fq,x);  NewGF(Fq,y);  NewGF(Fq,z);
  GFSetZ(z,getZ(P));
  if(IsZero(z)) return pointmodq(Emodq);
  GFSetZ(x,getX(P)); x/=z; 
  GFSetZ(y,getY(P)); y/=z; 
  //  cout<<"("<<x<<","<<y<<")"<<endl;
  return pointmodq(x,y,Emodq);
}

//#define DEBUG_ISO_TYPE 2

void curvemodqbasis::set_basis() 
{
  ffmodq(*this); // to initialize the class
  P1=pointmodq(*this);
  P2=P1;
  if(lazy_flag)
    {
      n2=1;
      one_generator(*this,n1,P1);
      return;
    }
  //#if defined(LiDIA_INTS) || defined(LiDIA_ALL)
  //  this->isomorphism_type(n1,n2,P1,P2);
  //#else
  my_isomorphism_type(*this,n1,n2,P1,P2);
  //#endif
  n=n1*n2;
#ifdef DEBUG_ISO_TYPE
  cout<<"Group structure of "<<(*this)<<" mod "<<::get_modulus(*this)<<": \n";
  if(n1>1) cout<<" gen 1 = "<<P1<<" (order "<<n1<<")\n";
  cout<<"Check: order is "<<order_point(P1)<<endl;
  if(n2>1) cout<<" gen 2 = "<<P2<<" (order "<<n2<<")"<<endl;
  cout<<"Check: order is "<<order_point(P2)<<endl;
  if(n2>1)
    {
      pointmodq Q1=(n1/n2)*P1; // order n2
      
      cout<<"Computing "<<n2<<"-Weil pairing of "<<Q1<<" and "<<P2
	  <<" mod "<<::get_modulus(*this)<<endl;

      long m = I2long(n2);
      gf_element mu = weil_pairing(Q1,P2,m);
      cout<<"Weil pairing of generators = "<<mu<<endl;
      gf_element mupow; power(mupow,mu,m);
      if (mupow==mu/mu)
	{
	  cout<<"OK, that's a "<<m<<"'th root of unity";
	  gf_element mupower = mu, one=mu/mu;
	  int m=1;
	  while(mupower!=one) {mupower*=mu; m++;}
	  cout<<" of exact order "<<m;
	  if(m==n2) cout<<" -OK"<<endl;
	  else cout<<" ???"<<endl;
	}
      else
	cout<<"WRONG, that's NOT a "<<m<<"'th root of unity"<<endl;
      
    }
#endif //  DEBUG_ISO_TYPE
  if(n1!=order_point(P1)) 
    {
      cerr<<"Error in isomorphism_type("<<(*this)<<") mod "<<::get_modulus((curvemodq)*this)
	  <<": n1 = "<<n1
	  <<" but point P1 = "<<P1<<" has order "<<order_point(P1)<<endl;
      n=n1=1; //  to prevent this reduction being used
    }
  if(n2!=order_point(P2)) 
    {
      cerr<<"Error in isomorphism_type("<<(*this)<<") mod "<<::get_modulus((curvemodq)*this)
	  <<": n2 = "<<n2
	  <<" but point P2 = "<<P2<<" has order "<<order_point(P2)<<endl;
      n=n2=1; //  to prevent this reduction being used
    }
}

pointmodq curvemodqbasis::get_gen(int i)
{
  if(i==1) return P1;
  if(i==2) return P2;
  return pointmodq(*this);
}

vector<pointmodq> curvemodqbasis::get_pbasis(int p)
{
  vector<pointmodq> ans;
  if((n%p)!=0) return ans;
#if 1 // defined(LiDIA_INTS) || defined(LiDIA_ALL)
  if((n1%p)==0) ans.push_back((n1/p)*P1);
  if((n2%p)==0) ans.push_back((n2/p)*P2);
#else
  ans = get_pbasis_via_divpol(p);
#endif
  return ans;
}

//#define DEBUG_PBASIS
vector<pointmodq> curvemodqbasis::get_pbasis_via_divpol(int p)
{
  vector<pointmodq> ans;
  if((n%p)!=0) return ans;
  FqPoly pdivpol = makepdivpol(*this, p);
#ifdef DEBUG_PBASIS  
  cout<<p<<"-division poly mod "<<get_modulus()<<" = "<<pdivpol<<endl;
#endif
  vector<gf_element> xi = roots(pdivpol);
#ifdef DEBUG_PBASIS  
  cout<<"roots of "<<p<<"-div pol mod "<<get_modulus()<<":  "<<xi<<endl;
#endif
  return get_pbasis_from_roots(p,xi);
}

vector<pointmodq> curvemodqbasis::get_pbasis_via_divpol(int p, const vector<bigint>& pdivpol)
{
  vector<pointmodq> ans;
  if((n%p)!=0) return ans;
  galois_field Fq = get_field(*this);
  NewFqPoly(Fq,pdivpolmodq);
  long i, deg = pdivpol.size()-1;
  SetDegree(pdivpolmodq,deg);
  for (i=0; i<=deg; i++) SetCoeff(pdivpolmodq,i,ZtoGF(Fq,pdivpol[i]));
#ifdef DEBUG_PBASIS  
  cout<<p<<"-division poly mod "<<get_modulus()<<" = "<<pdivpolmodq<<endl;
#endif
  vector<gf_element> xi = roots(pdivpolmodq);
#ifdef DEBUG_PBASIS  
  cout<<"roots of "<<p<<"-div pol mod "<<get_modulus()<<":  "<<xi<<endl;
#endif
  return get_pbasis_from_roots(p,xi);
}

//#define DEBUG_PBASIS
vector<pointmodq> curvemodqbasis::get_pbasis_from_roots(int p,  const vector<gf_element>& xi)
{
  vector<pointmodq> ans;
  if(xi.size()==0)
    {
#ifdef DEBUG_PBASIS  
      //      cout<<"no "<<p<<"-division points mod "<<get_modulus()<<endl;
#endif
      return ans;
    }
  unsigned int i;
  if(p==2)
    {
      for(i=0; (i<xi.size())&&(ans.size()<2); i++)
	{
	  pointmodq P(*this);
	  if(P.set_x_coordinate(xi[i])) ans.push_back(P);
	}
#ifdef DEBUG_PBASIS  
      cout<<"basis for 2-division points mod "<<get_modulus()<<": "<<ans<<endl;    
#endif
      return ans;
    }
  // p is now odd  
  unsigned int p12 = (p-1)/2, p212 = (p*p-1)/2;
  if(xi.size()==p212) // might have full p-torsion...
    {
      pointmodq P(*this);
      if(P.set_x_coordinate(xi[0]))
	 // then we do have full p-torsion, else _none_ of the xi will give rational points
	{
	  // store x-coords of multiples of p
	  ans.push_back(P);
	  vector<gf_element> xjp;   
	  pointmodq Q=P;
	  for(i=0; i<p12; i++)
	    {
	      {xjp.push_back(Q.get_x()); Q+=P;}
	    }
	  // now look for a point with x-coord not in that list
	  for(i=1; (i<xi.size())&&(ans.size()==1); i++)
	    if(find(xjp.begin(),xjp.end(),xi[i])==xjp.end())
	      {
		P.set_x_coordinate(xi[i]);
		ans.push_back(P);
	      }
	}
    } // end of all xi Fq-rational case
  else // we have at least one x
    {
      for(i=0; (i<xi.size())&&(ans.size()==0); i++)
	{
	  pointmodq P(*this);
	  if(P.set_x_coordinate(xi[i])) ans.push_back(P);
	}      
    }
#ifdef DEBUG_PBASIS  
  cout<<"basis for "<<p<<"-division points mod "<<get_modulus()<<": "<<ans<<endl;
#endif
  return ans;
}

// Baby-step-giant-step, point order, group structure
// EC discrete log via baby-step-giant-step adapted from LiDIA

const long MAX_BG_STEPS = 3000000;
#ifdef DEBUG_ISO_TYPE
const int debug_iso_type=DEBUG_ISO_TYPE;
#else
const int debug_iso_type=0;
#endif

bigint my_bg_algorithm(const pointmodq& PP,
		    const pointmodq& QQ,
                    const bigint& lower,
		    const bigint& upper,
		    bool info)
{
  //  cout<<"In my_bg_algorithm() with P="<<PP<<", Q="<<QQ<<", bounds "<<lower<<","<<upper<<endl;
  const bigint zero = BIGINT(0);
  const bigint minus_one = BIGINT(-1); // return value on failure
  if (PP.is_zero() && !QQ.is_zero()) return minus_one;

  if ((is_zero(lower)) && QQ.is_zero()) return zero;

  if (PP.get_curve() != QQ.get_curve())
    {
      cerr<<"bg_algorithm: Points P and Q on different curves"<<endl;
      return minus_one;
    }

  if ((is_negative(lower)) || (upper<lower))
    {
      cerr<<"bg_algorithm: lower bound > upper bound"<<endl;
      return minus_one;
    }

  pointmodq P(PP), Q(QQ);
  pointmodq H(P.get_curve()), H2(P.get_curve()), H3(P.get_curve());

  long i;
  bigint number_baby, number_giant, j, h;

  if (info)
    cout<<"\nBabystep Giantstep algorithm: "<<flush;
  
  if (upper - lower < BIGINT(30))    // for very small intervals
    {
      if (info)
        cout<<"\nTesting "<<(upper - lower) <<" possibilities ... "
                 <<flush; 
      H=lower*P;
      h=lower;
      if (H == Q) return h;
      do
        {
          H+=P;
          h++;
          if (H == Q)  return h;
        }
      while (h <= upper); 
      return minus_one;
    }
  
  //**** otherwise we use the Babystep Giantstep idea **************
  
  h = 1 + sqrt((upper - lower));  // compute number of babysteps
  if (h > MAX_BG_STEPS) h = MAX_BG_STEPS;
  number_baby=h;
  
  map<bigint,long> HT;
  
  H2 = Q-lower*P;
  H = pointmodq(P.get_curve());

  //****** Babysteps, store [x(i * P),i]  *********************
  
  if (info)
    cout << " (#Babysteps = " << number_baby << flush;
  
  for (i = 1; i <= number_baby; i++) 
    {
      H+=P;  // H = i*P and H2 = Q-lower*P

      if (H == H2)   // i * P = Q - lower* P, solution = lower+i
        {
	  if (info) cout<<") "<<flush;
#ifdef DEBUG
          assert((lower + i) * P ==  Q);
#endif
          return (lower + i);
        } // H==H2 case
  
      if (!H.is_zero()) // store [x(H),i] in table
	HT[LiftGF(H.get_x())]=i;
    }

  // Now for all i up to number_baby we have a table of pairs [x(i*P),i]
  // and H  =  number_baby*P
  // and H2 =  Q-lower*P

  // We will subtract H from H2 repeatedly, so H2=Q-lower*P-j*H in the loop

  //****** Giantsteps ***************************************************/
  
  number_giant = 1+((upper - lower)/(number_baby));
  
  if (info)
    cout << ", #Giantsteps = " << number_giant << ") " << endl;
  
  bigint step_size = number_baby; 

  for (j = 0; j <= number_giant; j++)
    {
      // Here H2=Q-(lower+j*step_size)*P
      if (H2.is_zero()) // on the nail, no need to check table
        {
	  h = lower + j * step_size;
#ifdef DEBUG
          assert(h*P == Q);
#endif
	  if (h <= upper) return h; else return minus_one;
        }  
      
      // look in table to see if H2= i*P for a suitable i

      map<bigint,long>::iterator HTi = HT.find(LiftGF(H2.get_x()));
      if(HTi!=HT.end())
	{
	  i = HTi->second;
	  H3=i*P;
	  if (H3 == H2) 
	    {
	      h = lower + i + j * step_size;
#ifdef DEBUG
                  assert(h * P ==  Q);
#endif
		  if (h <= upper) return h; else  return minus_one;
	    }
	} // H2 is in table
      H2-=H;
    } // loop on j
  return minus_one;
}

bigint my_order_point(const pointmodq& P, const bigint& mult)
{
  vector<bigint> plist = pdivs(mult);
  unsigned int i; bigint m, p, ans = BIGINT(1);
  if(P.is_zero()) return ans;
  for(i=0; i<plist.size(); i++)
    {
      p = plist[i];
      m = mult;  
      divide_out(m,p);
      pointmodq Q = m*P;
      while(!Q.is_zero()) {Q=p*Q; ans*=p;}
    }
  return ans;
}

bigint my_order_point(const pointmodq& P, const bigint& lower, const bigint& upper)
{
  return my_order_point(P,my_bg_algorithm(P,pointmodq(P.get_curve()),lower,upper));
}

bigint my_order_point(const pointmodq& P)
{
  bigint q = get_field(P.get_curve()).characteristic();
  bigint lower, upper;
  set_hasse_bounds(q,lower,upper);
  return my_order_point(P,lower,upper);
}

// returns minimal m>0 s.t. m*Q is in <P> with m*Q=a*P; n is assumed
// to be the order of P.  Special case: if <Q> and <P> are disjoint,
// then m=order(Q) and a=0.  On input, m holds order(Q) if known, else 0
bigint linear_relation( pointmodq& P,  pointmodq& Q, bigint& a)
{
  static bigint zero = BIGINT(0);
  static bigint one = BIGINT(1);
  bigint n = order_point(P), m = order_point(Q), g,n1,m1,h;
  int debug_linear_relation=0;
  if(debug_linear_relation)
    cout<<"In linear_relation() with P = "<<P<<" of order "<<n
	<<" and Q = "<<Q<<" of order "<<m<<endl;
  g=gcd(n,m);
  if(debug_linear_relation) cout<<"gcd = "<<g<<endl;
  if(g==one) {a=zero; return g /* =1 */;}
  n1=n/g; m1=m/g;
  pointmodq P1=n1*P; // both of exact order g: 
  pointmodq Q1=m1*Q; // now see if Q1 is a mult of P1
  if(debug_linear_relation)
    cout<<"P1 = "<<P1<<" and Q1 = "<<Q1<<endl;
  h=g;           // holds h s.t. h*Q1 is in <P1>
  vector<bigint> dlist = posdivs(g);
  sort(dlist.begin(),dlist.end());
  a=-1;
  for(unsigned int i=0; (i<dlist.size())&&(a==-1); i++)
    {
      h = dlist[i];
      a = my_bg_algorithm(P1,h*Q1,zero,g-1);
      if(debug_linear_relation)	cout<<"h = "<<h<<"; a = "<<a<<endl;
    }
  a=a*n1;
  m=h*m1;
  // debugging:
  if(m*Q!=a*P) {cerr<<"Error in linear relation with P="<<P<<", n="<<n<<", Q="<<Q<<": returns a="<<a<<" and m="<<m<<endl;}
  return m;
}

void set_hasse_bounds(const bigint& q, bigint& l, bigint& u)
{
  static const bigint one=BIGINT(1);
  sqrt(u, q << 2);
  l = q + one - u;    // lower bound of Hasse interval
  if (is_negative(l))  l=one;
  u = q + one + u;    // upper bound of Hasse interval  
}

// Given positive integers m,n, replace them by divisors which are
// coprime and have the same lcm
//#define DEBUG 1
bigint tidy_lcm(bigint& m, bigint& n)
{
#ifdef DEBUG
  bigint m0=m, n0=n;
#endif
  bigint g=gcd(m,n);
  bigint l=m*n/g;  // = lcm(m,n)
  g=gcd(m,n/g);  // divisible by primes dividing n to a higher power than m
  while(g!=BIGINT(1)) {m/=g; g=gcd(m,g);}
  n=l/m;
#ifdef DEBUG
  if((m*n==lcm(m0,n0)) &&
     (gcd(m,n)==BIGINT(1)) &&
     (m0%m==BIGINT(0)) && 
     (n0%n==BIGINT(0))) 
    {
      cout<<"tidy_lcm("<<m0<<","<<n0<<") changes them to "<<m<<","<<n<<" and returns "<<l<<endl;
      return l;
    }
  cout<<"Error in tidy_lcm("<<m0<<","<<n0<<")"<<endl;
  return BIGINT(0);
#else
  return l;
#endif
}

// merge_points_1: given a point P of order ordP, and another point Q,
// replaces P with a point of order lcm(ordP,order(Q)) (and updates
// ordP)

void merge_points_1(pointmodq& P, bigint& ordP,  pointmodq& Q)
{
  // First easy case:  order(Q) divides order(P), do nothing

  if ((ordP * Q).is_zero())
    {
      if(debug_iso_type>1)  cout<<"Order(Q) divides order(P)" <<endl;
      return;
    }
  bigint ordQ = order_point(Q);
  if(debug_iso_type>1)  cout<<"Order(Q) = "<< ordQ <<endl;

  // Second easy case:  order(P) divides order(Q), swap P & Q

  if (ordQ%ordP==0)
    {
      if(debug_iso_type>1)  cout<<"Order(P) divides order(Q)" <<endl;
      P=Q; ordP=ordQ;
      return;
    }

  // General case:
  // Construct a point whose order is lcm(ordP, ordQ):

  bigint nP=ordP, nQ=ordQ;
  bigint nPQ=tidy_lcm(nP,nQ);
  // Now (1) nP*nQ = lcm(ordP,ordQ)
  //     (2) gcd(nP,nQ)=1
  //     (3) nP|ordP and nQ|ordQ
  // So (ordP/nP)*P has order nP, similarly for Q, and
  // these orders are coprime so the sum of the points
  // has order nP*nQ=lcm(ordP,ordQ) as required:

  P = (ordP/nP)*P + (ordQ/nQ)*Q; 
  ordP = nPQ;
  if(debug_iso_type) 
    {
      cout<<"Changed P = "<<P<<":\t order(P) = "<<nPQ<<endl;
      if(order_point(P)!=nPQ) cout<<"that's wrong!"<<endl;
    }
  set_order_point(P,nPQ);
}	

// merge_points_2: given independent points P1, P2 of orders n1,n2,
// and another point Q, EITHER (unusual) replaces P1 with a point of
// higher order and updates n2target and replaces P2 with 0, OR
// (usual) replaces P2 with a point of higher order still independent
// of P1.  Here n2target is such that we expect the final group order
// to be n1*n2target, based on an assumption that P1 does have maximal
// order n1, but this assumption must be revised if using Q reveals a
// point of order greater than P1

void merge_points_2(pointmodq& P1, bigint& n1, pointmodq& P2, bigint& n2, 
		    bigint& n2target, pointmodq& Q)
{
  // Case 1: Q cannot improve if its order divides n2:

  pointmodq Q1 = n2*Q;
  if(Q1.is_zero())
    {
      if(debug_iso_type>0)  cout<<"Order(Q) divides n2=" <<n2<<endl;
      return;
    }

  pointmodq Q2 = (n1/n2)*Q1; // = n1*Q
  if(!(Q2.is_zero()))
    {
  // Case 2: P1 needs updating and we discard P2
      if(debug_iso_type>0)  
	cout<<"Order(Q) does not divide n1="<<n1<<", updating P1" <<endl;
      bigint oldn1=n1;
      merge_points_1(P1,n1,Q);
      n2target=(n2target*oldn1)/n1;
      if(debug_iso_type>0)  
	cout<<"New P1 has order " <<n1<<", assuming group structure "
	    <<n1<<"*"<<n2target<<endl;
      if(n2>1) {P2 = pointmodq(P2.get_curve()); n2=1;}
      return;
    }

  // General Case 3:

  // We find a multiple a*P1 such that Q-a*P1 is killed by n2target so
  // we can apply the Weil Pairing of order n2target
  Q1 = n2target*Q;
  Q2 = n2target*P1;  // has exact order n1/n2target
  bigint a = my_bg_algorithm(Q2,Q1,BIGINT(0),n1/n2target);
  if(a==BIGINT(-1)) // dlog failed, n1 must be wrong
  {
      if(debug_iso_type) 
      {
	  cout<<"Dlog of "<<Q1<<" w.r.t. "<<Q2<<" (order "<<n1/n2target
	      <<") does not exist, so current n1 must be too small"<<endl;
      }
      return;
  }
  if(debug_iso_type) 
    {
      cout<<"Dlog of "<<Q1<<" w.r.t. "<<Q2<<" (order "<<n1/n2target<<") is "<<a<<endl;
      cout<<"Check: a*Q2-Q1 =  "<<a*Q2-Q1<<" ( should be zero)"<<endl;
    }
  Q = Q-a*P1;  // this is killed by n2target
  if(Q.is_zero()) // then Q is a multiple of P, so we gain nothing
    {
      if(debug_iso_type) cout<<"Q-a*P1 = "<<Q<<", no use"<<endl;
      return;
    }
  if(debug_iso_type) 
    {
      cout<<"Replacing Q by Q-a*P1 = "<<Q<<" where a = "<< a << endl;
      cout<<"whose order divides n2target; computing Weil pairing of order "
	  <<n2target<<endl;
    }
  // At this point we have not changed the subgroup generated by P1,Q
  // (and have not touched P2) but now Q has order dividing n2target
  // (as for P2).  We now use the Weil pairing of P1 and Q, which is
  // an n2target'th root of unity
  Q1 = (n1/n2target)*P1;
  if(debug_iso_type) 
    {
      cout<<"order((n1/n2target)*P1) = "<<Q1<<" is "<<order_point(Q1)<<endl;
      cout<<"order(Q) =                "<<Q<<" is "<<order_point(Q)<<endl;
    }
  gf_element zeta = weil_pairing(Q1,Q,I2long(n2target));
  if(debug_iso_type)  cout<<"zeta = "<< zeta <<endl;
  if(IsZero(zeta))
  {
      cerr<<"Error: weil_pairing returns 0!"<<endl;
      cerr<<"n1 = "<<n1<<endl;
      cerr<<"n2 = "<<n2<<endl;
      cerr<<"n2target = "<<n2target<<endl;
      cerr<<"order((n1/n2target)*P1) = "<<Q1<<" is "<<order_point(Q1)<<endl;
      cerr<<"order(Q) =                "<<Q<<" is "<<order_point(Q2)<<endl;
  }
  bigint m = order(zeta);
  if(debug_iso_type)  cout<<"order = "<< m <<endl;
  // Compare this with n2 to see if we have gained:
  bigint l = lcm(n2,m);
  if(l==n2) return; // no gain
  bigint ordQ = my_order_point(Q,n2target);
  Q1 = (ordQ/m)*Q;   // of order m
  if(l==m) // replace P2, n2 by Q1, m
    {
      P2 = Q1;
      n2 = m;
      return;
    }
  // Now P2,Q1 have orders n2,m & both are independent of P1,
  // so we combine them to get a point of order l=lcm(n2,m)
  // still independent of P1
  bigint n2d=n2, md=m;
  l = tidy_lcm(n2d,md);
  P2 = (n2/n2d)*P2 + (m/md)*Q1; // of order n2d*md=l
  n2 = l;
  if(debug_iso_type) 
    {
      cout<<"Changed P2 = "<<P2<<":\t order(P2) = "<<n2<<endl;
      if(order_point(P2)!=n2) cout<<"that's wrong!"<<endl;
    }
}

// returns list of integers n2 such that
// (1) lower <= n1*n2 <= upper
// (2) n2|gcd(n1,q-1)
vector<bigint> n2list(const bigint& n1, 
		      const bigint& lower, const bigint& upper, 
		      const bigint& q)
{
  bigint n2min = lower/n1, n2max = upper/n1, n2,  g = gcd(n1,q-1);
  if(n2min*n1<lower) n2min++;
  vector<bigint> ans;
  for(n2=n2min; n2<=n2max; n2++) if(div(n2,g)) ans.push_back(n2);
  return ans;
}

// find a point of "large" order
void one_generator(curvemodq& Cq, bigint& n1, pointmodq& P1)
{
  galois_field Fq = get_field(Cq);
  bigint q = Fq.characteristic();
  bigint upper, lower; // bounds on group order
  set_hasse_bounds(q,lower,upper);
  if(debug_iso_type)
    cout<<"Lower and upper bounds on group order: ["
	<<lower<<","<<upper<<"]"<<endl;

  P1 = Cq.random_point();
  if(debug_iso_type) cout<<"P1 = "<<P1<<":\t"<<flush;
  n1 = my_order_point(P1,lower,upper);
  if(debug_iso_type) cout<<"Order(P1) = "<< n1 <<endl;

  int n;
  for(n=1; ((n<=10)&&(2*n1<=upper));  n++)
    { 
      pointmodq Q = Cq.random_point();
      if(debug_iso_type>1)  cout<<"Q = "<<Q<<":\t"<<flush;
      merge_points_1(P1,n1,Q);
      if(debug_iso_type>1)  
	{
	  cout<<"now P1 = "<<P1<<":\tof order "<<n1<<endl;
	}
    }
}

// find full Z-basis
void my_isomorphism_type(curvemodq& Cq, 
			 bigint& n1, bigint& n2, pointmodq& P1, pointmodq& P2)
{
  galois_field Fq = get_field(Cq);
  bigint q = Fq.characteristic();
  bigint upper, lower; // bounds on group order
  set_hasse_bounds(q,lower,upper);
  if(debug_iso_type)
    cout<<"Lower and upper bounds on group order: ["
	<<lower<<","<<upper<<"]"<<endl;
  int group_order_known=0;
  if((q<100)||(q==181)||(q==331)||(q==547))
    {
      Cq.set_group_order_via_legendre();
      lower = upper = Cq.group_order();
      group_order_known=1;
      if(debug_iso_type)
	cout<<"Lower and upper bounds on group order adjusted to actual order "
	    <<lower<<" since prime field size <100 or =181, 331, 547"<<endl;
    }

  pointmodq P(Cq), Q(Cq), Q1(Cq);
  bigint ordP, ordP2, ordQ;
  
  P = Cq.random_point();
  if(debug_iso_type) cout<<"P = "<<P<<":\t"<<flush;
  if(group_order_known)  ordP = my_order_point(P,lower);
  else                   ordP = my_order_point(P,lower,upper);
  if(debug_iso_type) cout<<"Order(P) = "<< ordP <<endl;

  vector<bigint> quotlist=n2list(ordP,lower,upper,q);
  if(debug_iso_type) 
    cout<<"Possible n2 values if n1 = "<<ordP<<": "<<quotlist<<endl;

  int n;
  for(n=1; ((n<=10)&&(2*ordP<=upper)) || quotlist.size()!=1;  n++)
    { 
      Q = Cq.random_point();
      if(debug_iso_type>1)  cout<<"Q = "<<Q<<":\t"<<flush;
      merge_points_1(P,ordP,Q);
      quotlist = n2list(ordP,lower,upper,q);
      if(debug_iso_type>1)  
	{
	  cout<<"now P = "<<P<<":\tof order "<<ordP<<endl;
	  cout<<"possible n2 values: "<<quotlist<<endl;
	}
    }

  // At this stage, P is likely to be the first generator; if so,
  // there is a second independent generator of order quot which we
  // must now find.  If not, then while looking for the second
  // generator we will come across a point whose order is not a
  // divisor of ordP, at which point we update P, ordP and quot.

  P1 = P; n1 = ordP; 
  P2 = pointmodq(Cq); n2 = 1; 
  bigint quot=quotlist[0]; // the only value in the list
  if(quot==1) // then there is no ambiguity (thanks to the special
	      // cases dealt with where the grou order is computed in
	      // advance)
    {
      if(debug_iso_type)
	cout<<"group is cyclic, generated by P = "<<P<<endl;
      return;
    }
  if(debug_iso_type)
    {
      cout<<"Maximal order found after using "<<n<<" random points is "
	  <<ordP<<endl;
      cout<<"with n2 = "<<quot<<endl;
      cout<<"Assuming that P does have maximal order,\n";
      cout<<"group structure is "<<ordP<<"*"<<quot<<"="<<(quot*ordP)<<endl;
    }
  ffmodq dummy(Cq);  // to initialize the function field's static data
  if(debug_iso_type) 
    cout<<"Looking for a second generator of order "<<quot<<endl;
  
  if(even(quot))
    {
      // We find the 2-torsion explicitly: this is better than
      // using random points
      pointmodq T = (ordP/2)*P1; // the one we have already
      if(debug_iso_type)
	cout<<"Existing 2-torsion point "<<T<<endl;	
      NewFqPoly(Fq,f);  FqPolyAssignX(f);  f=f-T.get_x();
      FqPoly x2divpol = makepdivpol(Cq,2);
      if(debug_iso_type)
	cout<<"2-division poly = "<<x2divpol<<endl;
      divide(x2divpol,x2divpol,f);
      if(debug_iso_type)
	cout<<"reduced 2-division poly = "<<x2divpol<<endl;
      // Now x2divpos will have as roots the two other x-coords of
      // 2-division points _unless_ P1 is not of maximal order in
      // which case it may be irreducible.
      gf_element a1=PolyCoeff(x2divpol,1);
      gf_element d =a1*a1-ItoGF(Fq,16)*PolyCoeff(x2divpol,0), rd;
      if(sqrt(Fq,d,rd))
	{
	  P2.set_x_coordinate((sqrt(Fq,d)-a1)/ItoGF(Fq,8));
	  if(debug_iso_type)
	    cout<<"New 2-torsion point "<<P2<<endl;	
	  n2=2;
	}
      else
	{
	  if(debug_iso_type)
	    cout<<"No more 2-torsion points! P1 must not be of maximal order"
		<<endl;	 
	}
    }
  
  while(n2!=quot)
    { 
      Q = Cq.random_point();
      if(debug_iso_type) cout<<"Using Q = "<<Q<<endl;
      merge_points_2(P1,n1,P2,n2,quot,Q);
      if(debug_iso_type) cout<<"Group order now "<<(n1*n2)<<endl;     
    }
  if(debug_iso_type) 
    cout<<"Finished: group structure ("<<n1<<","<<n2<<") with "
	<<"generators "<<P1<<","<<P2<<endl;
  return;
}

// find full Z-basis
void my_isomorphism_type_new(curvemodq& Cq, 
			 bigint& n1, bigint& n2, pointmodq& P1, pointmodq& P2)
{
  galois_field Fq = get_field(Cq);
  bigint q = Fq.characteristic();
  bigint upper, lower; // bounds on group order
  set_hasse_bounds(q,lower,upper);
  if(debug_iso_type)
    cout<<"Lower and upper bounds on group order: ["
	<<lower<<","<<upper<<"]"<<endl;

  int group_order_known=0;
  if((q<100)||(q==181)||(q==331)||(q==547))
    {
      Cq.set_group_order_via_legendre();
      lower = upper = Cq.group_order();
      group_order_known=1;
      if(debug_iso_type)
	cout<<"Lower and upper bounds on group order adjusted to actual order "
	    <<lower<<" since prime field size <100 or =181, 331, 547"<<endl;
    }

  pointmodq P(Cq), Q(Cq), Q1(Cq), P3(Cq);
  P1=P; P2=P; n1=n2=BIGINT(1);
  bigint m,a, ordQ, oldn1;

  // At all stages the current subgroup is generated by P1, P2 with
  // orders n1,n2 which are disjoint.  We stop when n1*n2 >= lower

  int count=0;

  while(n1*n2<lower)
    {
      count++;
      Q = Cq.random_point();
      if(debug_iso_type) cout<<"Q = "<<Q<<":\t"<<flush;
      if(group_order_known)  ordQ = my_order_point(Q,lower);
      else                   ordQ = my_order_point(Q,lower,upper);
      ordQ = order_point(Q);
      if(debug_iso_type) cout<<"Order(Q) = "<< ordQ <<endl;
 
      Q1=n1*Q;

      if(Q1.is_zero()) // P1,n1 will not change but we may increase n2
	{
	  if(debug_iso_type) 
	    cout<<"Case 1: n2 may increase"<<endl;
	  //	  if(count>5)
	    {
	  m = linear_relation(P1,Q,a);
	  if(debug_iso_type) 
	    cout<<"linear relation gives m="<<m<<", a="<<a<<endl;
	  if(m>1)
	    {
	      Q=Q-(a/m)*P1; // has order m and is disjoint from P1
	      set_order_point(Q,m);
	      if(n2==1)
		{
		  P2=Q; n2=m;
		  if(debug_iso_type)
		    {
		      cout<<"Adding second generator "<<P2<<" of order "
			    <<n2<<endl
			    <<"Group order now "<<n1*n2<<"="<<n1<<"*"<<n2<<endl;
		    }
		}
	      else
		{
		  // Now we must merge P2 and Q:
		  a=n2; // holds old value
		  merge_points_1(P2,n2,Q);
		  if(debug_iso_type)
		    {
		      if(n2>a) 
			cout<<"Replacing second generator by "<<P2<<" of order "
			    <<n2<<", gaining index "<<n2/a<<endl
			    <<"Group order now "<<n1*n2<<"="<<n1<<"*"<<n2<<endl;
		    }
		}
	    }
	  }
	}
      else // Q1 nonzero: n1 will increase
	{
	  if(debug_iso_type) 
	    cout<<"Case 2: n1 may increase"<<endl;
	  oldn1=n1;
	  if(n2>1)
	    {
	      P3=(n1/n2)*P1;  // so P2,P3 are a basis for n2-torsion
	      set_order_point(P3,n2);
	      if(debug_iso_type)
		cout<<"storing generator "<<P3<<" of "<<n2<<"-torsion"<<endl;
	    }
	  merge_points_1(P1,n1,Q);
	  if(debug_iso_type)
	    {
	      cout<<"Replacing first  generator by "<<P1<<" of order "
		  <<n1<<", gaining index "<<n1/oldn1<<endl
		  <<"Group order now "<<n1*n2<<"="<<n1<<"*"<<n2<<endl;
	    }
	  // Now replace P2 by a point of order n2 s.t. it and
	  // (n1/n2)*P1 are still abasis for n2-torsion:
	  if(n2>1)
	    {
	      m = linear_relation(P1,P3,a);
	      if(debug_iso_type) 
		cout<<"linear relation gives m="<<m<<", a="<<a<<endl;
	      P3=P3-(a/m)*P1;
	      set_order_point(P3,m);
	      if(debug_iso_type) 
		cout<<"First  P2 component ="<<P3<<endl;
	      if(m==n2) {P2=P3;} else
		{
		  m = linear_relation(P1,P2,a);
		  if(debug_iso_type) 
		    cout<<"linear relation gives m="<<m<<", a="<<a<<endl;
		  P2=P2-(a/m)*P1;
		  set_order_point(P2,m);
		  if(debug_iso_type) 
		    cout<<"Second  P2 component ="<<P2<<endl;
		  merge_points_1(P2,m,P3);
		  if(debug_iso_type) 
		    cout<<"Combined P2 component ="<<P2<<endl;
		}
	    }
	}
      if(debug_iso_type)
	{
	  if(order_point(P1)!=n1)
	    cerr<<"Generator P1 = "<<P1<<" has order "<<order_point(P1)
		<<" and not "<<n1<<endl;
	  if(order_point(P2)!=n2)
	    cerr<<"Generator P2 = "<<P2<<" has order "<<order_point(P2)
		<<" and not "<<n2<<endl;
	  if(linear_relation(P1,P2,a)!=n2)
	    cerr<<"Generators not independent!"<<endl;
	  cout<<"Generators: P1 = "<<P1<<" of order "
	      <<n1<<", P2 = "<<P2<<" of order "<<n2<<endl;
	  cout<<"Group order now "<<n1*n2<<"="<<n1<<"*"<<n2<<endl;
	}
    }
}
