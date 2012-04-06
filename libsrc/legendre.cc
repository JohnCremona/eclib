// legendre.cc: implementations of functions for solving legendre equations
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
#include "mmatrix.h"
#include "conic.h"
#include "legendre.h"
#include "illl.h"

//#define DEBUG_MINV
//#define DEBUG_LEM1
//#define WATCH_REDUCTION  
//#define DEBUG_LEGENDRE
//#define DEBUG_LEGENDRE_PARAM
//#define CHECK_ALL
//#define DEBUG_LLL
//#define MORDELL_REDUCE // else use JC's reduction via quadratics
#ifndef TRACE_HOLZER
#define TRACE_HOLZER 0 // =1 for verbose reduction of solutions
#endif
//#define CHECK_CERTS  // define this to check all certificates produced
//#define CHECK_LATTICE // define this to check that solutions lie in correct lattice
//#define CHECK_INDEX
//#define CHECK_CLAIMS
//#define REDUCE_INTERMEDIATES // reduces intermediate solutions
//#define HOLZER_MEASURES // shows Holzer measure of solutions before/after reduction

//#define TL2(nnn) cout<<"Calling lem2() at point "<<nnn<<endl;
#define TL2(nnn) 
#ifdef DEBUG_LEGENDRE
//#if(1)
#define BACK(nnn,uuu,code) cout<<"Backtracking at point "<<nnn<<" with u = "<<uuu<<" and code "<<code<<endl;
#else
#define BACK(nnn,uuu,code) 
#endif

#ifdef CHECK_LATTICE
#define CHECK_LEG(a,b,c,n,p,q,x,y,z) check_leg(a,b,c,n,p,q,x,y,z)
#else
#define CHECK_LEG(a,b,c,n,p,q,x,y,z) check_leg(a,b,c,x,y,z)
#endif

//#define CHECK_LEM3
#define USE_SMALL_P

void lem3(const bigint& a, const bigint& b,
	  bigint& m1, bigint& m2, bigint& m3, bigint& c1, bigint& c2);


void minv(const bigint& a1, const bigint& a2, 
	  const bigint& b1, const bigint& b2, 
	  const bigint& c1, const bigint& c2, 
	  bigint& xmin, bigint& ymin)
{
  bigint x=a1, y=a2, xx=b1, yy=b2;
  bigint n1 = c1*sqr(x)  + c2*sqr(y);
  bigint dot= c1*x*xx    + c2*y*yy;
  bigint alpha; 
  nearest(alpha,dot, n1);  // nearest int to quotient dot/n1
  int reduced = is_zero(alpha);
  if(!reduced)
    {
      xx-=alpha*x;
      yy-=alpha*y;
    }
  bigint n2 = c1*sqr(xx) + c2*sqr(yy);
  reduced=(n2>=n1);
  while(!reduced)
    {
      swap(n1,n2);
      swap(x,xx);  
      swap(y,yy);
      
      dot = c1*x*xx + c2*y*yy;
      nearest(alpha,dot,n1);
      reduced = is_zero(alpha);
      if(!reduced)
	{
	  xx-=alpha*x;
	  yy-=alpha*y;
	  n2 = c1*sqr(xx) + c2*sqr(yy);
	}
      reduced=(n2>=n1);
    }
#ifdef CHECK_ALL
#ifdef DEBUG_MINV
  cout<<"minv called with a=["<<a1<<","<<a2<<"], b=["<<b1<<","<<b2
      <<"], c1="<<c1<<", c2="<<c2<<" returns ["<<x<<","<<y<<"]"<<endl;
#endif // DEBUG_MINV
#endif // CHECK_ALL
  xmin=x; ymin=y;
}

// (e+f*i) is (a+b*i) mod (c+d*i) in Z[i]
void GIreduce(const bigint& a, const bigint& b, const bigint& c, 
	      const bigint& d, bigint& e, bigint& f)
{
  bigint n = sqr(c)+sqr(d);
  bigint q1, q2;
  nearest(q1,a*c+b*d,n); 
  nearest(q2,b*c-a*d,n);
  e = a - c*q1 + d*q2;
  f = b - c*q2 - d*q1;
}

// (e+f*i) = gcd ( (a+b*i) , (c+d*i) ) in Z[i]
void GIgcd(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& d, bigint& e, bigint& f)
{
  int cont=1;
  if(is_zero(c)&&is_zero(d)) {e=a;f=b;cont=0;}
  bigint x1=a, y1=b, x2=c, y2=d, x3, y3, g;
  while(cont)
    {
      GIreduce(x1,y1,x2,y2,x3,y3);
      if(is_zero(x3)&&is_zero(y3)) {e=x2; f=y2; cont=0;}
      x1=x2; x2=x3; y1=y2; y2=y3;
    }
  if((e<=0)&&(f> 0)) {g= f; f=-e; e=g; return;}
  if((e< 0)&&(f<=0)) {::negate(e); ::negate(f); return;}
  if((e>=0)&&(f< 0)) {g=-f; f= e; e=g; return;}
  return;
}

// Solve Legendre's equation when ab=1:

void lem1plus(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& n, const bigint& p, const bigint& q, 
	  bigint& x, bigint& y, bigint& z)
{
  static bigint one, zero; one=1; zero=0;
  GIgcd(q,one,c,zero,x,y); z=1;
  if(b==-one) {y=-y;}  // So we have qx=-y, qy=x (mod c)
#ifdef DEBUG_LEM1
  if(!CHECK_LEG(a,b,c,n,p,q,x,y,z)) 
    {
      cout<<"Wrong solution in lem1plus!\n";
      show_all(a,b,c,n,p,q,x,y,z);
    }
#endif // DEBUG_LEM1
  return;
}

// Solve Legendre's equation when ab=-1:
// (Trivial solution is (1,1,0), but we want a solution in the correct lattice, 
//  satisfying y=qx (mod c) where q^2=1 (mod c))
//

void lem1minus(const bigint& a, const bigint& b, const bigint& c, 
	       const bigint& n, const bigint& p, const bigint& q, 
	       bigint& x, bigint& y, bigint& z)
{
#ifdef DEBUG_LEM1
  cout<<"lem1minus: ab=-1, c="<<c<<" and q="<<q<<endl;
#endif
// Easy cases first:
  if(q==1)
    {
      z=0; x=1; y=a; return;
    }
  bigint absc=abs(c);
  if((q==-1)||(q==absc-1))
    {
      z=0; x=1; y=-a; return;
    }
  bigint cplus=gcd(q-1,c);
  bigint cminus=gcd(q+1,c);
  bigint d=cplus*cminus;
  z = d/absc;  // = 1 or 2
  if(d==absc)        {x=(cminus-cplus)/2;}
  else if(d==2*absc) {x=cminus/2-cplus;}
  else cout<<"Error in lem1minus: c="<<c<<", cplus="<<cplus<<
             ", cminus="<<cminus<<endl;
  y=cminus-x;
  if(b*c>0) 
    {d=x; x=y; y=d;}  // swapping x and y when a=-1, b= 1, c>0
                      //                   or  a= 1, b=-1, c<0
  if(a<0) {x=-x;}
#ifdef DEBUG_LEM1
  if((a*sqr(x)+b*sqr(y)+c*sqr(z)==0)&&div(c,a*x-q*y)) {;}  else
    {
      cout<<"Bad solution in lem1minus\n";
      show_all(a,b,c,n,p,q,x,y,z);
    }
#endif // DEBUG_LEM1
  return;
}

// Solve Legendre's equation when |ab|=1 or |bc|=1 or |ca|=1:

void lem1(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& n, const bigint& p, const bigint& q, 
	  bigint& x, bigint& y, bigint& z)
{
#ifdef DEBUG_LEM1
  cout<<"lem1 called with (a,b,c)=("<<a<<","<<b<<","<<c<<")"<<endl;
#endif // DEBUG_LEM1
  const bigint& bc=b*c;
  if(bc==-1) 
    {
      lem1minus(b,c,a,p,q,n,y,z,x);
      return;
    }
  if(bc== 1) 
    {
      lem1plus(b,c,a,p,q,n,y,z,x);
      return;
    }
  
  const bigint& ab=a*b;
  if(ab==-1) 
    {
      lem1minus(a,b,c,n,p,q,x,y,z);
      return;
    }
  if(ab== 1) 
    {
      lem1plus(a,b,c,n,p,q,x,y,z);
      return;
    }
  
  const bigint& ca=c*a;
  if(ca==-1) 
    {
      lem1minus(c,a,b,q,n,p,z,x,y);
      return;
    }
  if(ca== 1) 
    {
      lem1plus(c,a,b,q,n,p,z,x,y);
      return;
    }
  cout<<"lem1 wrongly called with ";show_eqn(a,b,c);
}

// Solve Legendre's equation ax^2+by^2+cz^2=0 using Rusin's reduction, 
// without assuming a,b,c pairwise coprime
// returns 0 if not soluble
int level;

int legendre_solve(const bigint& a, const bigint& b, const bigint& c, 
		   bigint& x, bigint& y, bigint& z, int use_lll)
{
  vector<bigint> factorbase = vector_union(vector_union(pdivs(a),pdivs(b)),pdivs(c));
  return legendre_solve(a,b,c,factorbase,x,y,z,use_lll);
}

int legendre_solve(const bigint& a, const bigint& b, const bigint& c, 
		   const vector<bigint>& factorbase,
		   bigint& x, bigint& y, bigint& z, int use_lll)
{
#ifdef DEBUG_LEGENDRE  
  cout<<"Starting to solve Legendre equation with coeffs "
      <<a<<", "<<b<<", "<<c<<"\n";
  cout<<"Factor Base = "<<factorbase<<endl;
#endif // DEBUG_LEGENDRE
//
// Step 1: Check signs are not all equal
//
  int sa=sign(a), sb=sign(b), sc=sign(c);
  if((sa==0)||(sb==0)||(sc==0)) return 0;
  if((sa==sb)&&(sa==sc)&&(sb==sc)) return 0;

//
// Step 2: Reduce to pairwise coprime coefficients
//
  bigint a1=a, b1=b, c1=c, g, n,p,q, xfac, yfac, zfac;
  xfac=1; yfac=1; zfac=1;
  int ntry=3;
  while(ntry>0) {
    ntry--;
    if((g=gcd(b1,c1))>1) {b1/=g; c1/=g; a1*=g; xfac*=g; ntry=3;}
    if(ntry>0) 
      {
	ntry--;
	if((g=gcd(c1,a1))>1) {c1/=g; a1/=g; b1*=g; yfac*=g; ntry=3;}
      }
    if(ntry>0) 
      {
	ntry--;
	if((g=gcd(a1,b1))>1) {a1/=g; b1/=g; c1*=g; zfac*=g; ntry=3;}
      }
  }  
#ifdef DEBUG_LEGENDRE  
  cout<<"After gcd reduction, new coeffs are "
      <<a1<<", "<<b1<<", "<<c1<<"\n";
  cout<<"scale factors are "
      <<xfac<<", "<<yfac<<", "<<zfac<<"\n";
#endif // DEBUG_LEGENDRE
//
// Step 3: Reduce to square-free coefficients
//         NB Here is the only place where factorization is done! 
  bigint a2,a3,b2,b3,c2,c3;
  vector<bigint> apdivs=factorbase, bpdivs=factorbase, cpdivs=factorbase;
  sqfdecomp(a1,apdivs,a2,a3); yfac*=a3; zfac*=a3;
  sqfdecomp(b1,bpdivs,b2,b3); xfac*=b3; zfac*=b3; 
  sqfdecomp(c1,cpdivs,c2,c3); xfac*=c3; yfac*=c3; 
  cancel1(xfac,yfac,zfac);
#ifdef DEBUG_LEGENDRE  
  cout<<"After squarefree reduction, new coeffs are "
      <<a2<<", "<<b2<<", "<<c2<<"\n";
  cout<<"scale factors are "
      <<xfac<<", "<<yfac<<", "<<zfac<<"\n";
#endif // DEBUG_LEGENDRE
// Note: a1=a2*a3^2, and apdivs contains all primes dividing a1, but
// any which do not divide a2 are effectively ignored by make_certificate()
//
// Step 4: Make solubility certificate (or detect insolubility)
//
  int res=make_certificate(a2,apdivs,b2,bpdivs,c2,cpdivs,n,p,q);
#ifdef DEBUG_LEGENDRE  
  if(res) 
    {
      cout<<"Legendre equation with coeffs "<<a2<<", "<<b2<<", "<<c2
	  <<" is not soluble!\n";
      switch(res) {
      case 1: cout<<"No certificate modulo a\n"; break;
      case 2: cout<<"No certificate modulo b\n"; break;
      case 3: cout<<"No certificate modulo c\n"; break;
      }
    }
  else 
    {
      cout<<"Legendre equation with coeffs "<<a2<<", "<<b2<<", "<<c2 <<" is soluble\n";
      cout<<"Certificate: "<<n<<" "<<p<<" "<<q<<endl;
    }
#endif // DEBUG_LEGENDRE
  if(res) return 0;
//
// Step 5: Find the solution using one of two methods
//
  if(use_lll)
    {
      legendre_via_lll(a2,b2,c2,n,p,q,x,y,z);
    }
  else // use JC+Rusin's factorization-free reduction
    {
      bigint u;
#ifdef DEBUG_LEGENDRE  
      level=0;
#endif // DEBUG_LEGENDRE
      int res = legendre_solve_cert_1(a2,b2,c2,n,p,q,x,y,z,u);
#ifdef DEBUG_LEGENDRE
      cout<<"Result code from level "<<(level+1)<<" = "<<res<<endl;
#endif
      if(res)
	{
	  cout<<"Problem: at top level, legendre_solve_cert returns ";
	  cout<<"nonzero code "<<res<<" and u = "<<u<<endl;
	}
#ifdef HOLZER_MEASURES
      cout<<"Before reduction of solution "; 
      cout<<"Holzer measure = "<<holzer_measure(a2,b2,c2,x,y,z)<<endl;
#endif // HOLZER_MEASURES
#ifdef DEBUG_LEGENDRE  
#ifdef CHECK_ALL
      if(!CHECK_LEG(a2,b2,c2,n,p,q,x,y,z)) cout<<" wrong solution!\n"; 
#endif // CHECK_ALL
#endif // DEBUG_LEGENDRE
//
// Step 6: Reduce the solution using one of two methods
//         (Not done if LLL method was used)
// NB this will in general move the solution off the certificate lattice
      cancel1(x,y,z);
#ifdef MORDELL_REDUCE
      legendre_reduce(a2,b2,c2,x,y,z,TRACE_HOLZER);
#else
      new_legendre_reduce(a2,b2,c2,x,y,z,TRACE_HOLZER);
#endif // MORDELL_REDUCE
#ifdef HOLZER_MEASURES
      cout<<"After reduction: ";
      cout<<"Holzer measure = "<<holzer_measure(a2,b2,c2,x,y,z)<<endl;
#endif // HOLZER_MEASURES
    }
//
// Step 7: Scale up the solution (see steps 2 and 3)
//
  x*=xfac; y*=yfac; z*=zfac;
  if(is_negative(x)) ::negate(x);
  if(is_negative(y)) ::negate(y);
  if(is_negative(z)) ::negate(z);
  cancel1(x,y,z);
//
// Step 8 (optional): Check the solution
//
#ifdef CHECK_ALL
#ifdef DEBUG_LEGENDRE  
  cout<<"Checking solution "; show_xyz(x,y,z); cout<<endl;
#endif // DEBUG_LEGENDRE
  if(!check_leg(a,b,c,x,y,z)) cout<<" wrong solution!\n"; 
#endif // CHECK_ALL
  return 1;
} // end of legendre_solve()

// Solve Legendre's equation ax^2+by^2+cz^2=0 using Rusin's reduction, 
// given "certificate" (n,p,q) satisfying a|n^2+bc, b|p^2+ac, c|q^2+ab.
//
// All this does is sort the coeffs and pass to lem4()
//
void legendre_solve_cert(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& n, const bigint& p, const bigint& q, 
		    bigint& x, bigint& y, bigint& z)
{
#ifdef DEBUG_LEGENDRE  
    cout<<"Solving Legendre's equation with coefficients ";
    show_eqn(a,b,c);
    cout<<"\n  and  ";
    show_cert(n,p,q);
    cout<<" using factorization-free reduction\n";
#endif // DEBUG_LEGENDRE
  x=y=z=0; 
#ifdef CHECK_CERTS
  if(!checkin(a,b,c,n,p,q)) return;
#endif // CHECK_CERTS

// Check if any coeff is a square (up to sign):
    bigint absa=abs(a), absb=abs(b), absc=abs(c), aa, bb, cc;
//
    if(isqrt(absb,bb))  // then |b| = bb^2 with bb>1, so
      {               // we can reduce b with lem2
	if(bb>1)
	  {
	    TL2(1)
	    lem2b(a,b,c,n,p,q,bb,x,y,z); return;
	  }
      }
    
// Now |b| is not a square

    if(isqrt(absc,cc)) 
      {
	if(cc>1)    // then |c| = cc^2 and we can reduce c with lem2
	  {
	    TL2(2)
	      lem2c(a,b,c,n,p,q,cc,x,y,z); return;
	  }
      }    
// Now |b|, |c| are not squares

    if(isqrt(absa,aa)) 
      {
	if(aa>1)    // then |a| = aa^2 and we can reduce a with lem2
	  {
	    TL2(3)
	      lem2a(a,b,c,n,p,q,aa,x,y,z); return;
	  }
      }    
// Now |a|, |b|, |c| are all not squares

// Sort so that |a|>=|b|>=|c|: 

  int perm;
  if(absa>absb)
    {
      if(absc>absa) perm=4; else {if(absb>absc) perm=0; else perm=1;}
    }
  else // absb>absa
    {
      if(absa>absc) perm=2; else {if(absc>absb) perm=5; else perm=3;}
    }

  switch(perm) {
  case 0: lem4(a,b,c,n,p,q,x,y,z); break;
  case 1: lem4(a,c,b,-n,-q,-p,x,z,y); break;
  case 2: lem4(b,a,c,-p,-n,-q,y,x,z); break;
  case 3: lem4(b,c,a,p,q,n,y,z,x); break;
  case 4: lem4(c,a,b,q,n,p,z,x,y); break;
  case 5: lem4(c,b,a,-q,-p,-n,z,y,x); break;
  }
#ifdef CHECK_ALL
  CHECK_LEG(a,b,c,n,p,q,x,y,z);
#endif // CHECK_ALL
}

static int permtable[6][4] = {{0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1}};

// Solve Legendre's equation ax^2+by^2+cz^2=0 using Rusin's reduction, 
// given "certificate" (n,p,q) satisfying a|n^2+bc, b|p^2+ac, c|q^2+ab.
//
// All this does is sort the coeffs and pass to lem4_1()
//
int legendre_solve_cert_1(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& n, const bigint& p, const bigint& q, 
		    bigint& x, bigint& y, bigint& z, bigint& u)
{
#ifdef DEBUG_LEGENDRE  
  level++;
  cout<<"\nLevel "<<level<<"\nCoefficients ";
    show_eqn(a,b,c);    cout<<"\n";
    show_cert(n,p,q);   cout<<"\n";
#endif // DEBUG_LEGENDRE
  x=y=z=0; 
#ifdef CHECK_CERTS
  if(!checkin(a,b,c,n,p,q)) return -1;
#endif // CHECK_CERTS

// Check if any coeff is a square (up to sign):
    bigint absa=abs(a), absb=abs(b), absc=abs(c);
//
    if(isqrt(absb,u))  // then |b| = u^2 with u>1, so we return
      {
	if(u>1) return 2;
      }
    
// Now |b| is not a square

    if(isqrt(absc,u))    // then |c| = u^2 with u>1 and we return
      {
	if(u>1) return 3;
      }    
// Now |b|, |c| are not squares

    if(isqrt(absa,u))     // then |a| = u^2 with u>1 and we return
      {
	if(u>1) return 1;
      }    
// Now |a|, |b|, |c| are all not squares

// Sort so that |a|>=|b|>=|c|: 

  int perm, res=-1, newres;
  if(absa>absb)
    {
      if(absc>absa) perm=4; else {if(absb>absc) perm=0; else perm=1;}
    }
  else // absb>absa
    {
      if(absa>absc) perm=2; else {if(absc>absb) perm=5; else perm=3;}
    }
  switch(perm) {
  case 0: res = lem4_1(a,b,c,n,p,q,x,y,z,u); break;
  case 1: res = lem4_1(a,c,b,-n,-q,-p,x,z,y,u); break;
  case 2: res = lem4_1(b,a,c,-p,-n,-q,y,x,z,u); break;
  case 3: res = lem4_1(b,c,a,p,q,n,y,z,x,u); break;
  case 4: res = lem4_1(c,a,b,q,n,p,z,x,y,u); break;
  case 5: res = lem4_1(c,b,a,-q,-p,-n,z,y,x,u); break;
  }
  if(res!=-1) {newres=permtable[perm][res];}
  else newres=res;
#ifdef DEBUG_LEGENDRE  
  if(newres!=res)
    cout<<"Permutation "<<perm<<" changes result code from "<<res<<" to "<<newres<<endl;
  level--;
#endif
  return newres; 
}

//Ensure that input is valid
int checkin(const bigint& a,const bigint& b,const bigint& c,
	    const bigint& n,const bigint& p,const bigint& q)
{
  int sa=sign(a), sb=sign(b), sc=sign(c);
  if((sa==0)||(sb==0)||(sc==0))
    {cout<<"checkin() error: coefficients all zero!"<<endl; return 0;}
  if((sa==sb)&&(sa==sc)&&(sb==sc))
    {cout<<"Input error: coefficients have same sign!"<<endl; return 0;}

  if(gcd(a,b)>1) 
    {cout<<"Input error: a and b not coprime!"<<endl; return 0;}
  if(gcd(b,c)>1) 
    {cout<<"Input error: b and c not coprime!"<<endl; return 0;}
  if(gcd(c,a)>1) 
    {cout<<"Input error: c and a not coprime!"<<endl; return 0;}

  if(ndiv(a,sqr(n)+(b*c))) 
    {cout<<"Input error: bad certificate for a"<<endl; return 0;}
  if(ndiv(b,sqr(p)+(a*c))) 
    {cout<<"Input error: bad certificate for b"<<endl; return 0;}
  if(ndiv(c,sqr(q)+(a*b))) 
    {cout<<"Input error: bad certificate for c"<<endl; return 0;}
  return 1;
}

// Check that purported solution is OK
int check_leg(const bigint& a, const bigint& b, const bigint& c,
	      bigint& x, bigint& y, bigint& z)
{
#ifdef DEBUG_LEGENDRE  
  cout<<"Checking solution "; show_xyz(x,y,z); 
  cout<<" for (a,b,c) = ("<<a<<","<<b<<","<<c<<"): ";
#endif // DEBUG_LEGENDRE
  bigint rhs=a*sqr(x)+b*sqr(y)+c*sqr(z);
  if(is_zero(rhs)) 
    {
#ifdef DEBUG_LEGENDRE  
      cout<<"OK"<<endl; 
#endif // DEBUG_LEGENDRE
      return 1;
    }
  else
    {
#ifndef DEBUG_LEGENDRE  
      cout<<"Checking solution "; show_xyz(x,y,z); 
      cout<<" for (a,b,c) = ("<<a<<","<<b<<","<<c<<"): ";
#endif // DEBUG_LEGENDRE
      cout<<"wrong!!"<<endl; 
      return 0;
    }
}

// Check that purported solution is OK & in correct lattice
int check_leg(const bigint& a, const bigint& b, const bigint& c,
	      const bigint& n, const bigint& p, const bigint& q, 
	      bigint& x, bigint& y, bigint& z)
{
  if(check_leg(a,b,c,x,y,z))
    {
#ifdef DEBUG_LEGENDRE  
      cout<<"Checking lattice congruences of solution for ";
      show_eqn(a,b,c);
      cout  <<"\n with ";
      show_cert(n,p,q);
      cout<<endl;
#endif // DEBUG_LEGENDRE
      int ok=1;
      if(!div(a,b*y-n*z)) {ok=0;cout<<"Lattice congruence mod a fails!\n";}
      if(!div(b,c*z-p*x)) {ok=0;cout<<"Lattice congruence mod b fails!\n";}
      if(!div(c,a*x-q*y)) {ok=0;cout<<"Lattice congruence mod c fails!\n";}
#ifdef DEBUG_LEGENDRE  
      if(ok) cout<<"Lattice congruences OK\n";
#endif // DEBUG_LEGENDRE
      return ok;
    }
  else return 0;
}

void legendre_reduce(const bigint& a, const bigint& b, const bigint& c, 
		     bigint& x0, bigint& y0, bigint& z0, int verb)
     // Given a, b, c,  and ax^2+by^2+cz^2=0
     // reduces x, y, z in place using Mordell's method (page 48)
     // to achieve Holzer's bounds |z|<=sqrt(ab) etc.
  // (just permutes & passes to conic_mordell_reduce() 
  // which assumes a>0, b>0, c<0)
  // we may assume a, b, c do not all have the same sign
{
#ifdef WATCH_REDUCTION  
  verb=1;
#endif
  if(is_positive(a))
    {
      if(is_positive(b))
	{
	  conic_mordell_reduce(a,b,c,x0,y0,z0,verb);
	}
      else
	{
	  if(is_positive(c))
	    {
	      conic_mordell_reduce(a,c,b,x0,z0,y0,verb);
	    }
	  else
	    {
	      conic_mordell_reduce(-c,-b,-a,z0,y0,x0,verb);      
	    }
	}
    }
  else
    {
      if(is_positive(b))
	{
	  if(is_positive(c))
	    {
	      conic_mordell_reduce(b,c,a,y0,z0,x0,verb);
	    }
	  else
	    {
	      conic_mordell_reduce(-a,-c,-b,x0,z0,y0,verb);      
	    }
	}
      else
	{
	  conic_mordell_reduce(-a,-b,-c,x0,y0,z0,verb);      
	}
    }
}

// Finds a certificate or returns 0 if none exists:
int make_certificate(const bigint& a, const bigint& b, const bigint& c, 
		     bigint& n, bigint& p, bigint& q)
{
  if(!sqrt_mod_m(n,-b*c,abs(a))) return 1;
  if(!sqrt_mod_m(p,-a*c,abs(b))) return 2;
  if(!sqrt_mod_m(q,-a*b,abs(c))) return 3;
  return 0;
}

int make_certificate(const bigint& a, const vector<bigint>& apdivs, 
		     const bigint& b, const vector<bigint>& bpdivs, 
		     const bigint& c, const vector<bigint>& cpdivs, 
		     bigint& n, bigint& p, bigint& q)
{
  if(!sqrt_mod_m(n,-b*c,abs(a),apdivs)) return 1;
  if(!sqrt_mod_m(p,-a*c,abs(b),bpdivs)) return 2;
  if(!sqrt_mod_m(q,-a*b,abs(c),cpdivs)) return 3;
  return 0;
}

// Check to see if  b is congruent to  +- c  mod a  (assumed positive!)
//  if not, how much to divide  a  by to ensure that congruence holds?
bigint should(const bigint& a, const bigint& b, const bigint& c)
{
  bigint u=gcd(a,b-c);
  bigint v=gcd(a,b+c);
  if(u>v) return a/u; else return a/v;
}

// minv finds shortest vector [xmin,ymin] in lattice generated by 
// integer vectors a=[a1,a2], b=[b1,b2] 
// with respect to norm c1*x^2+c2*y^2 (c1, c2>0)
// Throughout: [x,y] has norm n1
//             [xx,yy] has norm n2 & inner prod = dot.

void legendre_via_lll(const bigint& a, const bigint& b, const bigint& c, 
		      const bigint& k1, const bigint& k2, const bigint& k3, 
		      bigint& x, bigint& y, bigint& z)
{
  int i;
  bigint g,u,v,adash,bdash,bc=b*c, alpha, beta, gamma; 
  g = bezout(b,c,u,v);
  if(g!=1) {cout<<"b and c not coprime!\n"; return;}
  g = bezout(a,bc,adash,bdash); 
  if(g!=1) {cout<<"a and b*c not coprime!\n"; return;}
  alpha = (c*bdash*k1)   % a;
  beta  = (u*adash*b*k3) % bc;
  gamma = (v*adash*c*k2) % bc;

  vec_m * vecs = new vec_m[4];
  for(i=0; i<=3; i++) vecs[i] = vec_m(3);
  vecs[0][1] = abs(a); vecs[0][2] = abs(b); vecs[0][3] = abs(c);
  vecs[1][1] = bc;               vecs[1][2] = 0;      vecs[1][3] = 0; 
  vecs[2][1] = a*beta;           vecs[2][2] = a;      vecs[2][3] = 0; 
  vecs[3][1] = alpha*beta+gamma; vecs[3][2] = alpha;  vecs[3][3] = 1; 
  
#ifdef DEBUG_LLL
  cout<<"Basis for lattice L:\n";
  cout<<vecs[1]<<"\n";
  cout<<vecs[2]<<"\n";
  cout<<vecs[3]<<"\n";
#endif // DEBUG_LLL
  
  // Now cut down to the sublattice of index 2:
  int oddn1 = odd(bc);
  int oddn2 = odd((sqr(a*beta)+a*b)/bc);
  int oddn3 = odd((a*sqr(alpha*beta+gamma)+b*sqr(alpha)+c)/(a*bc));
  static bigint two; two=2;

  if(oddn1)
    {
      if(oddn2) vecs[2]-=vecs[1];
      if(oddn3) vecs[3]-=vecs[1];
      vecs[1] *= two;
    }
  else if(oddn2)
    {
      if(oddn3) vecs[3]-=vecs[2];
      vecs[2] *= two;
    }
  else if(oddn3)
    {
      vecs[3] *= two;
    }
  else cout<<"Problem in legendre_via_lll: all vectors are even!\n";
  
  
#ifdef DEBUG_LLL
  cout<<"Basis for lattice L0 before reduction:\n";
  cout<<vecs[1]<<"\n";
  cout<<vecs[2]<<"\n";
  cout<<vecs[3]<<"\n";
#endif // DEBUG_LLL
  lll_reduce(3,vecs);
#ifdef DEBUG_LLL
  cout<<"Lattice basis after reduction:\n";
  cout<<vecs[1]<<"\n";
  cout<<vecs[2]<<"\n";
  cout<<vecs[3]<<"\n";
#endif // DEBUG_LLL
  
  vec_m xyz, b1=vecs[1], b2=vecs[2], b3=vecs[3]; 
  for(i=0; i<13; i++)
    {
      switch(i) {
      case 0: default: xyz=b1; break;
      case 1: xyz=b1-b2; break;
      case 2: xyz=b1+b2; break;
      case 3: xyz=b2; break;
      case 4: xyz=b3; break;
      case 5: xyz=b3+b1; break;
      case 6: xyz=b3-b1; break;
      case 7: xyz=b3+b2; break;
      case 8: xyz=b3-b2; break;
      case 9: xyz=b3+b2-b1; break;
      case 10: xyz=b3-b2+b1; break;
      case 11: xyz=b3-b2-b1; break;
      case 12: xyz=b3+b2+b1; break;
      }
      x=xyz[1]; y=xyz[2]; z=xyz[3];
      bigint fxyz = a*sqr(x)+b*sqr(y)+c*sqr(z);
      if(fxyz==0) 
	{
	  if(i) 
	    {
	      cout<<"Message from legendre_via_lll: \nsolution does "
		  <<"not come from first vector but from case "
		  <<i<<endl;
	      cout<<"LLL basis: \n";
	      cout<<b1<<"\n"<<b2<<"\n"<<b3<<"\n";
	      cout<<"weights: "<<vecs[0]<<endl;
	    }
	  delete [] vecs;
	  return;
	} 
    }
  delete [] vecs;
  cout<<"Problem in legendre_via_lll: no vector gives a solution!"<<endl;
  x=0; y=0; z=0;
}

//
// Given one solution, returns quadratics parametrizing all solutions,
// with discriminants -4bc, -4ac, -4ab.  Here a, b, c are assumed
// pairwise coprime but not square-free, and a>0, b>0, c<0.

void legendre_param(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& x0, const bigint& y0, const bigint& z0, 
		    quadratic& qx, quadratic& qy, quadratic& qz)
{
  bigint a1=a, y1=y0, z1=z0, u, v, e;
  bigint g=abs(gcd(y0,z0));
  static bigint zero; zero=0;
  int scale = (g>1);
  if(scale) {a1/=sqr(g); y1/=g; z1/=g;}
  bigint z12=sqr(z1);
  
  qx.set( a1*x0, 2*b*y1,  -b*x0);
  qy.set(-a1*y1, 2*a1*x0,  b*y1);
  qz.set( a1*z1, zero,     b*z1);

#ifdef DEBUG_LEGENDRE_PARAM
  if(!testparamsol(a1,0,c,-b,qx,qy,qz,0))
    cout<<"Parametric solution wrong at (1)\n";
  else 
    cout<<"Correct parametrization (1):\n"
	<<qx<<"\n"<<qy<<"\n"<<qz<<"\n";
#endif // DEBUG_LEGENDRE_PARAM
  e=bezout(y1,z12,u,v);
  e=(u*x0)%z12;
  unimod m;  // the mij are not used here...

  qx.x_shift(e,m);  qx.set_coeff(1,qx[1]/z1);  qx.set_coeff(2,qx[2]/z12);
  qy.x_shift(e,m);  qy.set_coeff(1,qy[1]/z1);  qy.set_coeff(2,qy[2]/z12);
  qz.x_shift(e,m);  qz.set_coeff(1,qz[1]/z1);  qz.set_coeff(2,qz[2]/z12);
#ifdef DEBUG_LEGENDRE_PARAM
  if(!testparamsol(a1,zero,c,-b,qx,qy,qz,0))
    cout<<"Parametric solution wrong at (2)\n";
  else 
    cout<<"Correct parametrization (2):\n"
	<<qx<<"\n"<<qy<<"\n"<<qz<<"\n";
#endif // DEBUG_LEGENDRE_PARAM

  m.reset();
  if(a*b>0)     {qz.reduce(m);  qx.transform(m);  qy.transform(m);}
  else
    { 
      if(a*c>0)	{qy.reduce(m);	qz.transform(m);  qx.transform(m);}
      else      {qx.reduce(m);  qy.transform(m);  qz.transform(m);}
    }
      
#ifdef DEBUG_LEGENDRE_PARAM
  if(!testparamsol(a1,zero,c,-b,qx,qy,qz,0))
    cout<<"Parametric solution wrong at (3)\n";
  else 
    cout<<"Correct parametrization (3):\n"
	<<qx<<"\n"<<qy<<"\n"<<qz<<"\n";
#endif // DEBUG_LEGENDRE_PARAM
  if(scale) 
    {
      qy.set(g*qy[0], g*qy[1], g*qy[2]);
      qz.set(g*qz[0], g*qz[1], g*qz[2]);
    }
#ifdef DEBUG_LEGENDRE_PARAM
  if(!testparamsol(a,zero,c,-b,qx,qy,qz,0))
    cout<<"Parametric solution wrong at (4)\n";
  else 
    cout<<"Correct parametrization (4):\n"
	<<qx<<"\n"<<qy<<"\n"<<qz<<"\n";
#endif // DEBUG_LEGENDRE_PARAM
}

void new_legendre_reduce(const bigint& a, const bigint& b, const bigint& c, 
			 bigint& x0, bigint& y0, bigint& z0, int verb)
     // Given a, b, c,  ax^2+by^2+cz^2=0
     // reduces x, y, z in place using quadratics
{
#ifdef WATCH_REDUCTION  
  verb=1;
#endif
  //  x0=abs(x0); y0=abs(y0); z0=abs(z0);
  if(verb) 
    {
      cout<<"Reducing solution "; show_xyz(x0,y0,z0);
      cout<<" for (a,b,c) = ("<<a<<","<<b<<","<<c<<")..."<<endl;
      //      cout<<"using quadratic parametrization\n";
    }
  int sa=sign(a), sb=sign(b), sc=sign(c);
  bigint test;
  int which, ok;
  if( ((sa<0)&&(sb>0)&&(sc>0)) || ((sa>0)&&(sb<0)&&(sc<0)) )
    {
      which=1; test=b*c;
      ok=(sqr(x0)<=test);
    }
  else 
    {
      if( ((sb<0)&&(sa>0)&&(sc>0)) || ((sb>0)&&(sa<0)&&(sc<0)) ) 
	{
	  which=2; test=a*c;
	  ok=(sqr(y0)<=test);
	}
      else 
	{
	  which=3; test=a*b;
	  ok=(sqr(z0)<=test);
	}
    }
  
  if(ok)
    {
      if(verb) cout<<"...nothing to do, already reduced\n";
      return;
    }
  
  quadratic qx, qy, qz;
  legendre_param(a,b,c,x0,y0,z0,qx,qy,qz);
  if(verb)
    {
      cout<<"Parametrizing quadratics are\n";
      cout<<qx<<endl;
      cout<<qy<<endl;
      cout<<qz<<endl;
    }
  bigint newx0 = abs(qx[0]);
  bigint newy0 = abs(qy[0]);
  bigint newz0 = abs(qz[0]);
  cancel(newx0,newy0,newz0);
  if(verb)
    {
      cout<<"New solution =    "; show_xyz(newx0,newy0,newz0); cout<<endl;
    }
  switch(which) {
  case 1: ok = (sqr(newx0)<=test); break;
  case 2: ok = (sqr(newy0)<=test); break;
  case 3: ok = (sqr(newz0)<=test); break;
  }
  if(!ok)
    {
      if(verb)
	{
	  cout<<"new_legendre_reduce() produces solution no smaller than old one!\n";
	  cout<<"Calling legendre_reduce()...\n";
	}
      legendre_reduce(a,b,c,newx0,newy0,newz0,verb);
      if(verb)
	{
	  cout<<"New solution =    "; show_xyz(newx0,newy0,newz0); cout<<endl;
	}
    }
  x0=newx0;
  y0=newy0;
  z0=newz0;
}

// Main descent: Either
// Find a square factor of  abc  and return,
// Or, find  bc=+-1 and call lem1
// Or compute a smaller (aa,bb,cc,n1,p1,q1) and call lem4 recursively
//
// This routine assumes |a|>=|b|>=|c| (else use legendre_solve_cert())
//
// NOTATION NOW AGREES WITH PAPER!
// (except aa=a', bb=b', cc=c')
//
int lem4x_1(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& k1, const bigint& k2, const bigint& k3, 
	   bigint& x, bigint& y, bigint& z, bigint& u);

int lem4_1(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& k1, const bigint& k2, const bigint& k3, 
	  bigint& x, bigint& y, bigint& z, bigint& u)
#ifdef DEBUG_LEGENDRE
{
  cout<<"lem4 at level "<<level<<" called with ";show_eqn_cert(a,b,c,k1,k2,k3);
  //  cout<<"At level "<<level<<", log|abc| = "<<log(abs(bigfloat(a*b*c)))<<endl;
  int res = lem4x_1(a,b,c,k1,k2,k3,x,y,z,u);
  cout<<"lem4 at level "<<level<<" returns     "; show_xyz(x,y,z); 
  cout<<endl;
  return res;
}

int lem4x_1(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& k1, const bigint& k2, const bigint& k3, 
	   bigint& x, bigint& y, bigint& z, bigint& u)
#endif // DEBUG_LEGENDRE
{
  bigint aa,bb,cc,gamma,alpha,d1,d2,d,e,w,w1,w2,t,bc,u2;
  bigint w1star,w2star,gammastar,gammafactor;
  bigint x1,y1,z1;
  bigint k1dash, k2dash, k3dash;
  static bigint zero, one; zero=0; one=1;
  int res;
//
// We have |a|>=|b|>=|c| and none are squares.  If |b|=|c|=1 call lem1
//
  if((abs(b)==1)) // &&(abs(c)==1))
    {
      lem1(a,b,c,k1,k2,k3,x,y,z);
      return 0;
    }

  w=(k1*invmod(c,a))%a;
  minv(one,w,zero,a,abs(b),abs(c),w1,w2);
  t=(b*sqr(w1)+c*sqr(w2))/a;
  bc=b*c;
  lem3(t,bc,aa,bb,cc,gamma,alpha);
#ifdef DEBUG_LEGENDRE
  cout<<"w1,w2="<<w1<<","<<w2<<", t="<<t<<", ";
  cout<<"gamma="<<gamma<<endl;
#endif
  //
  if(cc>0) { ::negate(cc); ::negate(aa); ::negate(bb);}
  // cout<<"aa,bb,cc="<<aa<<","<<bb<<","<<cc<<endl;
  // Now t = aa*cc*gamma^2,  (aa=n2*n3 later)
  //    bc = bb*cc*alpha^2.
#ifdef CHECK_CLAIMS
  if( t!=aa*cc*sqr(gamma)) cout<<"Identity 1 fails\n";
  if(bc!=bb*cc*sqr(alpha))  cout<<"Identity 2 fails\n";
#endif  

  // 12 (?) chances to find a square factor of abc follow...
  u=gcd(alpha,b);
  if(u>1) {BACK(4,u,2) return 2;}
  u=gcd(alpha,c);
  if(u>1) {BACK(5,u,3) return 3;}

  // From this point alpha=1, 
  //     t = aa*cc*gamma^2,
  //    bc = bb*cc.
#ifdef CHECK_CLAIMS
  if(alpha!=1) cout<<"alpha="<<alpha<<", not 1\n";
  if(bc!=bb*cc)  cout<<"Identity 3 fails\n";
#endif  

  d1=gcd(cc,c);
  e=gcd(d1,w1);  u=abs(d1/e);
  if(u>1) {BACK(6,u,3) return 3;}

  // From this point, d1|w1, d1|c
#ifdef CHECK_CLAIMS
  if(!div(d1,w1)) cout<<"d1 ndiv w1\n";
  if(!div(d1,c)) cout<<"d1 ndiv c\n";
#endif  

  d2=cc/d1;
  e=gcd(d2,w2);  u=abs(d2/e); 
  if(u>1) {BACK(7,u,2) return 2;}

  // From this point, d2|w2, d2|b, cc=d1*d2
#ifdef CHECK_CLAIMS
  if(!div(d2,w2)) cout<<"d2 ndiv w2\n";
  if(!div(d2,b)) cout<<"d2 ndiv b\n";
  if(cc!=d1*d2) cout<<"Identity 4 fails\n";
#endif  

  d=gcd(b,gamma);
  e=gcd(d,w2); u=abs(d/e);
  if(u>1) {BACK(8,u,2) return 2;}
  if(d>1) {BACK(9,d,2) u=d; return 2;}

  // From this point, (b,gamma)=1
#ifdef CHECK_CLAIMS
  if(gcd(b,gamma)!=1) cout<<"(b,gamma) not coprime\n";
#endif  

  d=gcd(c,gamma);
  e=gcd(d,w1); u=abs(d/e);
  if(u>1) {BACK(10,u,3) return 3;}
  if(d>1) {BACK(11,d,3) u=d; return 3;}

  // From this point, (c,gamma)=1
  // Also: (d1,w2)=(d2,w1)=1
#ifdef CHECK_CLAIMS
  if(gcd(c,gamma)!=1) cout<<"***(c,gamma) not coprime\n";
  if(gcd(d1,w2)!=1) cout<<"***(d1,w2) not coprime\n";
  if(gcd(d2,w1)!=1) cout<<"***(d2,w1) not coprime\n";
#endif  

  bigint n1, n2, n3, cc1, cc2;
  bigint m1=c/d1, m2=b/d2;
  // cout<<"a, gamma, bc = "<<a<<", "<<gamma<<", "<<bc<<endl;

  lem3(a,aa,n1,n2,n3,cc1,cc2);
// So a = n3*n1*cc1^2,
//   aa = n3*n2*cc2^2  
#ifdef CHECK_CLAIMS
  if(!(a==n3*n1*sqr(cc1))) cout<<"***Identity 5 fails\n";
  if(!(aa==n3*n2*sqr(cc2))) cout<<"***Identity 6 fails\n";
#endif
  if(cc1>1) // we have a square factor of a: reduce
    {
      BACK(12,cc1,1) u=cc1; return 1;
    }
  if(cc2>1) // we have a square factor of aa: divide it out
            // and adjust gamma
    {
//    cout<<"cc2="<<cc2<<endl;
      aa/=sqr(cc2); 
      gamma*=cc2;
//    cout<<"a, aa, gamma, bc = "<<a<<", "<<aa<<", " <<gamma<<", "<<bc<<endl;
    }
// Now a = n3*n1
//    aa = n3*n2  (pairwise coprime)
#ifdef CHECK_CLAIMS
  if(a!=n3*n1) cout<<"***Identity 7 fails\n";
  if(aa!=n3*n2) cout<<"***Identity 8 fails\n";
  if((gcd(n2,n3)!=1)||(gcd(n3,n1)!=1)||(gcd(n2,n1)!=1))
    cout<<"*** n1, n2, n3 not pairwise coprime!"<<endl;
#endif

  do // while(res)
    {
      w1star=w1; w2star=w2; gammastar=gamma;
      gammafactor=gcd(w1,gamma);    // =gcd(w2,gamma)
      if(gammafactor>1)
	{
#ifdef DEBUG_LEGENDRE
	  cout<<"**** gammafactor = "<<gammafactor<<" ***"<<endl;
#endif
	  w1star/=gammafactor; w2star/=gammafactor; 
	  gammastar/=gammafactor;
//Hence yet another chance to get a square factor of b or c:
	  e=gcd(d1,w1star);  u=abs(d1/e);
	  if(u>1) {BACK(13,u,3) return 3;}

// From this point, d1|w1star
#ifdef CHECK_CLAIMS
	  if(!div(d1,w1star)) cout<<"d1 ndiv w1star\n";
#endif  

	  e=gcd(d2,w2star);  u=abs(d2/e); 
	  if(u>1) {BACK(14,u,2) return 2;}

// From this point, d2|w2star
#ifdef CHECK_CLAIMS
	  if(!div(d2,w2star)) cout<<"d2 ndiv w2star\n";
#endif  
	}
  //
  // Now (w1star,gammastar)=(w2star,gammastar)=1
  // and (w1star,n2)=(w2star,n2)=1
  // and d1|w1star, d2|w2star.   
  //
#ifdef CHECK_CLAIMS
      //      cout<<"checking coprimality..."<<endl;
      if(gcd(w1,n2)!=1) cout<<"***(w1,n2) not coprime\n";
      if(gcd(w2,n2)!=1) cout<<"***(w2,n2) not coprime\n";
      if(gcd(w1star,n2)!=1) cout<<"***(w1star,n2) not coprime\n";
      if(gcd(w2star,n2)!=1) cout<<"***(w2star,n2) not coprime\n";
      d=gcd(w1star,w2star);
      if(!div(sqr(d),n2*n1*sqr(n3))) {cout<<"***Identity 9 fails\n";
      // cout<<"aa,bb,cc="<<aa<<","<<bb<<","<<cc<<endl;
      }
      if(gcd(d,n2)!=1) cout<<"***Identity 10 fails\n";
#endif
//
//  Make new certificate
//
      bigint inv_w2_mod_n2 = invmod(w2star, n2);
      bigint inv_w1_mod_d2 = invmod(w1star, d2);
      bigint inv_w2_mod_d1 = invmod(w2star, d1);
      bigint aaa=invmod(a*gammastar,bc);
      k1dash=chrem(-b*w1star*inv_w2_mod_n2, -k1 ,n2,n3);
      k2dash=chrem(k3*w1star*aaa, k2*w2star*aaa, m1, m2);
      k3dash=chrem( k2*aa*gammastar*inv_w1_mod_d2, 
                   -k3*aa*gammastar*inv_w2_mod_d1, d2, d1);
// cout<<"k2,k3="<<k2<<", "<<k3<<endl;
// cout<<"w1, w2="<<w1<<", "<<w2<<endl;
// cout<<"aaa="<<aaa<<endl;
// cout<<"m1, m2="<<m1<<", "<<m2<<endl;

  //Descend ...
#ifdef DEBUG_LEGENDRE
      cout<<"New equation constructed by lem4() at level "<<level<<": \n";
      show_eqn_cert(aa,bb,cc,k1dash,k2dash,k3dash);
#endif
#ifdef CHECK_ALL
      if(!checkin(aa,bb,cc,k1dash,k2dash,k3dash))
	cout<<"New certificate is wrong|\n";
#endif // CHECK_ALL
      res = legendre_solve_cert_1(aa,bb,cc,k1dash,k2dash,k3dash,x1,y1,z1,u);
#ifdef DEBUG_LEGENDRE
      cout<<"Result code from level "<<(level+1)<<" = "<<res;
      if(u>1) cout<<" (u="<<u<<")";
      cout<<endl<<endl;
#endif
      switch(res)
	{
	case 0: break; // Solution found!
	case 1: // Square factor of aa=n2*n3 found
	  d = gcd(u,n3); if(d>1) {u=d; return 1;}  // since n3|a
	  // now u^2|n2
	  u2=sqr(u); n2/=u2; aa/=u2; 
	  gamma*=u;
	  break; // we will try again...
	case 2: case 3:
	  d = gcd(u,b);  if(d>1) {u=d; return 2;}
	  d = gcd(u,c);  if(d>1) {u=d; return 3;}
	default: return -1; // should not happen
	}
  // Now res=0 if we have  solution, 
  //        =1 if we have adjusted aa and want to try again.
    }
  while (res);

#ifdef CHECK_ALL
  if(!CHECK_LEG(aa,bb,cc,k1dash,k2dash,k3dash,x1,y1,z1)) 
    {
      cout<<" wrong solution for new eqn!\n"; 
      show_all(aa,bb,cc,k1dash,k2dash,k3dash,x1,y1,z1);
    }
#endif // CHECK_ALL
#ifdef CHECK_INDEX
  bigint det_factor=gamma*aa;
  cout<<"Mapping solution back, reduced det = "<<det_factor<<" = "<<gamma<<"*"<<aa;
  if(abs(gamma)>1) cout<<"\t***!!!***";
  cout<<endl;
#endif

  x = -gamma*aa*x1;               //  x = -gamma*n3*x1;
  y = (m1*(w2/d2)*y1+w1*z1);  //  y = (m1*(w2/d2)*y1+w1*z1)/n2;
  z = (m2*(w1/d1)*y1-w2*z1);  //  z=  (m2*(w1/d1)*y1-w2*z1)/n2;
#ifdef CHECK_ALL
  if(!CHECK_LEG(a,b,c,k1,k2,k3,x,y,z)) 
    {
      cout<<" wrong solution for original eqn!\n"; 
      show_all(a,b,c,k1,k2,k3,x,y,z);
    }
#endif

  if(abs(n2)>1)
    {
#ifdef DEBUG_LEGENDRE
      int divide_ok=1;
      if(!divide_exact(x,n2,x1))
	{
	  cout<<"x="<<x<<" not divisible by n2="<<n2<<endl;
	  cout<<"denominator = "<<n2/gcd(x,n2)<<endl;
	  divide_ok=0;
	}
      if(!divide_exact(y,n2,y1))
	{
	  cout<<"y="<<y<<" not divisible by n2="<<n2<<endl;
	  cout<<"denominator = "<<n2/gcd(y,n2)<<endl;
	  divide_ok=0;
	}
      if(!divide_exact(z,n2,z1))
	{
	  cout<<"z="<<z<<" not divisible by n2="<<n2<<endl;
	  cout<<"denominator = "<<n2/gcd(y,n2)<<endl;
	  divide_ok=0;
	}
      if(divide_ok)
	{	
	  x=x1; y=y1; z=z1;
	  cout<<"Successfully divided out by factor n2="<<n2<<endl;
	  d=gcd(gcd(x,y),z);
	  if(d>1) cout<<"Solution not primitive now: gcd = "<<d<<endl;
	}
      else
	{
	  d=gcd(gcd(x,y),z);
	  x/=d; y/=d; z/=d;
	  cout<<"Dividing out by factor d = "<<d<<endl;
	  cout<<"rather than n2 = "<<n2;
	  cout<<" (denom(d/n2) = "<<n2/gcd(d,n2)<<endl;
	}
#else
      x/=n2; y/=n2; z/=n2;
#endif	  
    }


#ifdef HOLZER_MEASURES
  cout<<"Holzer measure = "<<holzer_measure(a,b,c,x,y,z)<<endl;
#endif

#ifdef REDUCE_INTERMEDIATES
  cout<<"Before reduction of solution "; show_xyz(x,y,z); cout<<endl;
#ifdef WATCH_REDUCTION  
  cout<<"Before reduction of solution "; show_xyz(x,y,z); cout<<endl;
#endif // WATCH_REDUCTION  
  cancel1(x,y,z);
#ifdef MORDELL_REDUCE
  legendre_reduce(a,b,c,x,y,z,TRACE_HOLZER);
#else
  new_legendre_reduce(a,b,c,x,y,z,TRACE_HOLZER);
#endif // MORDELL_REDUCE
#ifdef WATCH_REDUCTION
  cout<<"After reduction: ";show_xyz(x,y,z);cout<<endl;
#ifdef HOLZER_MEASURES
  cout<<"Holzer measure = "<<holzer_measure(a,b,c,x,y,z)<<endl;
#endif
#endif // WATCH_REDUCTION  
#endif // REDUCE_INTERMEDIATES

#ifdef DEBUG_LEGENDRE  
  bigint f1=should(abs(a),k1*z,b*y);
  bigint f2=should(abs(b),k2*x,c*z);
  bigint f3=should(abs(c),k3*y,a*x);
  if( (f1!=1) || (f2!=1) || (f3!=1) )
    cout<<"    Found factor (lattice-deviation); ["
	<<f1<<","<<f2<<","<<f3<<"]"<<endl;
#endif // DEBUG_LEGENDRE
  return 0;
}

bigfloat sqr(const bigfloat& x) {return x*x;}

bigfloat holzer_measure(const bigint& a, const bigint& b, const bigint& c, 
			const bigint& x, const bigint& y, const bigint& z)
// max{|a|x^2,|b|y^2,|c|z^2}/|abc|   ( < 1 for a Holzer-reduced solution)
{
  bigfloat ax2=I2bigfloat(abs(a)*sqr(x)),
    by2=I2bigfloat(abs(b)*sqr(y)),
    cz2=I2bigfloat(abs(c)*sqr(z));
  bigfloat ans=ax2;
  if(ans<by2) ans=by2;
  if(ans<cz2) ans=cz2;
  ans/=I2bigfloat(a*b*c);
  return ans<0? -ans : ans;
}

// Given a, b, lem3 returns m1 etc so that a=c1^2*m1*m3, b=c2^2*m2*m3 
// with m1, m2, m3 pairwise coprime.   At all  times these equations hold, 
// and at each step the product m1*m2*m3 is decreased by a factor d, 
// so the process terminates when the coprimality condition is satisfied. 

// New version: divides out all squares of primes < 20 from a, b at the start.
int nsmallp=8;
static long smallp[8]={2, 3, 5, 7, 11, 13, 17, 19};
static long smallpsq[8]={4, 9, 25, 49, 121, 169, 289, 361};

void lem3(const bigint& a, const bigint& b,
	  bigint& m1, bigint& m2, bigint& m3, bigint& c1, bigint& c2)
{
  m1=a; m2=b; m3=1; c1=1; c2=1;
  if((a==0)||(b==0)) return;  // shouldn't happen
  int i;  long p, psq, r; bigint d,q;
#ifdef USE_SMALL_P
  for(i=0; i<nsmallp; i++)
    {
      p=smallp[i]; psq=smallpsq[i];
      while(::divides(m1,psq,q,r)) {m1=q; c1*=p;}
      while(::divides(m2,psq,q,r)) {m2=q; c2*=p;}
    }
#endif
  int flag12=1, flag13=1, flag23=1;  // flagij=0 if we know mi,mj are coprime
  while(flag12||flag13||flag23)
    {
//    cout<<m1<<", "<<m2<<", "<<m3<<endl;
      if(flag12)
	{
	  d=abs(gcd(m1,m2)); flag12=0;
	  if(d>1) 
	    {
	      m1/=d; m2/=d; m3*=d; 
	      flag13=flag23=1;
	    }
	}
      if(flag13)
	{
	  d=abs(gcd(m1,m3)); flag13=0;
	  if(d>1) 
	    {
	      m1/=d; m2*=d; m3/=d; c1*=d; 
	      flag12=flag23=1;
	    }
	}
      if(flag23)
	{
	  d=abs(gcd(m2,m3)); flag23=0;
	  if(d>1) 
	    {
	      m1*=d; m2/=d; m3/=d; c2*=d; 
	      flag12=flag13=1;
	    }
	}
    }
#ifdef CHECK_LEM3
  if( (a==sqr(c1)*m1*m3) && (b==sqr(c2)*m2*m3) 
      && (gcd(m1,m3)==1) && (gcd(m2,m3)==1) && (gcd(m1,m2)==1) )
    {;}
  else
    {
      cout<<"Error in lem3("<<a<<","<<b<<"), returning\n"
	  <<"c1="<<c1<<", c2="<<c2<<", m1="<<m1<<", m2="<<m2<<", m3="
	  <<m3<<endl;
    }
#endif
}

//Peel off a known square factor  u  from coefficient  a
void lem2a(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   const bigint& u,
	   bigint& x, bigint& y, bigint& z)
{
#ifdef DEBUG_LEGENDRE
  cout<<"lem2a called with (a,b,c)=(" <<a<<","<<b<<","<<c<<"), " 
      <<" and u = "<<u<<endl;
#endif // DEBUG_LEGENDRE
  x=y=z=0;
  bigint u2=sqr(u), a1, r;
  if((!::divides(a,u2,a1,r))||((u2<=1)))
    {
      cout<<"lem2a wrongly called with (a,b,c)=("
	  <<a<<","<<b<<","<<c<<")"<<endl;
      cout<<" and u = "<<u<<endl;
      return;
    }
  bigint n1 = n%a1;
  bigint p1 = (p*invmod(u,b))%b;
  bigint q1 = (q*invmod(u,c))%c;
  legendre_solve_cert(a1,b,c,n1,p1,q1,x,y,z);
  y *= u; z *= u;
//
// NB the lattice congruence mod a is NOT preserved by this scaling
//
  cancel1(x,y,z);
#ifdef CHECK_ALL
  if(!CHECK_LEG(a,b,c,n,p,q,x,y,z)) {cout<<"Wrong solution in lem2a!\n";
				     show_all(a,b,c,n,p,q,x,y,z);}
#endif // CHECK_ALL

#ifdef REDUCE_INTERMEDIATES
#ifdef MORDELL_REDUCE
  legendre_reduce(a,b,c,x,y,z,TRACE_HOLZER);
#else
  new_legendre_reduce(a,b,c,x,y,z,TRACE_HOLZER);
#endif // MORDELL_REDUCE
#endif // REDUCE_INTERMEDIATES
}

//Peel off a known square factor  u  from coefficient  b
void lem2b(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   const bigint& u,
	   bigint& x, bigint& y, bigint& z)
{
#ifdef DEBUG_LEGENDRE
  cout<<"lem2b called with (a,b,c)=("<<a<<","<<b<<","<<c<<"), " 
      <<" and u = "<<u<<endl;
#endif // DEBUG_LEGENDRE
  lem2a(b,c,a,p,q,n,u,y,z,x);
}

//Peel off a known square factor  u  from coefficient  c
void lem2c(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   const bigint& u,
	   bigint& x, bigint& y, bigint& z)
{
#ifdef DEBUG_LEGENDRE
  cout<<"lem2c called with (a,b,c)=("<<a<<","<<b<<","<<c<<"), " 
      <<" and u = "<<u<<endl;
#endif // DEBUG_LEGENDRE
  lem2a(c,a,b,q,n,p,u,z,x,y);
}

// These versions of lem4 are now obsolete

// Main descent: Either
// Find a square factor of  abc  and call lem2,
// Or, find  bc=+-1 and call lem1
// Or compute a smaller (aa,bb,cc,n1,p1,q1) and call lem4 recursively
//
// This routine assumes |a|>=|b|>=|c| (else use legendre_solve_cert())
//
void lem4x(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   bigint& x, bigint& y, bigint& z);

void lem4(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& n, const bigint& p, const bigint& q, 
	  bigint& x, bigint& y, bigint& z)
#ifdef DEBUG_LEGENDRE
{
  cout<<"lem4 called with ";show_eqn_cert(a,b,c,n,p,q);
  lem4x(a,b,c,n,p,q,x,y,z);
  cout<<"lem4 returns "; show_xyz(x,y,z); cout<<endl;
}

void lem4x(const bigint& a, const bigint& b, const bigint& c, 
	   const bigint& n, const bigint& p, const bigint& q, 
	   bigint& x, bigint& y, bigint& z)
#endif // DEBUG_LEGENDRE
{
  bigint aa,bb,cc,alpha,beta,d1,d2,u,d,e,k,k1,k2,t,bc;
  static bigint one, zero; one=1; zero=0;
//
// We have |a|>=|b|>=|c| and none are squares.  If |b|=|c|=1 call lem1
//
  if(abs(b)==1)
    {
      lem1(a,b,c,n,p,q,x,y,z);
      return;
    }

  k=(n*invmod(c,a))%a;
  minv(one,k,zero,a,abs(b),abs(c),k1,k2);
  // cout<<"k1,k2="<<k1<<","<<k2<<endl;
  //
  t=(b*sqr(k1)+c*sqr(k2))/a;
  // cout<<"t="<<t<<endl;

  bc=b*c;
  lem3(t,bc,aa,bb,cc,alpha,beta);
  //
  if(cc>0) { ::negate(cc); ::negate(aa); ::negate(bb);}
  // cout<<"aa,bb,cc="<<aa<<","<<bb<<","<<cc<<endl;
  // Now t = aa*cc*alpha^2,
  //    bc = bb*cc*beta^2.
#ifdef CHECK_CLAIMS
  if( t!=aa*cc*sqr(alpha)) cout<<"Identity 1 fails\n";
  if(bc!=bb*cc*sqr(beta))  cout<<"Identity 2 fails\n";
#endif  

  // 12 (?) chances to find a square factor of abc follow...
  bigint  betab=gcd(beta,b);
  if(betab>1) {TL2(4) lem2b(a,b,c,n,p,q,betab,x,y,z); return;}
  bigint betac=gcd(beta,c);
  if(betac>1) {TL2(5) lem2c(a,b,c,n,p,q,betac,x,y,z); return;}

  // From this point beta=1, 
  //     t = aa*cc*alpha^2,
  //    bc = bb*cc.
#ifdef CHECK_CLAIMS
  if(beta!=1) cout<<"beta="<<beta<<", not 1\n";
  if(bc!=bb*cc)  cout<<"Identity 3 fails\n";
#endif  

  d1=gcd(cc,c);
  e=gcd(d1,k1);  u=abs(d1/e);
  if(u>1) {TL2(6) lem2c(a,b,c,n,p,q,u,x,y,z); return;}

  // From this point, d1|k1, d1|c
#ifdef CHECK_CLAIMS
  if(!div(d1,k1)) cout<<"d1 ndiv k1\n";
  if(!div(d1,c)) cout<<"d1 ndiv c\n";
#endif  

  d2=cc/d1;
  e=gcd(d2,k2);  u=abs(d2/e); 
  if(u>1) {TL2(7) lem2b(a,b,c,n,p,q,u,x,y,z); return;}

  // From this point, d2|k2, d2|b, cc=d1*d2
#ifdef CHECK_CLAIMS
  if(!div(d2,k2)) cout<<"d2 ndiv k2\n";
  if(!div(d2,b)) cout<<"d2 ndiv b\n";
  if(cc!=d1*d2) cout<<"Identity 4 fails\n";
#endif  

  d=gcd(b,alpha);
  e=gcd(d,k2); u=abs(d/e);
  if(u>1) {TL2(8) lem2b(a,b,c,n,p,q,u,x,y,z); return;}
  if(d>1) {TL2(9) lem2b(a,b,c,n,p,q,d,x,y,z); return;}

  // From this point, (b,alpha)=1
#ifdef CHECK_CLAIMS
  if(gcd(b,alpha)!=1) cout<<"(b,alpha) not coprime\n";
#endif  

  d=gcd(c,alpha);
  e=gcd(d,k1); u=abs(d/e);
  if(u>1) {TL2(10) lem2c(a,b,c,n,p,q,u,x,y,z); return;}
  if(d>1) {TL2(11) lem2c(a,b,c,n,p,q,d,x,y,z); return;}

  // From this point, (c,alpha)=1
  // Also: (d1,k2)=(d2,k1)=(k1,alpha)=(k2,alpha)=1
#ifdef CHECK_CLAIMS
  if(gcd(c,alpha)!=1) cout<<"***(c,alpha) not coprime\n";
  if(gcd(d1,k2)!=1) cout<<"***(d1,k2) not coprime\n";
  if(gcd(d2,k1)!=1) cout<<"***(d2,k1) not coprime\n";
#endif  

#define AGREE_PAPER
  bigint aa1, aa2, aa3, cc1, cc2;
  bigint m1=c/d1, m2=b/d2, aaa=invmod(a*alpha,bc);
  // cout<<"a, alpha, bc, aaa = "<<a<<", "<<alpha<<", "<<bc<<", "<<aaa<<endl;
#ifdef AGREE_PAPER // else Rusin's code
  lem3(a,aa,aa3,aa1,aa2,cc1,cc2);
// So a = aa2*aa3*cc1^2,
//   aa = aa2*aa1*cc2^2  
#ifdef CHECK_CLAIMS
  if(!(a==aa2*aa3*sqr(cc1))) cout<<"***Identity 5 fails\n";
  if(!(aa==aa2*aa1*sqr(cc2))) cout<<"***Identity 6 fails\n";
#endif
  if(cc1>1) // we have a square factor of a: reduce
    {
      TL2(12) lem2a(a,b,c,n,p,q,cc1,x,y,z); 
      return;
    }
  if(cc2>1) // we have a square factor of aa: divide it out
            // and adjust alpha
    {
// cout<<"cc2="<<cc2<<endl;
      aa/=sqr(cc2); 
      alpha*=cc2;
      aaa=invmod(a*alpha,bc);
// cout<<"a, aa, alpha, bc, aaa = "<<a<<", "<<aa<<", "<<alpha<<", "<<bc<<", "<<aaa<<endl;
    }
// Now a = aa2*aa3,   aa1=n2, aa2=n3, aa3=n1 in paper notation
//    aa = aa2*aa1  (pairwise coprime)
#ifdef CHECK_CLAIMS
  if(a!=aa2*aa3) cout<<"***Identity 7 fails\n";
  if(aa!=aa2*aa1) cout<<"***Identity 8 fails\n";
  if((gcd(aa1,aa2)!=1)||(gcd(aa2,aa3)!=1)||(gcd(aa1,aa3)!=1))
    cout<<"***aa1, aa2, aa3 not pairwise coprime!"<<endl;
#endif

  d=gcd(k1,alpha);
  if(d>1) // then d also divides k2 so we can divide it out
    {
      k1/=d; k2/=d; alpha/=d;
      aaa=invmod(a*alpha,bc);
    }
  d=gcd(k2,alpha);
  if(d>1) // then d also divides k1 so we can divide it out
    {
      k1/=d; k2/=d; alpha/=d;
      aaa=invmod(a*alpha,bc);
    }
  // Now (k1,alpha)=(k2,alpha)=1
#ifdef CHECK_CLAIMS
  if(gcd(k1,alpha)!=1) cout<<"***(k1,alpha) not coprime\n";
  if(gcd(k2,alpha)!=1) cout<<"***(k2,alpha) not coprime\n";
#endif  

  d=gcd(k1,k2);

#ifdef CHECK_CLAIMS
  if(!div(sqr(d),aa1*aa3*sqr(aa2))) {cout<<"***Identity 9 fails\n";
   // cout<<"aa,bb,cc="<<aa<<","<<bb<<","<<cc<<endl;
   // cout<<"k1,k2="<<k1<<","<<k2<<endl;
   // cout<<"aa1,aa2,aa3="<<aa1<<","<<aa2<<","<<aa3<<endl;
   // cout<<"alpha="<<alpha<<endl;				 
			       }
#endif
  if(d>1)
    {
      e=gcd(d,aa1);
    }
//  From here, (d,aa1)=1; moreover, (k1,aa1)=(k2,aa1)=1
//
#ifdef CHECK_CLAIMS
  if(gcd(k1,aa1)!=1) cout<<"***Identity 10 fails\n";
  if(gcd(k2,aa1)!=1) cout<<"***Identity 11 fails\n";
  if(gcd(d,aa1)!=1) cout<<"***Identity 12 fails\n";
   // cout<<"aa,bb,cc="<<aa<<","<<bb<<","<<cc<<endl;
   // cout<<"k1,k2="<<k1<<","<<k2<<endl;
   // cout<<"aa1,aa2,aa3="<<aa1<<","<<aa2<<","<<aa3<<endl;
   // cout<<"alpha="<<alpha<<endl;
#endif
//
//  Make new certificate
//
  bigint n1=chrem(-b*k1*invmod(k2, aa1), 
		  -n ,aa1,aa2);
// cout<<"p,q="<<p<<", "<<q<<endl;
// cout<<"k1, k2="<<k1<<", "<<k2<<endl;
// cout<<"aaa="<<aaa<<endl;
// cout<<"m1, m2="<<m1<<", "<<m2<<endl;

  bigint p1=chrem(q*k1*aaa, 
		  p*k2*aaa, m1, m2);
  bigint q1=chrem( p*aa*alpha*invmod(k1,d2), 
		  -q*aa*alpha*invmod(k2,d1), d2, d1);

#else // Rusin's code
//Special tests if k1>1;
  if(gcd(k1,aa)>1) 
    {
#ifdef DEBUG_LEGENDRE
      cout<<"    Special K..."<<endl;
#endif // DEBUG_LEGENDRE
      lem3(a,aa,aa3,aa1,aa2,cc1,cc2);
      if(cc1>1) 
	{
#ifdef DEBUG_LEGENDRE
	  cout<<"    Special K...cc1="<<cc1<<">1: square factor of a"<<endl;
#endif // DEBUG_LEGENDRE
	  TL2(16) lem2a(a,b,c,n,p,q,cc1,x,y,z); 
	  return;
	}
      if(cc2>1) 
	{
#ifdef DEBUG_LEGENDRE
	  cout<<"    Special K...cc2>1"<<endl;
#endif // DEBUG_LEGENDRE
	  aa/=sqr(cc2); 
	  alpha*=cc2;
	}
// The reductions after this point should never happen
      d=gcd(k1,alpha);
      if(d>1) {k1/=d;k2/=d;alpha/=d;}
      d=gcd(k1,k2);
      if(d>1) 
	{
	  e=gcd(d,aa3);
	  if(e>1) 
	    {
#ifdef DEBUG_LEGENDRE
	      cout<<"    Special K...first e>1"<<endl;
#endif // DEBUG_LEGENDRE
	      TL2(17) lem2a(a,b,c,n,p,q,e,x,y,z); 
	      return;
	    }
	  e=gcd(d,aa1);
	  if(e>1) 
	    {
#ifdef DEBUG_LEGENDRE
	      cout<<"    Special K...second e>1"<<endl;
#endif // DEBUG_LEGENDRE
	      aa1/=sqr(e); 
	      alpha*=e;
	    }
	}
#ifdef DEBUG_LEGENDRE
      cout<<" ...but no square factors of abc found"<<endl;
#endif // DEBUG_LEGENDRE
    }
  else 
    {
      aa1=aa; aa2=1;
    }

//
//  Make new certificate
//
  bigint n1=chrem(c*k2*invmod(k1, aa1),n ,aa1,aa2);
  bigint p1=chrem(q*k1*invmod(aaa,m1), -p*k2*invmod(aaa,m2), m1, m2);
  bigint q1=chrem(-p*aa*alpha*invmod(k1,d2), 
		  q*aa*alpha*invmod(k2,d1), d2, d1);
#endif // AGREE_PAPER

  //Descend & climb back;
#ifdef CHECK_ALL
  cout<<"New equation constructed by lem4(): \n";
  show_eqn_cert(aa,bb,cc,n1,p1,q1);
  if(!checkin(aa,bb,cc,n1,p1,q1))
    {
      cout<<"New certificate is wrong|\n";
    }
#endif // CHECK_ALL
  bigint x1,y1,z1;
  bigint det_factor=alpha*aa;
  legendre_solve_cert(aa,bb,cc,n1,p1,q1,x1,y1,z1);
#ifdef CHECK_ALL
  if(!CHECK_LEG(aa,bb,cc,n1,p1,q1,x1,y1,z1)) 
    {
      cout<<" wrong solution for new eqn!\n"; 
      show_all(aa,bb,cc,n1,p1,q1,x1,y1,z1);
    }
#endif // CHECK_ALL
#ifdef CHECK_INDEX
  cout<<"Mapping solution back, reduced det = "<<det_factor<<" = "<<alpha<<"*"<<aa;
  if(abs(alpha)>1) cout<<"\t***!!!***";
  cout<<endl;
#endif

#ifdef AGREE_PAPER
  x = -alpha*aa*x1;         //  x = -alpha*aa2*x1;
  y = (m1*(k2/d2)*y1+k1*z1);//  y = (m1*(k2/d2)*y1+k1*z1)/aa1;
  z=  (m2*(k1/d1)*y1-k2*z1);//  z=  (m2*(k1/d1)*y1-k2*z1)/aa1;
// divide_exact(y,aa1,y);
// divide_exact(z,aa1,z);
#ifdef CHECK_ALL
  if(!CHECK_LEG(a,b,c,n,p,q,x,y,z)) 
    {
      cout<<" wrong solution for original eqn!\n"; 
      show_all(a,b,c,n,p,q,x,y,z);
    }
#endif
  d=gcd(gcd(x,y),z);
  if(d>1) 
    {
#ifdef DEBUG_LEGENDRE
      cout<<"Dividing out by factor d="<<d<<", aa1="<<aa1<<endl;
#endif
      x/=d; y/=d; z/=d;
    }
#else // Rusin's code
  x=-det_factor*x1;
  y=m1*(k2/d2)*y1+k1*z1;
  z=m2*(k1/d1)*y1-k2*z1;

  d=gcd(y,z); x1=gcd(d,x);
  //option; try making x,y,z divisible by most of A1;
#ifdef TWEAK_K1
  y2=m1*(k2/d2)*y1-k1*z1;
  z2=-m2*(k1/d1)*y1-k2*z1;
  e=gcd(z2,y2);x2=gcd(e,x);
  if(x2>x1) {y=y2;z=z2;x1=x2;d=e;}
#endif // TWEAK_K1
  if(x1>1) {x/=x1; y/=x1; z/=x1; d/=x1;}
//if d>1, then d^2 | a ...
#endif // AGREE_PAPER

#ifdef WATCH_REDUCTION  
  cout<<"Before reduction of solution "; show_xyz(x,y,z); cout<<endl;
#endif // WATCH_REDUCTION  
#ifdef REDUCE_INTERMEDIATES
#ifdef MORDELL_REDUCE
  legendre_reduce(a,b,c,x,y,z,TRACE_HOLZER);
#else
  new_legendre_reduce(a,b,c,x,y,z,TRACE_HOLZER);
#endif // MORDELL_REDUCE
#endif // REDUCE_INTERMEDIATES
#ifdef WATCH_REDUCTION
  cout<<"After reduction: ";show_xyz(x,y,z);cout<<endl;
#endif // WATCH_REDUCTION  

#ifdef DEBUG_LEGENDRE  
  bigint f1=should(abs(a),n*z,b*y);
  bigint f2=should(abs(b),p*x,c*z);
  bigint f3=should(abs(c),q*y,a*x);
  if( (f1!=1) || (f2!=1) || (f3!=1) )
    cout<<"    Found factor (lattice-deviation); ["
	<<f1<<","<<f2<<","<<f3<<"]"<<endl;
#endif // DEBUG_LEGENDRE
}

