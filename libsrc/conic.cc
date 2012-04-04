// conic.cc: implementations of functions for solving conics
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
#include "conic.h"
#include "legendre.h"

//#define DEBUG_CONIC
//#define DEBUG_CONIC_2
#ifndef TRACE_FACTORIZATION
#define TRACE_FACTORIZATION 0
#endif
//#define DEBUG_REDUCE
#define REDUCE_INTERMEDIATES // reduces intermediate solutions
//#define MORDELL_REDUCE // else use JC's (faster!) reduction via quadratics

// CONIC_METHODS:
//
// 0:  simple recursion, no reduction
// 1:  recursion with algebraic reduction (Denis Simon)
// 2:  recursion with lattice reduction (best non-factn-free)
// 3:  obsolete
// 4:  Uses factorization-free reduction from legendre.cc (best)
// 5:  Uses LLL method from legendre.cc

int solve_conic(const quadratic& q, const bigint& d,
		bigint& x, bigint& y, bigint& z, int method)
{
  return solve_conic(q[0],q[1],q[2],d,x,y,z,method);
}

int solve_conic(const quadratic& q, const bigint& d,
		       const vector<bigint>& factorbase,
		       bigint& x, bigint& y, bigint& z, int method)
{
  return solve_conic(q[0],q[1],q[2],d,factorbase,x,y,z,method);
}

int solve_conic(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                bigint& x, bigint& y, bigint& z, int method)
{
  vector<bigint> factorbase = pdivs(2*d);
  //  cout<<"factorbase(1) = "<<factorbase<<endl;
  if(is_zero(b))
    {
      factorbase=vector_union(factorbase,pdivs(a));
      factorbase=vector_union(factorbase,pdivs(c));
    }
  else
    {
      bigint disc = b*b-4*a*c;
      factorbase=vector_union(factorbase,pdivs(a));
      factorbase=vector_union(factorbase,pdivs(disc));
    }
  //  cout<<"factorbase(2) = "<<factorbase<<endl;
  return solve_conic(a,b,c,d,factorbase,x,y,z,method);
}

int solve_conic(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
		const vector<bigint>& factorbase,
                bigint& x, bigint& y, bigint& z, int method)
     // Solves axx+bxz+czz=dyy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, c, d, bb-4ac non-zero
{
  if(method>3)
    {
      int use_lll=0;
      if(method==5) use_lll=1;
      if(is_zero(b)) return legendre_solve(a,-d,c,factorbase,x,y,z,use_lll);
      bigint aa=sqr(b)-4*a*c, bb=a*d;
      bigint one, zero; one=1; zero=0;
      if(!legendre_solve(one,-aa,-bb,factorbase,x,z,y,use_lll))  return 0;
#ifdef DEBUG_CONIC
      testsol(one,zero,-aa,bb,x,y,z,1);
#endif // DEBUG_CONIC
      x=x-b*z; y*=a; z*=(2*a);
#ifdef DEBUG_CONIC
      testsol(a,b,c,d,x,y,z,1);
#endif // DEBUG_CONIC
      cancel(x,y,z);
      return 1;
    }
  else // method = 0, 1, 2, 3
    {
      int verb=0;
#ifdef DEBUG_CONIC
      verb=1;
#endif // DEBUG_CONIC
//All this function does is complete the square (if nec.) and pass to the next
      bigint aa, bb, xx, yy, zz, a1,a2,b1,b2;
      vector<bigint> aplist, bplist, cplist, dplist;

      int nondiag=!is_zero(b);
      bb=a*d;  
      aa=-a*c; if(nondiag) aa=sqr(b)-4*aa;
      aplist=factorbase;
      bplist=factorbase;
      sqfdecomp(aa,aplist,a1,a2);  // aa=a1*a2^2
      sqfdecomp(bb,bplist,b1,b2);  // bb=b1*b2^2

      if(!testlocsol(a1,aplist,b1,bplist)) return 0;
      if(solve_conic_diag(a1,aplist,b1,bplist,x,y,z,method))
	{
	  conic_diag_reduce(a1,b1,x,y,z,verb);
	  x*=(a2*b2); y*=a2; z*=b2; 
	  if(nondiag) x-=b*z;
	  y *= a;
	  z *= a; if(nondiag) zz<<=1;
	  cancel(x,y,z);
	  return 1;
	}
      else // shouldn't happen as we tested solubility earlier
	{
	  cout << "Problem in solve_conic (parameters (a,b,c,d)=("<<a<<","<<b<<","<<c<<","<<d<<"))\n";
	  cout << "testlocsol() predicted solubility but no solution found!\n";
	  x=0; y=0; z=0; return 0;
	}
    }
}
  
int solve_conic_diag_nontriv(const bigint& a, const vector<bigint>& aplist,
			     const bigint& b, const vector<bigint>& bplist,
			     bigint& x, bigint& y, bigint& z,
			     int method);
     // Solves xx-azz=byy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, b non-zero square-free, their prime divisors in aplist, bplist 
     // Here |b| >= |a|, |b| >=2, a!=1

int solve_conic_diag(const bigint& a, const vector<bigint>& aplist,
		     const bigint& b, const vector<bigint>& bplist,
		     bigint& x, bigint& y, bigint& z,
		     int method)
     // Solves xx-azz=byy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, b non-zero square-free, their prime divisors in aplist, bplist 
     //
     // Here trivial cases are dealt with, non-trivial passed on
{
#ifdef DEBUG_CONIC
  cout << "In solve_conic_diag with a = " << a << ", b = " << b << endl;
#endif // DEBUG_CONIC
  if(is_one(b)) 
    {
      x=1; y=1; z=0; 
#ifdef DEBUG_CONIC
      cout << "...returns ";  show_xyz(x,y,z);
#endif // DEBUG_CONIC
      return 1;
    }
  if(is_one(a)) 
    {
      x=1; y=0; z=1; 
#ifdef DEBUG_CONIC
      cout << "...returns ";  show_xyz(x,y,z);
#endif // DEBUG_CONIC
      return 1;
    }

  if(abs(a)>abs(b)) 
    {
      int res = solve_conic_diag(b,bplist,a,aplist,x,z,y,method);
#ifdef DEBUG_CONIC
      cout << "...returns ";  show_xyz(x,y,z);
#endif // DEBUG_CONIC
      return res;
    }

  // Now |a|<=|b|, neither a nor b = 1

  if(is_one(-b)) 
    {
#ifdef DEBUG_CONIC
      cout << "...returns fail (a=b=-1)\n";
#endif // DEBUG_CONIC
      return 0; // since a must be -1 too and xx+zz=-yy insoluble
    }

  if(b==-a) 
    {
      x=0; y=1; z=1; 
#ifdef DEBUG_CONIC
      cout << "...returns ";  show_xyz(x,y,z);
#endif // DEBUG_CONIC
      return 1;
    }

  if(b==a) 
    {
      bigint m1; m1=-1;
      int res = solve_conic_diag(m1,pdivs(BIGINT(1)),a,aplist,y,x,z,method);
      x*=a;
#ifdef DEBUG_CONIC
      cout << "...returns ";  show_xyz(x,y,z);
#endif // DEBUG_CONIC
      return res;
    }

  // Now |b| >= 2 , |b|>|a| and a is NOT square
  return solve_conic_diag_nontriv(a,aplist,b,bplist,x,y,z,method);
}

int solve_conic_diag_nontriv(const bigint& a, const vector<bigint>& aplist,
			     const bigint& b, const vector<bigint>& bplist,
			     bigint& x, bigint& y, bigint& z,
			     int method)
     // Solves xx-azz=byy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, b non-zero square-free, their prime divisors in aplist, bplist 
     // Here |b| > |a|, |b| >=2, a!=1
{
  if(method<3) {
  bigint x0, t, t0, t1, x1, y1, z1;  
  vector<bigint> tplist;
  int res = modsqrt(a,bplist,x0); // Solves x^2=a mod b, returns success/fail
  if(!res) 
    {
#ifdef DEBUG_CONIC
      cout << "...returns fail (no modular sqrt)\n";
#endif // DEBUG_CONIC
      return 0;
    }
  
  t=b;
  bigint newt = (sqr(x0)-a)/b;
#ifdef DEBUG_CONIC
  cout<<"x0 = "<<x0<<", t = "<<newt<<endl;
#endif // DEBUG_CONIC
  bigint m11, m21, m33, temp, fac;
  m11=1; m33=1;
  switch(method)
    {
    case 0:  // no reduction at all
      t=newt;
      m11=x0; m21=1; m33=t;
      break; // end of method=0
    case 1: // Simon's reduction
#ifdef DEBUG_CONIC_2
      cout<<"t="<<t<<", newt="<<newt<<endl;
#endif // DEBUG_CONIC_2
      while(abs(newt)<abs(t))
	{
	  temp = x0*m11 +  a*m21;
	  m21  =    m11 + x0*m21;
	  m11  = temp;
	  x0=mod(-x0,newt);
	  t=newt;
	  newt=(sqr(x0)-a)/t;
#ifdef DEBUG_CONIC_2
	  cout<<"newt="<<newt<<endl;
#endif // DEBUG_CONIC_2
	  fac=gcd(gcd(m11,m21),m33);
	  m11/=fac; m21/=fac; m33/=fac;
	  m33=m33*t;
	}	
      break; // end of method=1
    case 2: // lattice reduction
#ifdef DEBUG_REDUCE
      cout<<"Before lattice reduction, x0="<<x0<<", t="<<newt<<endl;
      long count=0;
#endif // DEBUG_REDUCE
      bigint absa=abs(a);
      bigint u1, v1, u2, v2, w;
      u1=1; v2=1;
      bigint n1=absa+sqr(x0), dot=b*x0, n2;
      bigint alpha = roundover(dot, n1);  // nearest int to quotient
      int reduced = is_zero(alpha);
      if(!reduced)
	{
	  u2-=alpha*u1;
	  v2-=alpha*v1;
	}
      n2 = absa*sqr(u2) + sqr(u2*x0+v2*b);
      reduced=(n2>=n1);
#ifdef DEBUG_REDUCE
      //	  cout<<"Entering lattice reduction loop, u1="
      //              <<u1<<", v1="<<v1<<", n1="<<n1<<endl;
      //	  cout<<"Entering lattice reduction loop, u2="
      //              <<u2<<", v2="<<v2<<", n2="<<n2<<endl;
      //	  cout<<"(alpha was "<<alpha<<")\n";
#endif // DEBUG_REDUCE
      while(!reduced)
	{
	  swap(u1,u2);  // w=u1; u1=u2; u2=w;
	  swap(v1,v2);  // w=v1; v1=v2; v2=w;
	  swap(n1,n2);  // w=n1; n1=n2; n2=w;
	  dot=absa*u1*u2+(u1*x0+v1*b)*(u2*x0+v2*b);
	  alpha=roundover(dot,n1);
	  reduced = is_zero(alpha);
	  if(!reduced)
	    {
	      u2-=alpha*u1;
	      v2-=alpha*v1;
	      n2 = absa*sqr(u2) + sqr(u2*x0+v2*b);
	    }
	  reduced=(n2>=n1);
#ifdef DEBUG_REDUCE
	  count++;
	  //	      cout<<"In lattice reduction loop, u1="
	  //              <<u1<<", v1="<<v1<<", n1="<<n1<<endl;
	  //	      cout<<"In lattice reduction loop, u2="
	  //              <<u2<<", v2="<<v2<<", n2="<<n2<<endl;
#endif // DEBUG_REDUCE
	}
      m11 = u1*x0+v1*b;  
      m21 = u1;
      t   = (sqr(m11)-a*sqr(m21))/b;
      m33 = t;
#ifdef DEBUG_REDUCE
      cout<<"("<<count<<" steps in reduction)\n";
      cout<<"After  lattice reduction, x0="<<m11<<", z0="<<m21
	<<", t="<<t<<endl;
#endif //DEBUG_REDUCE
      break; // end of method=2
    } // end of switch(method)

  //Now use recursion unless t=square already
  if(is_one(t))
    {
      x=m11; y=1; z=m21;
#ifdef DEBUG_CONIC
      cout << "...returns ";  show_xyz(x,y,z);
#endif //DEBUG_CONIC
      return 1;
    }
  sqfdecomp(t,t0,t1,tplist,TRACE_FACTORIZATION);
  if(is_one(t0))
    {
      x=m11; y=t1; z=m21;
      bigint g = cancel1(x,y,z);
      int verb=0;
#ifdef DEBUG_CONIC
      verb=1;
      if(method==2)
	cout << "... (after cancelling "<<g<<"): ";  show_xyz(x,y,z);
#endif //DEBUG_CONIC
#ifdef REDUCE_INTERMEDIATES
      conic_diag_reduce(a,b,x,y,z,verb);
#endif // REDUCE_INTERMEDIATES
#ifdef DEBUG_CONIC
      cout << "...returns ";  show_xyz(x,y,z);
#endif //DEBUG_CONIC
      return 1;
    }
  m33/=t1;

  res = solve_conic_diag(a,aplist,t0,tplist,x1,y1,z1,method);
  if(!res) 
    {
#ifdef DEBUG_CONIC
      cout << "...returns fail\n";
#endif //DEBUG_CONIC
      return 0; // since a must be -1 too and xx+zz=-yy insoluble
    }
  x =  m11*x1 + a*m21*z1;
  z =  m21*x1 +   m11*z1;
  y =  m33*y1;

  bigint g = cancel1(x,y,z);
  int verb=0;
#ifdef DEBUG_CONIC
  verb=1;
  if(method==2)
    cout << "... (after cancelling "<<g<<"): ";  show_xyz(x,y,z);
#endif // DEBUG_CONIC
#ifdef REDUCE_INTERMEDIATES
  conic_diag_reduce(a,b,x,y,z,verb);
#endif // REDUCE_INTERMEDIATES
#ifdef DEBUG_CONIC
  cout << "...returns ";  show_xyz(x,y,z);
#endif //DEBUG_CONIC
  return 1;
}      // end of if(method<3)
  else // method==3
    {
      bigint x0, t, t0, t1, x1, y1, z1;  int res; vector<bigint> tplist;
      res = modsqrt(a,bplist,x0); // Solves x^2=a mod b, returns success/fail
      if(!res) 
	{
#ifdef DEBUG_CONIC
	  cout << "...returns fail (no modular sqrt)\n";
#endif // DEBUG_CONIC
	  return 0; // since a must be -1 too and xx+zz=-yy insoluble
	}

  // Always have v ^2 - a*u ^2  = b*c
  //             vv^2 - a*uu^2  = c*newc

      bigint u, v, temp, c, uu, vv, nor, newc;
      u=1; v=x0;
      c = (x0*x0-a)/b;
      vv = mod(-x0,c), uu=1;
      nor = vv*vv-a*uu*uu;
      //  if(!div(c,nor)) cout<<"Error: "<<nor<<" not divisible by c\n";
      newc = nor/c;
#ifdef DEBUG_CONIC_2
      cout<<"   c="<<c<<"\n";
      cout<<"newc="<<newc<<"\n";
      //  cout<<"(v,u)=("<<v<<","<<u<<")\n";
#endif //DEBUG_CONIC_2
      while(!is_zero(newc)&&(abs(newc)<abs(c)))  // update u, v
	{
	  temp = (u*vv+v*uu)/c;    // exact by defn of uu, vv
	  v    = (v*vv+a*u*uu)/c;  // exact since v^2-au^2=0(mod c)
	  u    = temp;
	  //      cout<<"(v,u)=("<<v<<","<<u<<")\n";
	  uu = mod(u,abs(newc));
	  vv = mod(-v,abs(newc));
	  //      cout<<"(vv,uu)=("<<vv<<","<<uu<<")\n";
	  c    = newc;
	  nor = vv*vv-a*uu*uu;
	  //      if(nor==0) 
	  //         cout<<"Error: nor=0 when uu="<<uu<<", vv="<<vv<<endl;
	  //      if(!div(c,nor)) 
	  //         cout<<"Error: "<<nor<<" not divisible by c\n";
	  newc = nor/c;
#ifdef DEBUG_CONIC_2
	  cout<<"newc="<<newc<<"\n";
#endif // DEBUG_CONIC_2
	}
  // now we have v^2-a*u^2=b*c with smallest c, and
  // use recursion to solve x^2-a*z^2=c*y^2

      sqfdecomp(c,t0,t1,tplist,TRACE_FACTORIZATION);
      res = solve_conic_diag(a,aplist,t0,tplist,x1,y1,z1,method);
      if(!res) 
	{
#ifdef DEBUG_CONIC
	  cout << "...returns fail\n";
#endif // DEBUG_CONIC_2
	  return 0; // since a must be -1 too and xx+zz=-yy insoluble
	}
      x =  v*x1 + a*u*z1;
      z =  u*x1 +   v*z1;
      y =  t0*t1*y1;

      bigint g = cancel1(x,y,z);
#ifdef DEBUG_CONIC
      cout << "...returns (after cancelling "<<g<<") ";  show_xyz(x,y,z);
#endif // DEBUG_CONIC_2
      return 1;
    } // end of method==3
}     // end of solve_conic_diag_nontriv()

bigint cancel1(bigint& x, bigint& y, bigint& z)
     // cancels common factors only
{
  bigint g=gcd(x,y);
  if(!is_one(g)) 
    {
      g=gcd(g,z);
      if(!is_one(g)) {x/=g; y/=g; z/=g;}
    }
  return g;
}

void cancel(bigint& x, bigint& y, bigint& z)
     // cancels common factors and leaves z>=0 or z=0 and x>=0
{
  cancel1(x,y,z);
  if(is_positive(z)) return;
  if(is_negative(z)) {::negate(x); ::negate(y); ::negate(z); return;}
  if(is_positive(y)) return;
  if(is_negative(y)) {::negate(x); ::negate(y); return;}
  if(is_negative(x)) {::negate(x);}
  return;
}

int testsol(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                const bigint& x, const bigint& y, const bigint& z, int verb)
{
  int triv = is_zero(x) && is_zero(y) && is_zero(z);
  if(triv)
    {
      if(verb) cout << "Trivial solution!\n";
      return 0;
    }
  bigint t = a*x*x+b*x*z+c*z*z-d*y*y;
  if(is_zero(t))
    {
      if(verb) cout << "Solution OK!\n";
      return 1;
    }
  else 
    {
      if(verb) cout << "Solution wrong!\n";
      return 0;
    }

}

int testlocsol(const bigint& a, 
	       const bigint& b, 
	       const bigint& c)
// tests if ax^2+by^2+cz^2=0 is soluble, where a, b, c are pairwise
// coprime and square-free
{
  vector<bigint> alist=pdivs(a);
  vector<bigint> blist=pdivs(b);
  vector<bigint> clist=pdivs(c);
  return testlocsol(a,alist,b,blist,c,clist);
} 

int testlocsol(const bigint& a, const vector<bigint>& alist, 
	       const bigint& b, const vector<bigint>& blist, 
	       const bigint& c, const vector<bigint>& clist) 
// tests if ax^2+by^2+cz^2=0 is soluble, where a, b, c are pairwise
// coprime and square-free, their prime factors being in alist etc.
{
  int as=sign(a), bs=sign(b), cs=sign(c);
  if((as==bs)&&(bs==cs)) 
    {
      //      cout<<"testlocsol("<<a<<","<<b<<","<<c<<") returning 0 because of signs\n";
      return 0;
    }  
  bigint p, two; two=2;
  bigint mab=-a*b; 
  vector<bigint>::const_iterator pr;
  pr=clist.begin(); 
  while(pr!=clist.end())
    {
      p=*pr++;
      if(p==two) continue;
      if(legendre(mab,p)!=1) 
	{
// 	  cout<<"testlocsol fails legendre(mab,p) with "
// 	      <<"(a,b,p)=("<<a<<","<<b<<","<<p<<")\n"; 
	  return 0;
	}
    }
  bigint mbc=-b*c; 
  pr=alist.begin();
  while(pr!=alist.end())
    {
      p=*pr++;
      if(p==two) continue;
      if(legendre(mbc,p)!=1)
	{
// 	  cout<<"testlocsol fails legendre(mbc,p) with "
// 	      <<"(b,c,p)=("<<b<<","<<c<<","<<p<<")\n"; 
	  return 0;
	}
    }
  bigint mca=-c*a;
  pr=blist.begin();
  while(pr!=blist.end())
    {
      p=*pr++;
      if(p==two) continue;
      if(legendre(mca,p)!=1)
	{
// 	  cout<<"testlocsol fails legendre(mca,p) with "
// 	      <<"(c,a,p)=("<<c<<","<<a<<","<<p<<")\n"; 
	  return 0;
	}
    }
  return 1;
} 

int testlocsol(const bigint& a, const vector<bigint>& alist, 
	       const bigint& b, const vector<bigint>& blist)
// tests if ax^2+by^2=z^2 is soluble, where a, b are
// square-free, their prime factors being in alist and blist.
{
  // Avoid any factorization and gcd computation using the primes given
  bigint p, a0, b0, c;
  a0=1; b0=1; c=-1;
  vector<bigint> a0list, b0list, clist;

  long sa=sign(a), sb=sign(b);
  if((sa<0)&&(sb<0)) 
    {
//       cout<<"testlocsol("<<a<<","<<b<<") returning 0 because of signs\n";
      return 0;  // nothing more to do as no real solution
    }
  if(sa<0) ::negate(a0);
  if(sb<0) ::negate(b0);
  
  vector<bigint>::const_iterator pr;
  pr=alist.begin(); 
  while(pr!=alist.end())
    {
      p=*pr++;
      if(div(p,b)) 
	{
	  c*=p; clist.push_back(p);
	}
      else
	{
	  a0*=p; a0list.push_back(p);
	}
    }
  pr=blist.begin(); 
  while(pr!=blist.end())
    {
      p=*pr++;
      if(!div(p,c)) {b0*=p; b0list.push_back(p);}
    }  
  
#if(0)
  if((a!=-a0*c)||(b!=-b0*c)||(abs(c)!=gcd(a,b)))
    {
      cout<<"Error: (a,b)=("<<a<<","<<b<<") gives a0="<<a0<<", b0="<<b0<<", c="<<c<<endl;
    }
  //  cout<<"Calling testlocsol(a,b,c) with a,b,c="<<a0<<","<<b0<<","<<c<<"\n";
  //  cout<<"alist="<<a0list<<endl;
  //  cout<<"blist="<<b0list<<endl;
  //  cout<<"clist="<<clist<<endl;
#endif

  // Now a=a0*c, b=b0*c, c=gcd(a,b), and a0list, b0list, clist hold their primes
  
  return testlocsol(a0,a0list,b0,b0list,c,clist);

}

void testmodsqrt()
{
  long i,m;
  int res, ok=1;
  bigint a,mm,x;
  cout << "Enter a modulus m: ";
  cin >> m; mm=m;
  vector<bigint> plist=pdivs(mm);
  int* flag = new int[m];
  for(i=0; i<m; i++) flag[i]=0;
  for(i=0; i<=m/2; i++) flag[(i*i)%m]=1;
  for(i=0; i<m; i++) 
    {
      a=i; //      cout<<"a = "<<a;

      res = modsqrt(a,plist,x);
#if(0)
      if(res)
	{
	  cout << "\tsqrt(a) mod m = " << x;
	  if((x*x-a)%mm==0) cout<<" --checks.";
	  else cout << "--WRONG!";
	}
      else 
	{
	  cout << "\tNo solution";
	}
#endif
      if(res!=flag[i]) 
	{
	  cout << "WRONG ANSWER for a="<<a<<endl;
	  ok=0;
	}
    }
  if(ok) cout << "All correct" << endl;
}

void testsqf()
{
  bigint a,m,m1,m2;
  vector<bigint> plist;
  while(cout << "Enter a nonzero integer m: ",  cin >> m, m!=0)
    {
      sqfdecomp(m,m1,m2,plist,TRACE_FACTORIZATION);
      cout << "sqfdecomp returns m1 = " << m1 << " and m2 = " << m2 << endl;
      cout << "(plist = " << plist << ")\n";
      sqfdecomp(m,plist,m1,m2);
      cout << "sqfdecomp returns m1 = " << m1 << " and m2 = " << m2 << endl;
      cout << "(plist = " << plist << ")\n";
    }
}

void testcancel()
{
  bigint x,y,z;
  cout << "Enter x, y, z to be cancelled: ";
  cin >> x >> y >> z;
  cout << "Before: (x:y:z) = ("<<x<<":"<<y<<":"<<z<<")\n";
  cancel(x,y,z);
  cout << "After:  (x:y:z) = ("<<x<<":"<<y<<":"<<z<<")\n";
}

int solve_conic_param(const quadratic& q, const bigint& d,
			     const vector<bigint>& factorbase,
			     quadratic& qx, quadratic& qy, quadratic& qz, 
			     int method, int verb)
{
  return solve_conic_param(q[0],q[1],q[2],d,factorbase,qx,qy,qz,method,verb);
}

int solve_conic_param(const quadratic& q, const bigint& d,
			     quadratic& qx, quadratic& qy, quadratic& qz, 
			     int method, int verb)
{
  return solve_conic_param(q[0],q[1],q[2],d,qx,qy,qz,method,verb);
}

int solve_conic_param(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                quadratic& qx, quadratic& qy, quadratic& qz, int method, int verb)
{
  vector<bigint> factorbase = pdivs(2*d);
  if(is_zero(b))
    {
      factorbase=vector_union(factorbase,pdivs(a));
      factorbase=vector_union(factorbase,pdivs(c));
    }
  else
    {
      bigint disc = b*b-4*a*c;
      factorbase=vector_union(factorbase,pdivs(a));
      factorbase=vector_union(factorbase,pdivs(disc));
    }
  return solve_conic_param(a,b,c,d,factorbase,qx,qy,qz,method,verb);
}

int solve_conic_param(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
		const vector<bigint>& factorbase,
                quadratic& qx, quadratic& qy, quadratic& qz, int method, int verb)
     // Solves axx+bxz+czz=dyy for (x,y,z) not (0,0,0) and returns 1
     // or returns 0 if not possible
     // Should have a, c, d, bb-4ac non-zero
     // qx,qy,qz are (preallocated) arrays of coeffs of parametrizing quadratics
     //   with leading coeffs qx[0],qy[0],qz[0] one solution
{
  if(verb) cout<<"In solve_conic_param() with (a,b,c,d)=("<<a<<","<<b<<","<<c<<","<<d<<"),\nfactorbase="<<factorbase<<"...\n";
  bigint disc = b*b-4*a*c, x0, y0, z0;

  bigint delta, zero, one; zero=0; one=1;
  if(isqrt(disc,delta)) // case where discrim is square
    {
      if(verb) cout<<"disc="<<disc<<" is square, root = "<<delta<<endl;
      x0=(delta-b)/2; y0=0; z0=a;
      if(verb) 
	{
	  cout << "Solution: (x:y:z) = ("<<x0<<":"<<y0<<":"<<z0<<")\n\n";
	  cout << "Finding parametric solution...\n";
	}
      qx.set(a*d*x0,   zero,    (delta+b)/2);
      qy.set(zero,     a*delta,  zero);
      qz.set(sqr(a)*d, zero,    -a);
      if(verb) 
	{
	  cout<<"original x coeffs: ["<<qx<<"\n";
	  cout<<"original y coeffs: ["<<qy<<"\n";
	  cout<<"original z coeffs: ["<<qz<<"\n";
	  testparamsol(a,b,c,d,qx,qy,qz,verb);
	}
      bigint a1=gcd(a,qx[2]);  bigint a2=a*a/a1;
      bigint t;
      if(divide_exact(qx[0],a2,t)) {qx.set_coeff(0,t);} else
	cout<<"Problem: "<<a2<<" ndiv qx[0] = "<<qx[0]<<endl;
      if(divide_exact(qx[1],a,t))  {qx.set_coeff(1,t);} else
	cout<<"Problem: "<<a <<" ndiv qx[1] = "<<qx[1]<<endl;
      if(divide_exact(qx[2],a1,t))  {qx.set_coeff(2,t);} else
	cout<<"Problem: "<<a1<<" ndiv qx[2] = "<<qx[2]<<endl;
      if(divide_exact(qy[0],a2,t))  {qy.set_coeff(0,t);} else
	cout<<"Problem: "<<a2<<" ndiv qy[0] = "<<qy[0]<<endl;
      if(divide_exact(qy[1],a,t))  {qy.set_coeff(1,t);} else
	cout<<"Problem: "<<a <<" ndiv qy[1] = "<<qy[1]<<endl;
      if(divide_exact(qy[2],a1,t))  {qy.set_coeff(2,t);} else
	cout<<"Problem: "<<a1<<" ndiv qy[2] = "<<qy[2]<<endl;
      if(divide_exact(qz[0],a2,t))  {qz.set_coeff(0,t);} else
	cout<<"Problem: "<<a2<<" ndiv qz[0] = "<<qz[0]<<endl;
      if(divide_exact(qz[1],a,t))  {qz.set_coeff(1,t);} else
	cout<<"Problem: "<<a <<" ndiv qz[1] = "<<qz[1]<<endl;
      if(divide_exact(qz[2],a1,t))  {qz.set_coeff(2,t);} else
	cout<<"Problem: "<<a1<<" ndiv qz[2] = "<<qz[2]<<endl;
      
      if(verb) 
	{
	  cout<<"reduced x coeffs: "<<qx<<"\n";
	  cout<<"reduced y coeffs: "<<qy<<"\n";
	  cout<<"reduced z coeffs: "<<qz<<"\n";
	  testparamsol(a,b,c,d,qx,qy,qz,verb);
	  cout<<"leaving solve_conic_param()\n";
	}
      return 1;
    }

//generic case, disc non-square

  int res = solve_conic(a,b,c,d,factorbase,x0,y0,z0,method);
  if(!res) return 0;

  if(verb) 
    {
      cout << "Solution: (x:y:z) = ("<<x0<<":"<<y0<<":"<<z0<<")\n\n";
      cout << "Finding parametric solution...\n";
    }

  qx.set(x0,  2*(b*x0+2*c*z0), x0*disc);
  qy.set(y0, zero,            -y0*disc);
  qz.set(z0, -2*(2*a*x0+b*z0), z0*disc);

  if(verb)
    {
      cout<<"original x coeffs: "<<qx<<"\n";
      cout<<"original y coeffs: "<<qy<<"\n";
      cout<<"original z coeffs: "<<qz<<"\n";
      testparamsol(a,b,c,d,qx,qy,qz,verb);
    }

  // The quadratics qx and qz have discriminants 16*c*d*y0^2, 16*a*d*y0^2,
  // and resultant 16*(b^2-4*a*c)*y0^4.  
  // We adjust them to remove factors 4*y0^2, 4*y0^2, 16*y0^4
  // (new version 27/7/98)

  bigint s,t;
  bigint g = bezout(x0,z0,s,t);
  bigint fac = 2*y0; bigint fac2=sqr(fac);
  bigint e = (x0*(2*t*a-s*b) + z0*(t*b-2*c*s)) % (2*(sqr(y0)));
  if(verb) 
    cout<<"Shifting by e = " << e << " and dividing out by factor "<<fac<<endl;
      
  qx.set_coeff(2,qx[2]+e*(qx[1]+e*qx[0]));
  qx.set_coeff(1,qx[1]+2*e*qx[0]);
  if(divide_exact(qx[1],fac,t))  {qx.set_coeff(1,t);} else
    cout<<"Problem: factor ndiv qx[1] = "<<qx[1]<<endl;
  if(divide_exact(qx[2],fac2,t)) {qx.set_coeff(2,t);} else
    cout<<"Problem: factor^2 ndiv qx[2] = "<<qx[2]<<endl;
  
  qy.set_coeff(2,qy[2]+e*(qy[1]+e*qy[0]));
  qy.set_coeff(1,qy[1]+2*e*qy[0]);
  if(divide_exact(qy[1],fac,t))   {qy.set_coeff(1,t);} else
    cout<<"Problem: factor ndiv qy[1] = "<<qy[1]<<endl;
  if(divide_exact(qy[2],fac2,t))  {qy.set_coeff(2,t);} else 
    cout<<"Problem: factor^2 ndiv qy[2] = "<<qy[2]<<endl;
  
  qz.set_coeff(2,qz[2]+e*(qz[1]+e*qz[0]));
  qz.set_coeff(1,qz[1]+2*e*qz[0]);
  if(divide_exact(qz[1],fac,t))   {qz.set_coeff(1,t);} else
    cout<<"Problem: factor ndiv qz[1] = "<<qz[1]<<endl;
  if(divide_exact(qz[2],fac2,t))  {qz.set_coeff(2,t);} else 
    cout<<"Problem: factor^2 ndiv qz[2] = "<<qz[2]<<endl;

  if(verb) cout<<"leaving solve_conic_param()\n";
  return 1;
}

int testparamsol(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
                const quadratic& qx, const quadratic& qy, const quadratic& qz, int verb)
     // Tests to see if a given parametrization is a solution
{
  int ok;
  bigint coeff;
  // x^4 coefficient:
  coeff = a*sqr(qx[0]) + b*qx[0]*qz[0] + c*sqr(qz[0]) - d*sqr(qy[0]);
  ok=is_zero(coeff);
  if(!ok) 
    {
      if(verb) cout<<"Coefficient of x^4 is wrong\n";
      return 0;
    }
  coeff = 2*a*qx[0]*qx[1] + b*(qx[0]*qz[1]+qx[1]*qz[0]) + 2*c*qz[0]*qz[1] 
       - 2*d*qy[0]*qy[1];
  ok=is_zero(coeff);
  if(!ok) 
    {
      if(verb) cout<<"Coefficient of x^3 is wrong\n";
      return 0;
    }
  coeff = a*(sqr(qx[1])+2*qx[0]*qx[2]) + 
    b*(qx[0]*qz[2]+qx[1]*qz[1]+qx[2]*qz[0]) + 
    c*(sqr(qz[1])+2*qz[0]*qz[2]) - d*(sqr(qy[1])+2*qy[0]*qy[2]);
  ok=is_zero(coeff);
  if(!ok) 
    {
      if(verb) cout<<"Coefficient of x^2 is wrong\n";
      return 0;
    }
  coeff = 2*a*qx[1]*qx[2] + b*(qx[1]*qz[2]+qx[2]*qz[1]) + 2*c*qz[1]*qz[2] 
    - 2*d*qy[1]*qy[2];
  ok=is_zero(coeff);
  if(!ok) 
    {
      if(verb) cout<<"Coefficient of x^1 is wrong\n";
      return 0;
    }
  coeff = a*sqr(qx[2]) + b*qx[2]*qz[2] + c*sqr(qz[2]) - d*sqr(qy[2]);
  ok=is_zero(coeff);
  if(!ok) 
    {
      if(verb) cout<<"Coefficient of x^0 is wrong\n";
      return 0;
    }
  if(verb) 
    {
      cout<<"Parametric solution is correct!\n";
      bigint discx = qx.disc();
      cout<<"x-disc = "<<discx<<"\n";
      bigint discy = qy.disc();
      cout<<"y-disc = "<<discy<<"\n";
      bigint discz = qz.disc();
      cout<<"z-disc = "<<discz<<"\n";
      bigint res = resultant(qx,qz);
      cout<<"resultant = "<<res<<"\n";
    }
  return 1;
}

void conic_mordell_reduce(const bigint& a, const bigint& b, const bigint& c, bigint& x0, bigint& y0, bigint& z0, int verb)
     // Given a>0, b>0, c<0, abc square-free and ax^2+by^2+cz^2=0
     // reduces x, y, z in place using Mordell's method (page 48)
     // to achieve Holzer's bounds |z|<=sqrt(ab) etc.
  // NOTE that we allow c not to be square-free here, 
  // so we may have g=gcd(x0,y0)>1 with g^2|c.
{
  // To check if Holzer's conditions are satisfied
  // it suffices to chack that |z| < sqrt(ab).
  int ok = (sqr(z0)<=(a*b));
  //  if(ok) return;
  bigint zero; zero=0;
  if(verb)
    {
      bigfloat xlim = sqrt(I2bigfloat(b)*I2bigfloat(-c));
      bigfloat ylim = sqrt(I2bigfloat(a)*I2bigfloat(-c));
      bigfloat zlim = sqrt(I2bigfloat(a)*I2bigfloat(b));
      cout<<"Using Mordell reduction to reduce the solution ";
      show_xyz(x0,y0,z0);
      cout<<" for (a,b,c) = ("<<a<<","<<b<<","<<c<<") ";
      cout<<" so that\n";
      cout<<"|x0| <= "<<xlim<<endl;
      cout<<"|y0| <= "<<ylim<<endl;
      cout<<"|z0| <= "<<zlim<<endl;
      cout<<"First check this is a solution...";
      testsol(a,zero,b,-c,x0,z0,y0,1);
    }

  int fail=0;  
  // will be set to 1 if a reduction step fails to reduce (for debugging)
  int steps=0;
  while(!(ok||fail))
    {
      if(verb) 
	cout<<"Holzer's conditions fail, reducing the solution using Mordell's method...\n";
      steps++;
      bigint d,g,gxy,xx,yy,zz,u,v,newx0,newy0,newz0;
      bigint x00=x0, y00=y0, c0=c;
      gxy=abs(gcd(x0,y0)); int fix=0;
      if(gxy>1) {fix=1; x00/=gxy; y00/=gxy; c0/=sqr(gxy);}
      int c_parity;
      if(odd(c0)) {c_parity=1; d=c0;} else {c_parity=0; d=c0/2;}
      g=bezout(y00,-x00,xx,yy);
      xx*=d; yy*=d;  // now y00*xx-x00*yy=d
      u = a*x00*xx+b*y00*yy; v = -c0*z0;
      if(v<0) {::negate(u); ::negate(v);}
      bigfloat zz_real; zz_real=0;
      if(verb>1)
	{
	  zz_real = I2bigfloat(u)/I2bigfloat(v);
	  cout<<"u = "<<u<<"\n";
	  cout<<"v = "<<v<<"\n";
	  cout<<"Z-real = "<<zz_real<<endl;
	}
      if(c_parity)
	{
	  int z_parity = odd(a*xx+b*yy);
	  if(z_parity==0) zz=2*roundover(u,2*v);
	  else
	    {
	      zz=roundover(u,v);
	      if(even(zz)) {if(v*zz>u) zz-=1; else zz+=1;}
	    }
	  if(verb>1)
	    {
	      cout<<"Z should be closest to "<<zz_real
		  <<" with parity "<<z_parity<<"\n";
	      cout<<"Z = "<<zz<<endl;
	    }
	}
      else
	{
	  zz = roundover(u,v);
	  if(verb>1)
	    {
	      cout<<"Z should be closest to "<<zz_real<<"\n";
	      cout<<"Z = "<<zz<<endl;
	    }
	}

      //compute new solution with reduced z0:

      u=a*sqr(xx)+b*sqr(yy)+c0*sqr(zz);
      v=2*(a*x00*xx+b*y00*yy+c0*z0*zz);
      if(c_parity)d*=2;
      if(!divide_exact(x00*u-xx*v,d,newx0)) 
	cout<<"Problem: new x0 not integral!\n";
      if(!divide_exact(y00*u-yy*v,d,newy0)) 
	cout<<"Problem: new y0 not integral!\n";
      if(!divide_exact(z0*u-zz*v,d,newz0)) 
	cout<<"Problem: new z0 not integral!\n";
      if(fix) {newx0*=gxy; newy0*=gxy;}
      cancel1(newx0,newy0,newz0);
      if((fail = (abs(newz0)>=abs(z0))))
	{
	  cout<<"Problem: new solution is NOT smaller than old one!\n";
	  cout<<"(a,b,c) = ("<<a<<","<<b<<","<<c<<")\n";
	  cout<<"Old solution = "; show_xyz(x0,y0,z0);
	  cout<<"\nNew solution = "; show_xyz(newx0,newy0,newz0);
	  cout<<"\nX = "<<xx<<"\n";
	  cout<<"Y = "<<yy<<"\n";
	  cout<<"Z = "<<zz<<"\n";
	}
      x0=newx0; y0=newy0; z0=newz0;
      
      if(verb)
	{
	  cout<<"Solution "<<steps<<" = ("<<x0<<":"<<y0<<":"<<z0<<")\n";
	  testsol(a,zero,b,-c,x0,z0,y0,0);
	}

  // Check if Holzer's conditions are satisfied yet.  

      ok = ((z0*z0)<=(a*b));
    }

  if((!fail) &&verb)
    {
      cout<<steps<<" reduction steps were needed.\n";
      cout<<"Reduced solution = (" <<x0<<":"<<y0<<":"<<z0<<")\n";
      if(verb) testsol(a,zero,b,-c,x0,z0,y0,verb);
      cout<<"Holzer's conditions are satisfied\n";
    }
}

void conic_diag_reduce(const bigint& a, const bigint& b, bigint& x, bigint& y, bigint& z, int verb)
  // As above but with a,b square-free only, and x,y,z satisfying
  // x^2-az^2=by^2.   
  // Calls conic_mordell_reduce() or new_legendre_reduce()
{
  bigint zero, one; zero=0; one=1;
  int debug=0;
#ifdef DEBUG_CONIC
  debug=1;
  if(verb) 
    {
      cout<<"\nAt start of conic_diag_reduce() with (a,b)=("<<a<<","<<b<<")\n";
      cout<<"Solution =";      show_xyz(x,y,z);  cout<<"\t";
      testsol(one,zero,-a,b,x,y,z,verb);
    }
#endif
  int sa=sign(a), sb=sign(b), icase;
  if(sa==-1) { icase=1;}   // a<0, b>0
  else
    {
      if(sb==1) {icase=2;} // a>0, b>0
      else  {icase=3;}     // a>0, b<0
    }

  bigint c = gcd(a,b);  y*=c; z*=c;
  bigint g = gcd(x,c), a0=a/c, b0=b/c;
  if(g>1) {x/=g; y/=g; z/=g;}
  if(debug) 
    {
      cout<<"Case "<<icase<<": testing before reduction...";
      switch(icase)
	{
	case 1: testsol(c, zero,-a0,b0,x,y,z,1); break;
	case 2: testsol(a0,zero, b0,c, z,x,y,1); break;
	case 3: testsol(c, zero,-b0,a0,x,z,y,1);
	}
    }
  switch(icase)
    {
#ifdef MORDELL_REDUCE
    case 1: conic_mordell_reduce(c,-a0,-b0,x,z,y,verb); break;
    case 2: conic_mordell_reduce(a0,b0,-c, z,y,x,verb);  break;
    case 3: conic_mordell_reduce(c,-b0,-a0,x,y,z,verb);
#else // use JC's quadratics method
    case 1: new_legendre_reduce(c,-a0,-b0,x,z,y,verb); break;
    case 2: new_legendre_reduce(a0,b0,-c, z,y,x,verb);  break;
    case 3: new_legendre_reduce(c,-b0,-a0,x,y,z,verb);
#endif
    }
  if(debug) 
    {
      cout<<"Case "<<icase<<": testing after  reduction...";
      switch(icase)
	{
	case 1: testsol(c, zero,-a0,b0,x,y,z,1); break;
	case 2: testsol(a0,zero, b0,c, z,x,y,1); break;
	case 3: testsol(c, zero,-b0,a0,x,z,y,1);
	}
    }
  x*=c;
  if(verb&&debug) 
    {
      cout<<"At end of conic_diag_reduce()..";
      testsol(one,zero,-a,b,x,y,z,verb);
    }
}

int testsol(const quadratic& q, const bigint& d,
		  const bigint& x, const bigint& y, const bigint& z, 
		  int verb)
{
  return testsol(q[0],q[1],q[2],d,x,y,z,verb);
}

int testparamsol(const quadratic& q, const bigint& d,
			const quadratic& qx, const quadratic& qy, const quadratic& qz, 
			int verb)
{
  return testparamsol(q[0],q[1],q[2],d,qx,qy,qz,verb);
}

// Output utilities:

void show_xyz(const bigint& x, const bigint& y, const bigint& z)
{
  cout << "(x:y:z) = ("<<x<<":"<<y<<":"<<z<<")";
}
void show_cert(const bigint& p, const bigint& q, const bigint& r)
{
  cout<<"Certificate: ("<<p<<", "<<q<<", "<<r<<")";
}
void show_eqn(const bigint& a, const bigint& b, const bigint& c)
{
  cout<<"(a,b,c) = ("<<a<<", "<<b<<", "<<c<<")";
}
void show_eqn_cert(const bigint& a, const bigint& b, const bigint& c,
		   const bigint& p, const bigint& q, const bigint& r)
{
  show_eqn(a,b,c);cout<<endl;
  show_cert(p,q,r);cout<<endl;
}
void show_all(const bigint& a, const bigint& b, const bigint& c, 
	      const bigint& p, const bigint& q, const bigint& r, 
	      bigint& x, bigint& y, bigint& z)
{
  show_eqn(a,b,c);cout<<endl;
  show_cert(p,q,r);cout<<endl;
  show_xyz(x,y,z);cout<<endl;
}
