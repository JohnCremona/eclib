// marith.h: declarations of integer arithmetic functions (multiprecision)
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
 
#ifndef _ECLIB_MARITH_H
#define _ECLIB_MARITH_H

#include "arith.h"

inline int is_zero(const ZZ& x) {return IsZero(x);}
inline int is_nonzero(const ZZ& x) {return !IsZero(x);}
inline int is_positive(const ZZ& x) {return sign(x)>0;}
inline int is_negative(const ZZ& x) {return sign(x)<0;}
inline int is_one(const ZZ& x) {return IsOne(x);}
inline int odd(const ZZ& a) {return IsOdd(a);}
inline int even(const ZZ& a) {return !IsOdd(a);}
inline void rshift(const ZZ& a, long i, ZZ& c) {RightShift(c,a,i);}
inline void lshift(const ZZ& a, long i, ZZ& c) {LeftShift(c,a,i);}
#ifdef setbit
#undef setbit
#endif
inline void setbit(ZZ& a, int e) {SetBit(a,e);}

// The following are not in NTL & need defining
inline void longasI(long& a, const ZZ& x) {a = I2long(x);}
inline void negate(ZZ& a) {a=-a;}
inline void sqrt(ZZ& a, const ZZ& b) {SqrRoot(a,b);}
inline ZZ sqrt(const ZZ& a) {ZZ b; sqrt(b,a); return b;}
inline void square(ZZ& a, const ZZ& b) {sqr(a,b);}
inline ZZ gcd(const ZZ& a, const ZZ& b) {return GCD(a,b);}
inline ZZ lcm(const ZZ& a, const ZZ& b)
 {if (IsZero(a) && IsZero(b)) return ZZ::zero(); else return a*(b/GCD(a,b));}

// In NTL add, sub, mul, div are defined with result in first place
inline void addx(const ZZ& a, const ZZ& b, ZZ& c)  {add(c,a,b);}
inline void subx(const ZZ& a, const ZZ& b, ZZ& c)  {sub(c,a,b);}
inline void divx(const ZZ& a, const ZZ& b, ZZ& c)  {div(c,a,b);}
inline void mulx(const ZZ& a, const ZZ& b, ZZ& c)  {mul(c,a,b);}
inline ZZ pow(const ZZ& a, long e)  {return power(a,e);}

//N.B. no power to ZZ exponent in NTL
inline long jacobi(const ZZ& a, const ZZ& p)  {return Jacobi(a,p);}
inline void sqrt_mod_p(ZZ & x, const ZZ & a, const ZZ & p)  
  {SqrRootMod(x,a,p); if(x>(p-x)) x= p-x;}
inline void power_mod(ZZ& ans, const ZZ& base, const ZZ& expo, const ZZ& m) 
 {PowerMod(ans,base,expo,m);}
inline void nearest(ZZ& c, const ZZ& a, const ZZ& b) 
 {ZZ a0=(a%b);  c = (a-a0)/b; if(2*a0>b) c+=1;}
inline ZZ roundover(const ZZ& a, const ZZ& b)
 {ZZ a0=(a%b); ZZ c = (a-a0)/b; if(2*a0>b) c+=1; return c;}

#define ZZ_mod_long(a,m) (a%m)

ZZ bezout(const ZZ& aa, const ZZ& bb, ZZ& xx, ZZ& yy);
int divides(const ZZ& a, const ZZ& b, ZZ& q, ZZ& r);
int divides(const ZZ& a, long b, ZZ& q, long& r);
// returns 1 iff remainder r==0

// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
ZZ rounded_division(const ZZ& a, const ZZ& b);

ZZ Iround(bigfloat x);
ZZ Ifloor(bigfloat x);
ZZ Iceil(bigfloat x);

ZZ posmod(const ZZ& a, const ZZ& b);  // a mod b in range 0--(b-1)
long posmod(const ZZ& a, long b);
int isqrt(const ZZ& a, ZZ& root);
int divide_exact(const ZZ& aa, const ZZ& bb, ZZ& c);
     // c = a/b with error message if remainder is non-zero
long divide_out(ZZ& a, const ZZ& d);
long divide_out(ZZ& a, long d);
// divides a by d as many times as possible returning number of times (but none if a=0!)

ZZ show(const ZZ& a);
vector<ZZ> show(const vector<ZZ>& a);


// extra_primes is a set holding extra big primes found or added manually.  
class extra_prime_class {
public:
  std::set<ZZ> the_primes;
  extra_prime_class() {;}
  ~extra_prime_class();
  void read_from_file(const string pfilename, int verb=0);
  void write_to_file(const string pfilename, int verb=0);
  int contains(const ZZ& p)
  {
    return the_primes.find(p)!=the_primes.end();
  }
  void add(const ZZ& p)  
  {
    if(p>maxprime()) the_primes.insert(p);
  }
  void show() 
  {
    cout << "Extra primes in list: ";
    copy(the_primes.begin(),the_primes.end(), ostream_iterator<ZZ>(cout, " "));
    cout << endl;
  }
};

extern extra_prime_class the_extra_primes;  // The one and only instance

void initprimes(const string pfilename, int verb=0);

//
// divisors
//
// The following uses trial division:
vector<ZZ> pdivs_trial(const ZZ& number, int trace=0);
// The following uses gp factorization externally if available:
vector<ZZ> pdivs_gp(const ZZ& number, int trace=0);
// The following uses pari library's factorization if available:
vector<ZZ> pdivs_pari(const ZZ& number, int trace=0);
// The following uses one of the above
vector<ZZ> pdivs(const ZZ& number, int trace=0);

vector<ZZ> posdivs(const ZZ& number, const vector<ZZ>& plist);
vector<ZZ> posdivs(const ZZ& number);

vector<ZZ> alldivs(const ZZ& number, const vector<ZZ>& plist);
vector<ZZ> alldivs(const ZZ& number);

vector<ZZ> sqdivs(const ZZ& number, const vector<ZZ>& plist);
vector<ZZ> sqdivs(const ZZ& number);

vector<ZZ> sqfreedivs(const ZZ& number, const vector<ZZ>& plist);
vector<ZZ> sqfreedivs(const ZZ& number);

void sqfdecomp(const ZZ& a, ZZ& a1, ZZ& a2, vector<ZZ>& plist, int trace_fact=0);
     // a must be non-zero, computes square-free a1 and a2>0 such that a=a1*a2^2
     // plist will hold prime factors of a

void sqfdecomp(const ZZ& a, vector<ZZ>& plist, ZZ& a1, ZZ& a2);
    // a must be non-zero, computes square-free a1 and a2>0 such that a=a1*a2^2
    // plist already holds prime factors of a

// test for squarefree
int is_squarefree(const ZZ& a);

// (positive) squarefree part of a nonzero integer
ZZ squarefree_part(const ZZ& a);

// squarefree product of two squarefree integers (with signs)
ZZ squarefree_product(const ZZ& a, const ZZ& b);

// Given a, b, lem3 returns m1 etc so that a=c1^2*m1*m12, b=c2^2*m2*m12 
// with m1, m2, m12 pairwise coprime.   At all  times these equations hold, 
// and at each step the product m1*m2*m12 is decreased by a factor d, 
// so the process terminates when the coprimality condition is satisfied. 
void rusin_lem3(const ZZ& a, const ZZ& b,
	  ZZ& m1, ZZ& m2, ZZ& m12, ZZ& c1, ZZ& c2);

// Solves x-a1(mod m1), x=a2(mod m2)
ZZ chrem(const ZZ& a1, const ZZ& a2, 
	     const ZZ& m1, const ZZ& m2);

//
// general purpose routines -- ZZ overloads
//

ZZ mod(const ZZ& a, const ZZ& b);     // a mod b in range +- half b
long mod(const ZZ& a, long b);

inline ZZ addmod(const ZZ& a, const ZZ& b, const ZZ& m)
{
  if (is_zero(a)) return b;
  if (is_zero(b)) return a;
  return mod(a+b,m);
}

inline ZZ xmm(const ZZ& a, const ZZ& b, const ZZ& m)
{
  if (a==1) return b;
  if (a==-1) return -b;
  if (b==1) return a;
  if (b==-1) return -a;
  return (a*b) % m;
}

long val(const ZZ& factor, const ZZ& number);
long val(long factor, const ZZ& number);
vector<int> valuations(const ZZ& n, const vector<ZZ>& primes);
vector<int> valuations(const ZZ& n, const vector<int>& primes);

int div(const ZZ& factor, const ZZ& number);
int div(long factor, const ZZ& number);
inline int ndiv(const ZZ& factor, const ZZ& number)
  { return !div(factor, number);}
inline int ndiv(long factor, const ZZ& number)
  { return !div(factor, number);}


long bezout(const ZZ& aa, long bb, ZZ& xx, ZZ& yy);
ZZ invmod(const ZZ& a, const ZZ& p);  // -a mod p
long invmod(const ZZ& a, long p);

int m1pow(const ZZ& a);
int chi2(const ZZ& a);
int chi4(const ZZ& a);
int hilbert2(const ZZ& a, const ZZ& b);
int hilbert2(const ZZ& a, long b);
int hilbert2(long a, const ZZ& b);
int legendre(const ZZ& a, const ZZ& b);
int legendre(const ZZ& a, long b);
int kronecker(const ZZ& d, const ZZ& n);
int kronecker(const ZZ& d, long n);
// See hilbert.h for hilbert symbol functions

long gcd(const ZZ& a, long b);

// Assuming a*d-b*c!=0, computes a reduced Z-basis for <(a,b),(c,d)>
void gauss_reduce(const ZZ& a0, const ZZ& b0, const ZZ& c0, const ZZ& d0,
                  ZZ& a, ZZ& b, ZZ& c, ZZ& d);

// Set a, b so that a/b=n (mod m) with |a|, |b| minimal; return success if a^2, b^2 <= m/2
int modrat(const ZZ& n, const ZZ& m, ZZ& a, ZZ& b);

int sqrt_mod_2_power(ZZ& x, const ZZ& a, int e);
int sqrt_mod_p_power(ZZ& x, const ZZ& a, const ZZ& p, int e);
int sqrt_mod_m(ZZ& x, const ZZ& a, const ZZ& m);
int sqrt_mod_m(ZZ& x, const ZZ& a, const ZZ& m, const vector<ZZ>& mpdivs);
// Second version of sqrt_mod_m requires mpdivs to hold a list
// of primes which contains all those which divide m; if it acontains
// extras that does not matter.

int modsqrt(const ZZ& a, const vector<ZZ>& bplist, ZZ& x);
     // Solves x^2=a mod b with b square-free, returns success/fail


// root-finding functions:

// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)
// and assign roots to a list of these
int nrootscubic(long b, long c, long d, long p, vector<long>& roots);

void ratapprox(bigfloat x, ZZ& a, ZZ& b, const ZZ& maxd);
void ratapprox(bigfloat x, long& a, long& b, long maxd=0);

int is_nth_power(const ZZ& x, int n);
ZZ prime_to_S_part(const ZZ& x,  const vector<ZZ>& S);
int is_S_unit(const ZZ& x,  const vector<ZZ>& S);

// // class to iterate through divisors of a factored positive integer

class divisor_iterator {
protected:
  int ok;            // flags that iteration not finished
  int np;            // number of primes
  int nd;            // number of divisors
  vector<ZZ> PP; // list of np primes
  vector<long> EE;   // list of np maximum exponents
  vector<long> ee;   // list of np current exponents
  vector<ZZ> NN; // list of np+1 partial products

public:
  // constructors
  divisor_iterator(const vector<ZZ>& P, const vector<long>& E);
  divisor_iterator(const ZZ& N);
  divisor_iterator();

  // increment if possible
  void increment();
  // check for end
  int is_ok() {return ok;}
  // reset
  void rewind()
    {
      ee.resize(np, 0);
      NN.resize(np+1, ZZ(1));
      ok=1;
    }
  // deliver current value
  ZZ value() {return NN[0];}
  // total number of divisors
  long ndivs() {return nd;}
  // report on current status
  void report();
};

// [n^e for 0 <= e <= maxexp]
vector<ZZ> powers(const ZZ& n, int maxexp);
// [n^e for e in exponents]
vector<ZZ> powers(const ZZ& n, const vector<int>& exponents);

// Compute N from its factorization (lists of primes and exponents) --
// (name taken from gp)
ZZ factorback(const vector<ZZ>&PP, const vector<int>& EE);

// Radical of N
ZZ radical(const ZZ& N);

// Maximum conductor for a given list of primes
ZZ MaxN(const vector<ZZ>&S);

// convert a list of longs to a list of ZZs:
vector<ZZ> ZZify(const vector<long>& L);

// multiply all integers in a list by a constant:
vector<ZZ> multiply_list(const ZZ& a, const vector<ZZ>& L);
inline vector<ZZ> multiply_list(long a, const vector<ZZ>& L)
{ return multiply_list(ZZ(a), L); }

// multiply all integers in L1 by all in L2:
vector<ZZ> multiply_lists(const vector<ZZ>& L1, const vector<ZZ>& L2);

// multiply all integers in L by p^e for e in exponents:
vector<ZZ> multiply_list_by_powers(const ZZ& p, const vector<int>& exponents, const vector<ZZ>& L);
inline vector<ZZ> multiply_list_by_powers(long p, const vector<int>& exponents, const vector<ZZ>& L)
{ return multiply_list_by_powers(ZZ(p), exponents, L); }

#endif
// end of file marith.h
