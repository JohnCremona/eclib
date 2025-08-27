// arith.h: declarations of arithmetic functions (single precision)
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
 
#ifndef _ECLIB_ARITH_H
#define _ECLIB_ARITH_H      1
                           //flags that this file has been included
#include "xmod.h"

/* Prime number class; adapted from Pari  */

typedef unsigned char *byteptr;

class primeclass {
  friend class primevar;
  byteptr pdiffptr;
  long NPRIMES, BIGGESTPRIME;
  byteptr p_aptr;  // points to "current" prime
  long p_ind;      // index  of "current" prime
  long p_val;      // value  of "current" prime

public:
  primeclass();  // will use 10^6 as default or read from file
  explicit primeclass(long maxnum);
  ~primeclass();
  void init(long maxnum); // called in constructor, or to make more primes
  long number(long n) ; // returns n'th prime (n=1 gives p=2)
  vector<long> getfirst(long n); // return primes 2..p_n as vector<long>
  void reset(void);
  int at_end(void);
  int advance(void);
  friend long nprimes(void);
  friend long maxprime(void);
};

extern primeclass the_primes;  // The one and only instance

inline long prime_number (long n)   /* returns n'th prime from global list */
  {return the_primes.number(n);}
inline vector<long> primes (long n)  /* returns list of first n primes */
  {return the_primes.getfirst(n);}
inline long nprimes(void) {return the_primes.NPRIMES;}
inline long maxprime(void) {return the_primes.BIGGESTPRIME;}
long prime_pi(long p); // returns i>=0 such that p is the i'th prime

class primevar {
public:
        long val;        /* current value */
        long ind;        /* current index */
private:
        byteptr ndiff;     /* pointer to next diff*/
        long maxindex;     /* max index */
public:
        primevar(long max=the_primes.NPRIMES, long i=1) 
          {maxindex=max; ind=i; val=the_primes.number(i); 
           ndiff=the_primes.pdiffptr+i;}
        void init(long max=the_primes.NPRIMES, long i=1) 
          {maxindex=max; ind=i; val=the_primes.number(i); 
           ndiff=the_primes.pdiffptr+i;}
        void operator++() {if ((ind++)<=maxindex) { val+=*ndiff++;}}
        void operator++(int) {if ((ind++)<=maxindex) { val+=*ndiff++;}}
        int ok() const {return ind<=maxindex;}
        int more() const {return ind<maxindex;}
        long value() const {return val;}
        long index() const {return ind;}
        operator long() const {return val;}
};


/* Usage of primevar: 
---To loop through first n primes:

   long p; 
   for(primevar pr(n); pr.ok(); pr++) {p=pr; ... ;}  // or:
   for(pr.init(n); pr.ok(); pr++) {p=pr; ... ;}  // iff pr is existing primevar

---To loop through all primes:

   primevar pr; //or for existing primevar, pr.init();
   long p;
   while(pr.ok()) {p=pr; pr=++; ...;}

*/

long primdiv(long);        /* returns first prime divisor */
vector<long> pdivs(long);      /* list of prime divisors */

vector<long> posdivs(long, const vector<long>& plist);  // all positive divisors
inline vector<long> posdivs(long n)
{
  return posdivs(n, pdivs(n));
}

vector<long> alldivs(long, const vector<long>& plist);  // absolutely all divisors
inline vector<long> alldivs(long n)
{
  return alldivs(n, pdivs(n));
}

vector<long> sqdivs(long, const vector<long>& plist);   //  divisors whose square divides
inline vector<long> sqdivs(long n)
{
  return sqdivs(n, pdivs(n));
}

vector<long> sqfreedivs(long, const vector<long>& plist); // square-free divisors
inline vector<long> sqfreedivs(long n)
{
  return sqfreedivs(n, pdivs(n));
}

// utilities for compatibility with bigint args
inline int odd(const int& a) {return a&1;}
inline int even(const int& a) {return !(a&1);}
inline int odd(const long& a) {return a&1;}
inline int even(const long& a) {return !(a&1);}
inline int is_zero(int n) {return n==0;}
inline int is_zero(long n) {return n==0;}
inline int is_nonzero(int n) {return n!=0;}
inline int is_nonzero(long n) {return n!=0;}
inline int is_one(int n) {return n==1;}
inline int is_one(long n) {return n==1;}
inline long I2long(long n) {return n;}

long mod(long a, long modb); /* modulus in range plus or minus half mod */
long mod(int a, long modb); /* modulus in range plus or minus half mod */
int mod(int a, int modb); /* modulus in range plus or minus half mod */
int mod(long a, int modb); /* modulus in range plus or minus half mod */
long posmod(long a, long modb); /* ordinary modulus, but definitely positive */
long posmod(int a, long modb); /* ordinary modulus, but definitely positive */
int posmod(int a, int modb); /* ordinary modulus, but definitely positive */
int posmod(long a, int modb); /* ordinary modulus, but definitely positive */

// gcc division truncates towards 0, while we need rounding, with a
// consistent behaviour for halves (they go up here).
//
// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
long rounded_division(long a, long b);

long gcd(long, long);
int gcd(int, int);
long lcm(long, long);
long bezout(long, long, long&, long&);
int bezout(int aa, int bb, int& xx, int& yy);

long invmod(long, long);

inline int xmm(int a, int b, int m)
{
  if (a==1) return b;
  if (a==-1) return -b;
  if (b==1) return a;
  if (b==-1) return -a;
  return (a*(int64_t)b) % m;
}

inline int addmod(int a, int b, int m)
{
  if (is_zero(a)) return b;
  if (is_zero(b)) return a;
  return mod(a+b,m);
}

inline long addmod(long a, long b, long m)
{
  if (is_zero(a)) return b;
  if (is_zero(b)) return a;
  return mod(a+b,m);
}

inline long xmm(long a, long b, long m)
{
  if (a==1) return b;
  if (a==-1) return -b;
  if (b==1) return a;
  if (b==-1) return -a;
  return (a*(int64_t)b) % m;
}

// Assuming a*d-b*c!=0, computes a reduced Z-basis for <(a,b),(c,d)>
void gauss_reduce(long a0, long b0, long c0, long d0,
                  long& a, long& b, long& c, long& d);

// Set a, b so that a/b=n (mod m) with |a|, |b| minimal; return success if a^2, b^2 <= m/2
int modrat(long n, long m, long& a, long& b);
int modrat(int n, int m, int& a, int& b);

long val(long factor, long number); // order of factor in number

inline int divides(long factor,long number) 
{
  return (factor==0) ? (number==0) : (number%factor==0);
}

inline int ndivides(long factor,long number)
{
  return (factor==0) ? (number!=0) : (number%factor!=0);
}

// a=b*q+r, return 1 iff r==0
int divrem(long a, long b, long& q, long& r);
int divrem(int a, int b, int& q, int& r);

inline int m1pow(long a) {return (a%2 ?  -1 : +1);}
inline int sign(long a)   {return (a==0? 0: (a>0? 1: -1));}

// set root to rounded sqrt(a) if a>=0, return 1 iff exact
int isqrt(long a, long& root);
// return rounded sqrt(a) (undefined for a<0)
long isqrt(const long a);

int chi2(long a);
int chi4(long a);
int hilbert2(long a, long b);
int legendre(long a, long b);
int kronecker(long d, long n);
int intlog2(long& n, long& e, int roundup);

int is_squarefree(long n);
int is_valid_conductor(long n);

long squarefree_part(long d);

// return list of integers from first to last inclusive
vector<long> range(long first, long last);

#endif
