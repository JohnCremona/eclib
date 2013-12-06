// arith.h: declarations of arithmetic functions (single precision)
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
 
#ifndef _ARITH_H
#define _ARITH_H      1
                           //flags that this file has been included
#include "interface.h"
#include <cstring> // for memset gcc >= 4.3
#include "xmod.h" // supercedes the macros which used to be here

/* Prime number class; adapted from Pari  */

typedef unsigned char *byteptr;

class primeclass {
  friend class primevar;
  byteptr pdiffptr;
  long NPRIMES, BIGGESTPRIME;

  void reset(void);
  int at_end(void);
  int advance(void); 
  byteptr p_aptr;  // points to "current" prime
  long p_ind;      // index  of "current" prime
  long p_val;      // value  of "current" prime

public:
  primeclass();  // will use 10^6 as default or read from file
  primeclass(long maxnum);
  ~primeclass();
  void init(long maxnum); // called in constructor, or to make more primes
  long number(long n) ; // returns n'th prime (n=1 gives p=2)
  vector<long> getfirst(long n); // return primes 2..p_n as vector<long>
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
  const vector<long>& plist = pdivs(n);
  return posdivs(n,plist);
}

vector<long> alldivs(long, const vector<long>& plist);  // absolutely all divisors
inline vector<long> alldivs(long n)
{
  const vector<long>& plist = pdivs(n);
  return alldivs(n,plist);
}

vector<long> sqdivs(long, const vector<long>& plist);   //  divisors whose square divides
inline vector<long> sqdivs(long n)
{
  const vector<long>& plist = pdivs(n);
  return sqdivs(n,plist);
}

vector<long> sqfreedivs(long, const vector<long>& plist); // square-free divisors
inline vector<long> sqfreedivs(long n)
{
  const vector<long>& plist = pdivs(n);
  return sqfreedivs(n,plist);
}


long mod(long a, long modb); /* modulus in range plus or minus half mod */
long posmod(long a, long modb); /* ordinary modulus, but definitely positive */

#ifndef NTL_INTS
inline int abs(int a) {return (a<0) ? (-a) : a;}
#endif
long gcd(long, long);
int gcd(int, int);
long lcm(long, long);
long bezout(long, long, long&, long&);
int intbezout(int aa, int bb, int& xx, int& yy);

long invmod(long, long);
int modrat(long, long, float, long&, long&);
int modrat(int, int, float, int&, int&);

inline int is_zero(long n) {return n==0;}

long val(long factor, long number); // order of factor in number

inline int divides(long factor,long number) 
{ 
  return (factor==0) ? (number==0) : (number%factor==0);
}

inline int ndivides(long factor,long number) 
{
  return (factor==0) ? (number!=0) : (number%factor!=0);
  //  return !::divides(factor,number);
}

inline long m1pow(long a) {return (a%2 ?  -1 : +1);}
inline int sign(long a)   {return (a==0? 0: (a>0? 1: -1));} 
inline int sign(double a) {return (a==0? 0: (a>0? 1: -1));} 

long chi2(long a);
long chi4(long a);
long hilbert2(long a, long b);
long legendre(long a, long b);
long kronecker(long d, long n);
int intlog2(long& n, long& e, int roundup);

int is_squarefree(long n);
int is_valid_conductor(long n);


#endif
