// marith.h: declarations of integer arithmetic functions (multiprecision)
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
 
#ifndef _ECLIB_MARITH_H
#define _ECLIB_MARITH_H

#include <eclib/arith.h>

bigint bezout(const bigint& aa, const bigint& bb, bigint& xx, bigint& yy);
int divides(const bigint& a, const bigint& b, bigint& q, bigint& r);
int divides(const bigint& a, long b, bigint& q, long& r);
// returns 1 iff remainder r==0

bigint Iround(bigfloat x);
bigint Ifloor(bigfloat x);
bigint Iceil(bigfloat x);

bigint posmod(const bigint& a, const bigint& b);  // a mod b in range 0--(b-1)
long posmod(const bigint& a, long b);
int isqrt(const bigint& a, bigint& root);
int divide_exact(const bigint& aa, const bigint& bb, bigint& c);
     // c = a/b with error message if remainder is non-zero
long divide_out(bigint& a, const bigint& d);
long divide_out(bigint& a, long d);
// divides a by d as many times as possible returning number of times (but none if a=0!)

bigint show(const bigint& a);
vector<bigint> show(const vector<bigint>& a);


// extra_primes is a set holding extra big primes found or added manually.  
class extra_prime_class {
public:
  std::set<bigint> the_primes;
  extra_prime_class() {;}
  ~extra_prime_class();
  void read_from_file(const string pfilename, int verb=0);
  void write_to_file(const string pfilename, int verb=0);
  int contains(const bigint& p)
  {
    return the_primes.find(p)!=the_primes.end();
  }
  void add(const bigint& p)  
  {
    if(p>maxprime()) the_primes.insert(p);
  }
  void show() 
  {
    cout << "Extra primes in list: ";
    copy(the_primes.begin(),the_primes.end(), ostream_iterator<bigint>(cout, " "));
    cout << endl;
  }
};

extern extra_prime_class the_extra_primes;  // The one and only instance

void initprimes(const string pfilename, int verb=0);

//
// divisors
//
// The following uses trial division:
vector<bigint> pdivs_trial(const bigint& number, int trace=0);
// The following uses gp factorization externally if available:
vector<bigint> pdivs_gp(const bigint& number, int trace=0);
// The following uses pari library's factorization if available:
vector<bigint> pdivs_pari(const bigint& number, int trace=0);
// The following uses one of the above
vector<bigint> pdivs(const bigint& number, int trace=0);

vector<bigint> posdivs(const bigint& number, const vector<bigint>& plist);
vector<bigint> posdivs(const bigint& number);

vector<bigint> alldivs(const bigint& number, const vector<bigint>& plist);
vector<bigint> alldivs(const bigint& number);

vector<bigint> sqdivs(const bigint& number, const vector<bigint>& plist);
vector<bigint> sqdivs(const bigint& number);

vector<bigint> sqfreedivs(const bigint& number, const vector<bigint>& plist);
vector<bigint> sqfreedivs(const bigint& number);

void sqfdecomp(const bigint& a, bigint& a1, bigint& a2, vector<bigint>& plist, int trace_fact=0);
     // a must be non-zero, computes square-free a1 and a2>0 such that a=a1*a2^2
     // plist will hold prime factors of a

void sqfdecomp(const bigint& a, vector<bigint>& plist, bigint& a1, bigint& a2);
    // a must be non-zero, computes square-free a1 and a2>0 such that a=a1*a2^2
    // plist already holds prime factors of a

// Given a, b, lem3 returns m1 etc so that a=c1^2*m1*m12, b=c2^2*m2*m12 
// with m1, m2, m12 pairwise coprime.   At all  times these equations hold, 
// and at each step the product m1*m2*m12 is decreased by a factor d, 
// so the process terminates when the coprimality condition is satisfied. 
void rusin_lem3(const bigint& a, const bigint& b,
	  bigint& m1, bigint& m2, bigint& m12, bigint& c1, bigint& c2);

// Solves x-a1(mod m1), x=a2(mod m2)
bigint chrem(const bigint& a1, const bigint& a2, 
	     const bigint& m1, const bigint& m2);

//
// general purpose routines -- bigint overloads
//

bigint mod(const bigint& a, const bigint& b);     // a mod b in range +- half b
long mod(const bigint& a, long b);

long val(const bigint& factor, const bigint& number);
long val(long factor, const bigint& number);

int div(const bigint& factor, const bigint& number);
int div(long factor, const bigint& number);
inline int ndiv(const bigint& factor, const bigint& number)
  { return !div(factor, number);}
inline int ndiv(long factor, const bigint& number)
  { return !div(factor, number);}


long bezout(const bigint& aa, long bb, bigint& xx, bigint& yy);
bigint invmod(const bigint& a, const bigint& p);  // -a mod p
long invmod(const bigint& a, long p);

int m1pow(const bigint& a);
int chi2(const bigint& a);
int chi4(const bigint& a);
int hilbert2(const bigint& a, const bigint& b);
int hilbert2(const bigint& a, long b);
int hilbert2(long a, const bigint& b);
int legendre(const bigint& a, const bigint& b);
int legendre(const bigint& a, long b);
int kronecker(const bigint& d, const bigint& n);
int kronecker(const bigint& d, long n);
// See hilbert.h for hilbert symbol functions

long gcd(const bigint& a, long b);

int modrat(const bigint& n, const bigint& m, const bigint& lim, 
           /* return values: */ bigint& a, bigint& b);

int sqrt_mod_2_power(bigint& x, const bigint& a, int e);
int sqrt_mod_p_power(bigint& x, const bigint& a, const bigint& p, int e);
int sqrt_mod_m(bigint& x, const bigint& a, const bigint& m);
int sqrt_mod_m(bigint& x, const bigint& a, const bigint& m, const vector<bigint>& mpdivs);
// Second version of sqrt_mod_m requires mpdivs to hold a list
// of primes which contains all those which divide m; if it acontains
// extras that does not matter.

int modsqrt(const bigint& a, const vector<bigint>& bplist, bigint& x);
     // Solves x^2=a mod b with b square-free, returns success/fail


// root-finding functions:

// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)
// roots are put in r which should be allocated of size 3
int nrootscubic(long b, long c, long d, long p, long* roots);

void ratapprox(bigfloat x, bigint& a, bigint& b, const bigint& maxd=BIGINT(0));
void ratapprox(bigfloat x, long& a, long& b, long maxd=0);

int is_nth_power(const bigint& x, int n);
bigint prime_to_S_part(const bigint& x,  const vector<bigint>& S);
int is_S_unit(const bigint& x,  const vector<bigint>& S);

#endif
// end of file marith.h
