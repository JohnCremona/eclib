// interface.h: used to provide common interface for LiDIA or NTL
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
 
//  The following macro switches can be used: it is best to define them
//  in Makefiles rather than changing this file.
//  1. NTL_INTS     Use NTL for bigints (no bigrationals)
//  2. NTL_ALL      Use NTL also for bigfloats (RR) and bigcomplexes (CC)
//  3. (neither)    Defaults to NTL_ALL

#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#ifdef NTL_ALL
#define NTL_INTS
#define MPFP
#endif

#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
//#include <numeric>
#include <ext/numeric>
#include <iterator>
using namespace std;
#include "templates.h"
#include <stdint.h>

#ifndef MININT
#define MININT numeric_limits<int>::min()
#endif
#ifndef MINLONG
#define MINLONG numeric_limits<long>::min()
#endif
#ifndef MAXINT
#define MAXINT numeric_limits<int>::max()
#endif
#ifndef MAXLONG
#define MAXLONG numeric_limits<long>::max()
#endif

// integers and rationals

// Some of the following were defined for compatibility with LiDIA, which is no longer supported

#include <NTL/ZZ.h>
#include <NTL/ZZXFactoring.h>
using namespace NTL;

#define bigint ZZ
//#define bigrational QQ  // not defined in NTL

#define BIGINT(val) to_ZZ(val)
inline bigint atoI(const char* s) {return to_ZZ(s);}

inline int odd(const int& a) {return a&1;}
inline int even(const int& a) {return !(a&1);}
inline int odd(const long& a) {return a&1;}
inline int even(const long& a) {return !(a&1);}

inline int is_zero(const bigint& x) {return IsZero(x);}
inline int is_positive(const bigint& x) {return sign(x)>0;}
inline int is_negative(const bigint& x) {return sign(x)<0;}
inline int is_one(const bigint& x) {return IsOne(x);}
inline int odd(const bigint& a) {return IsOdd(a);}
inline int even(const bigint& a) {return !IsOdd(a);}
inline void rshift(const bigint& a, long i, bigint& c) {RightShift(c,a,i);}
inline void lshift(const bigint& a, long i, bigint& c) {LeftShift(c,a,i);}
#ifdef setbit
#undef setbit
#endif
inline void setbit(bigint& a, int e) {SetBit(a,e);}
inline long lg(const bigint& x) {return NumBits(x)-1;}
inline int is_long(const bigint& a) {return (a<=MAXLONG)&&(a>=MINLONG);}
inline int is_int(const bigint& a) {return (a<=MAXINT)&&(a>=MININT);}

// The following are not in NTL & need defining
int I2int(const bigint& x);    // too long to inline
long I2long(const bigint& x);  // too long to inline
inline double I2double(const bigint& x) {return to_double(x);}
inline void longasI(long& a, const bigint& x) {a = I2long(x);}
inline void negate(bigint& a) {a=-a;}
inline void sqrt(bigint& a, const bigint& b) {SqrRoot(a,b);}
inline bigint sqrt(const bigint& a) {bigint b; sqrt(b,a); return b;}
inline void square(bigint& a, const bigint& b) {sqr(a,b);}
inline bigint gcd(const bigint& a, const bigint& b) {return GCD(a,b);}
inline bigint lcm(const bigint& a, const bigint& b)
 {if (IsZero(a) && IsZero(b)) return ZZ::zero(); else return a*(b/GCD(a,b));}

// In NTL add, sub, mul, div are defined with result in first place
inline void addx(const bigint& a, const bigint& b, bigint& c)  {add(c,a,b);}
inline void subx(const bigint& a, const bigint& b, bigint& c)  {sub(c,a,b);}
inline void divx(const bigint& a, const bigint& b, bigint& c)  {div(c,a,b);}
inline void mulx(const bigint& a, const bigint& b, bigint& c)  {mul(c,a,b);}
inline bigint pow(const bigint& a, long e)  {return power(a,e);}

//N.B. no power to bigint exponent in NTL
inline long jacobi(const bigint& a, const bigint& p)  {return Jacobi(a,p);}
inline void sqrt_mod_p(bigint & x, const bigint & a, const bigint & p)  
  {SqrRootMod(x,a,p); if(x>(p-x)) x= p-x;}
inline void power_mod(bigint& ans, const bigint& base, const bigint& expo, const bigint& m) 
 {PowerMod(ans,base,expo,m);}
inline void nearest(bigint& c, const bigint& a, const bigint& b) 
 {bigint a0=(a%b);  c = (a-a0)/b; if(2*a0>b) c+=1;}
inline bigint roundover(const bigint& a, const bigint& b)
 {bigint a0=(a%b); bigint c = (a-a0)/b; if(2*a0>b) c+=1; return c;}

#define bigint_mod_long(a,m) (a%m)


// Reals and Complexes

#ifdef NTL_ALL

#include <NTL/RR.h>
#define bigfloat RR
RR Pi();
RR Euler();
RR atan(const RR&);
RR asin(const RR&);
inline RR pow(const RR& a, int e)  {return power(a,e);}
inline RR pow(const RR& a, long e)  {return power(a,e);}
namespace NTL {
inline RR cosh(const RR& x) {return (exp(x)+exp(-x))/2;}
inline RR sinh(const RR& x) {return (exp(x)-exp(-x))/2;}
inline RR tan(const RR& x) {return sin(x)/cos(x);}
RR atan2(const RR&, const RR&);
inline int is_approx_zero(const RR& x)
//  {return abs(x)<power2_RR(2-RR::precision());}
{
  if (IsZero(x)) return 1;
  long n = x.exponent()+RR::precision()-1;
  if (n>=0) return 0;
  // cout<<"x="<<x<<", exponent="<<x.exponent()<<", mantissa="<<x.mantissa()<<", precision="<<RR::precision()<<endl;
  // cout<<"is_approx_zero() returns "<<(x.mantissa()<power2_ZZ(-n))<<endl;
  return abs(x.mantissa())<power2_ZZ(-n);}
} // namespace NTL
#include <complex>
typedef complex<RR> CC;
#define bigcomplex CC

inline void set_precision(long n)
  {RR::SetPrecision(long(n*3.33));RR::SetOutputPrecision(n);}
inline void set_bit_precision(long n)
  {RR::SetPrecision(n);}
inline void set_precision(const string prompt)
  {long n; cerr<<prompt<<": "; cin>>n; set_precision(n);}
inline long decimal_precision() {return long(RR::precision()*0.3);}
inline long bit_precision() {return RR::precision();}
inline int is_approx_zero(const bigcomplex& z)
  {return is_approx_zero(z.real())&&is_approx_zero(z.imag());}
inline RR to_bigfloat(const int& n) {return to_RR(n);}
inline RR to_bigfloat(const long& n) {return to_RR(n);}
inline RR to_bigfloat(const double& x) {return to_RR(x);}
inline RR I2bigfloat(const bigint& x) { return to_RR(x);}
inline int doublify(const bigfloat& x, double& d){ d=to_double(x); return 0;}
inline long longify(bigfloat x) {return to_long(x);}
inline int is_zero(bigfloat x) {return IsZero(x);}
inline int is_zero(bigcomplex z) {return IsZero(z.real()) && IsZero(z.imag());}
inline void Iasb(bigint& a, bigfloat x) {RoundToZZ(a,x);}
inline void Iasb(long& a, bigfloat x) {ZZ n; RoundToZZ(n,x); a=I2long(n);}
istream& operator>>(istream& is, CC& z);
inline CC pow(const CC& a, int e)  {return exp(to_RR(e)*log(a));}
inline CC pow(const CC& a, long e)  {return exp(to_RR(e)*log(a));}
inline CC pow(const CC& a, const RR& e)  {return exp(e*log(a));}


//////////////////////////////////////////////////////////////////
#else  // C doubles and libg++ Complexes
//////////////////////////////////////////////////////////////////

// reals

#define bigfloat double

inline long decimal_precision() {return 15;}
inline int is_zero(double x) {return fabs(x)<1e-15;}
inline int is_approx_zero(double x) {return fabs(x)<1e-10;}
inline void set_precision(long n) {cout.precision(n);}
inline void set_precision(const string prompt)  {cout.precision(15);}
#define Pi()    (3.1415926535897932384626433832795028841)
#define Euler() (0.57721566490153286060651209008240243104)
inline double round(double x) {return floor(x+0.5);}
inline void Iasb(long& a, double x) {a = (long)x;}
inline long longify(double x) {return (long)x;}
inline int doublify(const bigfloat& x, double& d) {d=x; return 0;}


inline double to_bigfloat(const int& n) {return double(n);}
inline double to_bigfloat(const long& n) {return double(n);}
inline double to_bigfloat(const double& x) {return x;}

inline bigfloat I2bigfloat(const bigint& x) {return I2double(x);}

// complexes

#include <complex>
#define bigcomplex complex<double>

inline int is_zero(const bigcomplex& z) 
 {return is_zero(z.real())&&is_zero(z.imag());}
inline int is_approx_zero(const bigcomplex& z)
 {return is_approx_zero(z.real())&&is_approx_zero(z.imag());}

//////////////////////////////////////////////////////////////////
#endif // NTL_ALL
//////////////////////////////////////////////////////////////////

#undef setbit

#endif // #define _INTERFACE_H_
