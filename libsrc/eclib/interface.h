// interface.h: used to provide common interface
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
 
//  The macro NO_MPFP can optionally be set.  It controls whether most
//   real and complex floating point functions are implemented using
//   doubles and complex doubles (if set) or using NTL bigfloats
//   (RR) and bigcomplexes (if not set, the default)

#ifndef _ECLIB_INTERFACE_H_
#define _ECLIB_INTERFACE_H_

#ifndef NO_MPFP // the default is for this *not* to be set
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
#ifdef _LIBCPP_VERSION
#include <numeric>
#else
#include <ext/numeric>
#endif
#include <iterator>

#include <eclib/templates.h>
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

const double LOG_10_2=0.30102999566398114L;

// integers and rationals

// Some of the following were defined for compatibility with LiDIA, which is no longer supported

#include <NTL/ZZ.h>
#include <NTL/ZZXFactoring.h>
using namespace NTL;

typedef ZZ bigint;


// Reals and Complexes

#ifdef MPFP

#include <NTL/RR.h>
typedef RR bigfloat;
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
{
  if (IsZero(x)) return 1;
  long n = x.exponent()+RR::precision()-1;
  if (n>=0) return 0;
  // cout<<"x="<<x<<", exponent="<<x.exponent()<<", mantissa="<<x.mantissa()<<", precision="<<RR::precision()<<endl;
  // cout<<"is_approx_zero() returns "<<(x.mantissa()<power2_ZZ(-n))<<endl;
  return abs(x.mantissa())<power2_ZZ(-n);
}
} // namespace NTL

#include "bigcomplex.h"

// NB Internal precision is always bit-precision.  We used to use
// decimal precision in the user interface with a conversion factor of
// 3.33 (approx log(10)/log(2)).  NTL has a default internal bit
// precision of 150 which can be changed using RR::SetPrecision(), and
// RR:SetOutputPrecision(d) sets the output to d decimal places
// (default 10).  See www.shoup.net/ntl/doc/tour-ex6.html and
// www.shoup.net/ntl/doc/RR.cpp.html.

// Set internal precision to n bits and output precision to (log_10(2)*n)-1 decimal places
inline void set_precision(long n)
{RR::SetPrecision(n); RR::SetOutputPrecision(long(LOG_10_2*n)-1);}

// Mostly for backward compatibility (used in saturate.cc) or for
// temporarily changing internal precision when no output is relevant:
inline void set_bit_precision(long n)
  {RR::SetPrecision(n);}

// Prompt user for precision
inline void set_precision(const string prompt)
  {long n; cerr<<prompt<<": "; cin>>n; set_precision(n);}

// read current precision converted to decimal (approximately)
inline long decimal_precision() {return long(RR::precision()*LOG_10_2);}

// read current bit precision
inline long bit_precision() {return RR::precision();}

inline int is_approx_zero(const bigcomplex& z)
  {return is_approx_zero(z.real())&&is_approx_zero(z.imag());}
inline RR to_bigfloat(const int& n) {return to_RR(n);}
inline RR to_bigfloat(const long& n) {return to_RR(n);}
inline RR to_bigfloat(const double& x) {return to_RR(x);}
inline RR I2bigfloat(const bigint& x) { return to_RR(x);}
inline double I2double(const bigint& x) {return to_double(x);}
inline int doublify(const bigfloat& x, double& d){ d=to_double(x); return 0;}
int longify(const bigfloat& x, long& a, int rounding=0);
inline int is_real_zero(bigfloat x) {return IsZero(x);}
inline int is_complex_zero(bigcomplex z) {return IsZero(z.real()) && IsZero(z.imag());}
inline void Iasb(bigint& a, bigfloat x) {RoundToZZ(a,x);}
inline void Iasb(long& a, bigfloat x) {ZZ n; RoundToZZ(n,x); a=I2long(n);}
istream& operator>>(istream& is, bigcomplex& z);
inline bigcomplex pow(const bigcomplex& a, int e)  {return (to_RR(e)*a.log()).exp();}
inline bigcomplex pow(const bigcomplex& a, long e)  {return (to_RR(e)*a.log()).exp();}
inline bigcomplex pow(const bigcomplex& a, const RR& e)  {return (e*a.log()).exp();}

//////////////////////////////////////////////////////////////////
#else  // C doubles and libg++ Complexes
//////////////////////////////////////////////////////////////////

// reals

typedef double bigfloat;

inline long decimal_precision() {return 14;}
inline long bit_precision() {return 52;}
inline int is_real_zero(double x) {return fabs(x)<1e-14;}
inline int is_approx_zero(double x) {return fabs(x)<1e-11;}
inline int sign(double x) {return (double(0) < x) - (x < double(0));}

// We cannot set internal bit precision in this mode, so we just set the output decimal precision
inline void set_precision(long n) {cout.precision(min(14,long(LOG_10_2*n)));}
inline void set_precision(const string prompt)  {cout.precision(14);}
#define Pi()    (double)(3.1415926535897932384626433832795028841)
#define Euler() (double)(0.57721566490153286060651209008240243104)

inline double round(double x) {return floor(x+0.5);}
inline double roundup(double x) {return ceil(x-0.5);}
inline void Iasb(long& a, double x) {a = (long)x;}
// return value is 1 for success, else 0
int longify(double x, long& a, int rounding=0);
inline int doublify(const bigfloat& x, double& d) {d=x; return 0;}
inline double inv(double x) {return 1/x;}
inline double sqr(double x) {return x*x;}
inline double power2_RR(long e) {return 1<<e;}
inline double power(double x, long n) {return pow(x,n);}

inline double to_bigfloat(const int& n) {return double(n);}
inline double to_bigfloat(const long& n) {return double(n);}
inline double to_bigfloat(const double& x) {return x;}
inline double I2double(const bigint& x) {return to_double(x);}
inline double I2bigfloat(const bigint& x) { return to_double(x);}

// complexes

#include <complex>
typedef complex<double> bigcomplex;

inline int is_complex_zero(const bigcomplex& z) 
 {return is_real_zero(z.real())&&is_real_zero(z.imag());}
inline int is_approx_zero(const bigcomplex& z)
 {return is_approx_zero(z.real())&&is_approx_zero(z.imag());}

//////////////////////////////////////////////////////////////////
#endif // MPFP
//////////////////////////////////////////////////////////////////

#undef setbit

// Utility to return a string from an environment variable with a
// default to use if the variable is not set or empty.

string getenv_with_default(string env_var, string def_val);

inline int is_long(const bigint& a) {return (a<=MAXLONG)&&(a>=MINLONG);}
inline int is_int(const bigint& a) {return (a<=MAXINT)&&(a>=MININT);}
int I2int(const bigint& x);    // too long to inline
long I2long(const bigint& x);  // too long to inline

#endif // #define _INTERFACE_H_
