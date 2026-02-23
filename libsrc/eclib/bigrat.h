// bigrat.h: declaration of bigrational number class
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
 
#if     !defined(_ECLIB_BIGRATIONAL_H)
#define _ECLIB_BIGRATIONAL_H      1       //flags that this file has been included

#include "rat.h"
#include "marith.h"

class bigrational {

public:
        // constructors
        bigrational() : n(0), d(1) {;}
        explicit bigrational(const ZZ& x) : n(x), d(1) {;}
        bigrational(const ZZ& x, const ZZ& y) : n(x), d(y) { cancel();}
        bigrational(const bigrational& q) :n(q.n), d(q.d) {;}
        explicit bigrational(const rational& q) :n(q.n), d(q.d) {;}
        void operator=(const bigrational& q) {n=q.n; d=q.d;}
        void operator=(const ZZ& a) {n=a; d=ZZ(1);}
        void operator=(const rational& q) {n=ZZ(q.n); d=ZZ(q.d);}

        // bigrational manipulations
        void cancel();                           // cancel *this in situ
        friend ZZ num(const bigrational&q) {return q.n;} // the numerator
        ZZ num() const {return n;}
        friend ZZ den(const bigrational&q) {return q.d;}  // the denominator
        ZZ den() const {return d;}
        friend bigrational recip(const bigrational&q){return bigrational(q.d, q.n);}  // reciprocal
        friend ZZ round(const bigrational&q){return q.n / q.d;}  // nearest integer

        // Binary Operator Functions
        bigrational operator+() const {return *this;}
        friend bigrational operator+(const bigrational&, const bigrational&);
        friend bigrational operator+(const ZZ&, const bigrational&);
        friend bigrational operator+(const bigrational&, const ZZ&);
        bigrational& operator+=(const bigrational&);
        bigrational& operator+=(const ZZ&);

        bigrational operator-() const {return bigrational(-n, d);}
        friend bigrational operator-(const bigrational&, const bigrational&);
        friend bigrational operator-(const ZZ&, const bigrational&);
        friend bigrational operator-(const bigrational&, const ZZ&);
        bigrational& operator-=(const bigrational&);
        bigrational& operator-=(const ZZ&);

        friend bigrational operator*(const bigrational&, const bigrational&);
        friend bigrational operator*(const bigrational&, const ZZ&);
        friend bigrational operator*(const ZZ&, const bigrational&);
        bigrational& operator*=(const bigrational&);
        bigrational& operator*=(const ZZ&);

        friend bigrational operator/(const bigrational&, const bigrational&);
        friend bigrational operator/(const bigrational&, const ZZ&);
        friend bigrational operator/(const ZZ&, const bigrational&);
        bigrational& operator/=(const bigrational&);
        bigrational& operator/=(const ZZ&);

        friend int operator==(const bigrational&, const bigrational&);
        friend int operator!=(const bigrational&, const bigrational&);
        friend int operator<(const bigrational&, const bigrational&);
        friend int operator>(const bigrational&, const bigrational&);

        friend ostream& operator<< (ostream&s, const bigrational&);
        friend istream& operator>> (istream& is, bigrational& r);

        friend ZZ floor(const bigrational& r);
        friend ZZ ceil(const bigrational& r);
        operator bigfloat() const {return I2bigfloat(n)/I2bigfloat(d);}  // conversion operator

        int is_zero() const {return ::is_zero(n);}
        int is_one() const {return ::is_one(n) && ::is_one(d);}
        int is_1728() const {return ::is_zero(n-1728*d);}

        int is_square(bigrational& r) const
        {
         ZZ x;
         int res = isqrt(n*d, x);
         if (res) // sqrt(n/d) = sqrt(n*d)/d
           r = bigrational(x,d);
         return res;
        }
        bigrational squarefree_part() const {return bigrational(::squarefree_part(n), ::squarefree_part(d));}
        friend inline bigrational squarefree_product(const bigrational& r, const bigrational& s)
        {return (r*s).squarefree_part();}
        // For a non-zero, computes square-free a1 and a2>0 such that a=a1*a2^2
        friend void sqfdecomp(const bigrational& a, bigrational& a1, bigrational& a2);

  // Implementation
private:
        ZZ n, d;
};


// Inline bigrational functions

inline void bigrational::cancel()                     // cancel *this in situ
{
  ZZ g(gcd(n,d));
  if (g>1) {n/=g; d/=g;}
  if (d<0) {n=-n; d=-d;}
}

inline bigrational& bigrational::operator+=(const bigrational& q2)
{
        n = n*q2.d+d*q2.n;
        d *= q2.d;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator+=(const ZZ& x)
{
        n += d*x;
        return *this;
}

inline bigrational& bigrational::operator-=(const bigrational& q2)
{
        n = n*q2.d-d*q2.n;
        d *= q2.d;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator-=(const ZZ& x)
{
        n -= d*x;
        return *this;
}

inline bigrational& bigrational::operator*=(const ZZ& x)
{
        n*=x;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator*=(const bigrational& r)
{
        n*=r.n;
        d*=r.d;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator/=(const ZZ& x)
{
        d*=x;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator/=(const bigrational& r)
{
        d*=r.n;
        n*=r.d;
        cancel();
        return *this;
}

// Definitions of non-member binary operator functions

inline bigrational operator+(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d + q2.n*q1.d, q1.d * q2.d);
}

inline bigrational operator+(const ZZ& x1, const bigrational& q2)
{
        return bigrational(x1*q2.d + q2.n, q2.d);
}

inline bigrational operator+(const bigrational& q1, const ZZ& x)
{
        return bigrational(q1.n + x*q1.d, q1.d);
}

inline bigrational operator-(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d - q2.n*q1.d, q1.d * q2.d);
}

inline bigrational operator-(const ZZ& x1, const bigrational& q2)
{
        return bigrational(x1*q2.d - q2.n, q2.d);
}

inline bigrational operator-(const bigrational& q1, const ZZ& x)
{
        return bigrational(q1.n - x*q1.d, q1.d);
}

inline bigrational operator*(const bigrational& q1, const ZZ& x)
{
        return bigrational(q1.n*x, q1.d);
}

inline bigrational operator*(const ZZ& x1, const bigrational& q2)
{
        return bigrational(q2.n*x1, q2.d);
}

inline bigrational operator*(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.n, q1.d*q2.d);
}

inline bigrational operator/(const bigrational& q1, const ZZ& x)
{
        return bigrational(q1.n, q1.d*x);
}

inline bigrational operator/(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d, q1.d*q2.n);
}

inline bigrational operator/(const ZZ& x1, const bigrational& q2)
{
  return bigrational(q2.d*x1, q2.n);
}

inline int operator==(const bigrational& q1, const bigrational& q2)
{
        return q1.n*q2.d == q2.n*q1.d;
}

inline int operator!=(const bigrational& q1, const bigrational& q2)
{
        return q1.n*q2.d != q2.n*q1.d;
}

inline int operator<(const bigrational& q1, const bigrational& q2)
{
        return q1.n*q2.d < q2.n*q1.d;
}

inline int operator>(const bigrational& q1, const bigrational& q2)
{
        return q1.n*q2.d > q2.n*q1.d;
}

inline ostream& operator<<(ostream& s, const bigrational& q)
{
  if(q.d==0) s<<"oo";
  else
    {
      s << q.n;
      if (q.d!=1) {s << "/" << q.d;}
    }
   return s;
}


inline istream& operator>> (istream& is, bigrational& r)
{
  ZZ n,d(1);
  is>>n>>ws;
  if(!is.eof())
    {
      char c;
      is.get(c);
      if(c=='/')
	{
	  is>>d;
	}
      else
	{
	  is.putback(c);
	}
    }
  r=bigrational(n,d);
  return is;
}
// NB gcd(n,d)=1 and d>0:

inline ZZ floor(const bigrational& r)
{
  return (r.n-(r.n%r.d))/r.d;
}

inline ZZ ceil(const bigrational& r)
{
  if(is_one(r.d)) return r.n;
  return 1 + (r.n-(r.n%r.d))/r.d;
}

inline int is_S_integral(const bigrational& x,  const vector<ZZ>& S)
{
  return is_S_unit(den(x), S);
}

inline int is_S_unit(const bigrational& x,  const vector<ZZ>& S)
{
  return is_S_unit(den(x), S) && is_S_unit(num(x), S);
}

inline bigrational prime_to_S_part(const bigrational& x,  const vector<ZZ>& S)
{
  return bigrational(prime_to_S_part(num(x), S), prime_to_S_part(den(x), S));
}

inline void sqfdecomp(const bigrational& a, bigrational& a1, bigrational& a2)
{
  if (a.is_zero())
    return;
  vector<ZZ> plist;
  sqfdecomp(a.n, a1.n, a2.n, plist);
  sqfdecomp(a.d, a1.d, a2.d, plist);
}
#endif
