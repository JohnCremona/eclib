// rat.h: declaration of rational number class
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
 
#if     !defined(_ECLIB_RATIONAL_H)
#define _ECLIB_RATIONAL_H      1       //flags that this file has been included

#include <eclib/arith.h>

class rational {
  friend class bigrational;
public:
        // constructors
        rational(long num_val=0, long den_val=1);
        rational(const rational& q);
        void operator=(const rational& q);

        // rational manipulations
        void cancel();                           // cancel *this in situ
        friend long num(const rational&);        // the numerator
        friend long den(const rational&);        // the denominator
        friend rational recip(const rational&);  // reciprocal
        friend long round(const rational&);      // nearest integer

        // Binary Operator Functions
        friend rational operator+(const rational&, const rational&);
        friend rational operator+(long, const rational&);
        friend rational operator+(const rational&, long);
        friend rational operator-(const rational&, const rational&);
        friend rational operator-(long, const rational&);
        friend rational operator-(const rational&, long);
        friend rational operator*(const rational&, const rational&);
        friend rational operator*(const rational&, long);
        friend rational operator*(long, const rational&);
        friend rational operator/(const rational&, const rational&);
        friend rational operator/(const rational&, long);
        friend rational operator/(long, const rational&);
        friend int operator==(const rational&, const rational&);
        friend int operator!=(const rational&, const rational&);
        friend ostream& operator<< (ostream&s, const rational&);
        friend istream& operator>> (istream& is, rational& r);
        rational& operator+=(const rational&);
        rational& operator+=(long);
        rational& operator-=(const rational&);
        rational& operator-=(long);
        rational& operator*=(const rational&);
        rational& operator*=(long);
        rational& operator/=(const rational&);
        rational& operator/=(long);
        rational operator+() const;
        rational operator-() const;
        friend long floor(const rational& r);
        friend long ceil(const rational& r);
        operator double();  // conversion operator
#ifdef MPFP
        operator bigfloat();  // conversion operator
#endif
// Implementation
private:
        long n, d;
};


// Inline rational functions

inline void rational::cancel()                     // cancel *this in situ
{
 long g = ::gcd(n,d);
 if (g>1) {n/=g; d/=g;}
 if (d<0) {n=-n; d=-d;}
}

inline rational::rational(long num_val, long den_val)
{
  n=num_val; d=den_val;
  (*this).cancel(); 
}

inline rational::rational(const rational& q) :n(q.n), d(q.d) {;}
inline void rational::operator=(const rational& q) {n=q.n; d=q.d;}

inline rational rational::operator+() const
{
        return *this;
}

inline rational rational::operator-() const
{
        return rational(-n, d);
}


// Definitions of compound-assignment operator member functions

inline rational& rational::operator+=(const rational& q2)
{
        n = n*q2.d+d*q2.n;
        d *= q2.d;
        (*this).cancel();
        return *this;
}

inline rational& rational::operator+=(long num_val2)
{
        n += d*num_val2;
        return *this;
}

inline rational& rational::operator-=(const rational& q2)
{
        n = n*q2.d-d*q2.n;
        d *= q2.d;
        (*this).cancel();
        return *this;
}

inline rational& rational::operator-=(long num_val2)
{
        n -= d*num_val2;
        return *this;
}

inline rational& rational::operator*=(const rational& q2)
{
        n *= q2.n;
        d *= q2.d;
        (*this).cancel();
        return *this;
}

inline rational& rational::operator*=(long num_val2)
{
        n*=num_val2;
        (*this).cancel();
        return *this;
}

inline rational& rational::operator/=(long num_val2)
{
        d*=num_val2;
        (*this).cancel();
        return *this;
}

inline rational::operator double() {return double(n)/double(d);}
#ifdef MPFP
inline rational::operator bigfloat() {return to_bigfloat(n)/to_bigfloat(d);}
#endif

// Definitions of non-member rational functions

inline long num(const rational& q)
{
        return q.n;
}

inline long den(const rational& q)
{
        return q.d;
}

inline rational recip(const rational& q)
{
        return rational(q.d, q.n);
}

inline long round(const rational& q)
{
        return q.n / q.d;    //provisional -- should fix rounding direction.
}

// Definitions of non-member binary operator functions

inline rational operator+(const rational& q1, const rational& q2)
{
        return rational(q1.n*q2.d + q2.n*q1.d, q1.d * q2.d);
}

inline rational operator+(long num_val1, const rational& q2)
{
        return rational(num_val1*q2.d + q2.n, q2.d);
}

inline rational operator+(const rational& q1, long num_val2)
{
        return rational(q1.n + num_val2*q1.d, q1.d);
}

inline rational operator-(const rational& q1, const rational& q2)
{
        return rational(q1.n*q2.d - q2.n*q1.d, q1.d * q2.d);
}

inline rational operator-(long num_val1, const rational& q2)
{
        return rational(num_val1*q2.d - q2.n, q2.d);
}

inline rational operator-(const rational& q1, long num_val2)
{
        return rational(q1.n - num_val2*q1.d, q1.d);
}

inline rational operator*(const rational& q1, long num_val2)
{
        return rational(q1.n*num_val2, q1.d);
}

inline rational operator*(long num_val1, const rational& q2)
{
        return rational(q2.n*num_val1, q2.d);
}

inline rational operator*(const rational& q1, const rational& q2)
{
        return rational(q1.n*q2.n, q1.d*q2.d);
}

inline rational operator/(const rational& q1, long num_val2)
{
        return rational(q1.n, q1.d*num_val2);
}

inline rational operator/(const rational& q1, const rational& q2)
{
        return rational(q1.n*q2.d, q1.d*q2.n);
}

inline rational operator/(long num_val1, const rational& q2)
{
  return rational(q2.d*num_val1, q2.n);
}

inline int operator==(const rational& q1, const rational& q2)
{
        return q1.n*q2.d == q2.n*q1.d;
}

inline int operator!=(const rational& q1, const rational& q2)
{
        return q1.n*q2.d != q2.n*q1.d;
}

inline ostream& operator<<(ostream& s, const rational& q)
{
  if(q.d==0) s<<"oo";
  else
    {
      s << q.n;
      if (q.d!=1) {s << "/" << q.d;}
    }
   return s;
}

inline istream& operator>> (istream& is, rational& r)
{
  long n,d=1;
  is>>n;
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
  r=rational(n,d);
  return is;
}

// NB gcd(n,d)=1 and d>0:

inline long floor(const rational& r)
{
  return (r.n-(r.n%r.d))/r.d;
}

inline long ceil(const rational& r)
{
  if(r.d==1) return r.n;
  return 1 + (r.n-(r.n%r.d))/r.d;
}

#endif
