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

#include <eclib/rat.h>
#include <eclib/marith.h>

class bigrational {

public:
        // constructors
        bigrational() : n(0), d(1) {;}
        explicit bigrational(const bigint& num_val) : n(num_val), d(1) {;}
        bigrational(const bigint& num_val, const bigint& den_val);
        bigrational(const bigrational& q);
        explicit bigrational(const rational& q);
        void operator=(const bigrational& q);
        void operator=(const bigint& a);
        void operator=(const rational& q);

        // bigrational manipulations
        void cancel();                           // cancel *this in situ
        friend bigint num(const bigrational&);        // the numerator
        friend bigint den(const bigrational&);        // the denominator
        friend bigrational recip(const bigrational&);  // reciprocal
        friend bigint round(const bigrational&);      // nearest integer

        // Binary Operator Functions
        friend bigrational operator+(const bigrational&, const bigrational&);
        friend bigrational operator+(const bigint&, const bigrational&);
        friend bigrational operator+(const bigrational&, const bigint&);
        friend bigrational operator-(const bigrational&, const bigrational&);
        friend bigrational operator-(const bigint&, const bigrational&);
        friend bigrational operator-(const bigrational&, const bigint&);
        friend bigrational operator*(const bigrational&, const bigrational&);
        friend bigrational operator*(const bigrational&, const bigint&);
        friend bigrational operator*(const bigint&, const bigrational&);
        friend bigrational operator/(const bigrational&, const bigrational&);
        friend bigrational operator/(const bigrational&, const bigint&);
        friend bigrational operator/(const bigint&, const bigrational&);
        friend int operator==(const bigrational&, const bigrational&);
        friend int operator!=(const bigrational&, const bigrational&);
        friend int operator<(const bigrational&, const bigrational&);
        friend int operator>(const bigrational&, const bigrational&);
        friend ostream& operator<< (ostream&s, const bigrational&);
        friend istream& operator>> (istream& is, bigrational& r);
        bigrational& operator+=(const bigrational&);
        bigrational& operator+=(const bigint&);
        bigrational& operator-=(const bigrational&);
        bigrational& operator-=(const bigint&);
        bigrational& operator*=(const bigrational&);
        bigrational& operator*=(const bigint&);
        bigrational& operator/=(const bigrational&);
        bigrational& operator/=(const bigint&);
        bigrational operator+();
        bigrational operator-();
        friend bigint floor(const bigrational& r);
        friend bigint ceil(const bigrational& r);
        operator bigfloat();  // conversion operator

        int is_zero() const {return ::is_zero(n);}
        int is_1728() const {return ::is_zero(n-1728*d);}

  // Implementation
private:
        bigint n, d;
};


// Inline bigrational functions

inline void bigrational::cancel()                     // cancel *this in situ
{
  bigint g(gcd(n,d));
  if (g>1) {n/=g; d/=g;}
  if (d<0) {n=-n; d=-d;}
}

inline bigrational::bigrational(const bigint& num_val, const bigint& den_val)
  : n(num_val), d(den_val)
{
  cancel(); 
}

inline bigrational::bigrational(const bigrational& q) :n(q.n), d(q.d) {;}
inline bigrational::bigrational(const rational& q) :n(q.n), d(q.d) {;}
inline void bigrational::operator=(const bigrational& q) {n=q.n; d=q.d;}
inline void bigrational::operator=(const bigint& a) {n=a; d=bigint(1);}
inline void bigrational::operator=(const rational& q) {n=bigint(q.n); d=bigint(q.d);}

inline bigrational bigrational::operator+()
{
        return *this;
}

inline bigrational bigrational::operator-()
{
        return bigrational(-n, d);
}


// Definitions of compound-assignment operator member functions

inline bigrational& bigrational::operator+=(const bigrational& q2)
{
        n = n*q2.d+d*q2.n;
        d *= q2.d;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator+=(const bigint& num_val2)
{
        n += d*num_val2;
        return *this;
}

inline bigrational& bigrational::operator-=(const bigrational& q2)
{
        n = n*q2.d-d*q2.n;
        d *= q2.d;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator-=(const bigint& num_val2)
{
        n -= d*num_val2;
        return *this;
}

inline bigrational& bigrational::operator*=(const bigint& num_val2)
{
        n*=num_val2;
        cancel();
        return *this;
}

inline bigrational& bigrational::operator/=(const bigint& num_val2)
{
        d*=num_val2;
        cancel();
        return *this;
}

inline bigrational::operator bigfloat() {return I2bigfloat(n)/I2bigfloat(d);}

// Definitions of non-member bigrational functions

inline bigint num(const bigrational& q)
{
        return q.n;
}

inline bigint den(const bigrational& q)
{
        return q.d;
}

inline bigrational recip(const bigrational& q)
{
        return bigrational(q.d, q.n);
}

inline bigint round(const bigrational& q)
{
        return q.n / q.d;    //provisional -- should fix rounding direction.
}

// Definitions of non-member binary operator functions

inline bigrational operator+(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d + q2.n*q1.d, q1.d * q2.d);
}

inline bigrational operator+(const bigint& num_val1, const bigrational& q2)
{
        return bigrational(num_val1*q2.d + q2.n, q2.d);
}

inline bigrational operator+(const bigrational& q1, const bigint& num_val2)
{
        return bigrational(q1.n + num_val2*q1.d, q1.d);
}

inline bigrational operator-(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d - q2.n*q1.d, q1.d * q2.d);
}

inline bigrational operator-(const bigint& num_val1, const bigrational& q2)
{
        return bigrational(num_val1*q2.d - q2.n, q2.d);
}

inline bigrational operator-(const bigrational& q1, const bigint& num_val2)
{
        return bigrational(q1.n - num_val2*q1.d, q1.d);
}

inline bigrational operator*(const bigrational& q1, const bigint& num_val2)
{
        return bigrational(q1.n*num_val2, q1.d);
}

inline bigrational operator*(const bigint& num_val1, const bigrational& q2)
{
        return bigrational(q2.n*num_val1, q2.d);
}

inline bigrational operator*(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.n, q1.d*q2.d);
}

inline bigrational operator/(const bigrational& q1, const bigint& num_val2)
{
        return bigrational(q1.n, q1.d*num_val2);
}

inline bigrational operator/(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d, q1.d*q2.n);
}

inline bigrational operator/(const bigint& num_val1, const bigrational& q2)
{
  return bigrational(q2.d*num_val1, q2.n);
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
  bigint n,d(1);
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
  r=bigrational(n,d);
  return is;
}
// NB gcd(n,d)=1 and d>0:

inline bigint floor(const bigrational& r)
{
  return (r.n-(r.n%r.d))/r.d;
}

inline bigint ceil(const bigrational& r)
{
  if(is_one(r.d)) return r.n;
  return 1 + (r.n-(r.n%r.d))/r.d;
}

inline int is_S_integral(const bigrational& x,  const vector<bigint>& S)
{
  return is_S_unit(den(x), S);
}

inline int is_S_unit(const bigrational& x,  const vector<bigint>& S)
{
  return is_S_unit(den(x), S) && is_S_unit(num(x), S);
}

inline bigrational prime_to_S_part(const bigrational& x,  const vector<bigint>& S)
{
  return bigrational(prime_to_S_part(num(x), S), prime_to_S_part(den(x), S));
}

#endif
