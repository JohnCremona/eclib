// bigrat.h: declaration of bigrational number class
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
 
#if     !defined(_BIGRATIONAL_H)
#define _BIGRATIONAL_H      1       //flags that this file has been included

#include "marith.h"
#include "rat.h"

class bigrational {

public:
        // constructors
        bigrational(bigint num_val=BIGINT(0), bigint den_val=BIGINT(1));
        bigrational(const bigrational& q);
        bigrational(const rational& q);
        void operator=(const bigrational& q);
        void operator=(const rational& q);

        // bigrational manipulations
        void cancel();                           // cancel *this in situ
        friend bigint num(const bigrational&);        // the numerator
        friend bigint den(const bigrational&);        // the denominator
        friend bigrational recip(const bigrational&);  // reciprocal
        friend bigint round(const bigrational&);      // nearest integer

        // Binary Operator Functions
        friend bigrational operator+(const bigrational&, const bigrational&);
        friend bigrational operator+(bigint, const bigrational&);
        friend bigrational operator+(const bigrational&, bigint);
        friend bigrational operator-(const bigrational&, const bigrational&);
        friend bigrational operator-(bigint, const bigrational&);
        friend bigrational operator-(const bigrational&, bigint);
        friend bigrational operator*(const bigrational&, const bigrational&);
        friend bigrational operator*(const bigrational&, bigint);
        friend bigrational operator*(bigint, const bigrational&);
        friend bigrational operator/(const bigrational&, const bigrational&);
        friend bigrational operator/(const bigrational&, bigint);
        friend bigrational operator/(bigint, const bigrational&);
        friend int operator==(const bigrational&, const bigrational&);
        friend int operator!=(const bigrational&, const bigrational&);
        friend ostream& operator<< (ostream&s, const bigrational&);
        friend istream& operator>> (istream& is, bigrational& r);
        bigrational& operator+=(const bigrational&);
        bigrational& operator+=(bigint);
        bigrational& operator-=(const bigrational&);
        bigrational& operator-=(bigint);
        bigrational& operator*=(const bigrational&);
        bigrational& operator*=(bigint);
        bigrational& operator/=(const bigrational&);
        bigrational& operator/=(bigint);
        bigrational operator+();
        bigrational operator-();
        operator bigfloat();  // conversion operator

// Implementation
private:
        bigint n, d;
};


// Inline bigrational functions

inline void bigrational::cancel()                     // cancel *this in situ
{
 bigint g = gcd(n,d);
 if (g>1) {n/=g; d/=g;}
 if (d<0) {n=-n; d=-d;}
}

inline bigrational::bigrational(bigint num_val, bigint den_val)
{
  n=num_val; d=den_val;
  (*this).cancel(); 
}

inline bigrational::bigrational(const bigrational& q) :n(q.n), d(q.d) {;}
inline bigrational::bigrational(const rational& q) :n(BIGINT(q.n)), d(BIGINT(q.d)) {;}
inline void bigrational::operator=(const bigrational& q) {n=q.n; d=q.d;}
inline void bigrational::operator=(const rational& q) {n=BIGINT(q.n); d=BIGINT(q.d);}

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
        (*this).cancel();
        return *this;
}

inline bigrational& bigrational::operator+=(bigint num_val2)
{
        n += d*num_val2;
        return *this;
}

inline bigrational& bigrational::operator-=(const bigrational& q2)
{
        n = n*q2.d-d*q2.n;
        d *= q2.d;
        (*this).cancel();
        return *this;
}

inline bigrational& bigrational::operator-=(bigint num_val2)
{
        n -= d*num_val2;
        return *this;
}

inline bigrational& bigrational::operator*=(bigint num_val2)
{
        n*=num_val2;
        (*this).cancel();
        return *this;
}

inline bigrational& bigrational::operator/=(bigint num_val2)
{
        d*=num_val2;
        (*this).cancel();
        return *this;
}

//inline bigrational::operator bigfloat() {return double(n)/double(d);}

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

inline bigrational operator+(bigint num_val1, const bigrational& q2)
{
        return bigrational(num_val1*q2.d + q2.n, q2.d);
}

inline bigrational operator+(const bigrational& q1, bigint num_val2)
{
        return bigrational(q1.n + num_val2*q1.d, q1.d);
}

inline bigrational operator-(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d - q2.n*q1.d, q1.d * q2.d);
}

inline bigrational operator-(bigint num_val1, const bigrational& q2)
{
        return bigrational(num_val1*q2.d - q2.n, q2.d);
}

inline bigrational operator-(const bigrational& q1, bigint num_val2)
{
        return bigrational(q1.n - num_val2*q1.d, q1.d);
}

inline bigrational operator*(const bigrational& q1, bigint num_val2)
{
        return bigrational(q1.n*num_val2, q1.d);
}

inline bigrational operator*(bigint num_val1, const bigrational& q2)
{
        return bigrational(q2.n*num_val1, q2.d);
}

inline bigrational operator*(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.n, q1.d*q2.d);
}

inline bigrational operator/(const bigrational& q1, bigint num_val2)
{
        return bigrational(q1.n, q1.d*num_val2);
}

inline bigrational operator/(const bigrational& q1, const bigrational& q2)
{
        return bigrational(q1.n*q2.d, q1.d*q2.n);
}

inline bigrational operator/(bigint num_val1, const bigrational& q2)
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
  char c;
  bigint n,d=BIGINT(1);
  is>>n>>ws;
  if(!is.eof()) 
    {
      is>>c;
      if(c=='/')
	{
	  is>>d;
	}
      else
	{
	  if(c!='\n')
	    {
	      is.putback(c);
	    }
	}
    }
  r=bigrational(n,d);
  return is;
}

#endif
