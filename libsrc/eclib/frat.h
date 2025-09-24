// FRAT.h: interface to FLINT fmpq rational numbers

#if     !defined(_FRAT_H)
#define _FRAT_H      1       //flags that this file has been included

#include "int.h"
#include <flint/fmpq.h>

class RAT {

private:
  fmpq_t q;

public:
  // constructors
  RAT() {
    fmpq_init(q); // sets to 0
  }
  explicit RAT(const INT& n) {
    fmpq_init(q); // sets to 0
    if (n.sign())
      {
        fmpz_t d;
        fmpz_init_set_si(d, 1);
        fmpq_set_fmpz_frac(q, n.z, d);
        fmpz_clear(d);
      }
  }
  RAT(const INT& n, const INT& d) {
    fmpq_init(q); // sets to 0
    fmpq_set_fmpz_frac(q, n.z, d.z);
  }
  explicit RAT(long nu, long de=1) {
    fmpq_init(q); // sets to 0
    fmpq_set_si(q, nu, de);
  }
  RAT(const RAT& x) {
    fmpq_init(q); // sets to 0
    fmpq_set(q, x.q);
  }
  ~RAT() {
    fmpq_clear(q);
  }

  RAT& operator=(const RAT& x) {
    fmpq_init(q); // sets to 0
    fmpq_set(q, x.q);
    return *this;
  }
  RAT operator=(long n) {
    RAT x(n);
    return x;
  }

  // RAT manipulations
  void cancel() {
    fmpq_canonicalise(q); // cancel *this in situ
  }
  INT num() const       // the numerator
  {
    INT n;
    fmpz_init_set(n.z, fmpq_numref(q));
    return n;
  }
  INT den() const        // the denominator
  {
    INT n;
    fmpz_init_set(n.z, fmpq_denref(q));
    return n;
  }
  RAT recip() const  // the reciprocal
  {
    RAT x;
    fmpq_inv(x.q, q);
    return x;
  }

  int sign() const
  {
    return fmpq_sgn(q);
  }

  INT round() const;      // nearest integer

  // Binary Operator Functions
  RAT operator+(const RAT& b) const {
    RAT c;
    fmpq_add(c.q, q, b.q);
    return c;
  }
  RAT operator+(const INT& b) const {
    RAT c;
    fmpq_add_fmpz(c.q, q, b.z);
    return c;
  }
  RAT operator+(long b) const {
    RAT c;
    fmpq_add_si(c.q, q, b);
    return c;
  }
  RAT operator+(int b) const {
    RAT c;
    fmpq_add_si(c.q, q, b);
    return c;
  }
  friend inline RAT operator+(const INT& a, const RAT& b) {return b+a;}
  friend inline RAT operator+(long a, const RAT& b) {return b+a;}
  friend inline RAT operator+( int a, const RAT& b) {return b+a;}
  RAT operator-(const RAT& b) const {
    RAT c;
    fmpq_sub(c.q, q, b.q);
    return c;
  }
  RAT operator-(INT b) const {
    RAT c;
    fmpq_sub_fmpz(c.q, q, b.z);
    return c;
  }
  RAT operator-(long b) const {
    RAT c;
    fmpq_sub_si(c.q, q, b);
    return c;
  }
  RAT operator-(int b) const {
    RAT c;
    fmpq_sub_si(c.q, q, b);
    return c;
  }
  friend inline RAT operator-(const INT& a, const RAT& b) {return -(b-a);}
  friend inline RAT operator-(long a, const RAT& b) {return -(b-a);}
  friend inline RAT operator-(int a, const RAT& b) {return -(b-a);}

  RAT operator*(const RAT& b) const {
    RAT c;
    fmpq_mul(c.q, q, b.q);
    return c;
  }
  RAT operator*(const INT& b) const {
    RAT c;
    fmpq_mul_fmpz(c.q, q, b.z);
    return c;
  }
  RAT operator*(long b) const {
    RAT c;
    fmpq_mul_si(c.q, q, b);
    return c;
  }
  RAT operator*(int b) const {
    RAT c;
    fmpq_mul_si(c.q, q, b);
    return c;
  }
  friend inline RAT operator*(const INT& a, const RAT& b) {return b*a;}
  friend inline RAT operator*(const long a, const RAT& b) {return b*a;}
  friend inline RAT operator*(const int a, const RAT& b) {return b*a;}

  RAT operator/(const RAT& b) const {
    RAT c;
    fmpq_div(c.q, q, b.q);
    return c;
  }
  RAT operator/(const INT& b) const {
    RAT c;
    fmpq_div_fmpz(c.q, q, b.z);
    return c;
  }
  RAT operator/(long b) const {
    RAT c;
    fmpz_t z;
    fmpz_init_set_si(z, b);
    fmpq_div_fmpz(c.q, q, z);
    fmpz_clear(z);
    return c;
  }
  RAT operator/(int b) const {
    RAT c;
    fmpz_t z;
    fmpz_init_set_si(z, b);
    fmpq_div_fmpz(c.q, q, z);
    fmpz_clear(z);
    return c;
  }
  friend inline RAT operator/(const INT& a, const RAT& b) {
    RAT c = b/a;
    return c.recip();
  }
  int operator==(const RAT& b) const {
    return fmpq_equal(q, b.q);
  }
  int operator==(INT b) {
    return fmpq_equal_fmpz(q, b.z);
  }
  int operator==(long b) {
    return fmpq_equal_si(q, b);
  }
  friend inline int operator==(INT a, RAT b) {return b==a;}
  friend inline int operator==(long a, RAT b) {return b==a;}
  friend inline int operator!=(RAT a, RAT b) {return !(a==b);}
  friend inline int operator!=(RAT a, INT b) {return !(a==b);}
  friend inline int operator!=(INT a, RAT b) {return b!=a;}
  friend inline int operator!=(long a, const RAT& b) {return b!=INT(a);}

  int operator<(const RAT& b) const {return fmpq_cmp(q, b.q)<0;}
  int operator<=(const RAT& b) const {return fmpq_cmp(q, b.q)<=0;}
  int operator<(const INT& b) const {return fmpq_cmp_fmpz(q, b.z)<0;}
  int operator<=(const INT& b) const {return fmpq_cmp_fmpz(q, b.z)<=0;}
  int operator<(long b) const {return fmpq_cmp_si(q, b)<0;}
  int operator<=(long b) const {return fmpq_cmp_si(q, b)<=0;}
  int operator>(const RAT& b) const {return fmpq_cmp(q, b.q)>0;}
  int operator>=(const RAT& b) const {return fmpq_cmp(q, b.q)>=0;}
  int operator>(const INT& b) const {return fmpq_cmp_fmpz(q, b.z)>0;}
  int operator>=(const INT& b) const {return fmpq_cmp_fmpz(q, b.z)>=0;}
  int operator>(long b) const {return fmpq_cmp_si(q, b)>0;}
  int operator>=(long b) const {return fmpq_cmp_si(q, b)>=0;}
  friend inline int operator<(const INT& a, const RAT& b) {return b>a;}
  friend inline int operator<=(const INT& a, const RAT& b) {return b>=a;}
  friend inline int operator<(long a, const RAT& b) {return b>a;}
  friend inline int operator<=(long a, const RAT& b) {return b>=a;}
  friend inline int operator>(const INT& a, const RAT& b) {return b<a;}
  friend inline int operator>=(const INT& a, const RAT& b) {return b<=a;}
  friend inline int operator>(long a, const RAT& b) {return b<a;}
  friend inline int operator>=(long a, const RAT& b) {return b<=a;}

  friend inline RAT max(const RAT& a, const RAT& b);
  friend inline RAT min(const RAT& a, const RAT& b);

  friend ostream& operator<< (ostream&, const RAT&);
  friend istream& operator>> (istream&, RAT&);

  void operator+=(const RAT& b) {fmpq_add(q, q, b.q);}
  void operator+=(const INT& b) {fmpq_add_fmpz(q, q, b.z);}
  void operator+=(int b) {fmpq_add_si(q, q, b);}
  void operator-=(const RAT& b) {fmpq_sub(q, q, b.q);}
  void operator-=(const INT& b) {fmpq_sub_fmpz(q, q, b.z);}
  void operator-=(int b) {fmpq_sub_si(q, q, b);}
  void operator*=(const RAT& b) {fmpq_mul(q, q, b.q);}
  void operator*=(const INT& b) {fmpq_mul_fmpz(q, q, b.z);}
  void operator*=(int b) {fmpq_mul_si(q, q, b);}
  void operator/=(const RAT& b) {fmpq_div(q, q, b.q);}
  void operator/=(const INT& b) {fmpq_div_fmpz(q, q, b.z);}
  void operator/=(int b) {
    fmpz_t bb;
    fmpz_init_set_si(bb,b);
    fmpq_div_fmpz(q, q, bb);
    fmpz_clear(bb);
  }

  RAT operator+() const {
    RAT x(*this);
    return x;
  }
  RAT operator-() const {
    RAT x(*this);
    fmpq_neg(x.q, x.q);
    return x;
  }
  RAT abs() const {
    RAT x(*this);
    fmpq_abs(x.q, x.q);
    return x;
  }
  INT floor() const;
  INT ceil() const;
  int is_square() const {
    return (num()*den()).is_square();
  }
  int is_square(RAT& a) const {
    INT nd = num()*den(), rnd;
    int b=nd.is_square(rnd);
    if (b) a = RAT(rnd, den());
    return b;
  }
};


// Definitions of non-member RAT functions


inline INT RAT::round() const {
  return rounded_division(num(), den(), 0); // 0 means halves go up
}

inline ostream& operator<<(ostream& s, const RAT& q)
{
  INT n(q.num()), d(q.den());
  if (is_zero(d))
    s<<"oo";
  else
    {
      s << n;
      if (d!=1)
        s << "/" << d;
    }
  return s;
}

inline istream& operator>> (istream& is, RAT& r)
{
  std::string n;
  is>>n;
  fmpq_set_str(r.q, n.c_str(), 10);
  fmpq_canonicalise(r.q);
  return is;
}

inline INT RAT::floor() const
{
  INT q, r;
  divrem(num(), den(), q, r);
  return q;
}

inline INT RAT::ceil() const
{
  INT q, r;
  divrem(num(), den(), q, r);
  return (is_zero(r)? q : q+1);
}

inline int sign(const RAT& r)
{
  return r.sign();
}

inline RAT max(const RAT& a, const RAT& b) {return (a>=b? a : b);}
inline RAT min(const RAT& a, const RAT& b) {return (a<=b? a : b);}


#endif
