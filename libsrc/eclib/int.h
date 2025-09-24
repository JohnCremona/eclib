// FILE INT.H: wrapper for FLINT's fmpz type

#if     !defined(_INT_H)
#define _INT_H      1       //flags that this file has been included

#include <iostream>
#include <vector>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>

class INT {
private:
  fmpz_t z;
public:
  // Constructors:
  INT() {fmpz_init(z);} // sets to 0
  explicit INT(int a) {fmpz_init_set_si(z, a);}
  explicit INT(long a) {fmpz_init_set_si(z, a); }
  INT(const INT& a) {fmpz_init_set(z, a.z);}
  explicit INT(fmpz_t a) {fmpz_init_set(z, a);}
  // Destructor:
  ~INT() {fmpz_clear(z);}
  INT& operator=(int a) {fmpz_set_si(z, a); return *this;}
  INT& operator=(long a) {fmpz_set_si(z, a); return *this;}
  INT& operator=(const INT& a) {fmpz_set(z, a.z); return *this;}
  int operator==(const INT& a) const {return fmpz_equal(z, a.z);}
  int operator==(int a) const {return fmpz_equal_si(z, a);}
  int operator==(long a) const {return fmpz_equal_si(z, a);}
  int operator!=(const INT& a) const {return !fmpz_equal(z, a.z);}
  int operator!=(int a) const {return !fmpz_equal_si(z, a);}
  int operator!=(long a) const {return !fmpz_equal_si(z, a);}
  int sign() const {return fmpz_sgn(z);}
  INT operator-() const {INT b; fmpz_neg(b.z,z); return b;}
  INT abs() const {INT b; fmpz_abs(b.z,z); return b;}
  int is_long() const {return fmpz_fits_si(z);}
  INT operator+(int a) const {INT b; fmpz_add_si(b.z,z,a); return b;}
  INT operator+(long a) const {INT b; fmpz_add_si(b.z,z,a); return b;}
  INT operator+(const INT& a) const {INT b; fmpz_add(b.z,z,a.z); return b;}
  INT operator-(int a) const {INT b; fmpz_sub_si(b.z,z,a); return b;}
  INT operator-(long a) const {INT b; fmpz_sub_si(b.z,z,a); return b;}
  INT operator-(const INT& a) const {INT b; fmpz_sub(b.z,z,a.z); return b;}
  INT operator*(int a) const {INT b; fmpz_mul_si(b.z,z,a); return b;}
  INT operator*(long a) const {INT b; fmpz_mul_si(b.z,z,a); return b;}
  INT operator*(const INT& a) const {INT b; fmpz_mul(b.z,z,a.z); return b;}
  INT operator/(const INT& a) const {INT b; fmpz_divexact(b.z,z,a.z); return b;}
  INT operator/(int a) const {INT b; fmpz_divexact_si(b.z,z,a); return b;}
  INT operator/(long a) const {INT b; fmpz_divexact_si(b.z,z,a); return b;}
  INT operator%(const INT& a) const {INT b; fmpz_mod(b.z,z,a.z); return b;}
  int operator%(int a) const {INT b; fmpz_mod_ui(b.z,z,a); return (int)I2long(b);}  // a must be >0
  long operator%(long a) const {INT b; fmpz_mod_ui(b.z,z,a); return I2long(b);} // a must be >0
  INT operator^(int e) const {INT b; fmpz_pow_ui(b.z,z,e); return b;}
  INT operator^(long e) const {INT b; fmpz_pow_ui(b.z,z,e); return b;}
  INT operator<<(long e) const {INT b; fmpz_mul_2exp(b.z, z, e); return b;} // mult by 2^e
  void operator +=(const INT& a) {fmpz_add(z,z,a.z);}
  void operator +=(int a) {fmpz_add_si(z,z,a);}
  void operator +=(long a) {fmpz_add_si(z,z,a);}
  void operator -=(const INT& a) {fmpz_sub(z,z,a.z);}
  void operator -=(int a) {fmpz_sub_si(z,z,a);}
  void operator -=(long a) {fmpz_sub_si(z,z,a);}
  void operator *=(const INT& a) {fmpz_mul(z,z,a.z);}
  void operator *=(int a) {fmpz_mul_si(z,z,a);}
  void operator *=(long a) {fmpz_mul_si(z,z,a);}
  void operator /=(const INT& a) {fmpz_divexact(z,z,a.z);}
  void operator /=(int a) {fmpz_divexact_si(z,z,a);}
  void operator /=(long a) {fmpz_divexact_si(z,z,a);}
  void operator<<=(long e) {fmpz_mul_2exp(z, z, e);} // mult by 2^e in place
  int is_zero() const {return fmpz_is_zero(z);}
  int is_nonzero() const {return !fmpz_is_zero(z);}
  int is_one() const {return fmpz_is_one(z);}
  int is_even() const {return fmpz_is_even(z);}
  int is_odd() const {return fmpz_is_odd(z);}
  int is_square() const {return fmpz_is_square(z);}
  int is_square(INT& a) const {if (!fmpz_is_square(z)) return 0; else {fmpz_sqrt(a.z,z); return 1;}}
  INT isqrt() const {INT b; if (fmpz_sgn(z)>0) fmpz_sqrt(b.z,z); return b;}
  long valuation(const INT& p) {INT q;   return fmpz_remove(q.z, z, p.z);}
  friend std::ostream& operator<<(std::ostream& s, const INT& a);
  friend std::istream& operator>>(std::istream& s, INT& x);
  friend INT fmma(const INT& a, const INT& b, const INT& c, const INT& d); // a*b+c*d
  friend INT fmms(const INT& a, const INT& b, const INT& c, const INT& d); // a*b-c*d
  friend INT mod(const INT& a, const INT& b); // a mod b in range (-b/2,b/2]
  friend INT gcd(const INT& a, const INT& b);
  friend INT bezout(const INT& a, const INT& b, INT& x, INT& y);
  friend int compare(const INT& a, const INT& b);
  friend int compare(const INT& a, int b);
  friend int compare(const INT& a, long b);
  friend int divrem(const INT& a, const INT& b, INT& quo, INT& rem);
  friend std::vector<INT> pdivs(const INT& a);
  friend std::vector<INT> sqdivs(const INT& a);
  friend int legendre(const INT& a, const INT& p) {return kronecker(a,p);}
  friend int kronecker(const INT& a, const INT& n);
  friend INT invmod(const INT&a, const INT& p);
  friend long invmod(const INT&a, long p);
  friend long I2long(const INT& a);
// return e and set q, where a=q*f^e, e maximal
  friend long divide_out(const INT&a, const INT& f, INT& q);
  friend void sqrt_mod_p(INT& b, const INT& a, const INT& p);
  friend void swap(INT& a, INT& b);

  friend class RAT;
  friend class REAL;
  friend void make_mat( fmpz_mat_t A, const std::vector<std::vector<INT>>& M);
};

inline int sign(const INT& a) {return a.sign();}
inline INT operator+(int a, const INT& b) {return b+a;}
inline INT operator+(long a, const INT& b) {return b+a;}
inline INT operator-(int a, const INT& b) {return -(b-a);}
inline INT operator-(long a, const INT& b) {return -(b-a);}
inline INT operator*(int a, const INT& b) {return b*a;}
inline INT operator*(long a, const INT& b) {return b*a;}
inline INT fmma(const INT& a, const INT& b, const INT& c, const INT& d) {INT e; fmpz_fmma(e.z, a.z, b.z, c.z, d.z); return e;}
inline INT fmms(const INT& a, const INT& b, const INT& c, const INT& d) {INT e; fmpz_fmms(e.z, a.z, b.z, c.z, d.z); return e;}
inline INT abs(const INT& a) {return a.abs();}
inline INT posmod(const INT& a, const INT& b) {return a%b;} // in range [0,|b|)
inline INT posmod(const INT& a, long b) {return INT(a%b);} // in range [0,|b|)
inline INT mod(const INT& a, const INT& b) {INT c; fmpz_smod(c.z,a.z,b.z); return c;} // in range (-|b|/2,|b|/2-1]
inline INT gcd(const INT& a, const INT& b) {INT c; fmpz_gcd(c.z,a.z,b.z); return c;}
inline INT bezout(const INT& a, const INT& b, INT& x, INT& y) {INT c; fmpz_xgcd_canonical_bezout(c.z, x.z, y.z, a.z, b.z); return c;}
inline int compare(const INT& a, const INT& b) {return fmpz_cmp(a.z,b.z);}
inline int compare(const INT& a, int b) {return fmpz_cmp_si(a.z,b);}
inline int compare(const INT& a, long b) {return fmpz_cmp_si(a.z,b);}
inline int operator==(int a, const INT& b) {return b==a;}
inline int operator==(long a, const INT& b) {return b==a;}
inline int operator<(const INT& a, const INT& b) {return compare(a,b)<0;}
inline int operator<(const INT& a, long b) {return compare(a,b)<0;}
inline int operator>(const INT& a, const INT& b) {return compare(a,b)>0;}
inline int operator>(const INT& a, long b) {return compare(a,b)>0;}
inline int operator<(long b, const INT& a) {return compare(a,b)>0;}
inline int operator<=(long b, const INT& a) {return compare(a,b)>=0;}
inline int operator>(long b, const INT& a) {return compare(a,b)<0;}
inline int operator>=(long b, const INT& a) {return compare(a,b)<=0;}
inline int operator<=(const INT& a, const INT& b) {return compare(a,b)<=0;}
inline int operator<=(const INT& a, long b) {return compare(a,b)<=0;}
inline int operator>=(const INT& a, const INT& b) {return compare(a,b)>=0;}
inline int operator>=(const INT& a, long b) {return compare(a,b)>=0;}
inline int divides(const INT& a, const INT& b) {return (b%a).is_zero();}
inline int divides(int a, const INT& b) {return (b%a) ==0;}
inline int divides(long a, const INT& b) {return (b%a) ==0;}
inline int kronecker(const INT& a, const INT& p) {return fmpz_kronecker(a.z, p.z);}
inline int legendre(const INT& a, int p) {return legendre(a,INT(p));}
inline int legendre(const INT& a, long p) {return legendre(a,INT(p));}
// return e and set q, where a=q*f^e, e maximal
inline long divide_out(const INT&a, const INT& f, INT& q) {return fmpz_remove(q.z, a.z, f.z);}
// return valuation of a at f
inline long val(const INT&f, const INT& a) {INT q; return divide_out(a, f, q);}
inline void sqrt_mod_p(INT& b, const INT& a, const INT& p) {fmpz_sqrtmod(b.z, a.z, p.z);}
inline void swap(INT& a, INT& b) {fmpz_swap(a.z, b.z);}
inline INT pow(const INT& a, int e) {return a^e;}
inline INT pow(const INT& a, long e) {return a^e;}
inline int isqrt(const INT& a, INT& root) {root=a.isqrt(); return root*root==a;}
inline INT isqrt(const INT& a) {return a.isqrt();}
inline INT max(const INT& a, const INT& b) {return (a>=b? a : b);}
inline INT min(const INT& a, const INT& b) {return (a<=b? a : b);}
inline int is_zero(const INT& a) {return a.is_zero();}
inline int is_nonzero(const INT& a) {return a.is_nonzero();}
inline int is_one(const INT& a) {return a.is_one();}

// functions implemented in int.cc:

// Return q=round(a/b) with halves going down if round_down=1 (default) else up
INT rounded_division(const INT& a, const INT& b, int round_down=1);

std::vector<INT> pdivs(const INT& a);
std::vector<INT> sqdivs(const INT& a);
INT sqrt_mod_p(const INT& a, const INT& p);

#endif
