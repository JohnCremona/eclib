// bigcomplex.h: complex class built on NTL's RR
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

#ifndef _ECLIB_BIGCOMPLEX_H_
#define _ECLIB_BIGCOMPLEX_H_

class bigcomplex{
public:
  explicit bigcomplex() : re(RR()), im(RR()) {;}
  explicit bigcomplex(const RR& r) : re(r), im(RR()) {;}
  explicit bigcomplex(const RR& r, const RR& i) : re(r), im(i) {;}
  bigcomplex(const bigcomplex&) = default;

  // standard class methods
  RR real() const {return re;};
  RR imag() const {return im;};
  RR norm() const {return sqr(re)+sqr(im);};
  RR arg() const { return atan2(im, re);};
  RR abs() const {return ::NTL::sqrt(norm());};
  bigcomplex conj() const {return bigcomplex(re,-im);};
  bigcomplex timesI() const {return bigcomplex(-im,re);};
  int IsReal() const {return ::NTL::IsZero(im);}
  int IsImaginary() const {return ::NTL::IsZero(re);}
  int IsZero() const {return ::NTL::IsZero(re) && ::NTL::IsZero(im);}


  bigcomplex operator= (const RR& r) {re=r; im=RR(); return *this;}
  bigcomplex operator+=(const RR& r) {re+=r; return *this;};
  bigcomplex operator+(const RR& r) const {bigcomplex z(*this); z.re+=r; return z;};
  bigcomplex operator-=(const RR& r) {re -= r; return *this;};
  bigcomplex operator-(const RR& r) const {bigcomplex z(*this); z.re-=r; return z;};
  bigcomplex operator*=(const RR& r) {re*=r; im*=r; return *this;};
  bigcomplex operator*(const RR& r) const {bigcomplex z(*this); z.re*=r; z.im*=r; return z;};
  bigcomplex operator/=(const RR& r) {re/=r; im/=r; return *this;};
  bigcomplex operator/(const RR& r) const {bigcomplex z(*this); z.re/=r; z.im/=r; return z;};

  bigcomplex operator=(const bigcomplex& z) {re=z.re; im=z.im; return *this;};
  bigcomplex operator+=(const bigcomplex& z) {re+=z.re; im+=z.im; return *this;};
  bigcomplex operator+(const bigcomplex& z) const {return bigcomplex(re+z.re, im+z.im);};
  bigcomplex operator+() {return *this;};
  bigcomplex operator-=(const bigcomplex& z) {re-=z.re; im-=z.im; return *this;};
  bigcomplex operator-(const bigcomplex& z) const {return bigcomplex(re-z.re, im-z.im);};
  bigcomplex operator-() const {return bigcomplex(-re,-im);};
  bigcomplex operator*=(const bigcomplex& z) {RR x = re; re=x*z.re-im*z.im; im = x*z.im+im*z.re; return *this;};
  bigcomplex operator*(const bigcomplex& z) const {bigcomplex w(*this); w*=z; return w;};
  bigcomplex operator/=(const bigcomplex& z) {RR r=z.norm(), x = re; re=(x*z.re+im*z.im)/r; im = (-x*z.im+im*z.re)/r; return *this;};
  bigcomplex operator/(const bigcomplex& z) const {bigcomplex w(*this); w/=z; return w;};

  // standard functions as class methods:
  bigcomplex sqrt() const {RR r = abs(); RR u = ::NTL::sqrt((r+re)/2), v=::NTL::sqrt((r-re)/2); if (im<0) v=-v; return bigcomplex(u,v);};
  bigcomplex cos() const {bigcomplex z(this->timesI()); return (z+z.conj())/to_RR(2);};
  bigcomplex sin() const {bigcomplex z(this->timesI()); return (z-z.conj())/bigcomplex(RR(),to_RR(2));};
  bigcomplex exp() const {RR e = ::NTL::exp(re);  return bigcomplex(e * ::NTL::cos(im), e * ::NTL::sin(im));};
  bigcomplex log() const {return bigcomplex(::NTL::log(abs()), arg());}


  operator string() const
  {
    ostringstream s;
    s << "(" << re << "," << im << ")";
    return s.str();
  }

  string pretty_string() const
  {
    ostringstream s;
    if (IsReal())
      {
        s<<re;
        return s.str();
      }
    if (IsImaginary())
      {
        if (IsOne(im))
          s << "";
        else
          if (IsOne(-im))
            s << "-";
          else s << im;
      }
    else
      {
        s<<re;
        if(im>0)
          s<<"+";
        else
          s<<"-";
        RR abs_im(::NTL::abs(im));
        if (abs_im>1)
          s<<abs_im;
      }
    s<<"i";
    return s.str();
  }

  friend ostream& operator<<(ostream& s, const bigcomplex& z)
  {
    s << string(z);
    return s;
  }

  // Implementation
private:
  RR re, im;
};


// inline functions
inline RR real(const bigcomplex& z) {return z.real();};
inline RR imag(const bigcomplex& z) {return z.imag();};
inline RR norm(const bigcomplex& z) {return z.norm();};
inline RR arg(const bigcomplex& z) { return z.arg();};
inline RR abs(const bigcomplex& z) {return z.abs();};
inline bigcomplex conj(const bigcomplex& z) {return z.conj();};
inline int IsReal(const bigcomplex& z) {return z.IsReal();}
inline int IsImaginary(const bigcomplex& z) {return z.IsImaginary();}
inline int IsZero(const bigcomplex& z) {return z.IsZero();}

// operators with left operand RR
inline bigcomplex operator+(const RR& r, const bigcomplex& z) {return bigcomplex(r)+z;};
inline bigcomplex operator-(const RR& r, const bigcomplex& z) {return bigcomplex(r)-z;};
inline bigcomplex operator*(const RR& r, const bigcomplex& z) {return bigcomplex(r)*z;};
inline bigcomplex operator/(const RR& r, const bigcomplex& z) {return bigcomplex(r)/z;};

// standard method functions as functions:
inline bigcomplex sqrt(const bigcomplex& z) {bigcomplex w(z); return w.sqrt();};
inline bigcomplex cos(const bigcomplex& z) {bigcomplex w(z); return w.cos();};
inline bigcomplex sin(const bigcomplex& z) {bigcomplex w(z); return w.sin();};
inline bigcomplex exp(const bigcomplex& z) {bigcomplex w(z); return w.exp();};
inline bigcomplex log(const bigcomplex& z) {bigcomplex w(z); return w.log();};


inline void swap(bigcomplex& z1, bigcomplex& z2) {bigcomplex z3(z1); z1=z2; z2=z3;}

#endif
