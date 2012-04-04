// unimod.h: declarations of classes unimod and scaled_unimod
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
 
#ifndef _UNIMOD_H
#define _UNIMOD_H      1      //flags that this file has been included

#include "interface.h"

class unimod {
  friend class quadratic;
  friend class cubic;
  friend class scaled_unimod;
protected:
  bigint m11, m12, m21, m22;
public:
  unimod()   {m11=1;m12=0;m21=0;m22=1;}
  unimod(const bigint& a11, const bigint& a12, 
	 const bigint& a21, const bigint& a22) 
    :m11(a11), m12(a12), m21(a21), m22(a22) {;}
  unimod(const unimod& m) :m11(m.m11), m12(m.m12), m21(m.m21), m22(m.m22) {;}
  void reset()
    {m11=1; m12=0; m21=0; m22=1;}
  void set(const bigint& a11, const bigint& a12, 
	   const bigint& a21, const bigint& a22) 
    {m11=a11; m12=a12; m21=a21; m22=a22;}
  void operator=(const unimod& m) {m11=m.m11; m12=m.m12; m21=m.m21; m22=m.m22;}
  bigint operator()(int i, int j) const
    {
      if(i==1) {if(j==1) return m11; else return m12;}
      else     {if(j==1) return m21; else return m22;}
    }
  void operator*=(const unimod& a)
    { 
      set(m11*a.m11 + m12*a.m21,  m11*a.m12 + m12*a.m22,
	  m21*a.m11 + m22*a.m21,  m21*a.m12 + m22*a.m22);
    }
  friend unimod operator*(const unimod& a, const unimod& b);
  bigint scale_factor() const {bigint ans; ans=1; return ans;}
  bigint det() const {return (m11*m22-m12*m21);}
  void x_shift(const bigint& a) {m12 += (a*m11); m22 += (a*m21);}
  void y_shift(const bigint& a) {m11 += (a*m12); m21 += (a*m22);}
  unimod inverse() // return inverse matrix assuming det=1
    {  
      unimod ans(m22,-m12,-m21,m11); return ans;
    }
  void invert() // multiples by [0,1; -1,0]
    {  bigint temp = -m11; m11 = m12; m12 = temp;
              temp = -m21; m21 = m22; m22 = temp;
    }
  void output(ostream& os=cout) const
    {
      os<<"["<<m11<<","<<m12<<";"<<m21<<","<<m22<<"]";
    }
  friend ostream& operator<<(ostream& os, const unimod& m);
};

class scaled_unimod : public unimod {
private: 
  bigint m00;
public:
  scaled_unimod() :unimod() {m00=1;}
  scaled_unimod(const bigint& a00, 
		const bigint& a11, const bigint& a12, 
		const bigint& a21, const bigint& a22) 
    :unimod(a11,a12,a21,a22), m00(a00) {;}
  scaled_unimod(const scaled_unimod& m) :unimod(m), m00(m.m00) {;}
  scaled_unimod(const unimod& m) :unimod(m) {m00=1;}
  void reset()  {m00=1; m11=1; m12=0; m21=0; m22=1;}
  void set(const bigint& a00, 
	   const bigint& a11, const bigint& a12, 
	   const bigint& a21, const bigint& a22) 
    {m00=a00; m11=a11; m12=a12; m21=a21; m22=a22;}
  void set(const bigint& a11, const bigint& a12, 
	   const bigint& a21, const bigint& a22) 
    {m00=1; m11=a11; m12=a12; m21=a21; m22=a22;}
  void operator=(const scaled_unimod& m) 
    {m00=m.m00; m11=m.m11; m12=m.m12; m21=m.m21; m22=m.m22;}
  void operator=(const unimod& m) 
    {m00=1; m11=m.m11; m12=m.m12; m21=m.m21; m22=m.m22;}
  bigint scale_factor() const {return m00;}
  void operator*=(const scaled_unimod& a)
    { 
      set(m00*a.m00,
	  m11*a.m11 + m12*a.m21,  m11*a.m12 + m12*a.m22,
	  m21*a.m11 + m22*a.m21,  m21*a.m12 + m22*a.m22);
    }
  void operator*=(const unimod& a)
    { 
      set(m00,
	  m11*a.m11 + m12*a.m21,  m11*a.m12 + m12*a.m22,
	  m21*a.m11 + m22*a.m21,  m21*a.m12 + m22*a.m22);
    }
  friend scaled_unimod operator*(const scaled_unimod& a, const scaled_unimod& b);
  void x_scale(const bigint& c) {m11*=c; m21*=c;}
  void y_scale(const bigint& c) {m12*=c; m22*=c;}
  void u_scale(const bigint& c) {m00*=c;}
  void output(ostream& os=cout) const
    {
      unimod::output(os);
      os<<" / "<<m00;
    }
  friend ostream& operator<<(ostream& os, const scaled_unimod& m);
};


inline ostream& operator<<(ostream& os, const unimod& m)
{
  return os<<"["<<m.m11<<","<<m.m12<<";"<<m.m21<<","<<m.m22<<"]";
}

inline ostream& operator<<(ostream& os, const scaled_unimod& m)
{
  return os<<"["<<m.m11<<","<<m.m12<<";"<<m.m21<<","<<m.m22<<"]"<<" / "<<m.m00;
}

unimod operator*(const unimod& a, const unimod& b);
scaled_unimod operator*(const scaled_unimod& a, const scaled_unimod& b);

#endif
