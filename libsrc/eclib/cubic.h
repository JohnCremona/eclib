// cubic.h: integer cubic class for unimodular transforms and reduction.
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
 //
// Notation: g(x,z) is replaced by g(m11*x+m12*z,m21*x+m22*z)
//
// Stored as bigint* arrays g of size 4 representing 
//         g[0]*x^3+g[1]*x^2+g[2]*x+g[3]

class unimod;

class cubic {
  friend class unimod;
private: 
  bigint * coeffs;  // will always have length 4
  // init just allocates memory
  void init();
public:
  void set(long a, long b, long c, long d) 
    {coeffs[0]=a; coeffs[1]=b; coeffs[2]=c; coeffs[3]=d;}
  void set(const  bigint& a, const bigint& b, const bigint& c, const bigint& d) 
    {coeffs[0]=a; coeffs[1]=b; coeffs[2]=c; coeffs[3]=d;}
  void set(const  cubic& q)
    {coeffs[0]=q.coeffs[0]; coeffs[1]=q.coeffs[1]; 
     coeffs[2]=q.coeffs[2]; coeffs[3]=q.coeffs[3];}
  cubic() 
    {init(); set(0,0,0,0);}
  ~cubic();
  cubic(const  bigint& a, const bigint& b, const bigint& c, const bigint& d) 
    {init(); set(a,b,c,d);}
  cubic(long a, long b, long c, long d) 
    {init(); set(a,b,c,d);}
  cubic(const  cubic& q)  
    {init(); set(q);}
  void operator=(const cubic& g) {set(g);}
  int operator==(const cubic& g)
  {return ((coeffs[0]==g.coeffs[0]) && (coeffs[1]==g.coeffs[1]) &&
           (coeffs[2]==g.coeffs[2]) && (coeffs[3]==g.coeffs[3]));}
  inline bigint coeff(int i) 
    {if((i>=0)&&(i<=3)) return coeffs[i]; else return coeffs[0];}
  inline bigint operator[](int i) const
    {if((i>=0)&&(i<=3)) return coeffs[i]; else return coeffs[0];}
  inline bigint a(void) const {return coeffs[0];}
  inline bigint b(void) const {return coeffs[1];}
  inline bigint c(void) const {return coeffs[2];}
  inline bigint d(void) const {return coeffs[3];}
  inline void set_coeff(int i, const bigint& a)
    {if((i>=0)&&(i<=3)) coeffs[i]=a;}
  inline bigint eval(const bigint& x, const bigint& z) const
    { bigint x2=sqr(x), z2=sqr(z);
      return a()*x*x2 + b()*x2*z + c()*x*z2 + d()*z*z2;}
  inline bigint eval(const bigint& x) const
    { bigint x2=sqr(x);
      return a()*x*x2 + b()*x2 + c()*x + d();}
  inline bigint disc() const
    { bigint b2=sqr(b()), c2=sqr(c()), ac=a()*c(), bd=b()*d();
      return -27*sqr(a()*d()) + 18*ac*bd - 4*ac*c2 -4*bd*b2 + b2*c2;
    }
  inline void output(ostream& os=cout) const
    {
      os<<"["<<a()<<","<<b()<<","<<c()<<","<<d()<<"]";
    }
  friend inline ostream& operator<<(ostream& os, const cubic& g);
  void transform(const unimod& m);
  // In the next 4 functions, m already holds a unimod and is updated:
  void x_shift(const bigint& e, unimod& m);
  void y_shift(const bigint& e, unimod& m);
  void invert(unimod& m); // apply [0,-1;1,0]
  void negate(unimod& m); // apply [-1,0;0,-1]
  void seminegate(unimod& m); // apply [1,0;0,-1] (det=-1)
  void reduce(unimod& m);

// Mathews quantities for use when disc<0:
  bigint mat_c1() const
    { return d()*(d()-b())+a()*(c()-a());}
  bigint mat_c2() const
    {  return a()*d() - (a()+b())*(a()+b()+c());}
  bigint mat_c3() const
    {  return a()*d() + (a()-b())*(a()-b()+c());}
  bigint p_semi() const
    {  return sqr(b())-3*a()*c(); }
  bigint q_semi() const
    {  return b()*c()-9*a()*d(); }
  bigint r_semi() const
    {  return sqr(c())-3*b()*d(); }
  bigint u_semi() const
    {  return 2*b()*sqr(b()) + 27*sqr(a())*d() - 9*a()*b()*c();}
  bigint j_c1() const;
  bigint j_c2() const;
  bigint j_c3() const;
  bigint j_c4() const;

  bigcomplex hess_root() const;
  bigfloat real_root() const;  // requires disc<0
  int is_hessian_reduced(); // for positive discriminant only
  void hess_reduce(unimod& m);
  void mathews_reduce(unimod& m);
  int is_jc_reduced(); // for negative discriminant only
  void jc_reduce(unimod& m);
  // Just shifts x, returns the shift amount:
  bigint shift_reduce();
  vector<bigrational> rational_roots() const;
  int is_reducible() const
  {return ((a()==0) || (rational_roots().size()>0));}
  int is_irreducible() const
  {return ((a()!=0) && (rational_roots().size()==0));}
};


inline ostream& operator<<(ostream& os, const cubic& g)
{
  return os<<"["<<g.a()<<","<<g.b()<<","<<g.c() <<","<<g.d()<<"]";
}

// Functions for listing all reduced cubics
// verbose=1 shows original cubics found before final reduction and elimination of duplicates
// verbose=2 also shows details of triple loop

// All reduced cubics with a single discriminant (positive or negative):
// Set include_reducibles=0 to omit reducible cubics and any with a=0
// Set gl2=1 to get GL(2,Z)-inequivalent cubics (default is SL(2,Z))
vector<cubic> reduced_cubics(const bigint& disc, int include_reducibles=1, int gl2=0, int verbose=0);

// All reduced cubics with discriminant in range (0,maxdisc] if maxdisc>0 or [maxdisc,0) if maxdisc<0
// (not yet implemented)
vector<cubic> reduced_cubics_range(const bigint& maxdisc, int verbose=0);
