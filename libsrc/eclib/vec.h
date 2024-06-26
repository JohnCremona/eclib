// vec.h: declaration of integer vector classes
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
 
// Not to be included directly by user: use vector.h which defines
// _ECLIB_VECTOR_H and includes this twice
//

class svec;
class subspace;

class vec;
mat restrict_mat(const mat& m, const subspace& s, int cr=0);
mat prestrict(const mat& m, const subspace& s, const scalar& pr, int cr=0);

class vec {
  friend class svec;
  friend class smat;
  friend class smat_elim;
  friend class mat;
  friend class subspace;
public:
  // constructors
  explicit vec(long n=0);
  explicit vec(const vector<scalar>& arr);
  vec(const vec&);                       // copy constructor
  // member functions & operators
  void init(long n=0);                 // (re)-initializes
  vec& operator=(const vec&);         // assignment
  scalar& operator[](long i);            // the i'th component
  scalar operator[](long i) const;       // the i'th component
  vec& operator+=(const vec&);
  void add_row(const mat&m, int i);
  void addmodp(const vec&, const scalar&);
  vec& operator-=(const vec&);
  void sub_row(const mat&m, int i);
  vec& operator*=(const scalar&);
  vec& operator/=(const scalar&);
  vec slice(long i, long j=-1) const;  // returns subvec with indices i..j of 1..i if j=-1
  vec operator[](const vec_i&) const;   // subscript composition
  vec operator[](const vec_l&) const;   // subscript composition
  void set(long i, const scalar& x);                  // sets v[i]=x
  void add(long i, const scalar& x);                  // v[i]+=x
  void add_modp(long i, const scalar& x, const scalar& p);                  // v[i]+=x mod p
  void red_modp(const scalar& p);              // reduce mod p in place
  scalar sub(long i) const;                    // same as v[i] (no ref)
  const vector<scalar> get_entries()const {return entries;}
  static vec iota(long n);              // (1,2,...,n)
  // non-member (friend) functions and operators
  friend long dim(const vec&);                  // the dimension
  friend scalar operator*(const vec&, const vec&);   // dot product
  friend scalar operator*(const svec&, const vec&);
  friend scalar content(const vec&);
  friend scalar maxabs(const vec&);
  friend vec operator*(const mat& m, const vec& v);
  friend int operator==(const vec&, const vec&);
  friend int operator!=(const vec&, const vec&);
  friend int trivial(const vec&);                  // is v all 0
  friend int member(const scalar& a, const vec& v);//tests if a=v[i] for some i
  friend vec reverse(const vec& order);
  // add/sub row i of mat to v (implemented in mat.cc)
  friend void add_row_to_vec(vec& v, const mat& m, long i);
  friend void sub_row_to_vec(vec& v, const mat& m, long i);
  friend ostream& operator<< (ostream&s, const vec&);
  friend istream& operator>> (istream&s, vec&);
  friend void swapvec(vec& v, vec& w);
  friend mat restrict_mat(const mat& m, const subspace& s, int cr);
  friend mat prestrict(const mat& m, const subspace& s, const scalar& pr, int cr);
  friend scalar dotmodp(const vec& v1, const vec& v2, const scalar& pr);

  // Implementation
private:
  vector<scalar> entries;
};

// Declaration of non-member, non-friend functions

scalar content(const vec&);
vec operator+(const vec&);                   // unary
vec operator-(const vec&);                   // unary
vec operator+(const vec&, const vec&);
vec addmodp(const vec&, const vec&, const scalar&);
vec operator-(const vec&, const vec&);
inline vec operator*(const scalar&, const vec&);       // componentwise
vec operator/(const vec&, const scalar&);       // componentwise
void makeprimitive(vec& v);
void elim(const vec& a, vec& b, long pos);
void elim1(const vec& a, vec& b, long pos);
void elim2(const vec& a, vec& b, long pos, const scalar& lastpivot);
vec express(const vec& v, const vec& v1, const vec& v2);
int lift(const vec& v, const scalar& pr, vec& ans);  //lifts a mod-p vector to a rational
                                   //and scales to a primitive vec in Z. Returns success flag

// inline function definitions

inline long dim(const vec& v) {return v.entries.size();}

inline int operator!=(const vec& v, const vec& w) { return !(v==w);}

inline vec operator+(const vec& v) { return v;}

inline vec operator-(const vec& v) { return scalar(-1)*v;}

inline vec operator+(const vec& v1, const vec& v2)
{ vec w(v1); w+=v2; return w;}

inline vec addmodp(const vec& v1, const vec& v2, const scalar& pr)
{ vec w(v1); w.addmodp(v2,pr); return w;}

inline vec reduce_modp(const vec& v, const scalar& p)
{
  vec w(v); w.red_modp(p); return w;
}

inline vec operator-(const vec& v1, const vec& v2)
{ vec w(v1); w-=v2; return w;}

inline vec operator*(const scalar& scal, const vec& v)
{ vec w(v); w*=scal; return w;}

inline vec operator/(const vec& v, const scalar& scal)
{ vec w(v); w/=scal; return w;}

inline void makeprimitive(vec& v)
{ scalar g=content(v); if (g>1) v/=g;}

inline void elim(const vec& a, vec& b, long pos)
{ (b*=a[pos])-=(b[pos]*a);}

inline void elim1(const vec& a, vec& b, long pos)
{ (b*=a[pos])-=(b[pos]*a); makeprimitive(b);}

inline void elim2(const vec& a, vec& b, long pos, const scalar& lastpivot)
{ ((b*=a[pos])-=(b[pos]*a))/=lastpivot;}
