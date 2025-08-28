// vector.h: manage declarations of integer vector classes
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
 
#if     !defined(_ECLIB_VECTOR_H)
#define _ECLIB_VECTOR_H      1       //flags that this file has been included

#include "marith.h"

template<class T> class Zvec;
template<class T> class sZvec;
template<class T> class sZmat;
template<class T> class sZmat_elim;
template<class T> class Zmat;
template<class T> class subZspace;
template<class T> class ssubZspace;
template<class T> class form_finderT;

typedef Zvec<int> vec_i;
typedef Zmat<int> mat_i;
typedef subZspace<int> subspace_i;
typedef sZvec<int> svec_i;
typedef sZmat<int> smat_i;
typedef sZmat_elim<int> smat_i_elim;
typedef ssubZspace<int> ssubspace_i;
typedef form_finderT<int> form_finder_i;

typedef Zvec<long> vec_l;
typedef Zmat<long> mat_l;
typedef subZspace<long> subspace_l;
typedef sZvec<long> svec_l;
typedef sZmat<long> smat_l;
typedef sZmat_elim<long> smat_l_elim;
typedef ssubZspace<long> ssubspace_l;
typedef form_finderT<long> form_finder_l;

typedef Zvec<bigint> vec_m;
typedef Zmat<bigint> mat_m;
typedef subZspace<bigint> subspace_m;
typedef sZvec<bigint> svec_m;
typedef sZmat<bigint> smat_m;
typedef sZmat_elim<bigint> smat_m_elim;
typedef ssubZspace<bigint> ssubspace_m;
typedef form_finderT<bigint> form_finder_m;

template<class T> T default_modulus();

template<class T> int dim(const Zvec<T>&);
template<class T> T operator*(const Zvec<T>&, const Zvec<T>&);
template<class T> T operator*(const sZvec<T>&, const Zvec<T>&);
template<class T> T content(const Zvec<T>&);
template<class T> T maxabs(const Zvec<T>&);
template<class T> Zvec<T> operator*(const Zmat<T>& m, const Zvec<T>& v);
template<class T> int operator==(const Zvec<T>&, const Zvec<T>&);
template<class T> int operator!=(const Zvec<T>&, const Zvec<T>&);
template<class T> int trivial(const Zvec<T>&);
template<class T> int member(const T& a, const Zvec<T>& v);//tests if a=v[i] for some i
template<class T> Zvec<T> reverse(const Zvec<T>& order);
template<class T> void add_row_to_vec(Zvec<T>& v, const Zmat<T>& m, long i);
template<class T> void sub_row_to_vec(Zvec<T>& v, const Zmat<T>& m, long i);
template<class T> ostream& operator<< (ostream&s, const Zvec<T>&);
template<class T> istream& operator>> (istream&s, Zvec<T>&);
template<class T> void swapvec(Zvec<T>& v, Zvec<T>& w);
template<class T> Zmat<T> operator*(const Zmat<T>&, const Zmat<T>&);
template<class T> Zvec<T> operator*(const Zmat<T>&, const Zvec<T>&);
template<class T> int operator==(const Zmat<T>&, const Zmat<T>&);
template<class T> istream& operator>> (istream&s, Zmat<T>&);
template<class T> Zmat<T> colcat(const Zmat<T>& a, const Zmat<T>& b);
template<class T> Zmat<T> rowcat(const Zmat<T>& a, const Zmat<T>& b);
template<class T> Zmat<T> directsum(const Zmat<T>& a, const Zmat<T>& b);
template<class T> void elimrows(Zmat<T>& m, long r1, long r2, long pos);
template<class T> void elimrows1(Zmat<T>& m, long r1, long r2, long pos);
template<class T> void elimrows2(Zmat<T>& m, long r1, long r2, long pos, const T& last);
template<class T> Zmat<T> echelon0(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                      long& rk, long& ny, T& d);
template<class T> void elimp(Zmat<T>& m, long r1, long r2, long pos, const T& pr);
template<class T> void elimp1(Zmat<T>& m, long r1, long r2, long pos, const T& pr);
template<class T> Zmat<T> echelonp(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                      long& rk, long& ny, T& d, const T& pr);
template<class T> Zmat<T> echmodp(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                     long& rk, long& ny, const T& pr);
template<class T> Zmat<T> echmodp_uptri(const Zmat<T>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                     long& rk, long& ny, const T& pr);
template<class T> Zmat<T> ref_via_flint(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                           long& rk, long& ny, const T& pr);
template<class T> Zmat<T> ref_via_ntl(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const T& pr);
template<class T> Zmat<T> rref(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const T& pr);
template<class T> long rank_via_ntl(const Zmat<T>& M, const T& pr);
template<class T> T det_via_ntl(const Zmat<T>& M, const T& pr);
template<class T> subZspace<T> combine(const subZspace<T>& s1, const subZspace<T>& s2);
template<class T> ssubZspace<T> combine(const ssubZspace<T>& s1, const ssubZspace<T>& s2);
template<class T> sZmat<T> restrict_mat(const sZmat<T>& m, const ssubZspace<T>& s);
template<class T> Zmat<T> restrict_mat(const Zmat<T>& m, const subZspace<T>& s, int cr=0);
template<class T> sZmat<T> restrict_mat(const sZmat<T>& m, const subZspace<T>& s);
template<class T> int liftmat(const Zmat<T>& mm, const T& pr, Zmat<T>& m, T& dd);
template<class T> int lift(const subZspace<T>& s, const T& pr, subZspace<T>& ans);
template<class T> subZspace<T> pcombine(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr);
template<class T> Zmat<T> matmulmodp(const Zmat<T>&, const Zmat<T>&, const T& pr);
template<class T> long population(const Zmat<T>& m);
template<class T> T maxabs(const Zmat<T>& m);
template<class T> double sparsity(const Zmat<T>& m);
template<class T> Zmat<T> prestrict(const Zmat<T>& m, const subZspace<T>& s, const T& pr, int cr=0);
template<class T> T dotmodp(const Zvec<T>& v1, const Zvec<T>& v2, const T& pr);
template<class T> int dim(const subZspace<T>& s);
template<class T> T denom(const subZspace<T>& s);
template<class T> Zvec<int> pivots(const subZspace<T>& s);
template<class T> Zmat<T> basis(const subZspace<T>& s);
template<class T> subZspace<T> combine(const subZspace<T>& s1, const subZspace<T>& s2);
template<class T> subZspace<T> pcombine(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr);
template<class T> int lift(const subZspace<T>& s, const T& pr, subZspace<T>& ans);
template<class T> int eqmodp(const sZvec<T>&, const sZvec<T>&, const T& p);
template<class T> ostream& operator<< (ostream&s, const sZvec<T>&);
template<class T> T operator*(const sZvec<T>&, const sZvec<T>&);
template<class T> T operator*(const sZvec<T>&, const Zvec<T>&);
template<class T> T operator*(const Zvec<T>& v, const sZvec<T>& sv);
template<class T> T dotmodp(const sZvec<T>&, const sZvec<T>&, const T& pr);
template<class T> T dotmodp(const sZvec<T>&, const Zvec<T>&, const T& pr);
template<class T> T dotmodp(const Zvec<T>& v, const sZvec<T>& sv, const T& pr);
template<class T> sZvec<T> operator+(const sZvec<T>& v1, const sZvec<T>& v2);
template<class T> sZvec<T> operator-(const sZvec<T>& v1, const sZvec<T>& v2);
template<class T> int operator==(const sZvec<T>& v1, const sZvec<T>& v2);
template<class T> int operator!=(const sZvec<T>& v1, const sZvec<T>& v2);
template<class T> int operator==(const sZvec<T>& v1, const Zvec<T>& v2);
template<class T> sZmat<T> transpose(const sZmat<T>&);
template<class T> sZmat<T> operator* ( const sZmat<T>&, const sZmat<T>&);
template<class T> T content(const sZvec<T>& v);
template<class T> T make_primitive(sZvec<T>& v);
template<class T> sZvec<T> operator* ( const sZmat<T>& A, const sZvec<T>& v );
template<class T> sZvec<T> operator* ( const sZvec<T>& v, const sZmat<T>& A );
template<class T> sZvec<T> mult_mod_p( const sZmat<T>& A, const sZvec<T>& v, const T& p  );
template<class T> sZvec<T> mult_mod_p( const sZvec<T>& v, const sZmat<T>& A, const T& p  );
template<class T> sZmat<T> mult_mod_p ( const sZmat<T>&, const sZmat<T>&, const T&);
template<class T> inline vector<int> dim(const sZmat<T>& A)
{return vector<int>{A.nro, A.nco};}
template<class T> Zvec<T> operator*  (const sZmat<T>& m, const Zvec<T>& v);
template<class T> sZvec<T> operator* ( const sZmat<T>& A, const sZvec<T>& v );
template<class T> sZvec<T> operator* ( const sZvec<T>& v, const sZmat<T>& A );
template<class T> sZvec<T> mult_mod_p( const sZmat<T>& A, const sZvec<T>& v, const T& p  );
template<class T> Zvec<T> mult_mod_p( const sZmat<T>& A, const Zvec<T>& v, const T& p  );
template<class T> sZvec<T> mult_mod_p( const sZvec<T>& v, const sZmat<T>& A, const T& p  );
template<class T> sZmat<T> operator* ( const sZmat<T>& A, const sZmat<T>& B );
template<class T> sZmat<T> mult_mod_p ( const sZmat<T>& A, const sZmat<T>& B, const T& p );
template<class T> sZmat<T> mult_mod_p_flint ( const sZmat<T>& A, const sZmat<T>& B, const T& p );
template<class T> T maxabs( const sZmat<T>& A);
template<class T> sZmat<T> transpose(const sZmat<T>&);
template<class T> int operator==(const sZmat<T>&, const sZmat<T>&);
template<class T> int eqmodp(const sZmat<T>&, const sZmat<T>&, const T& p);
template<class T> ostream& operator<< (ostream&s, const sZmat<T>&);
template<class T> istream& operator>> (istream&s, sZmat<T>&);
template<class T> int get_population (const sZmat<T>& );
template<class T> inline double density (const sZmat<T>& m)
{return (((double)(get_population(m)))/m.nro)/m.nco;}
template<class T> void random_fill_in( sZmat<T>&, int, int );
template<class T> int liftmat(const sZmat<T>& mm, T pr, sZmat<T>& m, T& dd);
template<class T> int liftmats_chinese(const sZmat<T>& mm1, T pr1, const sZmat<T>& mm2, T pr2,
                              sZmat<T>& m, T& dd);
template<class T> int dim(const ssubZspace<T>& s)     {return s.bas().ncols();}
template<class T> sZmat<T> basis(const ssubZspace<T>& s)  {return s.bas();}

template<class T>
class Zvec {
  friend class sZvec<T>;
  friend class sZmat<T>;
  friend class sZmat_elim<T>;
  friend class Zmat<T>;
  friend class subZspace<T>;
public:
  // constructors
  explicit Zvec(long n=0);
  explicit Zvec(const vector<T>& arr);
  Zvec(const Zvec<T>&);                       // copy constructor
  // member functions & operators
  void init(long n=0);                 // (re)-initializes
  Zvec& operator=(const Zvec<T>&);         // assignment
  T& operator[](long i);            // the i'th component
  T operator[](long i) const;       // the i'th component
  Zvec& operator+=(const Zvec<T>&);
  void add_row(const Zmat<T>&m, int i);
  void addmodp(const Zvec<T>&, const T&);
  Zvec& operator-=(const Zvec<T>&);
  void sub_row(const Zmat<T>&m, int i);
  Zvec& operator*=(const T&);
  Zvec& operator/=(const T&);
  Zvec slice(long i, long j=-1) const;  // returns subZvec with indices i..j,  or 1..i if j=-1
  Zvec operator[](const Zvec<int>&) const;   // subscript composition
  Zvec operator[](const Zvec<long>&) const;   // subscript composition
  void set(long i, const T& x);                  // sets v[i]=x
  void add(long i, const T& x);                  // v[i]+=x
  void add_modp(long i, const T& x, const T& p);                  // v[i]+=x mod p
  void reduce_mod_p(const T& p);              // reduce mod p in place
  T sub(long i) const;                    // same as v[i] (no ref)
  const vector<T> get_entries()const {return entries;}
  static Zvec iota(long n);              // (1,2,...,n)
  // non-member (friend) functions and operators
  friend int dim<>(const Zvec<T>&);                  // the dimension
  friend T operator*<>(const Zvec<T>&, const Zvec<T>&);   // dot product
  friend T operator*<>(const sZvec<T>&, const Zvec<T>&);
  friend T content<>(const Zvec<T>&);
  friend T maxabs<>(const Zvec<T>&);
  friend Zvec operator*<>(const Zmat<T>& m, const Zvec<T>& v);
  friend int operator==<>(const Zvec<T>&, const Zvec<T>&);
  friend int operator!=<>(const Zvec<T>&, const Zvec<T>&);
  friend int trivial<>(const Zvec<T>&);                  // is v all 0
  friend int member<>(const T& a, const Zvec<T>& v);//tests if a=v[i] for some i
  friend Zvec reverse<>(const Zvec<T>& order);
  // add/sub row i of Zmat to v (implemented in mat.cc)
  friend void add_row_to_vec<>(Zvec<T>& v, const Zmat<T>& m, long i);
  friend void sub_row_to_vec<>(Zvec<T>& v, const Zmat<T>& m, long i);
  friend ostream& operator<< <>(ostream&s, const Zvec<T>&);
  friend istream& operator>> <>(istream&s, Zvec<T>&);
  friend void swapvec<>(Zvec<T>& v, Zvec<T>& w);

  // Implementation
private:
  vector<T> entries;
};

// Declaration of non-member, non-friend functions

template<class T> T content(const Zvec<T>&);
template<class T> Zvec<T> operator+(const Zvec<T>&);                   // unary
template<class T> Zvec<T> operator-(const Zvec<T>&);                   // unary
template<class T> Zvec<T> operator+(const Zvec<T>&, const Zvec<T>&);
template<class T> Zvec<T> addmodp(const Zvec<T>&, const Zvec<T>&, const T&);
template<class T> Zvec<T> operator-(const Zvec<T>&, const Zvec<T>&);
template<class T> Zvec<T> operator*(const T&, const Zvec<T>&);       // componentwise
template<class T> Zvec<T> operator/(const Zvec<T>&, const T&);       // componentwise
template<class T> void make_primitive(Zvec<T>& v);
template<class T> void elim(const Zvec<T>& a, Zvec<T>& b, long pos);
template<class T> void elim1(const Zvec<T>& a, Zvec<T>& b, long pos);
template<class T> void elim2(const Zvec<T>& a, Zvec<T>& b, long pos, const T& lastpivot);
template<class T> Zvec<T> express(const Zvec<T>& v, const Zvec<T>& v1, const Zvec<T>& v2);
template<class T> int lift(const Zvec<T>& v, const T& pr, Zvec<T>& ans);  //lifts a mod-p vector to a rational
                                   //and scales to a primitive vec in Z. Returns success flag

// inline function definitions

template<class T> inline int dim(const Zvec<T>& v) {return v.entries.size();}
template<class T> inline int operator!=(const Zvec<T>& v, const Zvec<T>& w) { return !(v==w);}
template<class T> inline Zvec<T> operator+(const Zvec<T>& v) { return v;}
template<class T> inline Zvec<T> operator-(const Zvec<T>& v) { return T(-1)*v;}
template<class T> inline Zvec<T> operator+(const Zvec<T>& v1, const Zvec<T>& v2)
{ Zvec<T> w(v1); w+=v2; return w;}
template<class T> inline Zvec<T> addmodp(const Zvec<T>& v1, const Zvec<T>& v2, const T& pr)
{ Zvec<T> w(v1); w.addmodp(v2,pr); return w;}
template<class T> inline Zvec<T> reduce_mod_p(const Zvec<T>& v, const T& p)
{ Zvec<T> w(v); w.reduce_mod_p(p); return w;}
template<class T> inline Zvec<T> operator-(const Zvec<T>& v1, const Zvec<T>& v2)
{ Zvec<T> w(v1); w-=v2; return w;}
template<class T> inline Zvec<T> operator*(const T& scal, const Zvec<T>& v)
{ Zvec<T> w(v); w*=scal; return w;}
template<class T> inline Zvec<T> operator/(const Zvec<T>& v, const T& scal)
{ Zvec<T> w(v); w/=scal; return w;}
template<class T> inline void make_primitive(Zvec<T>& v)
{ T g=content(v); if (g>1) v/=g;}
template<class T> inline void elim(const Zvec<T>& a, Zvec<T>& b, long pos)
{ (b*=a[pos])-=(b[pos]*a);}
template<class T> inline void elim1(const Zvec<T>& a, Zvec<T>& b, long pos)
{ (b*=a[pos])-=(b[pos]*a); make_primitive(b);}
template<class T> inline void elim2(const Zvec<T>& a, Zvec<T>& b, long pos, const T& lastpivot)
{ ((b*=a[pos])-=(b[pos]*a))/=lastpivot;}

// conversions between vectors of different scalar types

vec_m to_vec_m(const vec_i& v);
vec_m to_vec_m(const vec_l& v);
inline vec_m to_vec_m(const vec_m& v) {return v;}

vec_l to_vec_l(const vec_i& v);
inline vec_l to_vec_l(const vec_l& v) {return v;}
vec_l to_vec_l(const vec_m& v);

inline vec_i to_vec_i(const vec_i& v) {return v;}
vec_i to_vec_i(const vec_l& v);
vec_i to_vec_i(const vec_m& v);

#endif
