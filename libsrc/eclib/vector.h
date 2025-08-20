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

#include <eclib/marith.h>

template<class T> class vecT;
template<class T> class svecT;
template<class T> class smatT;
template<class T> class smatT_elim;
template<class T> class matT;
template<class T> class subspaceT;
template<class T> class ssubspaceT;

typedef vecT<int> vec_i;
typedef matT<int> mat_i;
typedef subspaceT<int> subspace_i;
typedef svecT<int> svec_i;
typedef smatT<int> smat_i;
typedef smatT_elim<int> smat_i_elim;
typedef ssubspaceT<int> ssubspace_i;

typedef vecT<long> vec_l;
typedef matT<long> mat_l;
typedef subspaceT<long> subspace_l;
typedef svecT<long> svec_l;
typedef smatT<long> smat_l;
typedef smatT_elim<long> smat_l_elim;
typedef ssubspaceT<long> ssubspace_l;

typedef vecT<bigint> vec_m;
typedef matT<bigint> mat_m;
typedef subspaceT<bigint> subspace_m;
typedef svecT<bigint> svec_m;
typedef smatT<bigint> smat_m;
typedef smatT_elim<bigint> smat_m_elim;
typedef ssubspaceT<bigint> ssubspace_m;

template<class T> int dim(const vecT<T>&);                  // the dimension
template<class T> T operator*(const vecT<T>&, const vecT<T>&);   // dot product
template<class T> T operator*(const svecT<T>&, const vecT<T>&);
template<class T> T content(const vecT<T>&);
template<class T> T maxabs(const vecT<T>&);
template<class T> vecT<T> operator*(const matT<T>& m, const vecT<T>& v);
template<class T> int operator==(const vecT<T>&, const vecT<T>&);
template<class T> int operator!=(const vecT<T>&, const vecT<T>&);
template<class T> int trivial(const vecT<T>&);                  // is v all 0
template<class T> int member(const T& a, const vecT<T>& v);//tests if a=v[i] for some i
template<class T> vecT<T> reverse(const vecT<T>& order);
template<class T> void add_row_to_vec(vecT<T>& v, const matT<T>& m, long i);
template<class T> void sub_row_to_vec(vecT<T>& v, const matT<T>& m, long i);
template<class T> ostream& operator<< (ostream&s, const vecT<T>&);
template<class T> istream& operator>> (istream&s, vecT<T>&);
template<class T> void swapvec(vecT<T>& v, vecT<T>& w);
template<class T> matT<T> operator*(const matT<T>&, const matT<T>&);
template<class T> vecT<T> operator*(const matT<T>&, const vecT<T>&);
template<class T> int operator==(const matT<T>&, const matT<T>&);
template<class T> istream& operator>> (istream&s, matT<T>&);
template<class T> matT<T> colcat(const matT<T>& a, const matT<T>& b);
template<class T> matT<T> rowcat(const matT<T>& a, const matT<T>& b);
template<class T> matT<T> directsum(const matT<T>& a, const matT<T>& b);
template<class T> void elimrows(matT<T>& m, long r1, long r2, long pos); //plain elimination, no clearing
template<class T> void elimrows1(matT<T>& m, long r1, long r2, long pos); //elimination + clearing
template<class T> void elimrows2(matT<T>& m, long r1, long r2, long pos, const T& last); //elimination + divide by last pivot
template<class T> matT<T> echelon0(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                      long& rk, long& ny, T& d);
template<class T> void elimp(matT<T>& m, long r1, long r2, long pos, const T& pr);
template<class T> void elimp1(matT<T>& m, long r1, long r2, long pos, const T& pr);
template<class T> matT<T> echelonp(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                      long& rk, long& ny, T& d, const T& pr);
template<class T> matT<T> echmodp(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                     long& rk, long& ny, const T& pr);
template<class T> matT<T> echmodp_uptri(const matT<T>& m, vecT<int>& pcols, vecT<int>& npcols,
                     long& rk, long& ny, const T& pr);
template<class T> matT<T> ref_via_flint(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                           long& rk, long& ny, const T& pr);
template<class T> matT<T> ref_via_ntl(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                         long& rk, long& ny, const T& pr);
template<class T> matT<T> rref(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                               long& rk, long& ny, const T& pr);
template<class T> long rank_via_ntl(const matT<T>& M, const T& pr);
template<class T> T det_via_ntl(const matT<T>& M, const T& pr);
template<class T> subspaceT<T> combine(const subspaceT<T>& s1, const subspaceT<T>& s2);
template<class T> matT<T> restrict_mat(const matT<T>& m, const subspaceT<T>& s, int cr=0);
template<class T> int liftmat(const matT<T>& mm, const T& pr, matT<T>& m, T& dd);
template<class T> int lift(const subspaceT<T>& s, const T& pr, subspaceT<T>& ans);
template<class T> subspaceT<T> pcombine(const subspaceT<T>& s1, const subspaceT<T>& s2, const T& pr);
template<class T> matT<T> matmulmodp(const matT<T>&, const matT<T>&, const T& pr);
template<class T> long population(const matT<T>& m); // #nonzero entries
template<class T> T maxabs(const matT<T>& m); // max entry
template<class T> double sparsity(const matT<T>& m); // #nonzero entries/#entries
template<class T> matT<T> prestrict(const matT<T>& m, const subspaceT<T>& s, const T& pr, int cr=0);
template<class T> T dotmodp(const vecT<T>& v1, const vecT<T>& v2, const T& pr);
template<class T> int dim(const subspaceT<T>& s);
template<class T> T denom(const subspaceT<T>& s);
template<class T> vecT<int> pivots(const subspaceT<T>& s);
template<class T> matT<T> basis(const subspaceT<T>& s);
template<class T> subspaceT<T> combine(const subspaceT<T>& s1, const subspaceT<T>& s2);
template<class T> subspaceT<T> pcombine(const subspaceT<T>& s1, const subspaceT<T>& s2, const T& pr);
template<class T> int lift(const subspaceT<T>& s, const T& pr, subspaceT<T>& ans);
template<class T> int eqmodp(const svecT<T>&, const svecT<T>&, const T& p);
template<class T> ostream& operator<< (ostream&s, const svecT<T>&);
template<class T> T operator*(const svecT<T>&, const svecT<T>&);
template<class T> T operator*(const svecT<T>&, const vecT<T>&);
template<class T> T operator*(const vecT<T>& v, const svecT<T>& sv);
template<class T> T dotmodp(const svecT<T>&, const svecT<T>&, const T& pr);
template<class T> T dotmodp(const svecT<T>&, const vecT<T>&, const T& pr);
template<class T> T dotmodp(const vecT<T>& v, const svecT<T>& sv, const T& pr);
template<class T> svecT<T> operator+(const svecT<T>& v1, const svecT<T>& v2);
template<class T> svecT<T> operator-(const svecT<T>& v1, const svecT<T>& v2);
template<class T> int operator==(const svecT<T>& v1, const svecT<T>& v2);
template<class T> int operator!=(const svecT<T>& v1, const svecT<T>& v2);
template<class T> int operator==(const svecT<T>& v1, const vecT<T>& v2);
template<class T> smatT<T> transpose(const smatT<T>&);
template<class T> smatT<T> operator* ( const smatT<T>&, const smatT<T>&);
template<class T> T content(const svecT<T>& v);
template<class T> T make_primitive(svecT<T>& v);
template<class T> svecT<T> operator* ( const smatT<T>& A, const svecT<T>& v );
template<class T> svecT<T> operator* ( const svecT<T>& v, const smatT<T>& A );
template<class T> svecT<T> mult_mod_p( const smatT<T>& A, const svecT<T>& v, const T& p  );
template<class T> svecT<T> mult_mod_p( const svecT<T>& v, const smatT<T>& A, const T& p  );
template<class T> smatT<T> mult_mod_p ( const smatT<T>&, const smatT<T>&, const T&);
template<class T> inline vector<int> dim(const smatT<T>& A)
{return vector<int>{A.nro, A.nco};}
template<class T> vecT<T> operator*  (const smatT<T>& m, const vecT<T>& v);
template<class T> svecT<T> operator* ( const smatT<T>& A, const svecT<T>& v );
template<class T> svecT<T> operator* ( const svecT<T>& v, const smatT<T>& A );
template<class T> svecT<T> mult_mod_p( const smatT<T>& A, const svecT<T>& v, const T& p  );
template<class T> vecT<T> mult_mod_p( const smatT<T>& A, const vecT<T>& v, const T& p  );
template<class T> svecT<T> mult_mod_p( const svecT<T>& v, const smatT<T>& A, const T& p  );
template<class T> smatT<T> operator* ( const smatT<T>& A, const smatT<T>& B );
template<class T> smatT<T> mult_mod_p ( const smatT<T>& A, const smatT<T>& B, const T& p );
template<class T> smatT<T> mult_mod_p_flint ( const smatT<T>& A, const smatT<T>& B, const T& p );
template<class T> T maxabs( const smatT<T>& A);
template<class T> smatT<T> transpose(const smatT<T>&);
template<class T> int operator==(const smatT<T>&, const smatT<T>&);
template<class T> int eqmodp(const smatT<T>&, const smatT<T>&, const T& p);
template<class T> ostream& operator<< (ostream&s, const smatT<T>&);
template<class T> istream& operator>> (istream&s, smatT<T>&);
template<class T> int get_population (const smatT<T>& );
template<class T> inline double density (const smatT<T>& m)
{return (((double)(get_population(m)))/m.nro)/m.nco;}
template<class T> void random_fill_in( smatT<T>&, int, int );
template<class T> int liftmat(const smatT<T>& mm, T pr, smatT<T>& m, T& dd);
template<class T> int liftmats_chinese(const smatT<T>& mm1, T pr1, const smatT<T>& mm2, T pr2,
                              smatT<T>& m, T& dd);
template<class T> int dim(const ssubspaceT<T>& s)     {return s.bas().ncols();}
template<class T> smatT<T> basis(const ssubspaceT<T>& s)  {return s.bas();}
template<class T> ssubspaceT<T> combine(const ssubspaceT<T>& s1, const ssubspaceT<T>& s2);
template<class T> smatT<T> restrict_mat(const smatT<T>& m, const ssubspaceT<T>& s);

template<class T>
class vecT {
  friend class svecT<T>;
  friend class smatT<T>;
  friend class smatT_elim<T>;
  friend class matT<T>;
  friend class subspaceT<T>;
public:
  // constructors
  explicit vecT(long n=0);
  explicit vecT(const vector<T>& arr);
  vecT(const vecT<T>&);                       // copy constructor
  // member functions & operators
  void init(long n=0);                 // (re)-initializes
  vecT& operator=(const vecT<T>&);         // assignment
  T& operator[](long i);            // the i'th component
  T operator[](long i) const;       // the i'th component
  vecT& operator+=(const vecT<T>&);
  void add_row(const matT<T>&m, int i);
  void addmodp(const vecT<T>&, const T&);
  vecT& operator-=(const vecT<T>&);
  void sub_row(const matT<T>&m, int i);
  vecT& operator*=(const T&);
  vecT& operator/=(const T&);
  vecT slice(long i, long j=-1) const;  // returns subvecT with indices i..j,  or 1..i if j=-1
  vecT operator[](const vecT<int>&) const;   // subscript composition
  vecT operator[](const vecT<long>&) const;   // subscript composition
  void set(long i, const T& x);                  // sets v[i]=x
  void add(long i, const T& x);                  // v[i]+=x
  void add_modp(long i, const T& x, const T& p);                  // v[i]+=x mod p
  void reduce_mod_p(const T& p);              // reduce mod p in place
  T sub(long i) const;                    // same as v[i] (no ref)
  const vector<T> get_entries()const {return entries;}
  static vecT iota(long n);              // (1,2,...,n)
  // non-member (friend) functions and operators
  friend int dim<>(const vecT<T>&);                  // the dimension
  friend T operator*<>(const vecT<T>&, const vecT<T>&);   // dot product
  friend T operator*<>(const svecT<T>&, const vecT<T>&);
  friend T content<>(const vecT<T>&);
  friend T maxabs<>(const vecT<T>&);
  friend vecT operator*<>(const matT<T>& m, const vecT<T>& v);
  friend int operator==<>(const vecT<T>&, const vecT<T>&);
  friend int operator!=<>(const vecT<T>&, const vecT<T>&);
  friend int trivial<>(const vecT<T>&);                  // is v all 0
  friend int member<>(const T& a, const vecT<T>& v);//tests if a=v[i] for some i
  friend vecT reverse<>(const vecT<T>& order);
  // add/sub row i of matT to v (implemented in mat.cc)
  friend void add_row_to_vec<>(vecT<T>& v, const matT<T>& m, long i);
  friend void sub_row_to_vec<>(vecT<T>& v, const matT<T>& m, long i);
  friend ostream& operator<< <>(ostream&s, const vecT<T>&);
  friend istream& operator>> <>(istream&s, vecT<T>&);
  friend void swapvec<>(vecT<T>& v, vecT<T>& w);

  // Implementation
private:
  vector<T> entries;
};

// Declaration of non-member, non-friend functions

template<class T> T content(const vecT<T>&);
template<class T> vecT<T> operator+(const vecT<T>&);                   // unary
template<class T> vecT<T> operator-(const vecT<T>&);                   // unary
template<class T> vecT<T> operator+(const vecT<T>&, const vecT<T>&);
template<class T> vecT<T> addmodp(const vecT<T>&, const vecT<T>&, const T&);
template<class T> vecT<T> operator-(const vecT<T>&, const vecT<T>&);
template<class T> vecT<T> operator*(const T&, const vecT<T>&);       // componentwise
template<class T> vecT<T> operator/(const vecT<T>&, const T&);       // componentwise
template<class T> void make_primitive(vecT<T>& v);
template<class T> void elim(const vecT<T>& a, vecT<T>& b, long pos);
template<class T> void elim1(const vecT<T>& a, vecT<T>& b, long pos);
template<class T> void elim2(const vecT<T>& a, vecT<T>& b, long pos, const T& lastpivot);
template<class T> vecT<T> express(const vecT<T>& v, const vecT<T>& v1, const vecT<T>& v2);
template<class T> int lift(const vecT<T>& v, const T& pr, vecT<T>& ans);  //lifts a mod-p vector to a rational
                                   //and scales to a primitive vec in Z. Returns success flag

// inline function definitions

template<class T> inline int dim(const vecT<T>& v) {return v.entries.size();}
template<class T> inline int operator!=(const vecT<T>& v, const vecT<T>& w) { return !(v==w);}
template<class T> inline vecT<T> operator+(const vecT<T>& v) { return v;}
template<class T> inline vecT<T> operator-(const vecT<T>& v) { return T(-1)*v;}
template<class T> inline vecT<T> operator+(const vecT<T>& v1, const vecT<T>& v2)
{ vecT<T> w(v1); w+=v2; return w;}
template<class T> inline vecT<T> addmodp(const vecT<T>& v1, const vecT<T>& v2, const T& pr)
{ vecT<T> w(v1); w.addmodp(v2,pr); return w;}
template<class T> inline vecT<T> reduce_mod_p(const vecT<T>& v, const T& p)
{ vecT<T> w(v); w.reduce_mod_p(p); return w;}
template<class T> inline vecT<T> operator-(const vecT<T>& v1, const vecT<T>& v2)
{ vecT<T> w(v1); w-=v2; return w;}
template<class T> inline vecT<T> operator*(const T& scal, const vecT<T>& v)
{ vecT<T> w(v); w*=scal; return w;}
template<class T> inline vecT<T> operator/(const vecT<T>& v, const T& scal)
{ vecT<T> w(v); w/=scal; return w;}
template<class T> inline void make_primitive(vecT<T>& v)
{ T g=content(v); if (g>1) v/=g;}
template<class T> inline void elim(const vecT<T>& a, vecT<T>& b, long pos)
{ (b*=a[pos])-=(b[pos]*a);}
template<class T> inline void elim1(const vecT<T>& a, vecT<T>& b, long pos)
{ (b*=a[pos])-=(b[pos]*a); make_primitive(b);}
template<class T> inline void elim2(const vecT<T>& a, vecT<T>& b, long pos, const T& lastpivot)
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
