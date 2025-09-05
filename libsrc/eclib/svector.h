// svector.h: manage declarations for sparse integer vector classes
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
 
#if     !defined(_ECLIB_SVECTOR_H)
#define _ECLIB_SVECTOR_H      1       //flags that this file has been included

#include "vector.h"

template<class T> int dim(const sZvec<T>& v);
template<class T> T operator*(const Zvec<T>& v, const sZvec<T>& sv);
template<class T> T dotmodp(const Zvec<T>& v, const sZvec<T>& sv, const T& pr);
template<class T> int operator!=(const sZvec<T>& v1, const Zvec<T>& v2);
template<class T> int operator==(const Zvec<T>& v1, const sZvec<T>& v2);
template<class T> int operator!=(const Zvec<T>& v1, const sZvec<T>& v2);


template<class T>
class sZvec {

  friend class sZmat<T>;
  friend class sZmat_elim<T>;

protected:
  int d;              // dimension
  map<int,T> entries;  // (i,x) in the table means v[i]=x

public:
     // constructors

  sZvec<T> (int dim=0) :d(dim) {;}
  explicit sZvec<T> (const Zvec<T> &);                  // conversion constructor

  // member functions & operators

  void clear() {entries.clear();}
  Zvec<T> as_vec( ) const;
  T  elem(int i) const;   // returns value of i'th entry
  void set(int i, const T& a);   // sets i'th entry to a
  void add(int i, const T& a);   // adds a to i'th entry
  void sub(int i, const T& a);   // subtracts a from i'th entry
  void add_mod_p(int i, const T& a, const T& p);   // adds a to i'th entry, mod p
  void sub_mod_p(int i, const T& a, const T& p);   // subtracts a from i'th entry, mod p
  sZvec<T>& add_scalar_times(const sZvec<T>&, const T&);
  sZvec<T>& operator+= (const sZvec<T>& w);
  sZvec<T>& operator-= (const sZvec<T>& w);
  sZvec<T>& operator*= (const T&);
  sZvec<T>& operator/= (const T&);
  void reduce_mod_p(const T& p);
  sZvec<T>& mult_by_scalar_mod_p(const T&, const T& p);
  sZvec<T>& add_scalar_times_mod_p(const sZvec<T>&, const T&, const T& p);
  // Same as previous except returns two sets of indices: "ons" are
  // indices for which an entry is created, and "offs" are indices for
  // which an entry is deleted
  sZvec<T>& add_scalar_times_mod_p(const sZvec<T>&, const T&, std::set<int>& ons, std::set<int>& offs,
			       const T& p);
  // two hand-coded special cases:
  sZvec<T>& add_mod_p(const sZvec<T>& w, const T& p);
  sZvec<T>& sub_mod_p(const sZvec<T>& w, const T& p);

  int size() const {return entries.size();}

  // functions to enable iteration over the entries without direct
  // access to the entries map:
  typename map<int,T>::const_iterator begin() const {return entries.begin();}
  typename map<int,T>::const_iterator end() const {return entries.end();}
  typename map<int,T>::iterator begin() {return entries.begin();}
  typename map<int,T>::iterator end() {return entries.end();}
  void erase(int i);  // erases v[i]; error if not set
  int first_index() const {return entries.upper_bound(0)->first;}
  std::set<int> support() const;

  // non-member (friend) functions and operators

  friend int dim<>(const sZvec<T>& v);
  friend int eqmodp<>(const sZvec<T>&, const sZvec<T>&, const T& p);
  friend ostream& operator<<<> (ostream&s, const sZvec<T>&);
  friend T operator*<>(const sZvec<T>&, const sZvec<T>&);
  friend T operator*<>(const sZvec<T>&, const Zvec<T>&);
  friend T operator*<>(const Zvec<T>& v, const sZvec<T>& sv);
  friend T dotmodp<>(const sZvec<T>&, const sZvec<T>&, const T& pr);
  friend T dotmodp<>(const sZvec<T>&, const Zvec<T>&, const T& pr);
  friend T dotmodp<>(const Zvec<T>& v, const sZvec<T>& sv, const T& pr);
  friend sZvec<T> operator+<>(const sZvec<T>& v1, const sZvec<T>& v2);
  friend sZvec<T> operator-<>(const sZvec<T>& v1, const sZvec<T>& v2);
  friend int operator==<>(const sZvec<T>& v1, const sZvec<T>& v2);
  friend int operator!=<>(const sZvec<T>& v1, const sZvec<T>& v2);
  friend int operator==<>(const sZvec<T>& v1, const Zvec<T>& v2);
  friend int operator!=<>(const sZvec<T>& v1, const Zvec<T>& v2);
  friend int operator==<>(const Zvec<T>& v1, const sZvec<T>& v2);
  friend int operator!=<>(const Zvec<T>& v1, const sZvec<T>& v2);
  friend int trivial<>(const Zvec<T>&);
  friend sZmat<T> transpose<>(const sZmat<T>&);
  friend sZmat<T> operator*<> ( const sZmat<T>&, const sZmat<T>&);
  friend T content<>(const sZvec<T>& v);
  friend T make_primitive<>(sZvec<T>& v); // divides by & returns content
  friend sZvec<T> operator*<> ( const sZmat<T>& A, const sZvec<T>& v );
  friend sZvec<T> operator*<> ( const sZvec<T>& v, const sZmat<T>& A );
  friend sZvec<T> mult_mod_p<>( const sZmat<T>& A, const sZvec<T>& v, const T& p  );
  friend sZvec<T> mult_mod_p<>( const sZvec<T>& v, const sZmat<T>& A, const T& p  );
  friend sZmat<T> mult_mod_p<>( const sZmat<T>&, const sZmat<T>&, const T&);
};

// Declaration of non-friend functions

template<class T>
inline sZvec<T> operator+(const sZvec<T>& v) {return v;}      // unary +
template<class T>
inline sZvec<T> operator-(const sZvec<T>& v)                  // unary -
{sZvec<T> ans(v); ans*=T(-1); return ans;}
template<class T>
inline sZvec<T> operator+(const sZvec<T>& v1, const sZvec<T>& v2)
{
  if(v1.entries.size()<v2.entries.size())
    {
      sZvec<T> ans(v2); ans+=v1; return ans;
    }
  else
    {
      sZvec<T> ans(v1); ans+=v2; return ans;
    }
}
template<class T>
inline sZvec<T> operator-(const sZvec<T>& v1, const sZvec<T>& v2)
{
  return v1 + (-v2);
}

template<class T>
inline sZvec<T> operator*(const T& scal, const sZvec<T>& v)
{sZvec<T> ans(v); ans*=scal; return ans;}

template<class T>
inline sZvec<T> operator/(const sZvec<T>& v, const T& scal)
{sZvec<T> ans(v); ans/=scal; return ans;}

template<class T>
inline int operator==(const sZvec<T>& v1, const sZvec<T>& v2)
{
  return (v1.d==v2.d) && (v1.entries == v2.entries);
}

template<class T>
inline int operator!=(const sZvec<T>& v1, const sZvec<T>& v2)
{
  return !(v1==v2);
}

#endif
