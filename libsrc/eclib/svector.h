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

#include <eclib/vector.h>

template<class T> int dim(const svecT<T>& v);
template<class T> T operator*(const vecT<T>& v, const svecT<T>& sv);
template<class T> T dotmodp(const vecT<T>& v, const svecT<T>& sv, const T& pr);
template<class T> int operator!=(const svecT<T>& v1, const vecT<T>& v2);
template<class T> int operator==(const vecT<T>& v1, const svecT<T>& v2);
template<class T> int operator!=(const vecT<T>& v1, const svecT<T>& v2);


template<class T>
class svecT {

  friend class smatT<T>;
  friend class smatT_elim<T>;

protected:
  int d;              // dimension
  map<int,T> entries;  // (i,x) in the table means v[i]=x

public:
     // constructors

  svecT<T> (int dim=0) :d(dim) {;}
  explicit svecT<T> (const vecT<T> &);                  // conversion constructor

  // member functions & operators

  void clear() {entries.clear();}
  vecT<T> as_vec( ) const;
  T  elem(int i) const;   // returns value of i'th entry
  void set(int i, const T& a);   // sets i'th entry to a
  void add(int i, const T& a);   // adds a to i'th entry
  void sub(int i, const T& a);   // subtracts a from i'th entry
  void add_mod_p(int i, const T& a, const T& p);   // adds a to i'th entry, mod p
  void sub_mod_p(int i, const T& a, const T& p);   // subtracts a from i'th entry, mod p
  svecT<T>& add_scalar_times(const svecT<T>&, const T&);
  svecT<T>& operator+= (const svecT<T>& w);
  svecT<T>& operator-= (const svecT<T>& w);
  svecT<T>& operator*= (const T&);
  svecT<T>& operator/= (const T&);
  void reduce_mod_p(const T& p);
  svecT<T>& mult_by_scalar_mod_p(const T&, const T& p);
  svecT<T>& add_scalar_times_mod_p(const svecT<T>&, const T&, const T& p);
  // Same as previous except returns two sets of indices: "ons" are
  // indices for which an entry is created, and "offs" are indices for
  // which an entry is deleted
  svecT<T>& add_scalar_times_mod_p(const svecT<T>&, const T&, std::set<int>& ons, std::set<int>& offs,
			       const T& p);
  // two hand-coded special cases:
  svecT<T>& add_mod_p(const svecT<T>& w, const T& p);
  svecT<T>& sub_mod_p(const svecT<T>& w, const T& p);

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

  friend int dim<>(const svecT<T>& v);
  friend int eqmodp<>(const svecT<T>&, const svecT<T>&, const T& p);
  friend ostream& operator<<<> (ostream&s, const svecT<T>&);
  friend T operator*<>(const svecT<T>&, const svecT<T>&);
  friend T operator*<>(const svecT<T>&, const vecT<T>&);
  friend T operator*<>(const vecT<T>& v, const svecT<T>& sv);
  friend T dotmodp<>(const svecT<T>&, const svecT<T>&, const T& pr);
  friend T dotmodp<>(const svecT<T>&, const vecT<T>&, const T& pr);
  friend T dotmodp<>(const vecT<T>& v, const svecT<T>& sv, const T& pr);
  friend svecT<T> operator+<>(const svecT<T>& v1, const svecT<T>& v2);
  friend svecT<T> operator-<>(const svecT<T>& v1, const svecT<T>& v2);
  friend int operator==<>(const svecT<T>& v1, const svecT<T>& v2);
  friend int operator!=<>(const svecT<T>& v1, const svecT<T>& v2);
  friend int operator==<>(const svecT<T>& v1, const vecT<T>& v2);
  friend int operator!=<>(const svecT<T>& v1, const vecT<T>& v2);
  friend int operator==<>(const vecT<T>& v1, const svecT<T>& v2);
  friend int operator!=<>(const vecT<T>& v1, const svecT<T>& v2);
  friend smatT<T> transpose<>(const smatT<T>&);
  friend smatT<T> operator*<> ( const smatT<T>&, const smatT<T>&);
  friend T content<>(const svecT<T>& v);
  friend T make_primitive<>(svecT<T>& v); // divides by & returns content
  friend svecT<T> operator*<> ( const smatT<T>& A, const svecT<T>& v );
  friend svecT<T> operator*<> ( const svecT<T>& v, const smatT<T>& A );
  friend svecT<T> mult_mod_p<>( const smatT<T>& A, const svecT<T>& v, const T& p  );
  friend svecT<T> mult_mod_p<>( const svecT<T>& v, const smatT<T>& A, const T& p  );
  friend smatT<T> mult_mod_p<>( const smatT<T>&, const smatT<T>&, const T&);
};

// Declaration of non-friend functions

template<class T>
inline svecT<T> operator+(const svecT<T>& v) {return v;}      // unary +
template<class T>
inline svecT<T> operator-(const svecT<T>& v)                  // unary -
{svecT<T> ans(v); ans*=T(-1); return ans;}
template<class T>
inline svecT<T> operator+(const svecT<T>& v1, const svecT<T>& v2)
{
  if(v1.entries.size()<v2.entries.size())
    {
      svecT<T> ans(v2); ans+=v1; return ans;
    }
  else
    {
      svecT<T> ans(v1); ans+=v2; return ans;
    }
}
template<class T>
inline svecT<T> operator-(const svecT<T>& v1, const svecT<T>& v2)
{
  return v1 + (-v2);
}

template<class T>
inline svecT<T> operator*(const T& scal, const svecT<T>& v)
{svecT<T> ans(v); ans*=scal; return ans;}

template<class T>
inline svecT<T> operator/(const svecT<T>& v, const T& scal)
{svecT<T> ans(v); ans/=scal; return ans;}

template<class T>
inline int operator==(const svecT<T>& v1, const svecT<T>& v2)
{
  return (v1.d==v2.d) && (v1.entries == v2.entries);
}

template<class T>
inline int operator!=(const svecT<T>& v1, const svecT<T>& v2)
{
  return !(v1==v2);
}

#endif
