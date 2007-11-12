// svec.h: declarations for sparse integer vector class svec
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
 
class svec {

  friend class smat;
  friend class smat_elim;

protected:
  int d;              // dimension
  map<int,scalar> entries;  // (i,x) in the table means v[i]=x

public:
     // constructors

  svec (int dim=0) :d(dim) {;}
  svec (const svec& v)                // copy constructor
    :d(v.d),     entries(v.entries)
  {;  }
  svec (const vec &);                  // conversion constructor
  svec (int dim, scalar* a);          // conversion constructor
   
  // member functions & operators

  void clear() {entries.clear();}
  vec as_vec( ) const;
  scalar  elem(int i) const;   // returns value of i'th entry
  void set(int i, scalar a);   // sets i'th entry to a
  void add(int i, scalar a);   // adds a to i'th entry
  void sub(int i, scalar a);   // subtracts a from i'th entry
  void add_mod_p(int i, scalar a, const scalar& p=BIGPRIME);   // adds a to i'th entry, mod p
  void sub_mod_p(int i, scalar a, const scalar& p=BIGPRIME);   // subtracts a from i'th entry, mod p
  svec& add_scalar_times(const svec&, scalar);
  svec& operator+= (const svec& w);
  svec& operator-= (const svec& w);
  svec& operator*= (scalar);
  svec& operator/= (scalar);
  void reduce_mod_p(const scalar& p=BIGPRIME);
  svec& mult_by_scalar_mod_p(scalar, const scalar& p=BIGPRIME);
  svec& add_scalar_times_mod_p(const svec&, scalar, const scalar& p=BIGPRIME);
  // Same as previous except returns two sets of indices: "ons" are
  // indices for which an entry is created, and "offs" are indices for
  // which an entry is deleted
  svec& add_scalar_times_mod_p(const svec&, scalar, std::set<int>& ons, std::set<int>& offs, 
			       const scalar& p=BIGPRIME);

  int size() const {return entries.size();}
     
  // functions to enable iteration over the entries without direct
  // access to the entries map:
  map<int,scalar>::const_iterator begin() const {return entries.begin();}
  map<int,scalar>::const_iterator end() const {return entries.end();}
  map<int,scalar>::iterator begin() {return entries.begin();}
  map<int,scalar>::iterator end() {return entries.end();}
  void erase(int i);  // erases v[i]; error if not set
  int first_index() const {return entries.upper_bound(0)->first;}
  std::set<int> support() const;

  // non-member (friend) functions and operators

  friend inline int dim(const svec& v)  {return v.d;}
  // Equality mod p:
  friend int eqmodp(const svec&, const svec&, const scalar& p=BIGPRIME);
  friend ostream& operator<< (ostream&s, const svec&);
  friend scalar operator*(const svec&, const svec&); //dot product
  friend scalar operator*(const svec&, const vec&);
  friend scalar dotmodp(const svec&, const svec&, scalar pr);
  friend scalar dotmodp(const svec&, const vec&, scalar pr);
  friend inline svec operator+(const svec& v1, const svec& v2);
  friend inline svec operator-(const svec& v1, const svec& v2);
  friend inline int operator==(const svec& v1, const svec& v2);
  friend inline int operator!=(const svec& v1, const svec& v2);
  friend int operator==(const svec& v1, const vec& v2);
  friend inline int operator!=(const svec& v1, const vec& v2) {return !(v1==v2);}
  friend inline int operator==(const vec& v1, const svec& v2) {return v2==v1;}
  friend inline int operator!=(const vec& v1, const svec& v2) {return v2!=v1;}
  friend smat transpose(const smat&);
  friend smat operator* ( const smat&, const smat&);
  friend smat mult_mod_p ( const smat&, const smat&, const scalar&);
  friend smat sidmat(scalar);
  friend scalar content(const svec& v);
  friend scalar make_primitive(svec& v); // divides by & returns content
  friend svec operator* ( const smat& A, const svec& v );
};

// Declaration of non-friend functions

inline svec operator+(const svec& v) {return v;}      // unary +
inline svec operator-(const svec& v)                  // unary -
{svec ans(v); ans*=-1; return ans;}
inline svec operator+(const svec& v1, const svec& v2)
{
  if(v1.entries.size()<v2.entries.size())   
    {
      svec ans(v2); ans+=v1; return ans;
    }
  else 
    {
      svec ans(v1); ans+=v2; return ans;
    }
}
inline svec operator-(const svec& v1, const svec& v2)
{
  if(v1.entries.size()<v2.entries.size())   
    {
      svec ans(v2); ans*=-1; ans+=v1; return ans;
    }
  else 
    {
      svec ans(v1); ans-=v2; return ans;
    }
}
inline svec operator*(scalar scal, const svec& v)
{svec ans(v); ans*=scal; return ans;}
inline svec operator/(const svec& v, scalar scal)
{svec ans(v); ans/=scal; return ans;}


inline int operator==(const svec& v1, const svec& v2)
{
  return (v1.d==v2.d) && (v1.entries == v2.entries); // STL built-in
}

inline int operator!=(const svec& v1, const svec& v2)
{
  return (v1.d!=v2.d) || !(v1.entries == v2.entries); // STL built-in
}
