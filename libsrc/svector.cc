// svector.cc: implementation of class svec (sparse integer vectors)
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

#include <eclib/svector.h>
#include <eclib/marith.h>

#if(0)

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i
#define svec svec_i

#include "svec.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l
#define svec svec_l

#include "svec.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m
#define svec svec_m
#define smat smat_m
#define smat_elim smat_m_elim

#include "svec.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

#endif

///////////////////////////////////////////////////////////////////////////

// Instantiate svecT template classes for T=int, long, bigint

template class svecT<int>;
template class svecT<long>;
template class svecT<bigint>;

// Definitions of member operators and functions:

template<class T>
svecT<T>::svecT(const vecT<T>& v)
  :d(dim(v))
{
  for(int i = 1; i <= d; i++ )
    {
      T vi = v[i];
      if(is_nonzero(vi))
        entries[i]=vi;
    }
}

template<class T>
vecT<T> svecT<T>::as_vec( ) const
{
  vecT<T> v(d); // initializes to 0
  for (const auto& wi : entries)
    v[wi.first] = wi.second;
  return v;
}

template<class T>
T svecT<T>::elem(int i)  const   // returns i'th entry
{
  auto vi = entries.find(i);
  if(vi==entries.end()) return T(0);
  return vi->second;
}

template<class T>
void svecT<T>::set(int i, const T& a)
{
  if(is_nonzero(a)) entries[i]=a;
}

template<class T>
void svecT<T>::erase(int i)
{
  auto vi = entries.find(i);
  if(vi==entries.end())
    {
      cerr<<"Error in svec::erase(): cannot delete missing entry #"<<i
	  <<" from v = "<<(*this)<<endl;
    }
  else entries.erase(vi);
}

template<class T>
std::set<int> svecT<T>::support() const
{
  std::set<int> ans;
  for (const auto& vi : entries)
    ans.insert(vi.first);
  return ans;
}

template<class T>
void svecT<T>::add(int i, const T& a)   // adds a to i'th entry
{
  if(is_zero(a)) return;
  auto  vi = entries.find(i);
  if(vi==entries.end())
    entries[i]=a;
  else
    {
      T sum = (vi->second)+a;
      if(is_zero(sum)) entries.erase(vi);
      else (vi->second)=sum;
    }
}

template<class T>
void svecT<T>::sub(int i, const T& a)   // subtracts a from i'th entry
{
  if(is_zero(a)) return;
  add(i, -a);
}

template<class T>
void svecT<T>::add_mod_p(int i, const T& aa, const T& p)
{
  T a=mod(aa,p);
  if(is_zero(a)) return;
  auto vi = entries.find(i);
  if(vi==entries.end())
    entries[i]=a;
  else
    {
      T sum = mod((vi->second)+a,p);
      if(is_zero(sum)) entries.erase(vi);
      else (vi->second)=sum;
    }
}

template<class T>
void svecT<T>::sub_mod_p(int i, const T& a, const T& p)
{
  add_mod_p(i, -a, p);
}

template<class T>
svecT<T>& svecT<T>::operator+=(const svecT<T>& w)
{
  if (d!=w.d)
    {
      cerr << "Incompatible svecs in svec::operator+=()"<<endl;
      return *this;
    }
  auto wi=w.entries.begin();
  auto vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=wi->second;
	      wi++;
	    }
	}
      else
	{
	  if((vi->first)<(wi->first)) {vi++;}
	  else
	    if((wi->first)<(vi->first))
	      {
		entries[wi->first]=wi->second;
		wi++;
	      }
	    else
	      {
		T sum = (vi->second) + (wi->second);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first);}
		wi++;
	      }
	}
    }
  return *this;
}

template<class T>
svecT<T>& svecT<T>::operator-=(const svecT<T>& w)
{
  this -> operator+=(-w);
  return *this;
}

template<class T>
svecT<T>& svecT<T>::add_scalar_times(const svecT<T>& w, const T& a)
{
  if (d!=w.d)
    {
      cerr << "Incompatible svecs in svec::add_scalar_times()"<<endl;
      return *this;
    }
  if(is_zero(a)) return *this;
  auto wi=w.entries.begin();
  auto vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=a*(wi->second);
	      wi++;
	    }
	}
      else
	{
	  if((vi->first)<(wi->first)) {vi++;}
	  else
	    if((wi->first)<(vi->first))
	      {
		entries[wi->first]=a*(wi->second);
		wi++;
	      }
	    else
	      {
		T sum = (vi->second) + a* (wi->second);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first);}
		wi++;
	      }
	}
    }
  return *this;
}

template<class T>
svecT<T>& svecT<T>::operator*=(const T& scal)
{
  for ( auto& vi : entries)
    (vi.second)*=scal;
  return *this;
}

template<class T>
void svecT<T>::reduce_mod_p(const T& p)
{
  auto vi = entries.begin();
  while( vi != entries.end() )
    {
      T a = mod(vi->second,p);
      if(is_nonzero(a))
        {
          (vi->second)=a;
          vi++;
        }
      else
        {
          vi = entries.erase(vi);
        }
    }
}

template<class T>
svecT<T>& svecT<T>::mult_by_scalar_mod_p(const T& scal, const T& p)
{
  T s = mod(scal,p);
  if(!is_one(s))
    for( auto& vi : entries)
      (vi.second)=xmodmul(vi.second,s,p);
  return *this;
}

template<class T>
svecT<T>& svecT<T>::add_scalar_times_mod_p(const svecT<T>& w, const T& scal, const T& p)
{
  if (d!=w.d)
    {
      cerr << "Incompatible svecs in svec::add_scalar_times()"<<endl;
      return *this;
    }
  T a = mod(scal,p);
  if(is_zero(a)) return *this;
  auto wi=w.entries.begin();
  auto vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=xmodmul(a,(wi->second),p);
	      wi++;
	    }
	}
      else
	{
	  if((vi->first)<(wi->first)) {vi++;}
	  else
	    if((wi->first)<(vi->first))
	      {
		entries[wi->first]=xmodmul(a,(wi->second),p);
		wi++;
	      }
	    else
	      {
		T sum = xmod((vi->second) + xmodmul(a, (wi->second),p),p);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); }
		wi++;
	      }
	}
    }
  //  reduce_mod_p(p);
  return *this;
}

template<class T>
svecT<T>& svecT<T>::add_scalar_times_mod_p(const svecT<T>& w, const T& scal,
                                           std::set<int>& ons, std::set<int>& offs,
                                           const T& p)
{
  ons.clear();
  offs.clear();
  if (d!=w.d)
    {
      cerr << "Incompatible svecs in svec::add_scalar_times()"<<endl;
      return *this;
    }
  T a = mod(scal,p);
  if(is_zero(a)) return *this;
  auto wi=w.entries.begin();
  auto vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=xmodmul(a,(wi->second),p);
	      ons.insert(wi->first);
	      wi++;
	    }
	}
      else
	{
	  if((vi->first)<(wi->first)) {vi++;}
	  else
	    if((wi->first)<(vi->first))
	      {
		entries[wi->first]=xmodmul(a,(wi->second),p);
		ons.insert(wi->first);
		wi++;
	      }
	    else
	      {
		T sum = xmod((vi->second) + xmodmul(a, (wi->second),p),p);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); offs.insert(wi->first);}
		wi++;
	      }
	}
    }
  //  reduce_mod_p(p);
  return *this;
}

template<class T>
svecT<T>& svecT<T>::add_mod_p(const svecT<T>& w, const T& p)
{
  if (d!=w.d)
    {
      cerr << "Incompatible svecs in svec::add_scalar_times()"<<endl;
      return *this;
    }
  auto wi=w.entries.begin();
  auto vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=wi->second;
	      wi++;
	    }
	}
      else
	{
	  if((vi->first)<(wi->first)) {vi++;} 
	  else
	    if((wi->first)<(vi->first)) 
	      {
		entries[wi->first]=wi->second;
		wi++;
	      } 
	    else
	      {
		T sum = xmod((vi->second) + (wi->second),p);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); }
		wi++;
	      }
	}
    }
  return *this;
}

template<class T>
svecT<T>& svecT<T>::sub_mod_p(const svecT<T>& w, const T& p)
{
  add_mod_p(-w, p);
  return *this;
}

template<class T>
svecT<T>& svecT<T>::operator/=(const T& scal)
{
  if(is_zero(scal))
    cerr<<"Attempt to divide svec by 0"<<endl;
  for( auto& vi : entries)
    (vi.second)/=scal;
  return *this;
}

// Definitions of non-member, friend operators and functions

template<class T>
int eqmodp(const svecT<T>& v1, const svecT<T>& v2, const T& p)
{
  if(v1.d!=v2.d) return 0;
  if (std::any_of(v1.entries.begin(), v1.entries.end(),
                  [v2,p] (const pair<int,T>& vi) {return xmod((vi.second)-(v2.elem(vi.first)),p)!=0;}))
    return 0;
  if (std::any_of(v2.entries.begin(), v2.entries.end(),
                  [v1,p] (const pair<int,T>& vi) {return xmod((vi.second)-(v1.elem(vi.first)),p)!=0;}))
    return 0;
  return 1;
}

template<class T>
int operator==(const svecT<T>& v1, const vecT<T>& v2)
{
  if(v1.d!=dim(v2)) return 0;
  for(int i=1; i<=v1.d; i++) if(v2[i]!=v1.elem(i)) return 0;
  return 1;
}

template<class T>
ostream& operator << (ostream& s, const svecT<T>& v)
{
  s<<"[";
  for(auto vi=v.entries.begin(); vi!=v.entries.end(); vi++)
    {
      if(vi!=v.entries.begin()) s<<",";
      s << "("<<vi->first<<":"<<vi->second<<")";
    }
  s<<"]";
  return s;
}

template<class T>
T operator*(const svecT<T>& v, const svecT<T>& w) //dot prod
{
  T ans(0);
  if((v.entries.size()==0)||(w.entries.size()==0)) return ans;
  auto vi=v.entries.begin(), wi=w.entries.begin();
  while((vi!=v.entries.end())&&(wi!=w.entries.end()))
    {
      if((vi->first)<(wi->first)) {vi++;} else
	if((wi->first)<(vi->first)) {wi++;} else
	  {
	    ans+=(vi->second)*(wi->second);
	    vi++; wi++;
	  }
    }
  return ans;
}

template<class T>
T operator*(const svecT<T>& v, const vecT<T>& w) //dot prod
{
  T ans(0);
  std::for_each(v.entries.cbegin(), v.entries.cend(),
                [&ans, w](auto vi){ans += (vi.second)* w[vi.first];});
  return ans;
}

template<class T>
T dotmodp(const svecT<T>& v, const vecT<T>& w, const T& pr)
{
  T ans(0);
  std::for_each(v.entries.cbegin(), v.entries.cend(),
                [&ans, w, pr](auto vi){ans=mod(ans+xmodmul(vi.second,w[vi.first],pr),pr);});
  return ans;
}

template<class T>
T dotmodp(const svecT<T>& v, const svecT<T>& w, const T& pr)
{
  T ans(0);
  if((v.entries.size()==0)||(w.entries.size()==0)) return ans;
  auto vi=v.entries.begin(), wi=w.entries.begin();
  while((vi!=v.entries.end())&&(wi!=w.entries.end()))
    {
      if((vi->first)<(wi->first)) {vi++;} else
	if((wi->first)<(vi->first)) {wi++;} else
	  {
	    ans=mod(ans+xmodmul(vi->second,wi->second,pr),pr);
	    vi++; wi++;
	  }
    }
  return ans;
}

template<class T>
T content(const svecT<T>& v)
{
  T g(0);
  for( const auto & vi : v.entries)
    {
      g=gcd(g,vi.second);
      if (is_one(g))
        return g;
    }
  return g;
}

template<class T>
T make_primitive(svecT<T>& v) // divides by & returns content
{
  T c=content(v);
  v/=c;
  return c;
}

template<class T> int dim(const svecT<T>& v)   {return v.d;}
template<class T> T operator*(const vecT<T>& v, const svecT<T>& sv)  {return sv*v;}
template<class T> T dotmodp(const vecT<T>& v, const svecT<T>& sv, const T& pr) {return dotmodp(sv,v,pr);}
template<class T> int operator!=(const svecT<T>& v1, const vecT<T>& v2) {return !(v1==v2);}
template<class T> int operator==(const vecT<T>& v1, const svecT<T>& v2) {return v2==v1;}
template<class T> int operator!=(const vecT<T>& v1, const svecT<T>& v2) {return v2!=v1;};

// Instantiate svecT template functions for T=int
template int operator*<int>(const svecT<int>&, const svecT<int>&);
template int operator*<int>(const svecT<int>&, const vecT<int>&);
template int content<int>(const svecT<int>& v);
template int make_primitive<int>(svecT<int>& v);
template int dotmodp<int>(const svecT<int>& v, const vecT<int>& w, const int& pr);
template int dotmodp<int>(const svecT<int>& v, const svecT<int>& w, const int& pr);
template int dim<int>(const svecT<int>& v);
template int eqmodp<int>(const svecT<int>&, const svecT<int>&, const int& p);
template ostream& operator<< <int>(ostream&s, const svecT<int>&);
template int operator*<int>(const vecT<int>& v, const svecT<int>& sv);
template int dotmodp<int>(const vecT<int>& v, const svecT<int>& sv, const int& pr);
template svecT<int> operator+<int>(const svecT<int>& v1, const svecT<int>& v2);
template svecT<int> operator-<int>(const svecT<int>& v1, const svecT<int>& v2);
template int operator==<int>(const svecT<int>& v1, const svecT<int>& v2);
template int operator!=<int>(const svecT<int>& v1, const svecT<int>& v2);
template int operator==<int>(const svecT<int>& v1, const vecT<int>& v2);
template int operator!=<int>(const svecT<int>& v1, const vecT<int>& v2);
template int operator==<int>(const vecT<int>& v1, const svecT<int>& v2);
template int operator!=<int>(const vecT<int>& v1, const svecT<int>& v2);
template svecT<int> operator+<int>(const svecT<int>& v);
template svecT<int> operator-<int>(const svecT<int>& v);
// template svecT<int> operator+<int>(const svecT<int>& v1, const svecT<int>& v2);
// template svecT<int> operator-<int>(const svecT<int>& v1, const svecT<int>& v2);
template svecT<int> operator*<int>(const int& scal, const svecT<int>& v);
template svecT<int> operator/<int>(const svecT<int>& v, const int& scal);
// template int operator==<int>(const svecT<int>& v1, const svecT<int>& v2);
// template int operator!=<int>(const svecT<int>& v1, const svecT<int>& v2);

// Instantiate svecT template functions for T=long
template long operator*<long>(const svecT<long>&, const svecT<long>&);
template long operator*<long>(const svecT<long>&, const vecT<long>&);
template long content<long>(const svecT<long>& v);
template long make_primitive<long>(svecT<long>& v);
template long dotmodp<long>(const svecT<long>& v, const vecT<long>& w, const long& pr);
template long dotmodp<long>(const svecT<long>& v, const svecT<long>& w, const long& pr);
template int dim<long>(const svecT<long>& v);
template int eqmodp<long>(const svecT<long>&, const svecT<long>&, const long& p);
template ostream& operator<< <long>(ostream&s, const svecT<long>&);
template long operator*<long>(const vecT<long>& v, const svecT<long>& sv);
template long dotmodp<long>(const vecT<long>& v, const svecT<long>& sv, const long& pr);
template svecT<long> operator+<long>(const svecT<long>& v1, const svecT<long>& v2);
template svecT<long> operator-<long>(const svecT<long>& v1, const svecT<long>& v2);
template int operator==<long>(const svecT<long>& v1, const svecT<long>& v2);
template int operator!=<long>(const svecT<long>& v1, const svecT<long>& v2);
template int operator==<long>(const svecT<long>& v1, const vecT<long>& v2);
template int operator!=<long>(const svecT<long>& v1, const vecT<long>& v2);
template int operator==<long>(const vecT<long>& v1, const svecT<long>& v2);
template int operator!=<long>(const vecT<long>& v1, const svecT<long>& v2);
template svecT<long> operator+<long>(const svecT<long>& v);
template svecT<long> operator-<long>(const svecT<long>& v);
// template svecT<long> operator+<long>(const svecT<long>& v1, const svecT<long>& v2);
// template svecT<long> operator-<long>(const svecT<long>& v1, const svecT<long>& v2);
template svecT<long> operator*<long>(const long& scal, const svecT<long>& v);
template svecT<long> operator/<long>(const svecT<long>& v, const long& scal);
// template int operator==<long>(const svecT<long>& v1, const svecT<long>& v2);
// template int operator!=<long>(const svecT<long>& v1, const svecT<long>& v2);

// Instantiate svecT template functions for T=long
template bigint operator*<bigint>(const svecT<bigint>&, const svecT<bigint>&);
template bigint operator*<bigint>(const svecT<bigint>&, const vecT<bigint>&);
template bigint content<bigint>(const svecT<bigint>& v);
template bigint make_primitive<bigint>(svecT<bigint>& v);
template bigint dotmodp<bigint>(const svecT<bigint>& v, const vecT<bigint>& w, const bigint& pr);
template bigint dotmodp<bigint>(const svecT<bigint>& v, const svecT<bigint>& w, const bigint& pr);
template int dim<bigint>(const svecT<bigint>& v);
template int eqmodp<bigint>(const svecT<bigint>&, const svecT<bigint>&, const bigint& p);
template ostream& operator<< <bigint>(ostream&s, const svecT<bigint>&);
template bigint operator*<bigint>(const vecT<bigint>& v, const svecT<bigint>& sv);
template bigint dotmodp<bigint>(const vecT<bigint>& v, const svecT<bigint>& sv, const bigint& pr);
template svecT<bigint> operator+<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
template svecT<bigint> operator-<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
template int operator==<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
template int operator!=<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
template int operator==<bigint>(const svecT<bigint>& v1, const vecT<bigint>& v2);
template int operator!=<bigint>(const svecT<bigint>& v1, const vecT<bigint>& v2);
template int operator==<bigint>(const vecT<bigint>& v1, const svecT<bigint>& v2);
template int operator!=<bigint>(const vecT<bigint>& v1, const svecT<bigint>& v2);
template svecT<bigint> operator+<bigint>(const svecT<bigint>& v);
template svecT<bigint> operator-<bigint>(const svecT<bigint>& v);
// template svecT<bigint> operator+<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
// template svecT<bigint> operator-<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
template svecT<bigint> operator*<bigint>(const bigint& scal, const svecT<bigint>& v);
template svecT<bigint> operator/<bigint>(const svecT<bigint>& v, const bigint& scal);
// template int operator==<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
// template int operator!=<bigint>(const svecT<bigint>& v1, const svecT<bigint>& v2);
// template smatT<bigint> transpose<bigint>(const smatT<bigint>&);
// template smatT<bigint> operator* <bigint>( const smatT<bigint>&, const smatT<bigint>&);
// template svecT<bigint> operator* <bigint>( const smatT<bigint>& A, const svecT<bigint>& v );
// template svecT<bigint> operator* <bigint>( const svecT<bigint>& v, const smatT<bigint>& A );
// template svecT<bigint> mult_mod_p<bigint>( const smatT<bigint>& A, const svecT<bigint>& v, const bigint& p  );
// template svecT<bigint> mult_mod_p<bigint>( const svecT<bigint>& v, const smatT<bigint>& A, const bigint& p  );
// template smatT<bigint> mult_mod_p<bigint>( const smatT<bigint>&, const smatT<bigint>&, const bigint&);
