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

///////////////////////////////////////////////////////////////////////////

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
