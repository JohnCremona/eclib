// svector.cc: implementation of class svec (sparse integer vectors)
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

#include "eclib/linalg.h"

// Definitions of member operators and functions:

template<class T>
sZvec<T>::sZvec(const Zvec<T>& v)
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
Zvec<T> sZvec<T>::as_vec( ) const
{
  Zvec<T> v(d); // initializes to 0
  for (const auto& wi : entries)
    v[wi.first] = wi.second;
  return v;
}

template<class T>
T sZvec<T>::elem(int i)  const   // returns i'th entry
{
  auto vi = entries.find(i);
  return (vi==entries.end()? T(0) : vi->second);
}

template<class T>
void sZvec<T>::set(int i, const T& a)
{
  entries.erase(i);
  if(is_nonzero(a))
    entries.insert({i,a});
}

template<class T>
void sZvec<T>::erase(int i)
{
  entries.erase(i);
}

template<class T>
std::set<int> sZvec<T>::support() const
{
  std::set<int> ans;
  for (const auto& vi : entries)
    ans.insert(vi.first);
  return ans;
}

template<class T>
void sZvec<T>::add(int i, const T& a)   // adds a to i'th entry
{
  if(is_zero(a)) return;
  auto vi = entries.find(i);
  if(vi==entries.end())
    entries[i]=a;
  else
    set(i, (vi->second)+a);
}

template<class T>
void sZvec<T>::sub(int i, const T& a)   // subtracts a from i'th entry
{
  if(is_zero(a)) return;
  add(i, -a);
}

template<class T>
void sZvec<T>::add_mod_p(int i, const T& aa, const T& p)
{
  T a=mod(aa,p);
  if(is_zero(a)) return;
  auto vi = entries.find(i);
  if(vi==entries.end())
    entries[i]=a;
  else
    set(i, mod((vi->second)+a, p));
}

template<class T>
void sZvec<T>::sub_mod_p(int i, const T& a, const T& p)
{
  add_mod_p(i, -a, p);
}

template<class T>
sZvec<T> operator+(const sZvec<T>& v)
{
  return v;
}
template sZvec<int> operator+<int>(const sZvec<int>&);
template sZvec<long> operator+<long>(const sZvec<long>&);
template sZvec<ZZ> operator+<ZZ>(const sZvec<ZZ>&);
template sZvec<INT> operator+<INT>(const sZvec<INT>&);

template<class T>
sZvec<T>& sZvec<T>::operator+=(const sZvec<T>& w)
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
sZvec<T> operator+(const sZvec<T>& v1, const sZvec<T>& v2)
{
  if(v1.entries.size()<v2.entries.size())
    {
      sZvec<T> ans(v2);
      ans+=v1;
      return ans;
    }
  else
    {
      sZvec<T> ans(v1);
      ans+=v2;
      return ans;
    }
}
template sZvec<int> operator+<int>(const sZvec<int>&, const sZvec<int>&);
template sZvec<long> operator+<long>(const sZvec<long>&, const sZvec<long>&);
template sZvec<ZZ> operator+<ZZ>(const sZvec<ZZ>&, const sZvec<ZZ>&);
template sZvec<INT> operator+<INT>(const sZvec<INT>&, const sZvec<INT>&);

template<class T>
sZvec<T> operator-(const sZvec<T>& v)
{
  sZvec<T> ans(v);
  ans*=T(-1);
  return ans;
}
template sZvec<int> operator-<int>(const sZvec<int>&);
template sZvec<long> operator-<long>(const sZvec<long>&);
template sZvec<ZZ> operator-<ZZ>(const sZvec<ZZ>&);
template sZvec<INT> operator-<INT>(const sZvec<INT>&);

template<class T>
sZvec<T>& sZvec<T>::operator-=(const sZvec<T>& w)
{
  this -> operator+=(-w);
  return *this;
}

template<class T>
sZvec<T> operator-(const sZvec<T>& v1, const sZvec<T>& v2)
{
  return v1 + (-v2);
}
template sZvec<int> operator-<int>(const sZvec<int>&, const sZvec<int>&);
template sZvec<long> operator-<long>(const sZvec<long>&, const sZvec<long>&);
template sZvec<ZZ> operator-<ZZ>(const sZvec<ZZ>&, const sZvec<ZZ>&);
template sZvec<INT> operator-<INT>(const sZvec<INT>&, const sZvec<INT>&);

template<class T>
sZvec<T>& sZvec<T>::add_scalar_times(const sZvec<T>& w, const T& a)
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
sZvec<T>& sZvec<T>::operator*=(const T& scal)
{
  for ( auto& vi : entries)
    (vi.second)*=scal;
  return *this;
}

template<class T>
sZvec<T> operator*(const T& scal, const sZvec<T>& v)
{
  sZvec<T> ans(v);
  ans*=scal;
  return ans;
}
template sZvec<int> operator*<int>(const int&, const sZvec<int>&);
template sZvec<long> operator*<long>(const long&, const sZvec<long>&);
template sZvec<ZZ> operator*<ZZ>(const ZZ&, const sZvec<ZZ>&);
template sZvec<INT> operator*<INT>(const INT&, const sZvec<INT>&);

template<class T>
void sZvec<T>::reduce_mod_p(const T& p)
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
sZvec<T>& sZvec<T>::mult_by_scalar_mod_p(const T& scal, const T& p)
{
  T s = mod(scal,p);
  if(!is_one(s))
    for( auto& vi : entries)
      (vi.second)=xmodmul(vi.second,s,p);
  return *this;
}

template<class T>
sZvec<T>& sZvec<T>::add_scalar_times_mod_p(const sZvec<T>& w, const T& scal, const T& p)
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
sZvec<T>& sZvec<T>::add_scalar_times_mod_p(const sZvec<T>& w, const T& scal,
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
sZvec<T>& sZvec<T>::add_mod_p(const sZvec<T>& w, const T& p)
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
sZvec<T>& sZvec<T>::sub_mod_p(const sZvec<T>& w, const T& p)
{
  add_mod_p(-w, p);
  return *this;
}

template<class T>
sZvec<T>& sZvec<T>::operator/=(const T& scal)
{
  if(is_zero(scal))
    cerr<<"Attempt to divide svec by 0"<<endl;
  for( auto& vi : entries)
    (vi.second)/=scal;
  return *this;
}

// Definitions of non-member, friend operators and functions

template<class T>
sZvec<T> operator/(const sZvec<T>& v, const T& scal)
{
  sZvec<T> ans(v);
  ans/=scal;
  return ans;
}
template sZvec<int> operator/<int>(const sZvec<int>&, const int&);
template sZvec<long> operator/<long>(const sZvec<long>&, const long&);
template sZvec<ZZ> operator/<ZZ>(const sZvec<ZZ>&, const ZZ&);
template sZvec<INT> operator/<INT>(const sZvec<INT>&, const INT&);

template<class T>
int eqmodp(const sZvec<T>& v1, const sZvec<T>& v2, const T& p)
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
template int eqmodp<int>(const sZvec<int>&, const sZvec<int>&, const int&);
template int eqmodp<long>(const sZvec<long>&, const sZvec<long>&, const long&);
template int eqmodp<ZZ>(const sZvec<ZZ>&, const sZvec<ZZ>&, const ZZ&);
template int eqmodp<INT>(const sZvec<INT>&, const sZvec<INT>&, const INT&);

template<class T>
int operator==(const sZvec<T>& v1, const sZvec<T>& v2)
{
  return (v1.d==v2.d) && (v1.entries == v2.entries);
}
template int operator==<int>(const sZvec<int>&, const sZvec<int>&);
template int operator==<long>(const sZvec<long>&, const sZvec<long>&);
template int operator==<ZZ>(const sZvec<ZZ>&, const sZvec<ZZ>&);
template int operator==<INT>(const sZvec<INT>&, const sZvec<INT>&);

template<class T>
int operator!=(const sZvec<T>& v1, const sZvec<T>& v2)
{
  return !(v1==v2);
}
template int operator!=<int>(const sZvec<int>&, const sZvec<int>&);
template int operator!=<long>(const sZvec<long>&, const sZvec<long>&);
template int operator!=<ZZ>(const sZvec<ZZ>&, const sZvec<ZZ>&);
template int operator!=<INT>(const sZvec<INT>&, const sZvec<INT>&);

template<class T>
int operator==(const sZvec<T>& v1, const Zvec<T>& v2)
{
  if(v1.d!=dim(v2)) return 0;
  for(int i=1; i<=v1.d; i++) if(v2[i]!=v1.elem(i)) return 0;
  return 1;
}
template int operator==<int>(const sZvec<int>&, const Zvec<int>&);
template int operator==<long>(const sZvec<long>&, const Zvec<long>&);
template int operator==<ZZ>(const sZvec<ZZ>&, const Zvec<ZZ>&);
template int operator==<INT>(const sZvec<INT>&, const Zvec<INT>&);

template<class T>
ostream& operator << (ostream& s, const sZvec<T>& v)
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
template ostream& operator<< <int>(ostream&, const sZvec<int>&);
template ostream& operator<< <long>(ostream&, const sZvec<long>&);
template ostream& operator<< <ZZ>(ostream&, const sZvec<ZZ>&);
template ostream& operator<< <INT>(ostream&, const sZvec<INT>&);

template<class T>
T operator*(const sZvec<T>& v, const sZvec<T>& w) //dot prod
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
template int operator*<int>(const sZvec<int>&, const sZvec<int>&);
template long operator*<long>(const sZvec<long>&, const sZvec<long>&);
template ZZ operator*<ZZ>(const sZvec<ZZ>&, const sZvec<ZZ>&);
template INT operator*<INT>(const sZvec<INT>&, const sZvec<INT>&);

template<class T>
T operator*(const sZvec<T>& v, const Zvec<T>& w) //dot prod
{
  T ans(0);
  std::for_each(v.entries.cbegin(), v.entries.cend(),
                [&ans, w](auto vi){ans += (vi.second)* w[vi.first];});
  return ans;
}
template int operator*<int>(const sZvec<int>&, const Zvec<int>&);
template long operator*<long>(const sZvec<long>&, const Zvec<long>&);
template ZZ operator*<ZZ>(const sZvec<ZZ>&, const Zvec<ZZ>&);
template INT operator*<INT>(const sZvec<INT>&, const Zvec<INT>&);

template<class T> T operator*(const Zvec<T>& v, const sZvec<T>& sv)
{
  return sv*v;
}
template int operator*<int>(const Zvec<int>&, const sZvec<int>&);
template long operator*<long>(const Zvec<long>&, const sZvec<long>&);
template ZZ operator*<ZZ>(const Zvec<ZZ>&, const sZvec<ZZ>&);
template INT operator*<INT>(const Zvec<INT>&, const sZvec<INT>&);

template<class T>
int trivial(const sZvec<T>& v)
{
  return v.entries.size()==0;
}
template int trivial<int>(const sZvec<int>&);
template int trivial<long>(const sZvec<long>&);
template int trivial<ZZ>(const sZvec<ZZ>&);
template int trivial<INT>(const sZvec<INT>&);

template<class T>
T dotmodp(const sZvec<T>& v, const Zvec<T>& w, const T& pr)
{
  T ans(0);
  std::for_each(v.entries.cbegin(), v.entries.cend(),
                [&ans, w, pr](auto vi){ans=mod(ans+xmodmul(vi.second,w[vi.first],pr),pr);});
  return ans;
}
template int dotmodp<int>(const sZvec<int>&, const Zvec<int>&, const int&);
template long dotmodp<long>(const sZvec<long>&, const Zvec<long>&, const long&);
template ZZ dotmodp<ZZ>(const sZvec<ZZ>&, const Zvec<ZZ>&, const ZZ&);
template INT dotmodp<INT>(const sZvec<INT>&, const Zvec<INT>&, const INT&);

template<class T>
T dotmodp(const sZvec<T>& v, const sZvec<T>& w, const T& pr)
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
template int dotmodp<int>(const sZvec<int>&, const sZvec<int>&, const int&);
template long dotmodp<long>(const sZvec<long>&, const sZvec<long>&, const long&);
template ZZ dotmodp<ZZ>(const sZvec<ZZ>&, const sZvec<ZZ>&, const ZZ&);
template INT dotmodp<INT>(const sZvec<INT>&, const sZvec<INT>&, const INT&);

template<class T> T dotmodp(const Zvec<T>& v, const sZvec<T>& sv, const T& pr)
{
  return dotmodp(sv,v,pr);
}
template int dotmodp<int>(const Zvec<int>&, const sZvec<int>&, const int&);
template long dotmodp<long>(const Zvec<long>&, const sZvec<long>&, const long&);
template ZZ dotmodp<ZZ>(const Zvec<ZZ>&, const sZvec<ZZ>&, const ZZ&);
template INT dotmodp<INT>(const Zvec<INT>&, const sZvec<INT>&, const INT&);

template<class T>
T content(const sZvec<T>& v)
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
template int content<int>(const sZvec<int>&);
template long content<long>(const sZvec<long>&);
template ZZ content<ZZ>(const sZvec<ZZ>&);
template INT content<INT>(const sZvec<INT>&);

template<class T>
T make_primitive(sZvec<T>& v) // divides by & returns content
{
  T c=content(v);
  v/=c;
  return c;
}
template int make_primitive<int>(sZvec<int>&);
template long make_primitive<long>(sZvec<long>&);
template ZZ make_primitive<ZZ>(sZvec<ZZ>&);
template INT make_primitive<INT>(sZvec<INT>&);

template<class T> int dim(const sZvec<T>& v)
{
  return v.d;
}
template int dim<int>(const sZvec<int>&);
template int dim<long>(const sZvec<long>&);
template int dim<ZZ>(const sZvec<ZZ>&);
template int dim<INT>(const sZvec<INT>&);

template<class T> int operator!=(const sZvec<T>& v1, const Zvec<T>& v2)
{
  return !(v1==v2);
}
template int operator!=<int>(const sZvec<int>&, const Zvec<int>&);
template int operator!=<long>(const sZvec<long>&, const Zvec<long>&);
template int operator!=<ZZ>(const sZvec<ZZ>&, const Zvec<ZZ>&);
template int operator!=<INT>(const sZvec<INT>&, const Zvec<INT>&);

template<class T> int operator==(const Zvec<T>& v1, const sZvec<T>& v2)
{
  return v2==v1;
}
template int operator==<int>(const Zvec<int>&, const sZvec<int>&);
template int operator==<long>(const Zvec<long>&, const sZvec<long>&);
template int operator==<ZZ>(const Zvec<ZZ>&, const sZvec<ZZ>&);
template int operator==<INT>(const Zvec<INT>&, const sZvec<INT>&);

template<class T> int operator!=(const Zvec<T>& v1, const sZvec<T>& v2)
{
  return v2!=v1;
}
template int operator!=<int>(const Zvec<int>& v1, const sZvec<int>& v2);
template int operator!=<long>(const Zvec<long>& v1, const sZvec<long>& v2);
template int operator!=<ZZ>(const Zvec<ZZ>& v1, const sZvec<ZZ>& v2);
template int operator!=<INT>(const Zvec<INT>& v1, const sZvec<INT>& v2);

// Instantiate sZvec template classes for T=int, long, ZZ, INT

template class sZvec<int>;
template class sZvec<long>;
template class sZvec<ZZ>;
template class sZvec<INT>;


