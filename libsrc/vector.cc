// vector.cc: manage implementations of integer vector classes
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
 
#include <eclib/marith.h>
#include <eclib/vector.h>
#include "random.cc"

// Instantiate vecT template classes for T=int, long, bigint

template class vecT<int>;
template class vecT<long>;
template class vecT<bigint>;

// Definitions of member operators and functions:

template<class T>
vecT<T>::vecT(long n)
{
  entries.resize(n, T(0));
}

template<class T>
vecT<T>::vecT(const vector<T>& arr) :entries(arr) {}

template<class T>
vecT<T>::vecT(const vecT<T>& v) :entries(v.entries) {} // copy constructor

template<class T>
void vecT<T>::init(long n)     // (re)-initializes
{
  entries.resize(n, T(0));
}

template<class T>
vecT<T>& vecT<T>::operator=(const vecT<T>& v)    // assignment
{
 if (this==&v) return *this;
 entries = v.entries;
 return *this;
}

template<class T>
T& vecT<T>::operator[](long i)
{
  return entries.at(i-1);
}

template<class T>
T vecT<T>::operator[](long i) const
{
  return entries.at(i-1);
}

template<class T>
vecT<T>& vecT<T>::operator+=(const vecT<T>& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const T& wi, const T& vi) { return vi + wi;});
  return *this;
}

template<class T>
void vecT<T>::addmodp(const vecT<T>& w, const T& pr)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [pr](const T& wi, const T& vi) { return mod(wi+vi,pr);});
}

template<class T>
vecT<T>& vecT<T>::operator-=(const vecT<T>& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const T& wi, const T& vi) { return vi - wi;});
  return *this;
}

template<class T>
vecT<T>& vecT<T>::operator*=(const T& scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const T& vi) {return vi * scal;});
  return *this;
}

template<class T>
vecT<T>& vecT<T>::operator/=(const T& scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const T& vi) {return vi / scal;});
  return *this;
}

template<class T>
vecT<T> vecT<T>::slice(long first, long last) const       // returns subvector
{
 if (last==-1) {last=first; first=1;}
 vecT<T> ans(last-first+1);
 std::copy(entries.begin()+first-1, entries.begin()+last, ans.entries.begin());
 return ans;
}

template<class T>
vecT<T> vecT<T>::operator[](const vecT<int>& index) const  // returns v[index[j]]
{
  vecT<T> w(dim(index));
  const vector<int>& vi = index.get_entries();
  std::transform(vi.begin(), vi.end(), w.entries.begin(),
                 [this](const int& i) {return entries.at(i-1);});
  return w;
}

template<class T>
vecT<T> vecT<T>::operator[](const vecT<long>& index) const  // returns v[index[j]]
{
  vecT<T> w(dim(index));
  const vector<long>& vi = index.get_entries();
  std::transform(vi.begin(), vi.end(), w.entries.begin(),
                 [this](const int& i) {return entries.at(i-1);});
  return w;
}

template<class T>
T vecT<T>::sub(long i) const
{
  return entries.at(i-1);
}

template<class T>
void vecT<T>::set(long i, const T& x)
{
  entries.at(i-1) = x;
}

template<class T>
void vecT<T>::add(long i, const T& x)
{
  entries.at(i-1) += x;
}

template<class T>
void vecT<T>::add_modp(long i, const T& x, const T& p)
{
  entries.at(i-1) = mod(entries.at(i-1)+x,p);
}

template<class T>
void vecT<T>::reduce_mod_p(const T& p)
{
  if (p==0) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [p](const T& vi) {return mod(vi,p);});
}

template<class T>
vecT<T> vecT<T>::iota(long n)
{
  vecT<T> v(n);
  std::iota(v.entries.begin(), v.entries.end(), T(1));
  return v;
}

// Definitions of non-member, friend operators and functions

template<class T>
T operator*(const vecT<T>& v, const vecT<T>& w)
{
  return std::inner_product(v.entries.begin(), v.entries.end(), w.entries.begin(), T(0));
}

template<class T>
int operator==(const vecT<T>& v, const vecT<T>& w)
{
  return v.entries == w.entries;
}

template<class T>
int trivial(const vecT<T>& v)
{
  return std::all_of(v.entries.begin(), v.entries.end(), [](const T& vi) {return vi==0;});
}

template<class T>
ostream& operator<<(ostream& s, const vecT<T>& v)
{
  s << "[";
  long i=0;
  for ( const auto& vi : v.entries)
    {
      if(i++)
        s<<",";
      s<<vi;
    }
  s << "]";
  return s;
}

template<class T>
istream& operator>>(istream& s, vecT<T>& v)
{
  for (T& vi : v.entries)
    s>>vi;
  return s;
}

// Definition of non-friend operators and functions

template<class T>
T content(const vecT<T>& v)
{
  return v.entries.empty()?
    T(1) :
    std::accumulate(v.entries.begin(), v.entries.end(), T(0),
                    [](const T& x, const T& y) {return gcd(x,y);});
}

template<class T>
T maxabs(const vecT<T>& v)
{
  return v.entries.empty()?
    T(0) :
    std::accumulate(v.entries.begin(), v.entries.end(), T(0),
                    [](const T& x, const T& y) {return max(x,abs(y));});
}

template<class T>
void swapvec(vecT<T>& v, vecT<T>& w)
{
  std::swap(v.entries, w.entries);
}

template<class T>
int member(const T& a, const vecT<T>& v)
{
  return std::find(v.entries.begin(), v.entries.end(), a) != v.entries.end();
}

template<class T>
vecT<T> reverse(const vecT<T>& order)
{
  vecT<T> ans(order);
  std::reverse(ans.entries.begin(), ans.entries.end());
  return ans;
}

template<class T>
vecT<T> express(const vecT<T>& v, const vecT<T>& v1, const vecT<T>& v2)
{
   T v1v1 = v1 * v1;
   T v1v2 = v1 * v2;
   T v2v2 = v2 * v2;
   T vv1 = v * v1;
   T vv2 = v * v2;
   vecT<T> ans({vv1*v2v2 - vv2*v1v2,  vv2*v1v1 - vv1*v1v2, v1v1*v2v2 - v1v2*v1v2});
   make_primitive(ans);
   if (ans[3]*v!=ans[1]*v1+ans[2]*v2)
     cerr << "Error in express: v is not in <v1,v2>"<<endl;
   return ans;
}

//#define DEBUG_LIFT

// int lift(const vecT<T>& v, const T& pr, vecT<T>& w)
// {
//   w = v;
//   w.reduce_mod_p(pr);
// }

template<class T>
int lift(const vecT<T>& v, const T& pr, vecT<T>& ans)
{
  long i0, i, j, d = dim(v);
  T nu, de;
  T lim = sqrt(pr>>1)-1;
  T maxallowed = 10*lim;
#ifdef DEBUG_LIFT
  cout<<"Lifting vector v = "<<v<<" mod "<<pr<<" (lim = "<<lim<<")"<<endl;
#endif
 // NB We do *not* make cumulative rescalings, since it is possible
 // for an apparently successful modrat reconstruction to give an
 // incorrect denominator.  I have an example with pr=2^30-35 where
 // the correct denominator is 4666 and one entry of the correct
 // primitive scaled vector is 47493 (greater than lim = 23170) but
 // since 47493/4666 = 587037152 = -10193/21607 (mod pr), rational
 // reconstruction returned nu=-10193, de = 21607.  If we kept the
 // (unsuccessful) scaling by 21607, all subsequent numerators would
 // be multiplied by this and we would never succeed.

 // This code allows for some entries to be >lim, and works as long as
 // (1) there is a lift with all entries at most 10*lim, (2) at least
 // one entry has the correct denominator, which is equaivalent to
 // requiring that in the primitive rescaling, there is an entry
 // coprime to the first non-zero entry.

 ans = reduce_mod_p(v, pr); // starts as a copy, and will be rescaled in place
#ifdef DEBUG_LIFT
  cout<<"After reduce_mod_p: v = "<<ans<<endl;
#endif
 if (maxabs(ans) <= maxallowed)
   {
#ifdef DEBUG_LIFT
     cout<<"No scaling needed, lift is "<<ans<<endl;
#endif
     return 1;
   }
 T vi0, inv_vi0, vi, maxvi(0);
 for(i0=1; i0<=d; i0++)
   {
     // scale so that i0'th entry is 1 mod p, then reduce vector
     // entries mod p to lie in (-p/2,p/2), and find the maximum
     // entry:
     while((vi0=ans[i0])==0) {i0++;} // skip over any zero entries
     inv_vi0=invmod(vi0,pr);
#ifdef DEBUG_LIFT
     cout<<"Scaling by "<<inv_vi0<<" (inverse of "<<vi0<<")"<<endl;
#endif
     for (i=1; i<=d; i++)
       {
         ans[i]=vi=mod(xmodmul(inv_vi0,ans[i],pr),pr);
         maxvi=max(maxvi,abs(vi));
       }
#ifdef DEBUG_LIFT
     cout<<"Reduced v = "<<ans<<", with max entry "<<maxvi<<endl;
#endif
     if(maxvi<=maxallowed) // no scaling needed!
           {
             // Normalize so first nonzero entry is positive:
             for(i0=1; i0<=d; i0++)
               {
                 while(ans[i0]==0) {i0++;}
                 if(ans[i0]<0) ans=-ans;
                 return 1;
               }
             return 0; // should not happen: means v==0!
           }

     for(i=1; (i<=d); i++)
       {
         modrat(ans[i],pr,nu,de);
         de=abs(de);
         if (de==1) continue; // loop on i
         // scale by de & recompute max entry:
#ifdef DEBUG_LIFT
         cout<<"Scaling by d="<<de<<endl;
#endif
         maxvi = 0;
         for (j=1; j<=d; j++)
           {
             ans[j] = vi = mod(xmodmul(de,ans[j],pr),pr);
             maxvi=max(maxvi,abs(vi));
           }
#ifdef DEBUG_LIFT
         cout<<"Now v = "<<ans<<", with max entry "<<maxvi<<endl;
#endif
         if(maxvi<=maxallowed)
           {
             // Normalize so first nonzero entry is positive:
             for(i0=1; i0<=d; i0++)
               {
                 while(ans[i0]==0) {i0++;}
                 if(ans[i0]<0) ans=-ans;
                 return 1;
               }
             return 0; // should not happen: means v==0!
           }
       }
   }
 // Normalize so first nonzero entry is positive:
 for(i0=1; i0<=d; i0++)
   {
     while(ans[i0]==0) {i0++;}
     if(ans[i0]<0) ans=-ans;
     return (maxvi<=lim);
   }
 return 0;
}

template<class T>
T dotmodp(const vecT<T>& v1, const vecT<T>& v2, const T& pr)
{
  auto a = [pr] (const T& x, const T& y) {return mod(x+y,pr);};
  auto m = [pr] (const T& x, const T& y) {return xmodmul(x,y,pr);};
  return std::inner_product(v1.entries.begin(), v1.entries.end(), v2.entries.begin(), T(0), a, m);
}

vec_m to_vec_m(const vec_i& v)
{
  const vector<int> & vi = v.get_entries();
  vector<bigint> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), [](const int& x) {return bigint(x);});
  return vec_m(w);
}

vec_m to_vec_m(const vec_l& v)
{
  const vector<long> & vi = v.get_entries();
  vector<bigint> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), [](const long& x) {return bigint(x);});
  return vec_m(w);
}

vec_i to_vec_i(const vec_m& v)
{
  const vector<bigint> & vi = v.get_entries();
  auto toint = [](const bigint& a) {return is_int(a)? I2int(a) : int(0);};
  vector<int> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), toint);
  return vec_i(w);
}

vec_i to_vec_i(const vec_l& v)
{
  const vector<long> & vi = v.get_entries();
  auto toint = [](const long& a) {return int(a);};
  vector<int> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), toint);
  return vec_i(w);
}

vec_l to_vec_l(const vec_m& v)
{
  const vector<bigint> & vi = v.get_entries();
  auto tolong = [](const bigint& a) {return is_long(a)? I2long(a) : long(0);};
  vector<long> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), tolong);
  return vec_l(w);
}

vec_l to_vec_l(const vec_i& v)
{
  const vector<int> & vi = v.get_entries();
  auto tolong = [](const int& a) {return long(a);};
  vector<long> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), tolong);
  return vec_l(w);
}

// Instantiate vecT template functions for T=int
template int dim<int>(const vecT<int>&);
template int operator*<int>(const vecT<int>&, const vecT<int>&);
template int content<int>(const vecT<int>&);
template int maxabs<int>(const vecT<int>&);
template int operator==<int>(const vecT<int>&, const vecT<int>&);
template int operator!=<int>(const vecT<int>&, const vecT<int>&);
template int trivial<int>(const vecT<int>&);                  // is v all 0
template int member<int>(const int& a, const vecT<int>& v);//tests if a=v[i] for some i
template vecT<int> reverse<int>(const vecT<int>& order);
template ostream& operator<<<int> (ostream&s, const vecT<int>&);
template istream& operator>><int> (istream&s, vecT<int>&);
template void swapvec<int>(vecT<int>& v, vecT<int>& w);
template int lift(const vecT<int>& v, const int& pr, vecT<int>& ans);

// Instantiate vecT template functions for T=long
template int dim<long>(const vecT<long>&);
template long operator*<long>(const vecT<long>&, const vecT<long>&);
template long content<long>(const vecT<long>&);
template long maxabs<long>(const vecT<long>&);
template int operator==<long>(const vecT<long>&, const vecT<long>&);
template int operator!=<long>(const vecT<long>&, const vecT<long>&);
template int trivial<long>(const vecT<long>&);                  // is v all 0
template int member<long>(const long& a, const vecT<long>& v);//tests if a=v[i] for some i
template vecT<long> reverse<long>(const vecT<long>& order);
template ostream& operator<<<long> (ostream&s, const vecT<long>&);
template istream& operator>><long> (istream&s, vecT<long>&);
template void swapvec<long>(vecT<long>& v, vecT<long>& w);
template int lift(const vecT<long>& v, const long& pr, vecT<long>& ans);

// Instantiate vecT template functions for T=bigint
template int dim<bigint>(const vecT<bigint>&);
template bigint operator*<bigint>(const vecT<bigint>&, const vecT<bigint>&);
template bigint content<bigint>(const vecT<bigint>&);
template bigint maxabs<bigint>(const vecT<bigint>&);
template int operator==<bigint>(const vecT<bigint>&, const vecT<bigint>&);
template int operator!=<bigint>(const vecT<bigint>&, const vecT<bigint>&);
template int trivial<bigint>(const vecT<bigint>&);                  // is v all 0
template int member<bigint>(const bigint& a, const vecT<bigint>& v);//tests if a=v[i] for some i
template vecT<bigint> reverse<bigint>(const vecT<bigint>& order);
template ostream& operator<<<bigint> (ostream&s, const vecT<bigint>&);
template istream& operator>><bigint> (istream&s, vecT<bigint>&);
template void swapvec<bigint>(vecT<bigint>& v, vecT<bigint>& w);
template int lift(const vecT<bigint>& v, const bigint& pr, vecT<bigint>& ans);
