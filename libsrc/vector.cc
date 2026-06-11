// vector.cc: implementations of integer vector classes
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
 
#include <random>
#include "eclib/convert.h"
#include "eclib/linalg.h"

// Definitions of member operators and functions:

template<class T>
Zvec<T>::Zvec(long n)
{
  entries.resize(n, T(0));
}

template<class T>
Zvec<T>::Zvec(const vector<T>& arr) :entries(arr) {}

template<class T>
Zvec<T>::Zvec(const Zvec<T>& v) :entries(v.entries) {} // copy constructor

template<class T>
void Zvec<T>::init(long n)     // (re)-initializes
{
  entries.resize(n, T(0));
}

template<class T>
void Zvec<T>::clear()
{
  entries.assign(entries.size(), T(0));
}

template<class T>
Zvec<T>& Zvec<T>::operator=(const Zvec<T>& v)    // assignment
{
 if (this==&v) return *this;
 entries = v.entries;
 return *this;
}

template<class T>
T& Zvec<T>::operator[](long i)
{
  return entries.at(i-1);
}

template<class T>
T Zvec<T>::operator[](long i) const
{
  return entries.at(i-1);
}

template<class T>
Zvec<T>& Zvec<T>::operator+=(const Zvec<T>& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const T& wi, const T& vi) { return vi + wi;});
  return *this;
}

template<class T>
void Zvec<T>::addmodp(const Zvec<T>& w, const T& pr)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [pr](const T& wi, const T& vi) { return mod(wi+vi,pr);});
}

template<class T>
Zvec<T>& Zvec<T>::operator-=(const Zvec<T>& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const T& wi, const T& vi) { return vi - wi;});
  return *this;
}

template<class T>
Zvec<T>& Zvec<T>::operator*=(const T& scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const T& vi) {return vi * scal;});
  return *this;
}

template<class T>
Zvec<T>& Zvec<T>::operator/=(const T& scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const T& vi) {return vi / scal;});
  return *this;
}

template<class T>
Zvec<T> Zvec<T>::slice(long first, long last) const       // returns subvector
{
 if (last==-1) {last=first; first=1;}
 Zvec<T> ans(last-first+1);
 std::copy(entries.begin()+first-1, entries.begin()+last, ans.entries.begin());
 return ans;
}

template<class T>
Zvec<T> Zvec<T>::operator[](const Zvec<int>& index) const  // returns v[index[j]]
{
  Zvec<T> w(index.dim());
  const vector<int>& vi = index.get_entries();
  std::transform(vi.begin(), vi.end(), w.entries.begin(),
                 [this](const int& i) {return entries.at(i-1);});
  return w;
}

template<class T>
Zvec<T> Zvec<T>::operator[](const Zvec<long>& index) const  // returns v[index[j]]
{
  Zvec<T> w(index.dim());
  const vector<long>& vi = index.get_entries();
  std::transform(vi.begin(), vi.end(), w.entries.begin(),
                 [this](const int& i) {return entries.at(i-1);});
  return w;
}

template<class T>
T Zvec<T>::sub(long i) const
{
  return entries.at(i-1);
}

template<class T>
void Zvec<T>::set(long i, const T& x)
{
  entries.at(i-1) = x;
}

template<class T>
void Zvec<T>::add(long i, const T& x)
{
  entries.at(i-1) += x;
}

template<class T>
void Zvec<T>::add_modp(long i, const T& x, const T& p)
{
  entries.at(i-1) = mod(entries.at(i-1)+x,p);
}

template<class T>
void Zvec<T>::reduce_mod_p(const T& p)
{
  if (p==0) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [p](const T& vi) {return mod(vi,p);});
}

template<class T>
Zvec<T> Zvec<T>::iota(long n)
{
  Zvec<T> v(n);
  std::iota(v.entries.begin(), v.entries.end(), T(1));
  return v;
}

template<class T>
Zvec<T> Zvec<T>::unit_vector(long n, long i)
{
  Zvec<T> v(n);
  v.set(i, T(1));
  return v;
}

// Definitions of non-member operators and functions

template<class T> int dim(const Zvec<T>& v)
{
  return v.dim();
}
template int dim<int>(const Zvec<int>&);
template int dim<long>(const Zvec<long>&);
template int dim<ZZ>(const Zvec<ZZ>&);
template int dim<INT>(const Zvec<INT>&);


template<class T>
T operator*(const Zvec<T>& v, const Zvec<T>& w)
{
  return std::inner_product(v.entries.begin(), v.entries.end(), w.entries.begin(), T(0));
}
template int operator*<int>(const Zvec<int>&, const Zvec<int>&);
template long operator*<long>(const Zvec<long>&, const Zvec<long>&);
template ZZ operator*<ZZ>(const Zvec<ZZ>&, const Zvec<ZZ>&);
template INT operator*<INT>(const Zvec<INT>&, const Zvec<INT>&);

template<class T>
int operator==(const Zvec<T>& v, const Zvec<T>& w)
{
  return v.entries == w.entries;
}
template int operator==<int>(const Zvec<int>&, const Zvec<int>&);
template int operator==<long>(const Zvec<long>&, const Zvec<long>&);
template int operator==<ZZ>(const Zvec<ZZ>&, const Zvec<ZZ>&);
template int operator==<INT>(const Zvec<INT>&, const Zvec<INT>&);

template<class T>
int operator!=(const Zvec<T>& v, const Zvec<T>& w)
{
  return v.entries != w.entries;
}
template int operator!=<int>(const Zvec<int>&, const Zvec<int>&);
template int operator!=<long>(const Zvec<long>&, const Zvec<long>&);
template int operator!=<ZZ>(const Zvec<ZZ>&, const Zvec<ZZ>&);
template int operator!=<INT>(const Zvec<INT>&, const Zvec<INT>&);

template<class T> Zvec<T> operator+(const Zvec<T>& v)
{
  return v;
}
template Zvec<int> operator+<int>(const Zvec<int>&);
template Zvec<long> operator+<long>(const Zvec<long>&);
template Zvec<ZZ> operator+<ZZ>(const Zvec<ZZ>&);
template Zvec<INT> operator+<INT>(const Zvec<INT>&);

template<class T> Zvec<T> operator-(const Zvec<T>& v)
{
  return T(-1)*v;
}
template Zvec<int> operator-<int>(const Zvec<int>&);
template Zvec<long> operator-<long>(const Zvec<long>&);
template Zvec<ZZ> operator-<ZZ>(const Zvec<ZZ>&);
template Zvec<INT> operator-<INT>(const Zvec<INT>&);

template<class T> Zvec<T> operator+(const Zvec<T>& v1, const Zvec<T>& v2)
{
  Zvec<T> w(v1); w+=v2; return w;
}
template Zvec<int> operator+<int>(const Zvec<int>&, const Zvec<int>&);
template Zvec<long> operator+<long>(const Zvec<long>&, const Zvec<long>&);
template Zvec<ZZ> operator+<ZZ>(const Zvec<ZZ>&, const Zvec<ZZ>&);
template Zvec<INT> operator+<INT>(const Zvec<INT>&, const Zvec<INT>&);

template<class T> Zvec<T> addmodp(const Zvec<T>& v1, const Zvec<T>& v2, const T& pr)
{
  Zvec<T> w(v1); w.addmodp(v2,pr); return w;
}
template Zvec<int> addmodp<int>(const Zvec<int>&, const Zvec<int>&, const int&);
template Zvec<long> addmodp<long>(const Zvec<long>&, const Zvec<long>&, const long&);
template Zvec<ZZ> addmodp<ZZ>(const Zvec<ZZ>&, const Zvec<ZZ>&, const ZZ&);
template Zvec<INT> addmodp<INT>(const Zvec<INT>&, const Zvec<INT>&, const INT&);

template<class T> Zvec<T> reduce_mod_p(const Zvec<T>& v, const T& p)
{
  Zvec<T> w(v); w.reduce_mod_p(p); return w;
}
template Zvec<int> reduce_mod_p<int>(const Zvec<int>&, const int&);
template Zvec<long> reduce_mod_p<long>(const Zvec<long>&, const long&);
template Zvec<ZZ> reduce_mod_p<ZZ>(const Zvec<ZZ>&, const ZZ&);
template Zvec<INT> reduce_mod_p<INT>(const Zvec<INT>&, const INT&);

template<class T> Zvec<T> operator-(const Zvec<T>& v1, const Zvec<T>& v2)
{
  Zvec<T> w(v1); w-=v2; return w;
}
template Zvec<int> operator-<int>(const Zvec<int>&, const Zvec<int>&);
template Zvec<long> operator-<long>(const Zvec<long>&, const Zvec<long>&);
template Zvec<ZZ> operator-<ZZ>(const Zvec<ZZ>&, const Zvec<ZZ>&);
template Zvec<INT> operator-<INT>(const Zvec<INT>&, const Zvec<INT>&);

template<class T> Zvec<T> operator*(const T& scal, const Zvec<T>& v)
{
  Zvec<T> w(v); w*=scal; return w;
}
template Zvec<int> operator*<int>(const int&, const Zvec<int>&);
template Zvec<long> operator*<long>(const long&, const Zvec<long>&);
template Zvec<ZZ> operator*<ZZ>(const ZZ&, const Zvec<ZZ>&);
template Zvec<INT> operator*<INT>(const INT&, const Zvec<INT>&);

template<class T> Zvec<T> operator/(const Zvec<T>& v, const T& scal)
{
  Zvec<T> w(v); w/=scal; return w;
}
template Zvec<int> operator/<int>(const Zvec<int>&, const int&);
template Zvec<long> operator/<long>(const Zvec<long>&, const long&);
template Zvec<ZZ> operator/<ZZ>(const Zvec<ZZ>&, const ZZ&);
template Zvec<INT> operator/<INT>(const Zvec<INT>&, const INT&);

template<class T> void make_primitive(Zvec<T>& v)
{
  T g=content(v); if (g>1) v/=g;
}
template void make_primitive<int>(Zvec<int>&);
template void make_primitive<long>(Zvec<long>&);
template void make_primitive<ZZ>(Zvec<ZZ>&);
template void make_primitive<INT>(Zvec<INT>&);

template<class T>
int trivial(const Zvec<T>& v)
{
  return std::all_of(v.entries.begin(), v.entries.end(), [](const T& vi) {return vi==0;});
}
template int trivial<int>(const Zvec<int>&);
template int trivial<long>(const Zvec<long>&);
template int trivial<ZZ>(const Zvec<ZZ>&);                  // is v all 0
template int trivial<INT>(const Zvec<INT>&);                  // is v all 0

template<class T>
ostream& operator<<(ostream& s, const Zvec<T>& v)
{
  vec_out(s, v.entries, "[", "]", ",");
  return s;
}
template ostream& operator<<<int> (ostream&s, const Zvec<int>&);
template ostream& operator<<<long> (ostream&s, const Zvec<long>&);
template ostream& operator<<<ZZ> (ostream&s, const Zvec<ZZ>&);
template ostream& operator<<<INT> (ostream&s, const Zvec<INT>&);

// Two possible input modes:
// (1) entries separated by whitespace e.g. "1 2 3"
// (2) same as default output, e.g. "[1,2,3]"
// In both cases the size of v must be already set.
template<class T>
istream& operator>>(istream& s, Zvec<T>& v)
{
  // cout << "Reading into vector v = " << v << endl;
  char c;
  s >> ws; // skip whitespace
  s >> c;
  // cout << "First character read: " << c << endl;
  switch (c) {
  case '[':
    // cout << "That was [" << endl;
    // Read the entries + one char (separator commas or final ']')
    for (T& vi : v.entries)
      {
        s>>ws>>vi>>ws>>c;
        // cout << "Read entry " << vi << " and character " << c << endl;
      }
    break;
  default:
    // put back the first non-whitespace char
    // cout << "That was not [ so putting it back" << endl;
    s.unget();
    // read the entries
    for (T& vi : v.entries)
      {
        s>>ws>>vi;
        // cout << "Read entry " << vi << endl;
      }
  }
  // cout << "Finally v = " << v << endl;
  return s;
}
template istream& operator>><int> (istream&s, Zvec<int>&);
template istream& operator>><long> (istream&s, Zvec<long>&);
template istream& operator>><ZZ> (istream&s, Zvec<ZZ>&);
template istream& operator>><INT> (istream&s, Zvec<INT>&);

template<class T>
T content(const Zvec<T>& v)
{
  return v.entries.empty()?
    T(1) :
    std::accumulate(v.entries.begin(), v.entries.end(), T(0),
                    [](const T& x, const T& y) {return gcd(x,y);});
}
template int content<int>(const Zvec<int>&);
template long content<long>(const Zvec<long>&);
template ZZ content<ZZ>(const Zvec<ZZ>&);
template INT content<INT>(const Zvec<INT>&);

template<class T>
T maxabs(const Zvec<T>& v)
{
  return v.entries.empty()?
    T(0) :
    std::accumulate(v.entries.begin(), v.entries.end(), T(0),
                    [](const T& x, const T& y) {return max(x,abs(y));});
}
template int maxabs<int>(const Zvec<int>&);
template long maxabs<long>(const Zvec<long>&);
template ZZ maxabs<ZZ>(const Zvec<ZZ>&);
template INT maxabs<INT>(const Zvec<INT>&);

template<class T>
void swapvec(Zvec<T>& v, Zvec<T>& w)
{
  std::swap(v.entries, w.entries);
}
template void swapvec<int>(Zvec<int>& v, Zvec<int>& w);
template void swapvec<long>(Zvec<long>& v, Zvec<long>& w);
template void swapvec<ZZ>(Zvec<ZZ>& v, Zvec<ZZ>& w);
template void swapvec<INT>(Zvec<INT>& v, Zvec<INT>& w);

template<class T>
int member(const T& a, const Zvec<T>& v)
{
  return std::find(v.entries.begin(), v.entries.end(), a) != v.entries.end();
}
template int member<int>(const int& a, const Zvec<int>& v);//tests if a=v[i] for some i
template int member<long>(const long& a, const Zvec<long>& v);//tests if a=v[i] for some i
template int member<ZZ>(const ZZ& a, const Zvec<ZZ>& v);//tests if a=v[i] for some i
template int member<INT>(const INT& a, const Zvec<INT>& v);//tests if a=v[i] for some i

template<class T>
Zvec<T> reverse(const Zvec<T>& order)
{
  Zvec<T> ans(order);
  std::reverse(ans.entries.begin(), ans.entries.end());
  return ans;
}
template Zvec<int> reverse<int>(const Zvec<int>& order);
template Zvec<long> reverse<long>(const Zvec<long>& order);
template Zvec<ZZ> reverse<ZZ>(const Zvec<ZZ>& order);
template Zvec<INT> reverse<INT>(const Zvec<INT>& order);

//#define DEBUG_LIFT

template<class T>
int lift(const Zvec<T>& v, const T& pr, Zvec<T>& ans)
{
  long i0, i, j, d = v.dim();
  T nu, de;
  T lim = isqrt(pr>>1)-1;
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
template int lift(const Zvec<int>& v, const int& pr, Zvec<int>& ans);
template int lift(const Zvec<long>& v, const long& pr, Zvec<long>& ans);
template int lift(const Zvec<ZZ>& v, const ZZ& pr, Zvec<ZZ>& ans);
template int lift(const Zvec<INT>& v, const INT& pr, Zvec<INT>& ans);

template<class T>
T dotmodp(const Zvec<T>& v1, const Zvec<T>& v2, const T& pr)
{
  auto a = [pr] (const T& x, const T& y) {return mod(x+y,pr);};
  auto m = [pr] (const T& x, const T& y) {return xmodmul(x,y,pr);};
  return std::inner_product(v1.entries.begin(), v1.entries.end(), v2.entries.begin(), T(0), a, m);
}
template int dotmodp(const Zvec<int>&, const Zvec<int>&, const int&);
template long dotmodp(const Zvec<long>&, const Zvec<long>&, const long&);
template ZZ dotmodp(const Zvec<ZZ>&, const Zvec<ZZ>&, const ZZ&);
template INT dotmodp(const Zvec<INT>&, const Zvec<INT>&, const INT&);

template<class T>
vec_m to_vec_m(const Zvec<T>& V)
{
  const vector<T> & Vi = V.get_entries();
  vector<ZZ> w(Vi.size());
  std::transform(Vi.begin(), Vi.end(), w.begin(), [](const T& x) {return to_ZZ(x);});
  return vec_m(w);
}
template vec_m to_vec_m<int>(const Zvec<int>& M);
template vec_m to_vec_m<long>(const Zvec<long>& M);
template vec_m to_vec_m<ZZ>(const Zvec<ZZ>& M);
template vec_m to_vec_m<INT>(const Zvec<INT>& M);

template<class T>
vec_i to_vec_i(const Zvec<T>& V)
{
  const vector<T> & Vi = V.get_entries();
  vector<int> w(Vi.size());
  std::transform(Vi.begin(), Vi.end(), w.begin(), [](const T& x) {return to_int(x);});
  return vec_i(w);
}
template vec_i to_vec_i<int>(const Zvec<int>& M);
template vec_i to_vec_i<long>(const Zvec<long>& M);
template vec_i to_vec_i<ZZ>(const Zvec<ZZ>& M);

template<class T>
vec_l to_vec_l(const Zvec<T>& V)
{
  const vector<T> & Vi = V.get_entries();
  vector<long> w(Vi.size());
  std::transform(Vi.begin(), Vi.end(), w.begin(), [](const T& x) {return to_long(x);});
  return vec_l(w);
}
template vec_l to_vec_l<int>(const Zvec<int>& M);
template vec_l to_vec_l<long>(const Zvec<long>& M);
template vec_l to_vec_l<ZZ>(const Zvec<ZZ>& M);

template<class T>
vec_I to_vec_I(const Zvec<T>& V)
{
  const vector<T> & Vi = V.get_entries();
  vector<INT> w(Vi.size());
  std::transform(Vi.begin(), Vi.end(), w.begin(), [](const T& x) {return to_INT(x);});
  return vec_I(w);
}
template vec_I to_vec_I<int>(const Zvec<int>& M);
template vec_I to_vec_I<long>(const Zvec<long>& M);
template vec_I to_vec_I<ZZ>(const Zvec<ZZ>& M);

// Function to generate a random integer vector of a given size wit
// entries taken uniformly from [minv..maxv], optionally repeating
// until the vector is primitive

vector<int> random_vector(size_t size, int minv, int maxv, int primitive)
{
  if (size==1 && primitive==1 && minv<=1 and maxv>=1)
    return {1};
  // We use static in order to instantiate the random engine and the
  // distribution once only.  It may provoke some thread-safety
  // issues.
  std::uniform_int_distribution<int> distribution(minv, maxv);
  static std::default_random_engine generator;

  vector<int> v(size);
  int cont = 0;
  while (cont==0 || (primitive && cont>1))
    {
      std::generate(v.begin(), v.end(),
                    [&distribution]() { return distribution(generator); });
      cont = std::accumulate(v.begin(), v.end(), 0,
                             [](const int& x, const int& y) {return gcd(x,y);});
    }
  return v;
}


// Instantiate Zvec template classes for T=int, long, ZZ, INT

template class Zvec<int>;
template class Zvec<long>;
template class Zvec<ZZ>;
template class Zvec<INT>;

