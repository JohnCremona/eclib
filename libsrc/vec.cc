// vec.cc: implementation of integer vector classes
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
 
// Only to be included by vector.cc

// Definitions of member operators and functions:

vec::vec(long n)
{
  entries.resize(n, 0);
}

vec::vec(vector<scalar> arr) :entries(arr) {}

vec::vec(const vec& v) :entries(v.entries) {} // copy constructor

void vec::init(long n)                 // (re)-initializes
{
  entries.resize(n, 0);
}

vec& vec::operator=(const vec& v)                    // assignment
{
 if (this==&v) return *this;
 entries = v.entries;
 return *this;
}

scalar& vec::operator[](long i)
{
  return entries.at(i-1);
}

scalar vec::operator[](long i) const
{
  return entries.at(i-1);
}

vec& vec::operator+=(const vec& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const scalar& wi, const scalar& vi) { return vi + wi;});
  return *this;
}

void vec::addmodp(const vec& w, scalar pr)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [pr](const scalar& wi, const scalar& vi) { return xmod(wi+vi,pr);});
}

vec& vec::operator-=(const vec& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const scalar& wi, const scalar& vi) { return vi - wi;});
  return *this;
}

vec& vec::operator*=(scalar scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const scalar& vi) {return vi * scal;});
  return *this;
}

vec& vec::operator/=(scalar scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const scalar& vi) {return vi / scal;});
  return *this;
}

vec vec::slice(long first, long last) const       // returns subvector
{
 if (last==-1) {last=first; first=1;}
 vec ans(last-first+1);
 // cout<<"slicing vector "<<(*this)<<" (size "<<entries.size()<<") from "<<first<<" to "<<last<<endl;
 std::copy(entries.begin()+first-1, entries.begin()+last, ans.entries.begin());
 // cout<<"slice is "<<ans<<" (size "<<ans.entries.size()<<")"<<endl;
 return ans;
}

vec vec::operator[](const vec& index) const  // returns v[index[j]]
{
  vec w(index.entries.size());
  std::transform(index.entries.begin(), index.entries.end(), w.entries.begin(),
                 [this](const scalar& i) {return entries.at(i-1);});
  return w;
}

scalar vec::sub(long i) const
{
  return entries.at(i-1);
}

void vec::set(long i, scalar x)
{
  entries.at(i-1) = x;
}

void vec::add(long i, scalar x)
{
  entries.at(i-1) += x;
}

void vec::add_modp(long i, scalar x, scalar p)
{
  entries.at(i-1) = xmod(entries.at(i-1)+x,p);
}

// Definitions of non-member, friend operators and functions

vec reduce_modp(const vec& v, const scalar& p)
{
  if (p==0) return v;
  vec w(v.entries.size());
  std::transform(v.entries.begin(), v.entries.end(), w.entries.begin(),
                 [p](const scalar& vi) {return xmod(vi,p);});
  return w;
}

scalar operator*(const vec& v, const vec& w)
{
  return std::inner_product(v.entries.begin(), v.entries.end(), w.entries.begin(), 0);
}

int operator==(const vec& v, const vec& w)
{
  return v.entries == w.entries;
}

int trivial(const vec& v)
{
  return std::all_of(v.entries.begin(), v.entries.end(), [](const scalar& vi) {return vi==0;});
}

ostream& operator<<(ostream& s, const vec& v)
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

istream& operator>>(istream& s, vec& v)
{
  for (scalar& vi : v.entries)
    s>>vi;
  return s;
}

// Definition of non-friend operators and functions

vec iota(scalar n)
{
 vec v(n);
 std::iota(v.entries.begin(), v.entries.end(), 1);
 return v;
}

scalar content(const vec& v)
{
  return v.entries.empty()?
    1 :
    std::accumulate(v.entries.begin(), v.entries.end(), 0,
                    [](const scalar& x, const scalar& y) {return gcd(x,y);});
}

void swapvec(vec& v, vec& w)
{
  std::swap(v.entries, w.entries);
}

int member(scalar a, const vec& v)
{
  return std::find(v.entries.begin(), v.entries.end(), a) != v.entries.end();
}

vec reverse(const vec& order)
{
  vec ans(order);
  std::reverse(ans.entries.begin(), ans.entries.end());
  return ans;
}

vec express(const vec& v, const vec& v1, const vec& v2)
{
   scalar v1v1 = v1 * v1;
   scalar v1v2 = v1 * v2;
   scalar v2v2 = v2 * v2;
   scalar vv1 = v * v1;
   scalar vv2 = v * v2;
   vec ans({vv1*v2v2 - vv2*v1v2,  vv2*v1v1 - vv1*v1v2, v1v1*v2v2 - v1v2*v1v2});
   makeprimitive(ans);
   if (ans[3]*v!=ans[1]*v1+ans[2]*v2)
     cerr << "Error in express: v is not in <v1,v2>"<<endl;
   return ans;
}

//#define LIFT_DEBUG
int lift(const vec& v, scalar pr, vec& ans)
{
  long i0, i, j, d = dim(v);
 scalar nu, de;
 int succ;
 float lim = sqrt(pr/2.0)-1;
 scalar maxallowed = 10*int(lim);
#ifdef LIFT_DEBUG
 cout<<"Lifting vector v = "<<v<<endl;
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

 ans = v; // starts as a copy, and will be rescaled in place
 scalar vi0, inv_vi0, vi, maxvi=0;
 for(i0=1; i0<=d; i0++)
   {
     // scale so that i0'th entry is 1 mod p, then reduce vector
     // entries mod p to lie in (-p/2,p/2), and find the maximum
     // entry:
     while((vi0=mod(v[i0],pr))==0) {i0++;} // skip over any zero entries
     inv_vi0=invmod(vi0,pr);
     for (i=1; i<=d; i++)
       {
         ans[i]=vi=mod(xmodmul(inv_vi0,ans[i],pr),pr);
         maxvi=max(maxvi,abs(vi));
       }
#ifdef LIFT_DEBUG
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
         succ=modrat(ans[i],pr,lim,nu,de);     de=abs(de);
         if ((!succ)||(de==1)) continue; // loop on i
         // scale by de & recompute max entry:
#ifdef LIFT_DEBUG
         cout<<"Scaling by d="<<de<<endl;
#endif
         maxvi = 0;
         for (j=1; j<=d; j++) 
           {
             ans[j] = vi = mod(xmodmul(de,ans[j],pr),pr);
             maxvi=max(maxvi,abs(vi));
           }
#ifdef LIFT_DEBUG
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

scalar dotmodp(const vec& v1, const vec& v2, scalar pr)
{
  auto a = [pr] (const scalar& x, const scalar& y) {return mod(x+y,pr);};
  auto m = [pr] (const scalar& x, const scalar& y) {return mod(x*y,pr);};
  return std::inner_product(v1.entries.begin(), v1.entries.end(), v2.entries.begin(), 0, a, m);
}
