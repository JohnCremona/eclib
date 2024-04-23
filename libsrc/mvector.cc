// mvector.cc: implementations of multiprecision integer vector class
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
 
#include <eclib/mvector.h>
#include <eclib/marith.h>

// Definitions of member operators and functions:

vec_m::vec_m(long n)
{
  static const bigint zero(0);
  entries.resize(n, zero);
}

vec_m::vec_m(const vec_i& v)
{
  entries.resize(v.entries.size());
  std::transform(v.entries.begin(), v.entries.end(), entries.begin(), [](const int& vi) {return bigint(vi);});
}

vec_m::vec_m(const vec_l& v)
{
  entries.resize(v.entries.size());
  std::transform(v.entries.begin(), v.entries.end(), entries.begin(), [](const long& vi) {return bigint(vi);});
}

void vec_m::init(long n)                 // (re)-initializes
{
  static const bigint zero(0);
  entries.resize(n, zero);
}

vec_m& vec_m::operator=(const vec_m& v)                    // assignment
{
 if (this==&v) return *this;
 entries = v.entries;
 return *this;
}

bigint vec_m::operator[](long i) const
{
  return entries.at(i-1);
}

bigint& vec_m::operator[](long i)
{
  return entries.at(i-1);
}

vec_m& vec_m::operator+=(const vec_m& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const bigint& wi, const bigint& vi) { return vi + wi;});
  return *this;
}

void vec_m::addmodp(const vec_m& w, const bigint& pr)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [pr](const bigint& wi, const bigint& vi) { return mod(wi+vi,pr);});
}

vec_m& vec_m::operator-=(const vec_m& w)
{
  std::transform(w.entries.begin(), w.entries.end(), entries.begin(), entries.begin(),
                 [](const bigint& wi, const bigint& vi) { return vi - wi;});
  return *this;
}

vec_m& vec_m::operator*=(const bigint& scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const bigint& vi) {return vi * scal;});
  return *this;
}

vec_m& vec_m::operator/=(const bigint& scal)
{
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const bigint& vi) {return vi / scal;});
  return *this;
}

vec_m vec_m::slice(long first, long last) const       // returns subvec_m
{
 if (last==-1) {last=first; first=1;}
 vec_m ans(last-first+1);
 std::copy(entries.begin()+first-1, entries.begin()+last, ans.entries.begin());
 return ans;
}

vec_m vec_m::operator[](const vec_i& index) const //returns v[index[j]]
{
  vec_m w(index.entries.size());
  std::transform(index.entries.begin(), index.entries.end(), w.entries.begin(),
                 [this](const int& i) {return entries.at(i-1);});
 return w;
}

vec_m vec_m::operator[](const vec_l& index) const
{
  vec_m w(index.entries.size());
  std::transform(index.entries.begin(), index.entries.end(), w.entries.begin(),
                 [this](const long& i) {return entries.at(i-1);});
 return w;
}

bigint vec_m::sub(long i) const
{
  return entries.at(i-1);
}

void vec_m::set(long i, const bigint& x)
{
  entries.at(i-1) = x;
}

void vec_m::add(long i, const bigint& x)
{
  entries.at(i-1) += x;
}

vec_l vec_m::shorten(long x) const  //converts to a vector of longs
{
  vec_l ans(entries.size());
  auto tolong = [](const bigint& a) {return is_long(a)? I2long(a) : long(0);};
  std::transform(entries.begin(), entries.end(), ans.entries.begin(), tolong);
  return ans;
}

vec_i vec_m::shorten(int x) const  //converts to a vector of ints
{
  vec_i ans(entries.size());
  auto toint = [](const bigint& a) {return is_long(a)? I2long(a) : int(0);};
  std::transform(entries.begin(), entries.end(), ans.entries.begin(), toint);
  return ans;
}

// Definitions of non-member, friend operators and functions

bigint operator*(const vec_m& v, const vec_m& w)
{
  static const bigint zero(0);
  return std::inner_product(v.entries.begin(), v.entries.end(), w.entries.begin(), zero);
}

int operator==(const vec_m& v, const vec_m& w)
{
  return v.entries == w.entries;
}

int trivial(const vec_m& v)
{
  return std::all_of(v.entries.begin(), v.entries.end(), [](const bigint& vi) {return is_zero(vi);});
}

ostream& operator<<(ostream& s, const vec_m& v)
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

istream& operator>>(istream& s, vec_m& v) // v cannot be const
{
  for (bigint& vi : v.entries)
    s>>vi;
  return s;
}

bigint content(const vec_m& v)
{
  static const bigint zero(0);
  return v.entries.empty()?
    bigint(1) :
    std::accumulate(v.entries.begin(), v.entries.end(), zero,
                    [](const bigint& x, const bigint& y) {return gcd(x,y);});
}

void swapvec(vec_m& v, vec_m& w)
{
  std::swap(v.entries, w.entries);
}

int member(const bigint& a, const vec_m& v)
{
  return std::find(v.entries.begin(), v.entries.end(), a) != v.entries.end();
}

// Definition of non-friend operators and functions

vec_m express(const vec_m& v, const vec_m& v1, const vec_m& v2)
{
  static bigint one(1);
  bigint v1v1 = v1 * v1;
  bigint v1v2 = v1 * v2;
  bigint v2v2 = v2 * v2;
  bigint vv1 = v * v1;
  bigint vv2 = v * v2;
  vec_m ans({vv1*v2v2 - vv2*v1v2, vv2*v1v1 - vv1*v1v2, v1v1*v2v2 - v1v2*v1v2});
  makeprimitive(ans);
  if (ans[3]*v!=ans[1]*v1+ans[2]*v2)
    cerr << "Error in express: v is not in <v1,v2>"<<endl;
  return ans;
}

vec_m lift(const vec_m& v, const bigint& pr)
{
 long i, d = dim(v);
 vec_m nums(d);
 bigint lim = sqrt(pr>>1), nu, de, g;
 g=1;               // g = least common denom. after lifting via modrat
 vec_m denoms(d);
 for (i=1; i<=d; i++)
   {
     modrat(v[i],pr,lim,nu,de);
     nums[i]=nu; denoms[i]=de;
     g=lcm(g,de);
   }
 for (i=1; i<=d; i++) nums[i]*=(g/denoms[i]);
 return nums;
}

int liftok(vec_m& v, const bigint& pr)
{
 long i, d = dim(v);
 bigint g, nu, de, lim = sqrt(pr>>1);
 int success=1;
 for (i=1, g=1; i<=d; i++)
   {
     bigint& vi=v[i];
     int succ=modrat(vi,pr,lim,nu,de);
     success = success && succ;
// Can't say success=success&&modrat(...) as then after first fail it does
// not call modrat at all due to compiler optimisation
     g=lcm(g,de);
   }
 for (i=1; i<=d; i++)
   v[i] = mod(g*v[i],pr);
 return success;
}

bigint dotmodp(const vec_m& v1, const vec_m& v2, const bigint& pr)
{
  static const bigint zero(0);
  auto a = [pr] (const bigint& x, const bigint& y) {return mod(x+y,pr);};
  auto m = [pr] (const bigint& x, const bigint& y) {return mod(x*y,pr);};
  return std::inner_product(v1.entries.begin(), v1.entries.end(), v2.entries.begin(), zero, a, m);
}

long dim(const vec_m& v) {return v.entries.size();}

int operator!=(const vec_m& v, const vec_m& w) { return !(v==w);}

vec_m operator+(const vec_m& v) {return v;}

vec_m operator-(const vec_m& v) {vec_m ans(v); ans*=BIGINT(-1); return ans;}

vec_m operator+(const vec_m& v1, const vec_m& v2) {vec_m ans(v1); ans+=v2; return ans;}

vec_m addmodp(const vec_m& v1, const vec_m& v2, const bigint& pr)
{vec_m ans(v1); ans.addmodp(v2,pr); return ans;}

vec_m operator-(const vec_m& v1, const vec_m& v2)
{vec_m ans(v1); ans-=v2; return ans;}

vec_m operator*(const bigint& scal, const vec_m& v)
{vec_m ans(v); ans*=scal; return ans;}

// NB cannot rely on coercion for next function as compiler tends to 
// construct a zero vector from the long and use dot product!
vec_m operator*(long scal, const vec_m& v)
{vec_m ans(v); ans*=BIGINT(scal); return ans;}

vec_m operator/(const vec_m& v, const bigint& scal)
{vec_m ans(v); ans/=scal; return ans;}

void makeprimitive(vec_m& v)
{ bigint g = content(v); if(!(is_one(g)||is_zero(g))) v/=g;}

void elim(const vec_m& a, vec_m& b, long pos)
{ (b*=a[pos])-=(b[pos]*a);}

void elim1(const vec_m& a, vec_m& b, long pos)
{ (b*=a[pos])-=(b[pos]*a); makeprimitive(b);}

void elim2(const vec_m& a, vec_m& b, long pos, const bigint& lastpivot)
{ ((b*=a[pos])-=(b[pos]*a))/=lastpivot;}

