// mvector.cc: implementations of multiprecision integer vector class
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
 
#include "marith.h"
#include "mvector.h"

// Definitions of member operators and functions:

vec_m::~vec_m()
{
  delete[] entries;
}

vec_m::vec_m(long n)
{
 d=n;
 entries=new bigint[n]; 
 if (!entries) {cout<<"Out of memory!\n"; abort();}
 bigint *v=entries; while(n--) *v++=0;
}

vec_m::vec_m(long n, bigint* array)   //array must have at least n elements!
{
 d=n;
 entries=new bigint[n];  
 if (!entries) {cout<<"Out of memory!\n"; abort();} 
 bigint *v=entries; while(n--) *v++=*array++;
}

vec_m::vec_m(const vec_m& v)                       // copy constructor
{
  d=v.d;
  entries=new bigint[d];  
  if (!entries) {cout<<"Out of memory!\n"; abort();}
  bigint *v1=entries, *v2=v.entries; long n=d;
  while(n--) *v1++=*v2++;
}

vec_m::vec_m(const vec_i& v)
{
  d=v.d;
  entries=new bigint[d];  
  if (!entries) {cout<<"Out of memory!\n"; abort();}
  bigint *v1=entries; int *v2=v.entries; long n=d;
  while(n--) *v1++=*v2++;
}

vec_m::vec_m(const vec_l& v)
{
  d=v.d;
  entries=new bigint[d];  
  if (!entries) {cout<<"Out of memory!\n"; abort();}
  bigint *v1=entries; long *v2=v.entries; long n=d;
  while(n--) *v1++=*v2++;
}

void vec_m::init(long n)                 // (re)-initializes 
{
 if (d!=n) // no point in deleting if same size
   {
     delete[] entries;
     d = n;
     entries=new bigint[d];  
     if (!entries) {cout<<"Out of memory!\n"; abort();}
   }
 bigint *v=entries; while(n--) *v++=0;
}

vec_m& vec_m::operator=(const vec_m& v)                    // assignment
{
 if (this==&v) return *this;
 if (d!=v.d) // no point in deleting if new is same size
   {
     delete[] entries;
     d = v.d;
     entries=new bigint[d];  
     if (!entries) {cout<<"Out of memory!\n"; abort();}
   }
 bigint *v1=entries, *v2=v.entries; long n=d; 
 while(n--) *v1++=*v2++;
 return *this;
}

bigint& vec_m::operator[](long i) const
{
 if ((i>0) && (i<=d)) return entries[i-1];
 else {cout << "bad subscript in vec_m::operator[]\n"; abort(); return entries[0];}
}

vec_m& vec_m::operator+=(const vec_m& q2)
{
  bigint* vi=entries, *wi=q2.entries; long i=d;
  if (d==q2.d) {while(i--)(*vi++)+=(*wi++);}
  else {cout << "Incompatible vec_ms in vec_m::operator+=\n"; abort();}
  return *this;
}

void vec_m::addmodp(const vec_m& w, const bigint& pr)
{
  bigint* vi=entries, *wi=w.entries; long i=d;
  if (d==w.d) {while(i--) {*vi = mod((*wi++)+(*vi),pr);vi++;}}
  else {cout << "Incompatible vec_ms in vec_m::addmodp\n";abort();}
}

vec_m& vec_m::operator-=(const vec_m& q2)
{
  bigint* vi=entries; bigint* wi=q2.entries; long i=d;
  if (d==q2.d) {while(i--)(*vi++)-=(*wi++);}
  else {cout << "Incompatible vec_ms in vec_m::operator-=\n";abort();}
  return *this;
}

vec_m& vec_m::operator*=(const bigint& scal)
{
  bigint* vi=entries; long i=d;
  while (i--) (*vi++) *= scal;
  return *this;
}

vec_m& vec_m::operator/=(const bigint& scal)
{
  bigint* vi=entries; long i=d;
  while (i--) (*vi++) /= scal;
  return *this;
}

vec_m vec_m::slice(long first, long last) const       // returns subvec_m
{
 if (last==-1) {last=first; first=1;}
 long n = last-first+1;
 vec_m ans(n);
 bigint *entriesi=entries+(first-1), *ansi=ans.entries; long i=n;
 while (i--) *ansi++ = *entriesi++;
 return ans;
}

vec_m vec_m::operator[](const vec_i& index) const//returns v[index[j]]
{long dd = dim(index); vec_m w(dd);
 for(long i=1; i<=dd;i++) w[i] = entries[index[i]-1];
 return w;
}

vec_m vec_m::operator[](const vec_l& index) const
{long dd = dim(index); vec_m w(dd);
 for(long i=1; i<=dd;i++) w[i] = entries[index[i]-1];
 return w;
}

bigint vec_m::sub(long i) const
{
 if ((i>0) && (i<=d)) return entries[i-1];
 else {bigint zero; cout << "bad subscript in vec_m::sub\n"; abort(); return zero;}
}

void vec_m::set(long i, const bigint& x)
{
 if ((i>0) && (i<=d)) entries[i-1]=x;
 else {cout << "bad subscript in vec_m::set\n";abort();}
}

void vec_m::add(long i, const bigint& x)
{
 if ((i>0) && (i<=d)) entries[i-1]+=x;
 else {cout << "bad subscript in vec_m::add\n";abort();}
}

vec_l vec_m::shorten(long x) const  //converts to a vector of longs
{
  vec_l ans(d);
  bigint *v1=entries; long *v2=ans.entries; long n=d;
  while(n--)
    {
      bigint& veci = *v1++;
      if(is_long(veci)) 
        *v2 = I2long(veci);
      else 
	{
	  cout << "Problem shortening bigint " << veci << " to a long!" << endl;
	  abort();
	}
      v2++;
    }
  return ans;
}

vec_i vec_m::shorten(int x) const  //converts to a vector of ints
{
  vec_i ans(d);
  bigint *v1=entries; int *v2=ans.entries; long n=d;
  while(n--)
    {
      bigint& veci = *v1++;
      if(is_int(veci)) 
        *v2 = I2int(veci);
      else 
	{
	  cout << "Problem shortening bigint " << veci << " to an int!" << endl;
	  abort();
	}
      v2++;
    }
  return ans;
}

// Definitions of non-member, friend operators and functions

bigint operator*(const vec_m& v, const vec_m& w)
{
 long dim=v.d; bigint dot;
 bigint* vi=v.entries, *wi=w.entries;
 if (dim==w.d) 
   while (dim--) dot+= (*vi++)*(*wi++);
 else 
   {
     cout << "Unequal dimensions in dot product\n";
     abort();
   }
 return dot;
}

int operator==(const vec_m& v, const vec_m& w)
{
   long dim=v.d;
   long equal = (dim==w.d);
   bigint* vi=v.entries, *wi=w.entries;
   while ((dim--) && equal) equal = ((*vi++)==(*wi++));
   return equal;
}

int trivial(const vec_m& v)
{
   int ans=1, i=v.d;   bigint* vi=v.entries;
   while ((i--)&&ans) ans=(is_zero(*vi++));
   return ans;
}

ostream& operator<<(ostream& s, const vec_m& v)
{
   long i=v.d; bigint* vi=v.entries;
   s << "[";
   while (i--) {s<<(*vi++); if(i)s<<",";}
   s << "]";
   return s;
}

istream& operator>>(istream& s, vec_m& v)
{
 long i = v.d;
 bigint* vi = v.entries;
 while (i--) s >> (*vi++);
 return s;
}

bigint mvecgcd(const vec_m& v)
{
 long i=v.d; bigint g; bigint *vi=v.entries;
 while ((i--)&&(!is_one(g))) g=gcd(g,*vi++);
 return g;
}

void swapvec(vec_m& v, vec_m& w)
{bigint *temp; 
 if (v.d==w.d) {temp=v.entries; v.entries=w.entries; w.entries=temp;}
 else {cout << "Attempt to swap vec_ms of different lengths!\n";abort();}
}

int member(const bigint& a, const vec_m& v)
{int ans=0; long i=dim(v); bigint* vi=v.entries;
 while (i--&&!ans) ans=(a==(*vi++));
 return ans;
}
 
// Definition of non-friend operators and functions

vec_m express(const vec_m& v, const vec_m& v1, const vec_m& v2)
{
  vec_m ans(3);
  static bigint one; one=1;
   bigint v1v1 = v1 * v1;
   bigint v1v2 = v1 * v2;
   bigint v2v2 = v2 * v2;
   bigint vv1 = v * v1;
   bigint vv2 = v * v2;
   ans[1]= vv1*v2v2 - vv2*v1v2;
   ans[2]= vv2*v1v1 - vv1*v1v2;
   ans[3]= v1v1*v2v2 - v1v2*v1v2;
   bigint g = mvecgcd(ans);
   if (g>one) ans/=g;
   if (ans[3]*v!=ans[1]*v1+ans[2]*v2)
     {
       cout << "Error in express: v is not in <v1,v2>\n";
       abort();
     }
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
     bigint& vi=v[i];
     modrat(vi,pr,lim,nu,de);
     nums[i]=nu; denoms[i]=de;
     g=lcm(g,de);
   }
 for (i=1; i<=d; i++) nums[i]*=(g/denoms[i]);
 return nums;
}

int liftok(vec_m& v, const bigint& pr)
{
 long i, d = dim(v); bigint g, nu, de; int success=1, succ;
 bigint lim = sqrt(pr>>1);
 for (i=1, g=1; i<=d; i++) 
   {
     bigint& vi=v[i]; succ=modrat(vi,pr,lim,nu,de);
     success = success && succ;
// Can't say success=success&&modrat(...) as then after first fail it does
// not call modrat at all due to clever compiler!
     g=lcm(g,de);
   }
 for (i=1; i<=d; i++) v[i] = mod(g*v[i],pr);
 return success;
}

bigint dotmodp(const vec_m& v1, const vec_m& v2, const bigint& pr)
{
  bigint ans;
  for(long i=1; i<=dim(v1); i++) ans=mod(ans+mod(v1[i]*v2[i],pr),pr);
  return ans;
}

long dim(const vec_m& v) {return v.d;}

int operator!=(const vec_m& v, const vec_m& w) { return !(v==w);}

vec_m operator+(const vec_m& v) 
{return v;}

vec_m operator-(const vec_m& v) 
{vec_m ans(v); ans*=BIGINT(-1); return ans;}

vec_m operator+(const vec_m& v1, const vec_m& v2)
{vec_m ans(v1); ans+=v2; return ans;}

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
{ bigint g=mvecgcd(v); if(!(is_one(g)||is_zero(g))) v/=g;}

void elim(const vec_m& a, vec_m& b, long pos)
{ (b*=a[pos])-=(b[pos]*a);}

void elim1(const vec_m& a, vec_m& b, long pos)
{ (b*=a[pos])-=(b[pos]*a); makeprimitive(b);}

void elim2(const vec_m& a, vec_m& b, long pos, const bigint& lastpivot)
{ ((b*=a[pos])-=(b[pos]*a))/=lastpivot;}

