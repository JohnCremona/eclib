// svec.cc: implementation of sparse integer vector class svec
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
 
// ONLY to be included by svector.cc

// Definitions of member operators and functions:

svec::svec(const vec& v)
  :d(dim(v))
{
  for(int i = 1; i <= d; i++ )
    {
      scalar vi = v[i];
      if(is_nonzero(vi))
        entries[i]=vi;
    }
}

vec svec::as_vec( ) const
{
  vec v(d); // initializes to 0
  for (const auto& wi : entries)
    v[wi.first] = wi.second;
  return v;
}

scalar svec::elem(int i)  const   // returns i'th entry
{
  auto vi = entries.find(i);
  if(vi==entries.end()) return scalar(0);
  return vi->second;
}

void svec::set(int i, const scalar& a)
{
  if(is_nonzero(a)) entries[i]=a;
}

void svec::erase(int i)
{
  auto vi = entries.find(i);
  if(vi==entries.end())
    {
      cerr<<"Error in svec::erase(): cannot delete missing entry #"<<i
	  <<" from v = "<<(*this)<<endl;
    }
  else entries.erase(vi);
}

std::set<int> svec::support() const
{
  std::set<int> ans;
  for (const auto& vi : entries)
    ans.insert(vi.first);
  return ans;
}

void svec::add(int i, const scalar& a)   // adds a to i'th entry
{
  if(is_zero(a)) return;
  auto  vi = entries.find(i);
  if(vi==entries.end())
    entries[i]=a;
  else
    {
      scalar sum = (vi->second)+a;
      if(is_zero(sum)) entries.erase(vi);
      else (vi->second)=sum;
    }
}

void svec::sub(int i, const scalar& a)   // subtracts a from i'th entry
{
  if(is_zero(a)) return;
  add(i, -a);
}

void svec::add_mod_p(int i, const scalar& aa, const scalar& p)
{
  scalar a=mod(aa,p);
  if(is_zero(a)) return;
  auto vi = entries.find(i);
  if(vi==entries.end())
    entries[i]=a;
  else
    {
      scalar sum = mod((vi->second)+a,p);
      if(is_zero(sum)) entries.erase(vi);
      else (vi->second)=sum;
    }
}

void svec::sub_mod_p(int i, const scalar& a, const scalar& p)
{
  add_mod_p(i, -a, p);
}

svec& svec::operator+=(const svec& w)
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
		scalar sum = (vi->second) + (wi->second);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first);}
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::operator-=(const svec& w)
{
  this -> operator+=(-w);
  return *this;
}

svec& svec::add_scalar_times(const svec& w, const scalar& a)
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
		scalar sum = (vi->second) + a* (wi->second);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first);}
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::operator*=(const scalar& scal)
{
  for ( auto& vi : entries)
    (vi.second)*=scal;
  return *this;
}

void svec::reduce_mod_p(const scalar& p)
{
  auto vi = entries.begin();
  while( vi != entries.end() )
    {
      scalar a = mod(vi->second,p);
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

svec& svec::mult_by_scalar_mod_p(const scalar& scal, const scalar& p)
{
  scalar s = mod(scal,p);
  if(!is_one(s))
    for( auto& vi : entries)
      (vi.second)=xmodmul(vi.second,s,p);
  return *this;
}

svec& svec::add_scalar_times_mod_p(const svec& w, const scalar& scal, const scalar& p)
{
  if (d!=w.d)
    {
      cerr << "Incompatible svecs in svec::add_scalar_times()"<<endl;
      return *this;
    }
  scalar a = mod(scal,p);
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
		scalar sum = xmod((vi->second) + xmodmul(a, (wi->second),p),p);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); }
		wi++;
	      }
	}
    }
  //  reduce_mod_p(p);
  return *this;
}

svec& svec::add_scalar_times_mod_p(const svec& w, const scalar& scal, std::set<int>& ons, std::set<int>& offs,
				   const scalar& p)
{
  ons.clear();
  offs.clear();
  if (d!=w.d)
    {
      cerr << "Incompatible svecs in svec::add_scalar_times()"<<endl;
      return *this;
    }
  scalar a = mod(scal,p);
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
		scalar sum = xmod((vi->second) + xmodmul(a, (wi->second),p),p);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); offs.insert(wi->first);}
		wi++;
	      }
	}
    }
  //  reduce_mod_p(p);
  return *this;
}

svec& svec::add_mod_p(const svec& w, const scalar& p)
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
		scalar sum = xmod((vi->second) + (wi->second),p);
		if(is_nonzero(sum)) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); }
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::sub_mod_p(const svec& w, const scalar& p)
{
  add_mod_p(-w, p);
  return *this;
}

svec& svec::operator/=(const scalar& scal)
{
  if(is_zero(scal))
    cerr<<"Attempt to divide svec by 0"<<endl;
  for( auto& vi : entries)
    (vi.second)/=scal;
  return *this;
}

// Definitions of non-member, friend operators and functions

int eqmodp(const svec& v1, const svec& v2, const scalar& p)
{
  if(v1.d!=v2.d) return 0;
  if (std::any_of(v1.entries.begin(), v1.entries.end(),
                  [v2,p] (const pair<int,scalar>& vi) {return xmod((vi.second)-(v2.elem(vi.first)),p)!=0;}))
    return 0;
  if (std::any_of(v2.entries.begin(), v2.entries.end(),
                  [v1,p] (const pair<int,scalar>& vi) {return xmod((vi.second)-(v1.elem(vi.first)),p)!=0;}))
    return 0;
  return 1;
}

int operator==(const svec& v1, const vec& v2)
{
  if(v1.d!=dim(v2)) return 0;
  for(int i=1; i<=v1.d; i++) if(v2[i]!=v1.elem(i)) return 0;
  return 1;
}

ostream& operator << (ostream& s, const svec& v)
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

scalar operator*(const svec& v, const svec& w) //dot prod
{
  scalar ans(0);
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

scalar operator*(const svec& v, const vec& w) //dot prod
{
  scalar ans(0);
  std::for_each(v.entries.cbegin(), v.entries.cend(),
                [&ans, w](auto vi){ans += (vi.second)* w[vi.first];});
  return ans;
}

scalar dotmodp(const svec& v, const vec& w, const scalar& pr)
{
  scalar ans(0);
  std::for_each(v.entries.cbegin(), v.entries.cend(),
                [&ans, w, pr](auto vi){ans=mod(ans+xmodmul(vi.second,w[vi.first],pr),pr);});
  return ans;
}

scalar dotmodp(const svec& v, const svec& w, const scalar& pr)
{
  scalar ans(0);
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

scalar content(const svec& v)
{
  scalar g(0);
  for( const auto & vi : v.entries)
    {
      g=gcd(g,vi.second);
      if (is_one(g))
        return g;
    }
  return g;
}

scalar make_primitive(svec& v) // divides by & returns content
{
  scalar c=content(v);
  v/=c;
  return c;
}
