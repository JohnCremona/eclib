// svec.cc: implementation of sparse integer vector class svec
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
 
// ONLY to be included by svector.cc

inline scalar xmm(scalar a, scalar b, scalar m)
{
  //return xmodmul(a,b,m);
  //return (a*b) % m;
  return (a*(int64_t)b) % m;
  //return (scalar)(((long)a*(long)b) % (long)m);
}
inline scalar xmm0(scalar a, scalar b)
{
  //return xmodmul(a,b,m);
  //return (a*b) % m;
  return (a*(int64_t)b) % BIGPRIME;
  //return (scalar)(((long)a*(long)b) % (long)m);
}



// Definitions of member operators and functions:

svec::svec(const vec& v)
{
  d = dim(v); scalar vi;
  for(int i = 1; i <= d; i++ )
    if((vi=v[i])) // assignment not equality!
      entries[i]=vi;
}

svec::svec (int dim, scalar* a)     // conversion constructor
  :d(dim)
{
  int i=0; scalar* ai=a;
  for( ; i<d; i++, ai++)
    if(*ai) 
      entries[i]=*ai; 
}

vec svec::as_vec( ) const
{
  vec v(d); // initializes to 0
  for(map<int,scalar>::const_iterator wi=entries.begin(); 
      wi!=entries.end(); wi++)
    v[wi->first] = wi->second;
  return v;
}

scalar svec::elem(int i)  const   // returns i'th entry
{
  map<int,scalar>::const_iterator vi = entries.find(i);
  if(vi==entries.end()) return 0;
  return vi->second;
}

void svec::set(int i, scalar a) 
{
  if(a) entries[i]=a;
#if(0) // 
  else // a==0 so kill the i'entry if it is there
    { 
      map<int,scalar>::iterator vi = entries.find(i);   
      if(vi!=entries.end()) entries.erase(vi);
    }
#endif
}

void svec::erase(int i)
{
  map<int,scalar>::iterator vi = entries.find(i);   
  if(vi==entries.end()) 
    {
      cout<<"Error in svec::erase(): cannot delete missing entry #"<<i
	  <<" from v = "<<(*this)<<endl; 
      abort();
    } 
  else entries.erase(vi);
}

std::set<int> svec::support() const
{
  std::set<int> ans;
  map<int,scalar>::const_iterator vi;
  for(vi=entries.begin(); vi!=entries.end(); vi++) 
    ans.insert(vi->first);
  return ans;
}

void svec::add(int i, scalar a)   // adds a to i'th entry
{
  if(!a) return;
  map<int,scalar>::iterator vi = entries.find(i);   
  if(vi==entries.end()) 
    entries[i]=a;
  else
    {
      scalar sum = (vi->second)+a;
      if(sum==0) entries.erase(vi);
      else (vi->second)=sum;
    }
}

void svec::sub(int i, scalar a)   // subtracts a from i'th entry
{
  if(!a) return;
  map<int,scalar>::iterator vi = entries.find(i);   
  if(vi==entries.end()) 
    entries[i]=-a;
  else
    {
      scalar sum = (vi->second)-a;
      if(sum==0) entries.erase(vi);
      else (vi->second)=sum;
    }
}

void svec::add_mod_p(int i, scalar a, const scalar& p)
{
  a=xmod(a,p);
  if(!a) return;
  map<int,scalar>::iterator vi = entries.find(i);   
  if(vi==entries.end()) 
    entries[i]=a;
  else
    {
      scalar sum = xmod((vi->second)+a,p);
      if(sum==0) entries.erase(vi);
      else (vi->second)=sum;
    }
}
void svec::sub_mod_p(int i, scalar a, const scalar& p)
{
  a=xmod(-a,p);
  if(!a) return;
  map<int,scalar>::iterator vi = entries.find(i);   
  if(vi==entries.end()) 
    entries[i]=a;
  else
    {
      scalar sum = xmod((vi->second)+a,p);
      if(sum==0) entries.erase(vi);
      else (vi->second)=sum;
    }
}

svec& svec::operator+=(const svec& w)
{
  if (d!=w.d)
    {
      cout << "Incompatible svecs in svec::operator+=()\n";
      abort();
      return *this;
    }
  map<int,scalar>::const_iterator  wi=w.entries.begin();
  map<int,scalar>::iterator vi=entries.begin();
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
		if(sum) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first);}
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::operator-=(const svec& w)
{
  if (d!=w.d)
    {
      cout << "Incompatible svecs in svec::operator-=()\n";
      abort();
      return *this;
    }
  map<int,scalar>::const_iterator  wi=w.entries.begin();
  map<int,scalar>::iterator vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=-wi->second;
	      wi++;
	    } 	    
	}
      else
	{
	  if((vi->first)<(wi->first))  {vi++;} 
	  else
	    if((wi->first)<(vi->first)) 
	      {
		entries[wi->first]=-wi->second;
		wi++;
	      } 
	    else
	      {
		scalar diff = (vi->second) - (wi->second);
		if(diff) {vi->second = diff; vi++;}
		else {vi++; entries.erase(wi->first);}
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::add_scalar_times(const svec& w, scalar a)
{
  if (d!=w.d)
    {
      cout << "Incompatible svecs in svec::add_scalar_times()\n";
      abort();
      return *this;
    }
  if(a==0) return *this;
  map<int,scalar>::const_iterator  wi=w.entries.begin();
  map<int,scalar>::iterator vi=entries.begin();
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
		if(sum) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first);}
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::operator*=(scalar scal)
{
  //  if(scal==0) cout<<"Attempt to multiply svec by 0\n"<<endl;
  for( map<int,scalar>::iterator vi=entries.begin(); 
       vi != entries.end(); vi++)
    (vi->second)*=scal;
  return *this;
}

void svec::reduce_mod_p(const scalar& p)
{
  scalar a;
  map<int,scalar>::iterator vi, vi2;
  for( vi=entries.begin(); vi != entries.end(); )
    {
      a = mod(vi->second,p);
      if(a) {(vi->second)=a; vi++;}
      else {vi2=vi; vi++; entries.erase(vi2->first);}
    }
}

svec& svec::mult_by_scalar_mod_p(scalar scal, const scalar& p)
{
  //  if(xmod(scal,p)==0) cout<<"Attempt to multiply svec by 0\n"<<endl;
  if(scal!=1)
    for( map<int,scalar>::iterator vi=entries.begin(); 
	 vi != entries.end(); vi++)
      (vi->second)=xmm(vi->second,scal,p);
  return *this;
}

svec& svec::add_scalar_times_mod_p(const svec& w, scalar a, const scalar& p)
{
  if (d!=w.d)
    {
      cout << "Incompatible svecs in svec::add_scalar_times()\n";
      abort();
      return *this;
    }
  if(a==0) return *this;
  map<int,scalar>::const_iterator  wi=w.entries.begin();
  map<int,scalar>::iterator vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=xmm(a,(wi->second),p);
	      wi++;
	    } 	    
	}
      else
	{
	  if((vi->first)<(wi->first)) {vi++;} 
	  else
	    if((wi->first)<(vi->first)) 
	      {
		entries[wi->first]=xmm(a,(wi->second),p);
		wi++;
	      } 
	    else
	      {
		scalar sum = xmod((vi->second) + xmm(a, (wi->second),p),p);
		if(sum) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); }
		wi++;
	      }
	}
    }
  //  reduce_mod_p(p);
  return *this;
}

svec& svec::add_scalar_times_mod_p(const svec& w, scalar a, std::set<int>& ons, std::set<int>& offs, 
				   const scalar& p)
{
  ons.clear();
  offs.clear();
  if (d!=w.d)
    {
      cout << "Incompatible svecs in svec::add_scalar_times()\n";
      abort();
      return *this;
    }
  if(a==0) return *this;
  map<int,scalar>::const_iterator  wi=w.entries.begin();
  map<int,scalar>::iterator vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=xmm(a,(wi->second),p);
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
		entries[wi->first]=xmm(a,(wi->second),p);
		ons.insert(wi->first);
		wi++;
	      } 
	    else
	      {
		scalar sum = xmod((vi->second) + xmm(a, (wi->second),p),p);
		if(sum) {vi->second = sum; vi++;}
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
      cout << "Incompatible svecs in svec::add_scalar_times()\n";
      abort();
      return *this;
    }
  map<int,scalar>::const_iterator  wi=w.entries.begin();
  map<int,scalar>::iterator vi=entries.begin();
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
		if(sum) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); }
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::sub_mod_p(const svec& w, const scalar& p)
{
  if (d!=w.d)
    {
      cout << "Incompatible svecs in svec::add_scalar_times()\n";
      abort();
      return *this;
    }
  map<int,scalar>::const_iterator  wi=w.entries.begin();
  map<int,scalar>::iterator vi=entries.begin();
  while(wi!=w.entries.end())
    {
      if(vi==entries.end())
	{
	  while(wi!=w.entries.end())
	    {
	      entries[wi->first]=-(wi->second);
	      wi++;
	    } 	    
	}
      else
	{
	  if((vi->first)<(wi->first)) {vi++;} 
	  else
	    if((wi->first)<(vi->first)) 
	      {
		entries[wi->first]=-(wi->second);
		wi++;
	      } 
	    else
	      {
		scalar sum = xmod((vi->second) - (wi->second),p);
		if(sum) {vi->second = sum; vi++;}
		else {vi++; entries.erase(wi->first); }
		wi++;
	      }
	}
    }
  return *this;
}

svec& svec::operator/=(scalar scal)
{
  if(scal==0) 
    {
      cout<<"Attempt to divide svec by 0\n"<<endl;
      abort();
    }
  for( map<int,scalar>::iterator vi=entries.begin(); 
       vi != entries.end(); vi++)
    (vi->second)/=scal;
  return *this;
}

// Definitions of non-member, friend operators and functions

int eqmodp(const svec& v1, const svec& v2, const scalar& p)
{
  if(v1.d!=v2.d) return 0;
  map<int,scalar>::const_iterator vi;
  for(vi=v1.entries.begin(); vi!=v1.entries.end(); vi++)
    if(xmod((vi->second)-(v2.elem(vi->first)),p)!=0) return 0;
  for(vi=v2.entries.begin(); vi!=v2.entries.end(); vi++)
    if(xmod((vi->second)-(v1.elem(vi->first)),p)!=0) return 0;
  return true;
}

int operator==(const svec& v1, const vec& v2)
{
  if(v1.d!=dim(v2)) return 0;
  for(int i=1; i<=v1.d; i++) if(v2[i]!=v1.elem(i)) return 0;
  return 1;
}

ostream& operator << (ostream& s, const svec& v)
{
  map<int,scalar>::const_iterator vi;
  s<<"[";
  for(vi=v.entries.begin(); vi!=v.entries.end(); vi++)
    {
      if(vi!=v.entries.begin()) s<<",";
      s << "("<<vi->first<<":"<<vi->second<<")";
    }
  s<<"]";
  return s;
}

scalar operator*(const svec& v, const svec& w) //dot prod
{
  scalar ans=0;
  if((v.entries.size()==0)||(w.entries.size()==0)) return ans;
  map<int,scalar>::const_iterator vi=v.entries.begin(), wi=w.entries.begin();
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
  scalar ans=0;
  if((v.entries.size()==0)) return ans;
  map<int,scalar>::const_iterator vi;
  for(vi=v.entries.begin(); vi!=v.entries.end(); vi++)
    ans+=(vi->second)* w[vi->first];
  return ans;
}

scalar dotmodp(const svec& v, const vec& w, scalar pr)
{
  scalar ans=0;
  map<int,scalar>::const_iterator vi;
  for(vi=v.entries.begin(); vi!=v.entries.end(); vi++)
    ans=xmod(ans+xmm(vi->second,w[vi->first],pr),pr);
  return ans;
}

scalar dotmodp(const svec& v, const svec& w, scalar pr)
{
  scalar ans=0;
  if((v.entries.size()==0)||(w.entries.size()==0)) return ans;
  map<int,scalar>::const_iterator vi=v.entries.begin(), wi=w.entries.begin();
  while((vi!=v.entries.end())&&(wi!=w.entries.end()))
    {
      if((vi->first)<(wi->first)) {vi++;} else
	if((wi->first)<(vi->first)) {wi++;} else
	  {
	    ans=xmod(ans+xmm(vi->second,wi->second,pr),pr);
	    vi++; wi++;
	  }
    }
  return ans;
}

scalar content(const svec& v)
{
  scalar ans=0;
  map<int,scalar>::const_iterator vi;
  for(vi=v.entries.begin(); (ans!=1)&&(vi!=v.entries.end()); vi++)
    ans=gcd(ans,vi->second);
  return ans;
}

scalar make_primitive(svec& v) // divides by & returns content
{
  scalar c=content(v);
  v/=c;
  return c;
}
