// templates.h:  some utility functions for vector<T> classes
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
 
template <class T >
inline ostream& operator<<(ostream& os, const vector<T>& v)
{
  os <<"[ ";
  copy(v.begin(),v.end(), ostream_iterator<T>(cout, " "));
  os << "]";
  return os;
}

template <class T >
inline void vec_out(ostream& os, const vector<T>& v, unsigned int n=0)
{
  unsigned int m=v.size();  bool trunc=0;
  if((n>0)&&(m>n)) {m=n; trunc=1;}
  os <<"[ ";
  copy(v.begin(),v.begin()+m, ostream_iterator<T>(cout, " "));
  if(trunc) os << "...";
  os << "]";
}

template <class T >
inline ostream& operator<<(ostream& os, const std::set<T>& v)
{
  os <<"{ ";
  copy(v.begin(),v.end(), ostream_iterator<T>(cout, " "));
  os << "}";
  return os;
}

// NB The following will only work properly, giving a result with no
// repeats, if the input vectors are sorted!
template <class T >
vector<T> vector_union(const vector<T>& a, const vector<T>& b)
{
  //    cout<<"merging "<<a<<" and "<<b<<" using set_union"<<endl;
  vector<T> c;
  set_union(a.begin(),a.end(),b.begin(),b.end(),inserter(c,c.end()));
  //    cout<<"result is "<<c<<endl;
  return c;
}

// the following returns true if the entries of a and b are equal from
// index "from" to index "from+l-1":

template <class T >
int startswith(const vector<T>& a, const vector<T>& b, long l, long from=0)
{
  return equal(a.begin()+from,a.begin()+from+l,b.begin()+from);
}

