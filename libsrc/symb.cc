// FILE SYMB.CC: Implementations for symbols
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

#include "eclib/symb.h"

// Friend of class symb:
ostream& operator<< (ostream& s, const symb& sy)
{
   s << "(" << sy.c << ":" << sy.d << ")";
   return s;
}

//#define DEBUG_NORMALIZE
symb symb::normalize() const
{
#ifdef DEBUG_NORMALIZE
  cout<<"Normalizing symbol "<<(*this)<<endl;
#endif
 long n=N->N;
 long u=N->unitdiv(c);
#ifdef DEBUG_NORMALIZE
  cout<<"scaling by u =  "<<u<<endl;
#endif
 long cc=N->reduce(xmodmul(c,u,n));
 long dd=N->reduce(xmodmul(d,u,n))%(n/cc);
#ifdef DEBUG_NORMALIZE
  cout<<"new c =  "<<cc<<endl;
  cout<<"new d =  "<<dd<<endl;
#endif
 symb ans(cc,dd,N);
#ifdef DEBUG_NORMALIZE
  cout<<"Returning normalized symbol "<<ans;
  int ok = (ans==(*this)) && ::div(ans.cee(),n);
  if(ok) cout<<" ok"; else cout<<" wrong!";
  cout<<endl;
#endif
  return ans;
}

// Constructor for modsym, converting from symb:
modsym::modsym(const symb& s)
{
 long c,d,h,x,y;
 c = s.cee(); d = s.dee();
 h = bezout(c , d, x, y);
 a=rational(-x , d/h);
 b=rational( y , c/h);
}

// Friend of class modsym:
ostream& operator<< (ostream& s, const modsym& m)
{
   s << "{" << (m.a) << "," << (m.b) << "}";
   return s;
}

//Members of class symblist:

symblist::symblist(long n)
{
  maxnum=n;
  num=0;
  list.resize(n);
}

void symblist::add(const symb& s, long start)
{
 if (index(s,start)==-1)
 {
  if (num<maxnum)
    {
      list[num]=s;
      long c = s.cee(), d=posmod(s.dee(),s.level()/c);
      hashtable[pair<long,long>(c,d)]=num;
      num++;
      //      cout<<"Adding symbol "<<s<<" as special number "<<num<<endl;
    }
  else
    {
      cerr << "Error in symblist::add: attempt to add too many symbols to list!"<<endl;
    }
 }
}

long symblist::index(const symb& s, long start) const
{
  // cout<<"index of "<<s;
 symb ss = s.normalize();
 long c = ss.cee(), d=ss.dee();
 auto j = hashtable.find(pair<long,long>(c,d));
 if(j==hashtable.end())
   return -1;
 // cout<<" is "<<j->second<<endl;
 return j->second;
}


symb symblist::item(long n) const
{
 if ((n>num)||(n<0))
   {
     cerr<<"Error in symblist::item: index out of range!"<<endl;
     return symb();
   }
 else return list[n];
}

//Member functions for class symbdata:
symbdata::symbdata(long n) :moddata(n),specials(nsymb2)
{
  //   cout << "In constructor symbdata::symbdata.\n";
  //   cout << "nsymb2 = " << nsymb2 << "\n";
 if (nsymb2>0)
 { long ic,id,d; symb s;
//N.B. dlist include d=1 at 0 and d=mod at end, which we don't want here
   for (ic=1; (ic<ndivs-1)&&(specials.count()<nsymb2); ic++)
   { long c=dlist[ic];
     long start=specials.count();
     dstarts[ic]=start;
     for (id=1; (id<N-phi)&&(specials.count()<nsymb2); id++)
     { d = noninvlist[id];
       if (::gcd(d,c)==1)
       {  s = symb(c,d,this);
          specials.add(s,start);     //only adds it if not there already!
       }
     }     // end of d loop
    }      // end of c loop
   if (specials.count()<nsymb2)
     {
       cout << "Problem: makesymbols found only " << specials.count() << " symbols ";
       cout << "out of " << nsymb2 << endl;
     }
   //   cout << "Special symbols: "; specials.display();
 }
}

long symbdata::index2(long c, long d) const
{ long kd = code(d);
  // cout<<"index2("<<c<<":"<<d<<"):"<<endl;
  if (kd>0)                // d invertible, with inverse kd
    {
      //      cout<<"d = "<<d<<" has code kd = "<<kd<< " so reducing c*kd = "<<c*kd<<" mod "<<N<<" --> "<<xmodmul(c,kd,N)<<"-->"<<reduce(xmodmul(c,kd,N))<<endl;
      return reduce(xmodmul(c,kd,N));   // (c:d) = (c*kd:1)
    }
  else
  { long kc = code(c);
    if (kc>0)              // (c:d) = (1:kc*d) if c invertible
      {
        //        cout<<"c = "<<c<<" has code kc = "<<kc<< " so reducing kc*d = "<<kc*d<<" mod "<<N<<" --> "<<xmodmul(kc,d,N)<<"-->"<<reduce(xmodmul(kc,d,N))<<endl;
        return   N-code(xmodmul(kc,d,N));
      }
    else
    {
     long start = dstarts[noninvdlist[-kc]];
     symb s(c,d,this);
     //     cout<<"About to compute index of symbol "<<s<<" from start = "<<start<<" in list of specials, size "<<specials.count() << endl;
     long ind = specials.index(s,start);
     if(ind<0)
       {
	 cout<<"error in index(): symbol "<<s<<" not in list!"<<endl;
       }
     return nsymb1+ind;
    }
  }
}

symb symbdata::symbol(long i) const
{ if (i<N) return symb(i,1,this);
  else if (i<nsymb1) return symb(1,noninvlist[i-N],this);
 else return specials[i-nsymb1]; // specials.item[i-nsymb1];
}

void symbdata::display() const
{ moddata::display();
  cout << "Number of special symbols = " << nsymb2 << "\n";
  specials.display();
}

void symbdata::check(void) const
{
  int ok=1; symb s;
  for (long i=0; i<nsymb; i++)
    {
      long j = index(s=symbol(i));
      if (i!=j)
        {
          cout << i << "-->" << s << "-->" << j << "\n";
          ok=0;
        }
    }
  if (ok) cout << "symbols check OK!\n";
  else cout << "symbols check found errors!\n";
}

modsym jumpsymb(const symb& s1, const symb& s2)
{
  //Assuming s1==s2, returns closed modular symbol {g1(0),g2(0)} where gi<->si
  long c1=s1.cee(), c2=s2.cee(), d1=s1.dee(), d2=s2.dee();
  return modsym(rational(-invmod(c1,d1),d1),rational(-invmod(c2,d2),d2));
}
