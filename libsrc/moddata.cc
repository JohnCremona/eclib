// FILE MODDATA.CC: Implementation of member functions for class moddata
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

#include "eclib/moddata.h"
#include "eclib/arith.h"

level::level(long n, long neigs)
  : N(n), plist(pdivs(n)), dlist(posdivs(n)), nap(neigs)
{
  //cout<<"Creating a class level with n = " << n << endl;
  npdivs=plist.size();
  ndivs=dlist.size();
  primelist=plist;
  primevar pr; long p; p0=0;
  while(primelist.size()<(unsigned)nap)
    {
      p=pr;
      if (ndivides(p,N)) 
	{
	  if(p0==0) p0=p;
	  primelist.push_back(p);
	}
      ++pr;
    }
  sqfac=1;
  for(long ip=0; ip<npdivs; ip++) 
    {
      p = plist[ip];
      if(::divides(p*p,n)) sqfac*=p;
    }
  long rootn=(long)(sqrt((double)n)+0.1); // rounded down
  squarelevel=(n==rootn*rootn);
}

long bezout_x(long aa, long bb, long& xx)
{long a=aa, b=bb, x=0, oldx = 1;
 while (b!=0)
 { long q = a/b;
   long c    = a    - q*b; a    = b; b = c;
   long newx = oldx - q*x; oldx = x; x = newx;
  }
 if (a<0) {xx=-oldx; return(-a);}
 else     {xx= oldx; return( a);}
}

moddata::moddata(long n) :level(n)
{
  //   cout << "In constructor moddata::moddata.\n";
 long x,nd,nnoninv;
 phi=psi=N;
 for(long i=0; i<npdivs; i++)
   {  long p = plist[i];
      phi -= phi/p;
      psi += psi/p;
    }
 nsymb = psi;
 nsymb1 = 2*N-phi;
 nsymb2 = nsymb-nsymb1;
 invlist.resize(N);          //codes
 noninvlist.resize(N-phi);   //list of non-units
 noninvdlist.resize(N-phi);  //list of divisors for each nonunit
 gcdtable.resize(N);         //list of gcds mod N
 unitdivlist.resize(N);      //list of units s.t. u*res | N
 nnoninv=0;
 for (long i=0; i<N; i++)            //set up codes
 { long d = bezout_x(i,N,x);
   gcdtable[i]=d;
   if (d==1) {unitdivlist[i] = invlist[i] = reduce(x); }
   else
   {invlist[i]=-nnoninv;
    noninvlist[nnoninv]=i;
    noninvdlist[nnoninv]=-1;
    if (d<N)
    {
     for (nd=0; (nd<ndivs)&&(dlist[nd]!=d); nd++) ;
     noninvdlist[nnoninv]=nd;
    }
    nnoninv++;
    if(::gcd(x,N)!=1) // adjust x so coprime to N
      {
	long m=N/d, mm, mpower, mmold,u,v;
	mpower=mm=m; mmold=1;
	while(mm!=mmold) 
	  {
	    mpower=xmodmul(m,mpower,N); 
	    mmold=mm; 
	    mm=::gcd(mpower,N);
	  }

	bezout(mm,N/mm,u,v);
	// Must be careful of overflow!
	x = (x*v)%mm; 
	x = (x*(N/mm))%N;
	x = (x+(u*mm))%N;
      }
    unitdivlist[i]=x;
   }
 }
 if (ndivs>0) {dstarts.resize(ndivs);}
}

void moddata::display() const
{
 cout << "Level = " << N << "\n";
 cout << "Number of symbols = " << nsymb << "\n";
 cout << ndivs << " non-trivial divisors: " << dlist << endl;
 cout << npdivs << " prime divisors: " << plist << endl;
 cout << "invlist: " << invlist << endl;
 cout << "noninvlist: " << noninvlist << endl;
 cout << "noninvdlist: " << noninvdlist << endl;
 cout << "gcdtable: " << gcdtable << endl;
 cout << "unitdivlist: " << unitdivlist << endl;
}

string of_filename(long n, char c)
{
  stringstream s;
  s << getenv_with_default("OF_DIR","./newforms");
  s  << "/" << c << n;
  return s.str();
}

string nf_filename(long n, char c)
{
  stringstream s;
  s << getenv_with_default("NF_DIR","./newforms");
  s  << "/" << c << n;
  return s.str();
}

string small_nf_filename(long n, char c)
{
  stringstream s;
  s << getenv_with_default("SNF_DIR","./smallnf");
  s  << "/" << c << n;
  return s.str();
}
