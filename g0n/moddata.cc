// FILE MODDATA.CC: Implementation of member functions for class moddata
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

#include "moddata.h"

level::level(long n, long neigs)
{
  //cout<<"Creating a class level with n = " << n << endl;
  modulus=n; 
  plist=pdivs(n); npdivs=plist.size();
  dlist=posdivs(n); ndivs=dlist.size();
  nap=neigs;
  primelist=plist;
  primevar pr; long p; p0=0;
  while(primelist.size()<(unsigned)nap)
    {
      p=pr;
      if (ndiv(p,modulus)) 
	{
	  if(p0==0) p0=p;
	  primelist.push_back(p);
	}
      pr++;
    }
  sqfac=1;
  for(long ip=0; ip<npdivs; ip++) 
    {
      p = plist[ip];
      if(::div(p*p,n)) sqfac*=p;
    }
  long rootn=(long)(sqrt((double)n)+0.1); // rounded down
  squarelevel=(n==rootn*rootn);
}

long bezout_x(long aa, long bb, long& xx)
{long a,b,c,x,oldx,newx,q;
 oldx = 1; x = 0; a = aa; b = bb;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
  }
 if (a<0) {xx=-oldx; return(-a);}
 else     {xx= oldx; return( a);}
}
 
moddata::moddata(long n) :level(n)
{
  //   cout << "In constructor moddata::moddata.\n";
 long i,p,x,d,nd,nnoninv;
 phi=psi=modulus;
 for(i=0; i<npdivs; i++)
   {  p = plist[i];
      phi -= phi/p;
      psi += psi/p;
    }
 nsymb = psi;
 nsymb1 = 2*modulus-phi;
 nsymb2 = nsymb-nsymb1;
 invlist.resize(modulus);          //codes
 noninvlist.resize(modulus-phi);   //list of non-units
 noninvdlist.resize(modulus-phi);  //list of divisors for each nonunit
 gcdtable.resize(modulus);         //list of gcds mod N
 unitdivlist.resize(modulus);      //list of units s.t. u*res | N
 nnoninv=0;
 for (i=0; i<modulus; i++)            //set up codes
 { d = bezout_x(i,modulus,x);
   gcdtable[i]=d;
   if (d==1) {unitdivlist[i] = invlist[i] = reduce(x); }
   else
   {invlist[i]=-nnoninv;
    noninvlist[nnoninv]=i;
    noninvdlist[nnoninv]=-1;
    if (d<modulus)
    {
     for (nd=0; (nd<ndivs)&&(dlist[nd]!=d); nd++) ;
     noninvdlist[nnoninv]=nd;
    }
    nnoninv++;
    if(::gcd(x,modulus)!=1) // adjust x so coprime to N
      {
	long m=modulus/d, mm, mpower, mmold,u,v;
	mpower=mm=m; mmold=1;
	while(mm!=mmold) 
	  {
	    mpower=xmodmul(m,mpower,modulus); 
	    mmold=mm; 
	    mm=::gcd(mpower,modulus);
	  }

	bezout(mm,modulus/mm,u,v);
	// Must be careful of overflow!
	x = (x*v)%mm; 
	x = (x*(modulus/mm))%modulus;
	x = (x+(u*mm))%modulus;
      }
    unitdivlist[i]=x;
#if(0)
    //Check:  
    if(::gcd(x,modulus)!=1)
      {
	cout<<"Error:  unitdivlist["<<i<<"] = "<<x<<" is not coprime to "<<modulus<<endl;
	abort();
      }
    if(((i*x-d)%modulus)!=0)
      {
	cout<<"Error:  unitdivlist["<<i<<"] = "<<x<<" is wrong"<<endl;
	abort();
      }
#endif
   }
 }
 if (ndivs>0) {dstarts.reserve(ndivs);}
}

void moddata::display() const
{
 cout << "Level = " << modulus << "\n";
 cout << "Number of symbols = " << nsymb << "\n";
 cout << ndivs << " non-trivial divisors: " << dlist << endl;
 cout << npdivs << " prime divisors: " << plist << endl;
 cout << "invlist: " << invlist << endl;
 cout << "noninvlist: " << noninvlist << endl;
 cout << "noninvdlist: " << noninvdlist << endl;
 cout << "gcdtable: " << gcdtable << endl;
 cout << "unitdivlist: " << unitdivlist << endl;
}

