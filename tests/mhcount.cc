// FILE MHCOUNT.CC: Program to list/count Manin-Heilbronn matrices
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
//
#include <eclib/arith.h>

void outmat(int n, long a, long b, long c, long d)
{
  cout<<"mats["<<n<<"]=mat22("<<a<<","<<b<<","<<c<<","<<d<<");\n";
}

int main(void)
{
  int seemats=0; int np=10;
  cout << "See the matrices? "; cin>>seemats;
  cout << "How many primes? "; cin >> np;
  for(primevar pr(np); pr.ok(); pr++)
    {
      long p = (long)pr; if(p==2) continue;
      long p2 = (p-1)/2;
      int nmats=0;
      cout << "p = " << p << ";\t"; 
      if(seemats)
        {
          cout << "\n";
          outmat(nmats,1,0,0,p);
        }
      nmats++;
      if(seemats)
        outmat(nmats,p,0,0,1);
      nmats++;
      for(int s=1; s>-2; s-=2)
      for(long r=1; r<=p2; r++)
	{
//	  cout<<"r = " << s*r << ":" << endl;
	  long x1=p, x2=-s*r, y1=0, y2=1, a=-p, b=s*r;
	  if(seemats)
            outmat(nmats,x1,x2,y1,y2);
	  nmats++;
	  while(b!=0)
	    {
	      long c=mod(a,b);
              long q=(a-c)/b;
	      long x3=q*x2-x1;
              long y3=q*y2-y1;
	      a=-b; b=c; x1=x2; x2=x3; y1=y2; y2=y3;
	      if(seemats)
                outmat(nmats,x1,x2,y1,y2);
	      nmats++;
	    }
	}
      if(seemats)
        cout<<"\n";
      cout << "nmats = " << nmats << ";" << endl;
    }
}


