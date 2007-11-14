// mptest.cc -- test program for marith functions
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
 //

#include "marith.h"
#define MAXPRIME 1000000

int main()
{
  initprimes("PRIMES",1);

  cout<<"long prime factors of 420:\n";
  vector<long> lplist = pdivs(420);
  cout<<lplist<<endl;
  cout<<" with posdivs (long version)\n";
  vector<long> lqlist = posdivs(420, lplist);
  cout<<lqlist<<endl;
  bigint num; num=420;
  vector<bigint> a(3);
  a[0]=10; a[1]=20; a[2]=30;
  vector<bigint> b = a;
  cout << "Elements of a are: " << a[0]<<"\t"<<a[1]<<"\t"<<a[2]<<endl;
  cout << "Elements of b are: " << b[0]<<"\t"<<b[1]<<"\t"<<b[2]<<endl;
  cout << "a = " << a << endl << "end of output of a" << endl;
  cout << "b = " << b << endl << "end of output of b" << endl;
  b=a;
  cout<<"After b=a, b = " << endl << b << endl;

  cout<<"bigint versions of divisor functions:\n\n";
  cout<<"bigint prime factors of 420:\n";
  vector<bigint> iplist = pdivs(num);
  cout<<iplist<<endl;
  cout<<" with posdivs (bigint version) \n";
  vector<bigint> iqlist = posdivs(num, iplist);
  cout<<iqlist<<endl;

  cout<<"and alldivs (bigint version) \n";
  //deliberately reuse space
  iqlist = alldivs(num, iplist);
  cout<<iqlist<<endl;

  cout<<"making a copy of iqlist\n";
  iplist=iqlist;
  cout<<" reuse iqlist for sqdivs\n";
  iqlist = sqdivs(num);
  cout<<"iqlist (should be sqdivs) "<<iqlist<<endl;
  cout<<"iplist (should be alldivs) "<<iplist<<endl;

  cout<<"making an vector<bigint> of iplist (alldivs)\n";
  vector<bigint> irary(iplist);
  cout<<"stream output: "<<irary<<endl;
  cout<<"irary(7) is "<<irary[7]<<"\tirary[3] is "<<irary[3]<<endl;

  cout<<"testing find function\n";
  vector<bigint>::iterator vi;
  vi = find(irary.begin(),irary.end(),35);
  if(vi==irary.end()) cout<<"35 is not there\n";
  else cout<<"35 is there:  "<<*vi<<" is item number "<<(vi-irary.begin())<<endl;
  vi = find(irary.begin(),irary.end(),13);
  if(vi==irary.end()) cout<<"13 is not there\n";
  else cout<<"13 is there:  "<<*vi<<" is item number "<<(vi-irary.begin())<<endl;

  cout<<"\n\nTest of sqrt and isqrt\n";
  bigint astop; astop=999;
  bigint aaa,roota; int res; 
  while (cout << "\nEnter a positive bigint a (999 to stop): ", cin >> aaa, aaa!=astop) 
   {
     res = sign(aaa);
     cout << "a = " << aaa << ", sign(a) = " << res << "\n";
     roota = sqrt(aaa);
     cout << "a = " << aaa << ", sqrt(a) = " << roota << " (rounded down)\n";
     res = isqrt(aaa,roota);
     if(res) cout << "a is a square with exact square root " << roota << "\n";
     else cout << "a is not a square\n";
   }

  cout<<"\n\nTest of ressol (sqrt mod p)\n";
  bigint bb,p,r;
  while(cout << "Enter a prime p: ", cin>>p, cout<<p<<endl, is_positive(p))
    {
      int nbad=0; long x;
      for(x=1; x<100; x++)
	{
	  bb=x;
	  if(legendre(bb,p)==1)
	    {
	      bb%=p; if(bb<0) bb+=p;
	      ressol(r,bb,p);
	      cout << "sqrt("<<x<<" mod p) = " << r << "\t";
	      if(div(p,(r*r-x))) cout << "---OK\n";
	      else {cout << "---WRONG\n"; nbad++;}
	    }
	  else
	    {
	      cout << x<< " is not a quadratic residue mod p!\n";
	    }
	}
      if(nbad==0) cout << "First 100 OK"<<endl;
      else cout << nbad << " wrong ones out of first 100"<<endl;
    }

 bigint m;
 while (cout << "\nEnter a bigint m (0 to stop): ", cin >> m, sign(m)!=0) 
   {
     cout << "m = " << m << endl;
     vector<bigint> plist=pdivs(m);
     cout << "m has " << plist.size() << " prime divisors: " << plist << endl;
     cout << "with exponents: "; 
     for(vector<bigint>::const_iterator pr = plist.begin(); pr!=plist.end(); pr++)
       cout << *pr <<":"<<val(*pr,m) << "\t";
     cout<<endl;

     vector<bigint> dlist=alldivs(m,plist);
     cout << "m has " << dlist.size() << " divisors: " << dlist << endl;
     dlist = posdivs(m,plist);
     cout << "m has " << dlist.size() 
          << " positive divisors: " << dlist << endl;
     dlist = sqdivs(m,plist);
     cout << "m has " << dlist.size() 
          << " positive divisors whose square divides m: " << dlist << endl;
     dlist = sqfreedivs(m,plist);
     cout << "m has " << dlist.size() 
          << " positive square-free divisors: " << dlist << endl;
 }

 cout<<endl;

}  // end of main    
