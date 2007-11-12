// ptest.cc -- test program for arith functions
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
 
#include "arith.h"
#include <functional> 

int main()
{
 long n, p;

 cout<<"Initialized table of " << nprimes() << " primes, up to "<<maxprime() << endl;

 cout<<"Enter an index: ";  cin >> n;
 p=prime_number(n);  cout << "prime_number("<<n<<") = "<<p<<endl;

 cout<<"How many primes do you want to see (as a vector<long>)? ";  cin >> n;
 vector<long> plist=primes(n);
 cout << plist << endl;

 cout<<"Enter a number to see if it is that list: "; cin>>p;
 vector<long>::iterator pi=find(plist.begin(),plist.end(),p);
 if(pi==plist.end()) cout<<"NOT in the list"<<endl;
 else    cout<<p<<" is list item "<<(pi-plist.begin())<<" (counting from 0)"<<endl;

 vector<long> v(10);
 iota(v.begin(),v.end(),1);
 cout<<"iota(10): "<<v<<endl;

 n=2310*210*64*17;
 cout<<"n = "<<n<<endl;
 vector<long> exps;

 transform(plist.begin(),plist.end(),inserter(exps,exps.end()),
	   bind2nd(ptr_fun(val),n));
 cout<<"exps = "<<exps<<endl;

 vector<long> plist1=primes(10);
 cout<<"Comparing previous prime list with "<<plist1<<endl;;
 cout<<"First starts with second, 10 items from 0: "<<startswith(plist,plist1,10)<<endl;
 cout<<"First starts with second,  5 items from 5: "<<startswith(plist,plist1,5,5)<<endl;
 cout<<"Second starts with first, 10 items from 0: "<<startswith(plist1,plist,10)<<endl;
 cout<<"Second starts with first,  5 items from 5: "<<startswith(plist1,plist,5,5)<<endl;

 //initialize a prime iterator for n primes
 cout<<"How many primes do you want to see (one by one)? ";  cin >> n;
 for(primevar pr(n); pr.ok(); pr++)
   cout << "Prime number " << pr.index() << " = " << pr << endl;
 
 long m;
 while (cout << "\nEnter an integer m (0 to stop): ", cin >> m, m!=0) 
   {
     cout << "Smallest prime factor of " << m << " is " << primdiv(m) << endl;
     plist=pdivs(m);
     cout << "m has " << plist.size() << " prime divisors: " << plist << endl;
     cout << "with exponents: "; 
     for(vector<long>::const_iterator pr = plist.begin(); pr!=plist.end(); pr++)
       cout << *pr <<":"<<val(*pr,m) << "\t";
     cout<<endl;

     vector<long> dlist=alldivs(m,plist);
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

 long a,b;
 int cont=1;
 while(cont)
   {
     cout<<"Enter integers a b (0 0 to stop): ";
     cin>>a>>b;
     long g=gcd(a,b);
     if(g==0) cont=0;
     else
       {
	 cout<<"gcd = "<<g<<endl;
	 long l=lcm(a,b);
	 cout<<"lcm = "<<l<<endl;
       }
   }
}  /* main() */
