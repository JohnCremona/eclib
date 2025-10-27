// mptest.cc -- test program for marith functions
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

#include <eclib/marith.h>
#define MAXPRIME 1000000

// test function for divisor_iterator
int test_divisor_iterator(const ZZ& N)
{
  if (is_zero(N)) return 1;
  vector<ZZ> divs1 = posdivs(N);
  vector<ZZ> divs2;
  divisor_iterator divN(N);
  while (divN.is_ok())
    {
      divs2.push_back(divN.value());
      divN.increment();
    }
  std::sort(divs1.begin(), divs1.end());
  std::sort(divs2.begin(), divs2.end());
  if (divs1!=divs2)
    {
      cout << "divisors are "<< divs1 <<endl;
      cout << "new list is  "<< divs2 <<endl;
    }
  return divs1==divs2;
}

int main()
{
  initprimes("PRIMES",1);

  cout<<"long prime factors of 420:\n";
  vector<long> lplist = pdivs(420);
  cout<<lplist<<endl;
  cout<<" with posdivs (long version)\n";
  vector<long> lqlist = posdivs(420, lplist);
  cout<<lqlist<<endl;
  ZZ num; num=420;
  vector<ZZ> a(3);
  a[0]=10; a[1]=20; a[2]=30;
  vector<ZZ> b = a;
  cout << "Elements of a are: " << a[0]<<"\t"<<a[1]<<"\t"<<a[2]<<endl;
  cout << "Elements of b are: " << b[0]<<"\t"<<b[1]<<"\t"<<b[2]<<endl;
  cout << "a = " << a << endl << "end of output of a" << endl;
  cout << "b = " << b << endl << "end of output of b" << endl;
  b=a;
  cout<<"After b=a, b = " << endl << b << endl;

  cout<<"ZZ versions of divisor functions:\n\n";
  cout<<"ZZ prime factors of 420:\n";
  vector<ZZ> iplist = pdivs(num);
  cout<<iplist<<endl;
  cout<<" with posdivs (ZZ version) \n";
  vector<ZZ> iqlist = posdivs(num, iplist);
  cout<<iqlist<<endl;

  cout<<"and alldivs (ZZ version) \n";
  //deliberately reuse space
  iqlist = alldivs(num, iplist);
  cout<<iqlist<<endl;

  cout<<"making a copy of iqlist\n";
  iplist=iqlist;
  cout<<" reuse iqlist for sqdivs\n";
  iqlist = sqdivs(num);
  cout<<"iqlist (should be sqdivs) "<<iqlist<<endl;
  cout<<"iplist (should be alldivs) "<<iplist<<endl;

  cout<<"making an vector<ZZ> of iplist (alldivs)\n";
  vector<ZZ> irary(iplist);
  cout<<"stream output: "<<irary<<endl;
  cout<<"irary(7) is "<<irary[7]<<"\tirary[3] is "<<irary[3]<<endl;

  ZZ radN = radical(num*num), N1, N2;
  cout << "\nRadical of " << num << "^2 = " << num*num <<" is " << radN << endl;
  sqfdecomp(num, N1, N2, iplist);
  cout << "Square-free decomposition of " << num << " is N1*N2^2 with N1 = " << N1 << " and N2 = " << N2 << endl;

  cout<<"\ntesting find function for the list " << irary << ":\n";
  auto vi = find(irary.begin(),irary.end(),35);
  if(vi==irary.end()) cout<<"35 is not there\n";
  else cout<<"35 is there:  "<<*vi<<" is item number "<<(vi-irary.begin())<<endl;
  vi = find(irary.begin(),irary.end(),13);
  if(vi==irary.end()) cout<<"13 is not there\n";
  else cout<<"13 is there:  "<<*vi<<" is item number "<<(vi-irary.begin())<<endl;

  cout << "\nTesting list multiplication functions\n";
  vector<ZZ> powers_of_2 = powers(ZZ(2), 10);
  cout << "Powers of 2 up to 2^10: " << powers_of_2 << endl;
  //  vector<int> odds = {1,3,5,7,9};
  vector<ZZ> odd_powers_of_2 = powers(ZZ(2), {1,3,5,7,9});
  cout << "Odd powers of 2 up to 2^10: " << odd_powers_of_2 << endl;

  vector<ZZ> L = ZZify({1,2,3,4,5});
  cout << "L   = " << L << endl;
  vector<ZZ> Lx3 = multiply_list(3, L);
  cout << "3*L = " << Lx3 << endl;

  vector<ZZ> L2 = ZZify({1,10,100});
  cout << "L2   = " << L2 << endl;
  vector<ZZ> LxL2 = multiply_lists(L,L2);
  cout << "L*L2 = " << LxL2 << endl;

  vector<int> ee = {0,1,2,3};
  vector<ZZ> L2e = multiply_list_by_powers(2, ee, L);
  cout << "L*[2^e for e in "<<ee<<"]: "<< L2e << endl;

  cout<<"\n\nTest of sqrt and isqrt\n";
  ZZ astop; astop=999;
  ZZ aaa,roota;
  while (cout << "\nEnter a positive ZZ a (999 to stop): ", cin >> aaa, aaa!=astop) 
   {
     int res = sign(aaa);
     cout << "a = " << aaa << ", sign(a) = " << res << "\n";
     roota = sqrt(aaa);
     cout << "a = " << aaa << ", sqrt(a) = " << roota << " (rounded down)\n";
     res = isqrt(aaa,roota);
     if(res) cout << "a is a square with exact square root " << roota << "\n";
     else cout << "a is not a square\n";
   }

  cout<<"\n\nTest of sqrt mod p\n";
  ZZ bb,p,r;
  while(cout << "Enter a prime p: ", cin>>p, cout<<p<<endl, is_positive(p))
    {
      int nbad=0; long x;
      for(x=1; x<100; x++)
	{
	  bb=x;
	  if(legendre(bb,p)==1)
	    {
	      bb%=p; if(bb<0) bb+=p;
	      sqrt_mod_p(r,bb,p);
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

 ZZ m;
 while (cout << "\nEnter a ZZ m (0 to stop): ", cin >> m, sign(m)!=0) 
   {
     cout << "m = " << m << endl;
     vector<ZZ> plist=pdivs(m);
     cout << "m has " << plist.size() << " prime divisors: " << plist << endl;
     vector<int> vals = valuations(m, plist);
     cout << "with exponents " << vals << endl;
     ZZ m2 = factorback(plist, vals);
     cout << "factorback recovers " << m2 << " (should be " << m << ")" << endl;

     vector<ZZ> dlist=alldivs(m,plist);
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

 // test divisor iterator
 cout << "\nTest of divisor iterator class" <<endl;
 vector<ZZ>  Nlist = {ZZ(1), ZZ(65536), ZZ(900), ZZ(666666)};
 for (auto N: Nlist)
   {
     cout << N << ": ";
     if (test_divisor_iterator(N))
       cout << "OK";
     else
       cout << "wrong";
     cout<<endl;
   }
}  // end of main
