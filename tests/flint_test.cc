// flint_test.cc: test of interface to FLINT wraper classes INT and RAT
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

#include "eclib/xsplit.h"   // which includes linalg.h
#include "eclib/frat.h"

void testINT(INT a, int factor=1)
{
  cout<<"a = "<<a;
  cout<<", -a = "<<-a;
  cout<<", a+20 = "<<a+20;
  cout<<", 20+a = "<<20+a;
  cout<<", a-2 = "<<a-2;
  cout<<", 2-a = "<<2-a;
  cout<<", a.abs() = "<<a.abs();
  cout<<", sign(a) = "<<sign(a);
  cout<<endl;
  INT b = (a^8) - 1;
  cout<<"b = a^8-1 = "<<b;
  if (factor)
    cout<<", prime factors: "<<pdivs(b);
  else
    cout<<" (not being factored)" << endl;
  cout<<endl;
}

int main ()
{
  INT a(-13), b(3), c;
  cout << "a = "<<a<<" with sign "<<sign(a)<<endl;
  cout << "b = "<<b<<" with sign "<<b.sign()<<endl;
  cout << "c = "<<c<<" with sign "<<sign(c)<<endl;
  cout << "-a = "<<-a<<" with sign "<<sign(-a)<<endl;
  cout << "a.abs() = "<<a.abs()<<endl;
  cout << "abs(b) = "<<abs(b)<<endl;
  c = a+b;
  cout << "a+b = "<<c<<endl;
  cout << "a-b = "<<a-b<<endl;
  cout << "a*b = "<<a*b<<endl;
  cout << "2+a = "<<2+a<<endl;
  cout << "a+2 = "<<a+2<<endl;
  cout << "2-a = "<<2-a<<endl;
  cout << "a-2 = "<<a-2<<endl;
  cout << "2*a = "<<2*a<<endl;
  cout << "a*2 = "<<a*2<<endl;
  cout << "a%b = "<<a%b<<endl;
  cout << "mod(a,b) = "<<mod(a,b)<<endl;
  cout << "posmod(a,b) = "<<posmod(a,b)<<endl;
  cout << "gcd(a,b) = "<<gcd(a,b)<<endl;
  testINT(a);
  testINT(b);
  INT x, y;
  INT g = bezout(a,b,x,y);
  cout << "bezout(a,b,x,y) = "<<g<<" with x="<<x<<", y="<<y<<endl;
  int e=9;
  INT p = a^e;
  cout << "a^"<<e<<" = "<<p<<endl;
  long f = 3;
  cout << "b^"<<f<<" = "<<(b^f)<<endl;
  cout << "a==b? " << (a==b) << endl;
  cout << "a!=b? " << (a!=b) << endl;
  cout << "a==-13? " << (a==-13) << endl;
  cout << "a!=-13? " << (a!=-13) << endl;
  cout << "a==0? " << (a==0) << endl;
  cout << "a!=0? " << (a!=0) << endl;
  INT n = a*(a+1)*((a^2)+1);
  cout<<"Prime factors of "<<n<<" : "<<pdivs(n) << endl;
  cout<<"Divisors whose square divides "<<n<<" : "<<sqdivs(n) << endl;
  INT pr(7);
  for (long i=0; i<pr; i++)
    cout<<"("<<i<<"|"<<pr<<") = "<<legendre(INT(i),pr)<<endl;
  cout<<"Is "<<a<<" square? "<<a.is_square()<<endl;
  INT a2 = a*a;
  cout<<"Is "<<a2<<" square? "<<a2.is_square()<<endl;
  INT s = a2.isqrt();
  cout<<"(a square root is "<<s<<")"<<endl;

  for (int i=0; i<2; i++)
    {
      cout<<"Enter an integer: "<<flush;
      cin >> a;
      cout<<"   value entered: "<<a;
      if (a.is_long())
        {
          cout<<" - as a long int: "<<I2long(a)<<endl;
        }
      else
        {
          cout<<" - does not fit in a long int"<<endl;
        }
    }
  cout << "Enter three integers separated by whitespace: " << flush;
  cin >> a >> b >> c;
  cout << "Values entered: " << a << " " << b << " " << c << endl;
  cout << "Enter three integers separated by commas: " << flush;
  char ch;
  cin >> a >> ws >> ch >> b >> ch >> c;
  cout << "Values entered: " << a << " " << b << " " << c << endl;
  // Test ofconstruction from a string
  INT big("847538457305748064257802375802345784320765208455280452806202");
  cout << "INT constructed from string: " << big << endl;
  testINT(big, 0); // don't factor big^8-1 

  RAT half(1,2);
  cout << "one half = " << half <<", rounds to "<<half.round()<<", floor="<<half.floor()<<", ceil="<<half.ceil()<<", recip="<<half.recip()<<endl;
  RAT third(1,3);
  cout << "one third = " << third <<", rounds to "<<third.round()<<", floor="<<third.floor()<<", ceil="<<third.ceil()<<", recip="<<third.recip()<<endl;

  cout<<"Enter a rational: "<<flush;
  RAT q;
  cin >> q;
  cout<<"   value entered: "<<q<<endl;
}
