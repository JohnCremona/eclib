#include "eclib/xsplit.h"   // which includes linalg.h
#include "eclib/frat.h"

void testINT(INT a)
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
  cout<<"b = a^4-1 = "<<b;
  cout<<", prime factors: "<<pdivs(b);
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
  cout<<"Enter an integer: "<<flush;
  cin >> a;
  cout<<"   value entered: "<<a<<endl;
  if (a.is_long())
    {
      cout<<" as a long int: "<<I2long(a)<<endl;
    }
  else
    {
      cout<<" does not fit in a long int!"<<endl;
    }

  RAT half(1,2);
  cout << "one half = " << half <<", rounds to "<<half.round()<<", floor="<<half.floor()<<", ceil="<<half.ceil()<<", recip="<<half.recip()<<endl;
  RAT third(1,3);
  cout << "one third = " << third <<", rounds to "<<third.round()<<", floor="<<third.floor()<<", ceil="<<third.ceil()<<", recip="<<third.recip()<<endl;

  cout<<"Enter a rational: "<<flush;
  RAT q;
  cin >> q;
  cout<<"   value entered: "<<q<<endl;
}
