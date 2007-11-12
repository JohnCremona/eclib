// modtest.cc: test program for mod functions

#include "marith.h"

int main()
{
  bigint a,b,c,d,mm; 
  cout<<"uninitialized bigint  = "<<mm<<endl;
  long m,m2;
  while(cout << "Enter modulus m: ",  cin >> m, m!=0)
    {
      m2=2*abs(m);
      cout<<"long modulus = "<<m<<endl;
      for(a=-m2; a<m2; a+=1)
	{
	  b = mod(a,m);
	  c = posmod(a,m);
	  d = a%m;
	  cout<<"a = "<<a<<"\t mod(a,m) = "<<b<<"\t posmod(a,m) = "<<c;
	  cout<<"\t a%m = " << d << endl;
	}
      mm=m;
      cout<<"bigint modulus = "<<mm<<endl;
      for(a=-m2; a<m2; a+=1)
	{
	  b = mod(a,mm);
	  c = posmod(a,mm);
	  d = a%mm;
	  cout<<"a = "<<a<<"\t mod(a,m) = "<<b<<"\t posmod(a,m) = "<<c;
	  cout<<"\t a%m = " << d << endl;
	}
    }
}  /* main() */
