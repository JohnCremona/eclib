// leggen.cc: program to construct a list of soluble legendre equations

//                 ax^2+by^2-cz^2
//		   where a, b, c are all primes with 10^k digits
//		   for given k.  We take a, b to be the smallest two primes
//		   over 10^k and search for succesive c

#include "marith.h"
#include "quadratic.h"
#include "conic.h"

int main()
{
  initprimes("PRIMES");

  bigint a,b,c,x,y,z;
  int neq, ndig, n=0;
  cerr<<"How many equations? "; cin>>neq;
  cerr<<"How many digits?    "; cin>>ndig;
  power(a,bigint(10),ndig);
  a=next_prime(a);
  b=next_prime(a);
  c=b;
  cerr<<"a = "<<a<<"\nb = "<<b<<endl;
  while(n<neq) 
    {
      c=next_prime(c);
      if(testlocsol(a,b,-c))
	{
	  n++;
	  cout<<a<<" "<<b<<" "<<(-c)<<"\n";
	}
    }
  cout<<"0 0 0\n";
  abort();
}



