// leggen2.cc: Program to generate solvable legendre equations with prime coefficients.

#include "marith.h"

int main()
{
  int e, n, i;
  cerr<<"Enter number of digits: \n"; cin>>e;  if(e<0) abort();
  cerr<<"Enter number of triples required: \n"; cin>>n; if(n<1) abort();

  bigint a, b, c;
  a=next_prime(pow(bigint(10),e));
  b=next_prime(a);
  c=b;
  i=0;
  while(i<n)
    {
      c=next_prime(c);
      if(jacobi(-a*b%c,c)==1)
	if(jacobi(a*c%b,b)==1)
	  if(jacobi(b*c%a,a)==1)
	    {
	      i++;
	      cout<<a<<" "<<b<<" -"<<c<<endl;
	    }	  
    }
  cout<<"0 0 0"<<endl;
}
