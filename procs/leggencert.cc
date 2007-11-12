// leggencert.cc: program to generate solvable legendre equations with prime coefficients.

#include "marith.h"

int make_certificate(const bigint& a, const bigint& b, const bigint& c, 
		     bigint& n, bigint& p, bigint& q);

int main()
{
  int e, n, i;
  cerr<<"Enter number of digits: \n"; cin>>e;  if(e<0) abort();
  cerr<<"Enter number of triples required: \n"; cin>>n; if(n<1) abort();

  bigint a, b, c;
  bigint k1,k2,k3;
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
	      make_certificate(a,b,-c,k1,k2,k3);
	      cout<<k1<<" "<<k2<<" "<<k3<<endl;
	    }	  
    }
  cout<<"0 0 0"<<endl;
}

int make_certificate(const bigint& a, const bigint& b, const bigint& c, 
		     bigint& n, bigint& p, bigint& q)
{
  if(!sqrt_mod_m(n,-b*c,abs(a))) return 1;
  if(!sqrt_mod_m(p,-a*c,abs(b))) return 2;
  if(!sqrt_mod_m(q,-a*b,abs(c))) return 3;
  return 0;
}

