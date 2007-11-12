// makecert.cc: Program to make and output legendre equation certificates

#include "marith.h"
#include "quadratic.h"
#include "conic.h"
#include "legendre.h"

#ifndef CONIC_METHOD
#define CONIC_METHOD 4
#endif
//#define TEST_PARAM

int main()
{
  initprimes("PRIMES");

  bigint a,b,c,x,y,z;
  bigint k1,k2,k3;

  while(1) 
    {
      cerr << "Enter coefficients a b c: ";
      cin >> a >> b >> c;
      cout<<a<<" "<<b<<" "<<c<<endl; 
      if(a==0) abort();
      make_certificate(a,b,c,k1,k2,k3);
      cout<<k1<<" "<<k2<<" "<<k3<<endl;
    }
  abort();
}



