// tc.cc: test program for complex template class

#include <stream.h>
#include <complex.h>
#define Complex complex<double> 

int main()
{
  cout.precision(15);
  Complex a(0,1), b(1,0), c;
  cout<<"a = "<<a<<"\tb = "<<b<<endl;
  c=a+b;
  cout<<"c=a+b = "<<c<<endl;
  c=a-b;
  cout<<"c=a-b = "<<c<<endl;
  c=a*b;
  cout<<"c=a*b = "<<c<<endl;
  c=a/b;
  cout<<"c=a/b = "<<c<<endl;
  c=sqrt(a);
  cout<<"c=sqrt(a) = "<<c<<endl;
}
