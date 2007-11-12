// testl.cc: test program for LiDIA library

#include <LiDIA/bigint.h>
using namespace LiDIA;

#include <iostream>
using namespace std;

int main()
{
  int i;
  bigint a;
  while(1)
    {
      cout<<"Enter an integer: ";
      cin>>a;
      cout<<"a = "<<a<<endl;
      if(is_int(a)) cout<<"a is NOT an int!"<<endl;
      
      int ok = a.intify(i);
      if(ok) cout<<"as an int, a = "<<i<<endl;
      else cout<<"problem with intify, returns "<<i<<endl;
      
    }
}
